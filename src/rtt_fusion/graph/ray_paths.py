"""Ray paths through a whole workflow, for drawing (CAD: *Trace rays* where a branch passes blocks
of other kinds – a FourierGrating grating or stack, an hp-FEM grating – which a Raytatouille system
alone cannot trace through).

`ray_paths` computes the detectors node-wise with rays only (no fields, the plan R everywhere) and
records every Raytatouille trace on the way (`note`, called by the backend). Blocks of other kinds
move the rays without Raytatouille (a grating to its plane, the orders on from there): the pieces
are joined by ray – a ray keeps its start point in the pupil (`source_xy`) and its bundle the names
of the blocks passed (`path`), so the piece a trace continues is the one of the same coherence
group and pupil point whose bundle's path ends where the new one's begins (the longest such)."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from rtt_fusion.graph.model import Category, Graph


@dataclass(frozen=True)
class Piece:
    """The recorded positions of the rays of one trace: `points` (n, k, 3), NaN after a loss;
    `after`: the block names the bundle had passed when it was traced, `owner` the block that
    traced it (its rays then carry `after + (owner,)`); `power`: the mean |amplitude|² of the
    bundle's live rays (a weak diffraction order is drawn or not)."""

    after: tuple[str, ...]
    owner: str
    field: int
    wavelength: int
    wavelength_um: float
    source_xy: Any
    points: Any
    power: float = 1.0


@dataclass(frozen=True)
class RayPath:
    """Polylines to draw: `points` (n, k, 3) with NaN after a loss (the rays of one piece, its
    predecessor's last point in front), `path` the blocks passed (for colouring by light path),
    `power` that of its bundle (`Piece`)."""

    path: tuple[str, ...]
    field: int
    wavelength_um: float
    points: Any
    power: float = 1.0


_RECORDING: ContextVar[list[Piece] | None] = ContextVar("rtt_fusion_ray_paths", default=None)


@contextmanager
def recording() -> Iterator[list[Piece]]:
    """Records the traces of the backend (`note`) while active, in this thread / task."""
    pieces: list[Piece] = []
    token = _RECORDING.set(pieces)
    try:
        yield pieces
    finally:
        _RECORDING.reset(token)


def active() -> bool:
    return _RECORDING.get() is not None


def note(bundle: Any, points: Any, rows: Any, owner: str) -> None:
    """Called by the backend after a trace: the recorded positions `points` (m, k, 3) of the
    bundle's rays `rows` (indices into the bundle). One polarisation only (the backend traces
    both with the same rays)."""
    pieces = _RECORDING.get()
    if pieces is None or bundle.key.polarisation != 0:
        return
    rows = np.asarray(rows, dtype=int)
    amplitude = np.abs(np.asarray(bundle.amplitude)) ** 2
    live = np.asarray(bundle.alive, dtype=bool) & (amplitude > 0)
    pieces.append(
        Piece(
            tuple(bundle.path),
            owner,
            int(bundle.key.field_index),
            int(bundle.key.wavelength_index),
            float(bundle.wavelength_um),
            np.asarray(bundle.source_xy, dtype=float)[rows],
            np.asarray(points, dtype=float),
            float(amplitude[live].mean()) if live.any() else 0.0,
        )
    )


class _Noting:
    """An evaluator of a block that moves rays without Raytatouille (a grating, a stack, a block
    of another plugin): its output bundles are noted as one point per ray, where the block left
    them – the corner of the drawn path there. Every other attribute is the evaluator's."""

    def __init__(self, evaluator: Any) -> None:
        self._evaluator = evaluator

    def __getattr__(self, name: str) -> Any:
        return getattr(self._evaluator, name)

    def __call__(self, ctx: Any, node: Any, socket: str, inputs: Any) -> Any:
        from rtt_fusion.graph.light import RayBundle

        out = self._evaluator(ctx, node, socket, inputs)
        for b in out if isinstance(out, list) else []:
            if isinstance(b, RayBundle) and b.path:
                alive = np.asarray(b.alive, dtype=bool)
                points = np.where(alive[:, None], np.asarray(b.position, dtype=float), np.nan)
                rows = np.arange(len(alive))
                note(replace(b, path=b.path[:-1]), points[:, None, :], rows, b.path[-1])
        return out


def _last(points: Any) -> Any:
    """The last finite point of each polyline (NaN rows where none is)."""
    ok = np.all(np.isfinite(points), axis=2)
    k = np.where(ok.any(axis=1), ok.shape[1] - 1 - np.argmax(ok[:, ::-1], axis=1), -1)
    out = np.full((len(points), 3), np.nan)
    has = k >= 0
    out[has] = points[np.flatnonzero(has), k[has]]
    return out


def joined(pieces: list[Piece]) -> list[RayPath]:
    """The pieces as polylines, each ray of a piece starting at the last point of the same ray in
    its predecessor (`Piece`, module docstring); a ray without one starts where it was traced."""
    out: list[RayPath] = []
    ends = [_last(p.points) for p in pieces]
    for p in pieces:
        best, best_len = None, -1
        for j, q in enumerate(pieces):
            reach = (*q.after, q.owner) if q.owner else None
            if (
                reach is None
                or q is p
                or (q.field, q.wavelength) != (p.field, p.wavelength)
                or len(reach) > len(p.after)
                or p.after[: len(reach)] != reach
                or len(reach) <= best_len
            ):
                continue
            best, best_len = j, len(reach)
        points = p.points
        if best is not None:
            q = pieces[best]
            index = {tuple(np.round(xy, 12)): i for i, xy in enumerate(q.source_xy)}
            front = np.full((len(points), 1, 3), np.nan)
            for i, xy in enumerate(p.source_xy):
                k = index.get(tuple(np.round(xy, 12)))
                if k is not None:
                    front[i, 0] = ends[best][k]
            points = np.concatenate([front, points], axis=1)
        if points.shape[1] >= 2:
            out.append(RayPath((*p.after, p.owner), p.field, p.wavelength_um, points, p.power))
    return out


def ray_paths(
    graph: Graph,
    registry: Any,
    evaluators: Mapping[str, Any],
    shared: Mapping[str, Any],
    detectors: list[str] | None = None,
    min_power: float = 0.0,
) -> tuple[list[RayPath], list[str]]:
    """The rays of the workflow as polylines, from the sources to the detectors (all, or the ones
    named): the detectors computed node-wise with rays only (every link on R, no fields near a
    focus) and the Raytatouille traces on the way joined (`joined`). Returns the paths and the
    executor's log. `shared` as for a run (rings, pixels, materials, ...); fewer rings draw
    faster. `min_power`: paths whose bundle has less than this share of the power of the strongest
    one through the same block are left out (weak diffraction orders)."""
    from rtt_fusion.graph.execute import ExecutionError, Executor

    g = replace(graph, links=[replace(lk, transport="auto") for lk in graph.links])
    noting = {
        t: e if t.startswith("rtt.") or not callable(e) else _Noting(e)
        for t, e in evaluators.items()
    }
    ex = Executor(
        g,
        registry,
        noting,
        {**shared, "auto_fields": False, "rtt_execution": "nodes"},
        available=frozenset("R"),
    )
    if detectors is None:
        detectors = [
            n.id
            for n in g.nodes.values()
            if (t := registry.get(n.type_id)) is not None and t.category == Category.DETECTOR
        ]
    log: list[str] = []
    with recording() as pieces:
        for d in detectors:
            node = g.nodes[d]
            outs = [
                s.name
                for s in registry.get(node.type_id).sockets
                if s.direction.value == "out" and s.kind.value == "data"
            ]
            for socket in outs[:1]:
                try:
                    ex.evaluate(d, socket)
                except ExecutionError as e:
                    log.append(f"{d}: {e}")
    log += list(ex.context.log)
    paths = joined(pieces)
    strongest: dict[str, float] = {}  # per block: the amplitudes hold the ray density, too
    for p in paths:
        strongest[p.path[-1]] = max(strongest.get(p.path[-1], 0.0), p.power)
    return [p for p in paths if p.power >= min_power * strongest[p.path[-1]]], log


__all__ = ["Piece", "RayPath", "active", "joined", "note", "ray_paths", "recording"]
