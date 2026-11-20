"""The time loop of a biosensor (D13 step 5, pilot K8): `bindoptix.time_loop`.

BindOptiX computes the biology over time – transport, capture, binding – and gives the occupancy
B (states per spot over time, `data.occupancy`). The optics is computed once per design: the
meta block sets the receptor state of an emitter block (`emitter.dipole_layer`, parameter
`state`) to each state of B in turn and evaluates the observed output once per state (a camera
image I_s behind the emitter, cached by the executor – only the emitter and what lies behind it
is recomputed). Since emitters are mutually incoherent, the image of spot j at time t is exactly

    I(t, j) = Σ_s f_s(t, j) · I_s,

f the fractions of B (each row sums to 1). States without fluorophores (e.g. "cell" of a capture
run) are named in `dark_states` and contribute nothing; any other state must exist in Q. The
optics is the same for every spot (one emitter block, the camera around one spot); spots with
different optics need one loop each. The result (`TimeSeries`) keeps the per-state images and B;
`image(t, j)` assembles an image, `measures` gives the chosen measure ("power", "peak", "centre")
over time per spot, and `display.plot` draws it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from rtt_fusion._types import FloatArray
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.meta import MEASURES, _collect, measure_of
from rtt_fusion.graph.model import Category, Node, NodeType, ParamSpec, data_in, data_out
from rtt_fusion.graph.registry import NodeTypeRegistry


@dataclass
class TimeSeries:
    """Camera images over time from per-state images and the occupancy B."""

    times_s: FloatArray
    spots: tuple[int, ...]
    states: tuple[str, ...]
    fractions: FloatArray  # (nt, S, n_states)
    per_state: dict[str, Image] = field(default_factory=dict)  # dark states absent
    measure: str = "power"

    def image(self, t: int, spot: int = 0) -> Image:
        """The image of spot index `spot` (into `spots`) at time index t."""
        lit = [s for s in self.states if s in self.per_state]
        first = self.per_state[lit[0]]
        values = np.zeros_like(np.asarray(first.values, dtype=float))
        groups: list[Any] = []
        power = 0.0
        for k, s in enumerate(self.states):
            if s not in self.per_state:
                continue
            f = float(self.fractions[t, spot, k])
            values = values + f * np.asarray(self.per_state[s].values)
            power += f * self.per_state[s].ray_power
            groups += self.per_state[s].groups
        return Image(values, first.extent, power, list(dict.fromkeys(groups)))

    @property
    def measures(self) -> FloatArray:
        """The measure of every image, (nt, S); linear measures straight from the states."""
        nt, ns = self.fractions.shape[:2]
        out = np.zeros((nt, ns))
        for t in range(nt):
            for j in range(ns):
                out[t, j] = measure_of(self.image(t, j), self.measure)
        return out

    def to_curves(self) -> tuple[list[tuple[str, list[float], list[float]]], str, str]:
        m = self.measures
        t = [float(x) for x in self.times_s]
        curves = [(f"spot {s}", t, [float(v) for v in m[:, j]]) for j, s in enumerate(self.spots)]
        return curves, "time (s)", f"image {self.measure}"


def _emitter(ctx: Context, node: Node) -> str:
    name = str(node.params.get("emitter", "") or "")
    if name:
        if name not in ctx.graph.nodes:
            raise ExecutionError("exec.param", f"emitter: no block {name!r}", node.id)
        return name
    found = [n for n, v in ctx.graph.nodes.items() if v.type_id == "emitter.dipole_layer"]
    if len(found) != 1:
        raise ExecutionError(
            "exec.param", f"emitter: name the emitter block ({len(found)} in the graph)", node.id
        )
    return found[0]


def time_loop(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> TimeSeries:
    if ctx.evaluate is None:
        raise ExecutionError("exec.meta", "meta blocks need the executor", node.id)
    occ = inputs.get("occupancy") or None
    if isinstance(occ, list):
        occ = occ[0] if len(occ) == 1 else None
    if occ is None or not hasattr(occ, "fractions"):
        raise ExecutionError("exec.data", "occupancy: connect B (data.occupancy)", node.id)
    measure = str(node.params.get("measure", "power"))
    if measure not in MEASURES:
        raise ExecutionError("exec.param", f"measure {measure!r} ({', '.join(MEASURES)})", node.id)
    dark = {
        s.strip() for s in str(node.params.get("dark_states", "") or "").split(",") if s.strip()
    }
    block, out_socket = _collect(ctx, node)
    emitter = ctx.graph.nodes[_emitter(ctx, node)]
    saved = emitter.params.get("state")
    states = tuple(str(s) for s in occ.states)
    lit = [s for s in states if s not in dark]
    if not lit:
        raise ExecutionError("exec.param", "every state is dark", node.id)
    per_state: dict[str, Image] = {}
    try:
        for i, s in enumerate(lit):
            emitter.params["state"] = s
            ctx.report(i / len(lit), f"{node.name or node.id}: optics of state {s!r}")
            with ctx.span(i / len(lit), (i + 1) / len(lit)):
                result = ctx.evaluate(block, out_socket)
            if hasattr(result, "image") and not isinstance(result, Image):
                result = result.image  # a display's view
            if not isinstance(result, Image):
                raise ExecutionError(
                    "exec.meta", f"the observed output gives a {type(result).__name__}", node.id
                )
            per_state[s] = result
    finally:
        if saved is None:
            emitter.params.pop("state", None)
        else:
            emitter.params["state"] = saved
    shapes = {np.asarray(im.values).shape for im in per_state.values()}
    if len(shapes) != 1:
        raise ExecutionError("exec.meta", "the images of the states differ in size", node.id)
    ctx.log.append(
        f"{node.id}: optics once per state ({', '.join(lit)}), {len(occ.times_s)} times × "
        f"{len(occ.spots)} spots from B"
    )
    return TimeSeries(
        np.asarray(occ.times_s, dtype=float),
        tuple(int(s) for s in occ.spots),
        states,
        np.asarray(occ.fractions, dtype=float),
        per_state,
        measure,
    )


def _key_extra(ctx: Context, node: Node) -> list[Any]:
    try:
        observed: Any = _collect(ctx, node)
    except ExecutionError:
        observed = "observe=?"
    try:
        emitter = _emitter(ctx, node)
    except ExecutionError:
        emitter = "?"
    return [observed, f"emitter={emitter}"]


time_loop.key_extra = _key_extra  # type: ignore[attr-defined]
time_loop.representations = frozenset()  # type: ignore[attr-defined]

TYPES = (
    NodeType(
        "bindoptix.time_loop",
        Category.META,
        (
            data_in("occupancy", "occupancy"),
            data_in("observe", "result", "observe"),
            data_out("results", "result", "results"),
        ),
        (
            ParamSpec(
                "emitter",
                "text",
                "",
                description="emitter block whose state is varied (empty: the only one)",
            ),
            ParamSpec(
                "collect",
                "text",
                "",
                description="observed output block.socket if nothing is linked to observe",
            ),
            ParamSpec(
                "dark_states",
                "text",
                "",
                description="states without fluorophores, comma-separated (e.g. cell)",
            ),
            ParamSpec(
                "measure",
                "choice",
                "power",
                options=MEASURES,
                description="measure of the images over time",
            ),
        ),
        through=(("results", ("occupancy",)),),
        description="Time loop of a biosensor: optics once per receptor state, images over time "
        "from BindOptiX's occupancy B (Σ f_s(t)·I_s)",
    ),
)

EVALUATORS: dict[str, Evaluator] = {"bindoptix.time_loop": time_loop}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    for t in TYPES:
        registry.add(t)
    return dict(EVALUATORS)
