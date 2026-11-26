"""Ray paths for drawing (`graph.ray_paths`): the Raytatouille traces of a node-wise run, joined
across blocks of other kinds (a grating sample) by ray – the drawn rays go on where a FourierGrating
block sent them (CAD: *Trace rays* through a workflow with gratings)."""

from __future__ import annotations

import json
import warnings
from dataclasses import replace
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.special import jv

from rtt_fusion.graph import Graph, Node, builtin_registry, graph_from_raytatouille, grating
from rtt_fusion.graph.light import CoherenceKey, RayBundle
from rtt_fusion.graph.ray_paths import joined, note, ray_paths, recording

SYSTEM = Path(__file__).parent / "data" / "linnik_mini.rtt.json"


def _bundle(path: tuple[str, ...], xy: list[list[float]], field: int = 0) -> RayBundle:
    n = len(xy)
    return RayBundle(
        key=CoherenceKey("s", field, 0, 0),
        wavelength_um=0.6,
        position=np.zeros((n, 3)),
        direction=np.tile([0.0, 0.0, 1.0], (n, 1)),
        opl=np.zeros(n),
        amplitude=np.ones(n, dtype=complex),
        source_xy=np.asarray(xy, dtype=float),
        alive=np.ones(n, dtype=bool),
        path=path,
    )


def _line(*zs: float) -> Any:
    return np.array([[[0.0, 0.0, z] for z in zs]])


def test_pieces_are_joined_by_ray_and_by_the_block_they_continue() -> None:
    with recording() as pieces:
        note(_bundle((), [[0, 0]]), _line(0.0), [0], "br")  # the source's start point
        note(_bundle(("br",), [[0, 0]]), _line(0.0, 5.0), [0], "A")  # a Raytatouille block
        note(_bundle(("br",), [[0, 0]]), _line(0.0, 7.0), [0], "A2")  # another block, same input
        note(_bundle(("br", "A"), [[0, 0]]), _line(9.0), [0], "G[1]")  # a grating's order 1
        note(_bundle(("br", "A", "G[1]"), [[0, 0]]), _line(9.0, 20.0), [0], "B")
        note(_bundle(("br", "A", "G[1]"), [[1, 0]]), _line(9.0, 20.0), [0], "B")  # another ray
        note(_bundle(("br", "A", "G[1]"), [[0, 0]], field=1), _line(9.0, 30.0), [0], "C")
        other = replace(_bundle(("br",), [[0, 0]]), key=CoherenceKey("s", 0, 0, 1))
        note(other, _line(0.0, 1.0), [0], "A")
    assert len(pieces) == 7  # the second polarisation is the same rays: not noted
    made = joined(pieces)
    by = {p.path: p.points[0, :, 2] for p in made if p.path[-1] != "B"}
    # the grating's point joins the end of A (not A2, which the bundle did not pass), B that point
    np.testing.assert_array_equal(by[("br", "A", "G[1]")], [5.0, 9.0])
    np.testing.assert_array_equal(by[("br", "A")], [0.0, 0.0, 5.0])
    first, second = (p.points[0, :, 2] for p in made if p.path[-1] == "B")
    np.testing.assert_array_equal(first, [9.0, 9.0, 20.0])
    # a ray the predecessor does not have: no point in front; another field point: no predecessor
    assert np.isnan(second[0])
    np.testing.assert_array_equal(by[("br", "A", "G[1]", "C")], [9.0, 30.0])


def _graph(table: str | None) -> tuple[Graph, Any, dict[str, Any], dict[str, Any]]:
    from rtt_fusion.graph import backend_rtt, replace_block

    rt = pytest.importorskip("raytatouille")
    d = json.loads(SYSTEM.read_text(encoding="utf-8"))
    g: Graph = graph_from_raytatouille(d)
    reg = builtin_registry()
    ev = dict(backend_rtt.EVALUATORS)
    ev.update(grating.register(reg))
    if table is not None:
        replace_block(g, reg, "sample", Node("sample", "fg.grating", params={"table": table}))
    shared = {"system": rt.System.from_json(json.dumps(d)), "rings": 2, "pixels": [16, 8]}
    return g, reg, ev, shared


def _points(paths: list[Any]) -> set[tuple[float, ...]]:
    out = set()
    for p in paths:
        for q in p.points.reshape(-1, 3):
            if np.isfinite(q).all():
                out.add(tuple(np.round(q, 7)))
    return out


def test_rays_through_a_grating_block_are_drawn_as_through_the_mirror(tmp_path: Path) -> None:
    """A "mirror table" in place of the Linnik's sample mirror: the drawn rays have the same
    corners as with the mirror, and every ray at the camera goes back to the source unbroken."""
    from test_grating_rtt import write_table

    table = write_table(tmp_path / "mirror.json", {0: -1.0}, 10.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mirror, log_m = ray_paths(*_graph(None))
        split, log_s = ray_paths(*_graph(table))
    assert mirror and split, (log_m, log_s)
    assert _points(split) == _points(mirror)
    assert {p.path[-1] for p in split} >= {"sample[0]", "camera"}
    by_path = {p.path: p for p in split}
    for p in split:
        if p.path[-1] != "camera":
            continue
        # walk back to the source: each piece's rays start where its predecessor's ended
        chain, here = [p], p.path[:-1]
        while here:
            prev = by_path.get(here)
            if prev is None:
                here = here[:-1]
                continue
            chain.append(prev)
            here = prev.path[:-1]
        assert len(chain[-1].path) == 2  # the first block, from the source's start points
        reached = np.isfinite(p.points[:, -1]).all(axis=1)
        assert reached.any()
        for later, earlier in pairwise(chain):
            np.testing.assert_allclose(later.points[reached, 0], earlier.points[reached, -1])


def test_weak_orders_are_left_out(tmp_path: Path) -> None:
    from test_grating_rtt import write_table

    amplitudes = {m: -float(jv(m, 1.0)) for m in (-1, 0, 1)}  # J1²/J0² = 0.33
    table = write_table(tmp_path / "phase.json", amplitudes, 10.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        every, _ = ray_paths(*_graph(table))
        strong, _ = ray_paths(*_graph(table), min_power=0.5)

    def orders(paths: list[Any]) -> list[str]:
        at = [p.path for p in paths if p.path[-1] == "camera"]
        return sorted(o for path in at for o in path if o.startswith("sample["))

    assert orders(every) == ["sample[-1]", "sample[0]", "sample[1]"]
    assert orders(strong) == ["sample[0]"]
