"""The time loop of a biosensor (D13 step 5, pilot K8): `bindoptix.time_loop`.

References (W98), independent of the loop's assembly:
- (a) incoherent sum: the camera image of one emitter block holding the emitters of both states
  together (state "mix") is the detector's own incoherent sum over all of them; at
  f = (½, ½) the loop's image is half of it, at f = (1, 0) the image of the "free" emitters;
- (b) the optics is computed once per lit state, however many times and spots B has (the
  executor's call counts);
- (c) dark states contribute nothing and need no ensemble in Q; a lit state missing in Q is an
  error;
- (d) with BindOptiX: the chain biooptomd.layer → emitter.dipole_layer → optics.ideal_imaging →
  optics.detector, B from `data.occupancy`, measure over time = Σ_s f_s(t) P_s.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion._types import FloatArray
from rtt_fusion.adapters import biooptomd as bmd
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import ExecutionError, Executor
from rtt_fusion.graph.model import Category, NodeType, data_out

pytest.importorskip("fieldtracing")

IMG = [0.0, 0.0, 50.0]


@dataclass
class _Occupancy:
    """B as BindOptiX's Occupancy carries it (times, spots, states, fractions)."""

    times_s: FloatArray
    spots: tuple[int, ...]
    states: tuple[str, ...]
    fractions: FloatArray


def _series() -> dict[str, FloatArray]:
    """Two states with the same mean lateral position (the emitters' anchor)."""
    free = np.array(
        [[[-8.0, 0.0, 0.0, 1.0, 0.0, 0.0, 6.0], [8.0, 0.0, 0.0, 0.0, 1.0, 0.0, 6.0]]]
    )  # in-plane dipoles, 6 nm
    bound = np.array(
        [[[0.0, -8.0, 0.0, 0.0, 0.0, 1.0, 3.0], [0.0, 8.0, 0.0, 0.6, 0.0, 0.8, 3.0]]]
    )  # upright and tilted, 3 nm
    return {"free": free, "bound": bound, "mix": np.concatenate([free, bound], axis=1)}


def _graph(
    tmp: Path, occupancy: Any, dark: str = "", occupancy_file: str = ""
) -> tuple[Graph, Any, Any]:
    reg, ev = workflow.registry()
    reg.add(NodeType("test.occupancy", Category.DATA, (data_out("occupancy", "occupancy"),)))
    ev = {**ev, "test.occupancy": lambda ctx, node, socket, inputs: occupancy}
    bmd.write_layer(tmp, dipoles=bmd.DipoleSeries(_series()))
    g = Graph()
    g.add(Node("q", "biooptomd.layer", params={"directory": str(tmp)}))
    g.add(
        Node(
            "e",
            "emitter.dipole_layer",
            params={
                "wavelength_um": 0.6,
                "n_substrate": 1.52,
                "n_solution": 1.333,
                "na": 1.3,
                "samples": 10,
                "n_theta": 24,
                "n_phi": 24,
                "distance": 0.005,
            },
        )
    )
    g.add(
        Node(
            "m",
            "optics.ideal_imaging",
            params={"point": [0, 0, 1e-6], "image_point": IMG, "magnification": 40.0, "na": 1.3},
        )
    )
    g.add(
        Node(
            "d",
            "optics.detector",
            params={
                "point": IMG,
                "normal": [0, 0, -1],
                "half_size": [0.03, 0.03],
                "pixels": [32, 32],
            },
        )
    )
    if occupancy_file:
        g.add(Node("b", "data.occupancy", params={"file": occupancy_file}))
    else:
        g.add(Node("b", "test.occupancy"))
    g.add(Node("t", "bindoptix.time_loop", params={"dark_states": dark}))
    g.connect(("q", "dipoles"), ("e", "dipoles"))
    g.connect(("e", "out"), ("m", "in"))
    g.connect(("m", "out"), ("d", "in"))
    g.connect(("b", "occupancy"), ("t", "occupancy"))
    g.connect(("d", "image"), ("t", "observe"))
    for link in g.links:
        if link.source[1] == "out":
            g.set_transport(link, "S")
    return g, reg, ev


def _occupancy(fractions: list[list[list[float]]], states: tuple[str, ...]) -> _Occupancy:
    f = np.asarray(fractions, dtype=float)
    return _Occupancy(np.arange(f.shape[0]) * 60.0, tuple(range(f.shape[1])), states, f)


def test_w98_images_over_time_are_the_incoherent_sum(tmp_path: Path) -> None:
    occ = _occupancy([[[1.0, 0.0], [0.5, 0.5]], [[0.5, 0.5], [0.0, 1.0]]], ("free", "bound"))
    g, reg, ev = _graph(tmp_path, occ)
    ex = Executor(g, reg, ev, available=workflow.available())
    series = ex.evaluate("t", "results")
    # the reference: one emitter block with all emitters, one detector evaluation
    g.nodes["e"].params["state"] = "mix"
    mix = Executor(g, reg, ev, available=workflow.available()).evaluate("d", "image")
    g.nodes["e"].params["state"] = "free"
    free = Executor(g, reg, ev, available=workflow.available()).evaluate("d", "image")
    half = np.asarray(series.image(0, 1).values)
    assert np.abs(half - 0.5 * np.asarray(mix.values)).max() < 1e-9 * np.asarray(mix.values).max()
    assert np.abs(np.asarray(series.image(0, 0).values) - np.asarray(free.values)).max() < (
        1e-12 * np.asarray(free.values).max()
    )
    # (b) optics once per state, not per time or spot
    assert ex.calls[("e", "out")] == 2 and ex.calls[("d", "image")] == 2
    assert "state" not in g.nodes["e"].params or g.nodes["e"].params["state"] == "free"


def test_w98_dark_states_and_missing_states(tmp_path: Path) -> None:
    occ = _occupancy([[[0.2, 0.3, 0.5]], [[0.0, 0.4, 0.6]]], ("free", "bound", "cell"))
    g, reg, ev = _graph(tmp_path / "a", occ, dark="cell")
    series = Executor(g, reg, ev, available=workflow.available()).evaluate("t", "results")
    p = {s: series.per_state[s].power for s in ("free", "bound")}
    want = [0.2 * p["free"] + 0.3 * p["bound"], 0.4 * p["bound"]]
    assert series.measures[:, 0] == pytest.approx(want, rel=1e-12)
    curves, x_label, _ = series.to_curves()
    assert x_label == "time (s)" and curves[0][2] == pytest.approx(want, rel=1e-12)
    g, reg, ev = _graph(tmp_path / "b", occ)  # "cell" not dark and not in Q
    with pytest.raises(ExecutionError) as e:
        Executor(g, reg, ev, available=workflow.available()).evaluate("t", "results")
    assert e.value.code == "exec.data" and "cell" in str(e.value)


def test_w98_k8_chain_with_bindoptix_occupancy(tmp_path: Path) -> None:
    pytest.importorskip("bindoptix")
    import bindoptix.binding

    f = np.array([[[1.0, 0.0]], [[0.7, 0.3]], [[0.1, 0.9]]])
    occ = bindoptix.binding.Occupancy(np.array([0.0, 600.0, 1200.0]), (4,), ("free", "bound"), f)
    (tmp_path / "b.json").write_text(occ.to_json(), encoding="utf-8")
    g, reg, ev = _graph(tmp_path / "run", None, occupancy_file=str(tmp_path / "b.json"))
    series = Executor(g, reg, ev, available=workflow.available()).evaluate("t", "results")
    assert series.spots == (4,)
    p = np.array([series.per_state["free"].power, series.per_state["bound"].power])
    assert series.measures[:, 0] == pytest.approx(f[:, 0, :] @ p, rel=1e-12)
