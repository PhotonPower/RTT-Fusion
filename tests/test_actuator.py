"""Actuators (`fusion.actuator`): one value sets a parameter of the Raytatouille system and moves
the blocks the stage carries; the executor applies them before a run, a meta block after each
value it sets; the result is the same in one process and behind a provider."""

from __future__ import annotations

import json
import math
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import Graph, Node, builtin_registry, optics
from rtt_fusion.graph.actuator import Move, apply_actuators, parse_moves
from rtt_fusion.graph.execute import ExecutionError, Executor
from rtt_fusion.graph.remote import BlockProvider, register_remote


def _source(**rows: Any) -> Node:
    parameters = [{"name": k, **v} for k, v in rows.items()]
    return Node("source", "rtt.source", params={"rtt_system": {"parameters": parameters}})


def _actuator(**params: Any) -> Node:
    return Node("stage", "fusion.actuator", params=params)


def test_moves_are_parsed() -> None:
    assert parse_moves("G1.dx, G2.dx; E1.rz*0.5\nE2.dz * -1") == [
        Move("G1", "dx"),
        Move("G2", "dx"),
        Move("E1", "rz", 0.5),
        Move("E2", "dz", -1.0),
    ]
    assert parse_moves("") == []
    for bad in ("G1", "G1.dw", "G1.dx*", "G1.dx*a"):
        with pytest.raises(ExecutionError, match="moves"):
            parse_moves(bad, "stage")


def test_a_stage_sets_the_parameter_and_carries_the_blocks() -> None:
    g = Graph()
    g.add(_source(DX={"value": 0.0, "min": -0.4, "max": 0.4}, Z={"expression": "80 + DX"}))
    g.add(Node("G1", "fg.grating", params={"point": [0, 0, 80], "rtt": {"/root/children/2": {}}}))
    g.add(Node("E1", "fg.stack", params={"point": [0, 0, 79]}))
    g.add(Node("G1@2", "fg.grating", params={"point": [0, 0, 80], "rtt": {}}))  # a copy of G1
    g.add(_actuator(value=0.25, parameter="DX", moves="G1.dx, E1.dx*2, G1.rz"))
    model_rows = g.nodes["source"].params["rtt_system"]["parameters"]
    apply_actuators(g)
    rows = {r["name"]: r for r in g.nodes["source"].params["rtt_system"]["parameters"]}
    assert rows["DX"]["value"] == 0.25 and rows["Z"] == {"name": "Z", "expression": "80 + DX"}
    assert model_rows[0]["value"] == 0.0  # the rows of the model are not changed
    assert g.nodes["G1"].params["pose_by_stage"] == [0.25, 0, 0, 0, 0, 0.25]
    assert g.nodes["G1@2"].params["pose_by_stage"] == [0.25, 0, 0, 0, 0, 0.25]
    assert g.nodes["E1"].params["pose_by_stage"] == [0.5, 0, 0, 0, 0, 0]
    assert g.nodes["G1"].params["keep_pose"] is True  # carried, not placed by the system
    assert "keep_pose" not in g.nodes["E1"].params  # no replaced element: nothing to keep

    # zero, drive, degrees
    g.nodes["stage"].params.update({"value": 45.0, "drive": 5.0, "zero": 40.0, "unit": "deg"})
    g.nodes["stage"].params.update({"parameter": "", "moves": "G1.rz"})
    apply_actuators(g)
    assert g.nodes["G1"].params["pose_by_stage"] == pytest.approx([0, 0, 0, 0, 0, math.radians(10)])
    assert "pose_by_stage" not in g.nodes["E1"].params  # no longer carried

    # muted: nothing moved
    g.nodes["stage"].muted = True
    apply_actuators(g)
    assert all("pose_by_stage" not in n.params for n in g.nodes.values())


@pytest.mark.parametrize(
    ("params", "message"),
    [
        ({"value": 0.5, "parameter": "DX"}, "outside its range"),
        ({"value": 1.0, "parameter": "Z"}, "is an expression"),
        ({"value": 1.0, "parameter": "DY"}, "no parameter 'DY'.*free parameters: DX"),
        ({"value": 1.0, "moves": "G9.dx"}, "no block 'G9'"),
    ],
)
def test_what_a_stage_cannot_do(params: dict[str, Any], message: str) -> None:
    g = Graph()
    g.add(_source(DX={"value": 0.0, "min": -0.4, "max": 0.4}, Z={"expression": "80 + DX"}))
    g.add(_actuator(**params))
    with pytest.raises(ExecutionError, match=message):
        apply_actuators(g)


def _lens_graph(stage: dict[str, Any] | None, offset: list[float] | None = None) -> Graph:
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": 1.0, "rings": 8, "center": [0, 0, -1]},
        )
    )
    lens = {"point": [0, 0, 0], "focal_length": 40.0}
    if offset is not None:
        lens["pose_offset"] = offset
    g.add(Node("lens", "optics.thin_lens", params=lens))
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, 39.0],
                "normal": [0, 0, -1],
                "half_size": [0.2, 0.2],
                "pixels": [21, 21],
            },
        )
    )
    g.connect(("src", "out"), ("lens", "in"))
    g.connect(("lens", "out"), ("cam", "in"))
    if stage is not None:
        g.add(Node("stage", "fusion.actuator", params=stage))
    return g


def _registry() -> tuple[Any, dict[str, Any]]:
    from rtt_fusion.graph import actuator, meta

    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(meta.register(reg))
    ev.update(actuator.register(reg))
    return reg, ev


def test_the_executor_applies_a_stage_like_a_pose_offset() -> None:
    reg, ev = _registry()
    shared = {"pixels": [21, 21]}
    moved = Executor(_lens_graph({"value": 0.05, "moves": "lens.dx"}), reg, ev, shared)
    offset = Executor(_lens_graph(None, [0.05, 0, 0, 0, 0, 0]), reg, ev, shared)
    still = Executor(_lens_graph(None), reg, ev, shared)
    a, b, c = (e.evaluate("cam", "image").values for e in (moved, offset, still))
    np.testing.assert_array_equal(a, b)
    assert np.abs(a - c).max() > 0.1 * c.max()  # the lens moved the focus


def test_a_sweep_drives_a_stage() -> None:
    """meta.sweep → the stage's drive: each value moves the lens, the sweep restores it."""
    reg, ev = _registry()
    g = _lens_graph({"value": 0.0, "moves": "lens.dx"})
    g.add(Node("sweep", "meta.sweep", params={"measure": "centre", "values": [0.0, 0.05, 0.1]}))
    g.connect(("sweep", "drive"), ("stage", "drive"))
    g.connect(("cam", "image"), ("sweep", "observe"))
    ex = Executor(g, reg, ev, {"pixels": [21, 21]})
    measures = ex.evaluate("sweep", "results").measures
    assert measures[0] > 2 * measures[2], measures  # the focus leaves the centre
    assert g.nodes["lens"].params["pose_by_stage"] == [0, 0, 0, 0, 0, 0]  # back at the start
    assert g.nodes["stage"].params.get("drive", 0.0) == 0.0


def _michelson() -> dict[str, Any]:
    """The Michelson with the reference mirror at Z_REF = 100 + DZ (a free parameter)."""
    d: dict[str, Any] = json.loads(
        (Path(__file__).parent / "data" / "michelson.rtt.json").read_text("utf-8")
    )
    d["aperture"] = {"type": "epd", "value": 1.0}
    d["parameters"] = [
        {"name": "DZ", "value": 0.0, "variable": True, "min": -0.01, "max": 0.01},
        {"name": "Z_REF", "expression": "100 + DZ"},
    ]
    mirror = next(c for c in d["root"]["children"] if c.get("name") == "reference mirror")
    mirror["pose"]["position"][2] = {"param": "Z_REF"}
    return d


def _piezo_sweep(remote: bool, tmp_path: Path) -> list[float]:
    """A sweep of the stage that sets DZ; the rtt.* blocks in the engine or behind a provider."""
    pytest.importorskip("raytatouille")
    from rtt_fusion.graph import backend_rtt, model_from_raytatouille, unroll
    from rtt_fusion.graph.workflow import registry

    reg, extra = registry()
    ev = {**backend_rtt.EVALUATORS, **extra}
    if remote:
        types = [t for t in reg.types() if t.type_id in backend_rtt.EVALUATORS]
        provider = BlockProvider(types, dict(backend_rtt.EVALUATORS), lambda g, plain: {})
        catalogue = [{"provider": "raytracer", "type": t} for t in provider.catalogue()]

        def call(name: str, request: dict[str, Any]) -> dict[str, Any]:
            wire = json.loads(json.dumps(request, allow_nan=False))
            answer: dict[str, Any] = json.loads(json.dumps(provider.handle(wire), allow_nan=False))
            return answer

        ev.update(register_remote(reg, catalogue, call, tmp_path / "portdata"))
    g = unroll(model_from_raytatouille(_michelson()), reg)
    g.nodes["camera"].params["half_size"] = [0.3, 0.3]
    g.add(Node("piezo", "fusion.actuator", params={"parameter": "DZ"}))
    values = [k * 0.6328e-3 / 16 for k in range(9)]  # λ/8 path difference per step
    g.add(Node("ps", "meta.sweep", params={"measure": "centre", "values": values}))
    g.connect(("ps", "drive"), ("piezo", "drive"))
    g.connect(("camera", "image"), ("ps", "observe"))
    ex = Executor(g, reg, ev, {"rings": 10, "pixels": [16, 16]})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        measures: list[float] = ex.evaluate("ps", "results").measures
    rows = {r["name"]: r for r in g.nodes["source"].params["rtt_system"]["parameters"]}
    assert rows["DZ"]["value"] == 0.0  # the sweep put the stage back
    return measures


def test_a_stage_on_a_system_parameter_in_one_process_and_behind_a_provider(
    tmp_path: Path,
) -> None:
    """The Michelson's reference mirror on a piezo stage (system parameter DZ): the fringe
    signal goes through a period over λ/2 of travel, the same with the rtt.* blocks behind a
    provider (the changed parameter reaches it with every evaluation)."""
    local = _piezo_sweep(False, tmp_path)
    remote = _piezo_sweep(True, tmp_path)
    np.testing.assert_allclose(remote, local, rtol=0, atol=1e-12 * max(local))
    assert max(local) - min(local) > 0.5 * max(local)  # the signal changes
    assert local[0] == pytest.approx(local[8], rel=0.05)  # λ/2 of travel: one period
