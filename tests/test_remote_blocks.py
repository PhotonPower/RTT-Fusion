"""Blocks evaluated in another process (CAD ADR 0010, D2/D3, W2b): the provider's types travel as
plain JSON, a remote block gives the same result as in the engine's process, bit for bit."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import Graph, Node, builtin_registry, optics
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.model import NodeType
from rtt_fusion.graph.remote import BlockProvider, RemoteEvaluator, register_remote


def _graph() -> Graph:
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": 1.0, "rings": 8, "center": [0, 0, -1]},
        )
    )
    g.add(Node("lens", "optics.thin_lens", params={"point": [0, 0, 0], "focal_length": 40.0}))
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, 39.0],
                "normal": [0, 0, -1],
                "half_size": [1.5, 1.5],
                "pixels": [32, 32],
            },
        )
    )
    g.connect(("src", "out"), ("lens", "in"))
    g.connect(("lens", "out"), ("cam", "in"))
    return g


def _types_and_evaluators() -> tuple[Any, dict[str, Any]]:
    reg = builtin_registry()
    return reg, optics.register(reg)


def test_w69_types_travel_as_plain_json() -> None:
    reg, _ = _types_and_evaluators()
    for t in reg.types():
        back = NodeType.from_dict(json.loads(json.dumps(t.to_dict())))
        assert back == dataclasses.replace(t, port_frames=None)


def test_w69_a_remote_block_equals_the_local_one(tmp_path: Path) -> None:
    reg, local = _types_and_evaluators()
    image: Image = Executor(_graph(), reg, local).evaluate("cam", "image")
    # the provider: a plugin process with the thin lens
    provider = BlockProvider(
        [reg.get("optics.thin_lens")], {"optics.thin_lens": local["optics.thin_lens"]}
    )
    calls: list[dict[str, Any]] = []

    def call(name: str, request: dict[str, Any]) -> dict[str, Any]:
        assert name == "lenses"
        wire = json.loads(json.dumps(request, allow_nan=False))  # through CAD's channel
        calls.append(wire)
        answer: dict[str, Any] = json.loads(json.dumps(provider.evaluate(wire), allow_nan=False))
        return answer

    engine_reg, engine_local = _types_and_evaluators()
    catalogue = [{"provider": "lenses", "type": entry} for entry in provider.catalogue()]
    remote = register_remote(engine_reg, catalogue, call, tmp_path / "portdata")
    assert isinstance(remote["optics.thin_lens"], RemoteEvaluator)
    evaluators = {**engine_local, **remote}
    again: Image = Executor(
        _graph(), engine_reg, evaluators, {"rings": 8, "materials": object()}
    ).evaluate("cam", "image")
    assert len(calls) == 1 and calls[0]["shared"] == {"rings": 8}  # only plain settings travel
    assert again.values.tobytes() == image.values.tobytes() and again.ray_power == image.ray_power
    assert not any((tmp_path / "portdata").iterdir())  # the files are gone after the request


def test_w69_errors_of_the_provider(tmp_path: Path) -> None:
    reg, local = _types_and_evaluators()

    def broken(ctx: Any, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
        raise ValueError("no glass")

    provider = BlockProvider([reg.get("optics.thin_lens")], {"optics.thin_lens": broken})
    remote = RemoteEvaluator("lenses", lambda _p, r: provider.evaluate(r), frozenset("R"), tmp_path)
    with pytest.raises(Exception, match="lenses: ValueError: no glass"):
        Executor(_graph(), reg, {**local, "optics.thin_lens": remote}).evaluate("cam", "image")
    with pytest.raises(ValueError, match="no evaluator"):
        BlockProvider([reg.get("optics.thin_lens")], {})
    unknown = provider.evaluate(
        {
            "node": {"id": "x", "type": "optics.mirror"},
            "socket": "out",
            "inputs": {},
            "folders": {},
            "out": str(tmp_path / "o"),
        }
    )
    assert not unknown["ok"] and "no block type" in unknown["error"]


# ---------------------------------------------------------------------------- W71


@pytest.mark.parametrize("execution", ["nodes", "paths"])
def test_w71_raytatouille_blocks_in_another_process(tmp_path: Path, execution: str) -> None:
    """ADR 0010, D6/W3: every rtt.* block of the Linnik evaluated by a provider behind a JSON
    channel, the run prepared once (its graph: the system from the blocks' data, the branches);
    the camera image equals the one computed in the engine's process."""
    pytest.importorskip("raytatouille")
    import time

    from rtt_fusion.graph import backend_rtt, model_from_raytatouille, unroll
    from rtt_fusion.graph.workflow import registry

    d = json.loads((Path(__file__).parent / "data" / "linnik_mini.rtt.json").read_text("utf-8"))
    reg, extra = registry()
    shared = {"rings": 12, "pixels": [32, 32], "rtt_execution": execution}
    rtt_types = [t for t in reg.types() if t.type_id in backend_rtt.EVALUATORS]
    start = time.perf_counter()
    local: Image = Executor(
        unroll(model_from_raytatouille(d), reg), reg, {**backend_rtt.EVALUATORS, **extra}, shared
    ).evaluate("camera", "image")
    t_local = time.perf_counter() - start
    provider = BlockProvider(
        rtt_types, dict(backend_rtt.EVALUATORS), prepare=lambda graph, plain: {"prepared": True}
    )
    kinds: list[str] = []

    def call(name: str, request: dict[str, Any]) -> dict[str, Any]:
        wire = json.loads(json.dumps(request, allow_nan=False))
        kinds.append(str(wire.get("kind", "evaluate")))
        answer: dict[str, Any] = json.loads(json.dumps(provider.handle(wire), allow_nan=False))
        return answer

    engine_reg, engine_extra = registry()
    catalogue = [{"provider": "raytracer", "type": e} for e in provider.catalogue()]
    remote = register_remote(engine_reg, catalogue, call, tmp_path / "portdata")
    assert set(remote) == {t.type_id for t in rtt_types}
    start = time.perf_counter()
    again: Image = Executor(
        unroll(model_from_raytatouille(d), engine_reg),
        engine_reg,
        {**engine_extra, **remote},
        shared,
    ).evaluate("camera", "image")
    t_remote = time.perf_counter() - start
    assert kinds.count("prepare") == 1 and kinds.count("evaluate") >= 2  # prepared once per run
    assert local.values.max() > 0
    np.testing.assert_allclose(again.values, local.values, rtol=0, atol=1e-12 * local.values.max())
    n = kinds.count("evaluate")
    print(f"W71 {execution}: local {t_local:.3f} s, remote {t_remote:.3f} s, {n} requests")
