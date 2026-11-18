"""The conversion block fusion.convert: R <-> H in a graph, later representations reported."""

from __future__ import annotations

import numpy as np
import pytest

from typing import Any

from rtt_fusion.graph import Graph, Node, NodeTypeRegistry, builtin_registry, conversion, optics
from rtt_fusion.graph.execute import ExecutionError, Executor
from rtt_fusion.graph.light import RayBundle
from rtt_fusion.graph.workflow import registry

pytest.importorskip("fieldtracing")
from rtt_fusion.graph.field import FieldBundle  # noqa: E402


def _graph(to: str) -> tuple[Graph, NodeTypeRegistry, dict[str, Any]]:
    reg = builtin_registry()
    ev = {**optics.register(reg), **conversion.register(reg)}
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": 0.25, "rings": 12, "gaussian_waist": 0.1},
        )
    )
    g.add(Node("conv", "fusion.convert", params={"to": to}))
    g.connect(("src", "out"), ("conv", "in"))
    return g, reg, ev


def test_the_workflow_registry_has_an_evaluator() -> None:
    _, extra = registry()
    assert "fusion.convert" in extra


def test_rays_to_field_and_back() -> None:
    g, reg, ev = _graph("H")
    fields = Executor(g, reg, ev).evaluate("conv", "out")
    assert fields and all(isinstance(f, FieldBundle) for f in fields)
    g2, reg2, ev2 = _graph("R")
    rays = Executor(g2, reg2, ev2).evaluate("conv", "out")
    assert rays and all(isinstance(r, RayBundle) for r in rays)
    # R -> R passes the light unchanged
    src = Executor(g2, reg2, ev2).evaluate("src", "out")
    assert np.allclose(np.asarray(rays[0].alive), np.asarray(src[0].alive))


def test_later_representations_are_reported() -> None:
    g, reg, ev = _graph("S")
    with pytest.raises(ExecutionError, match="exec.convert_unsupported"):
        Executor(g, reg, ev).evaluate("conv", "out")
