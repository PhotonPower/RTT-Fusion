"""Meta blocks (R20): sweep, replicates, optimisation around the acyclic graph.

References, independent of the code under test – the two-beam interference of W5 on the axis of
the analytic Michelson: I(δ) = 2 I_Arm (1 + cos 2k₀δ), I_Arm = 1/4, δ the extra length of arm B:
- W24 sweep over δ: the centre follows I(δ) at every value; the source is computed once (cache:
  only the blocks behind the changed mirror run again).
- W25 replicates over δ equally spaced over one period λ/2: the mean image is 2 I_Arm on the axis
  (the cosine averages out); optimisation of δ within [λ/4, 3λ/4] for the maximum: δ = λ/2.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from rtt_fusion.graph import Node, meta
from rtt_fusion.graph.execute import ExecutionError, Executor
from rtt_fusion.graph.meta import MetaResult

sys.path.insert(0, str(Path(__file__).parent))
from test_execute import I_ARM, WAVELENGTH_MM, michelson

X0 = 20.0  # position of mirror B along x (mm), see test_execute.michelson


def with_meta(type_id: str, **params: object) -> Executor:
    ex = michelson()
    meta.register(ex.registry)
    ex.evaluators.update(meta.EVALUATORS)
    ex.graph.add(
        Node(
            "m",
            type_id,
            params={
                "parameter": "mB.point[0]",
                "collect": "cam.image",
                "measure": "centre",
                **params,
            },
        )
    )
    return ex


def expected(delta: float) -> float:
    return float(2 * I_ARM * (1 + np.cos(2 * 2 * np.pi / WAVELENGTH_MM * delta)))


def test_w24_sweep_follows_two_beam_interference() -> None:
    deltas = np.linspace(0.0, WAVELENGTH_MM / 2, 9)
    ex = with_meta("meta.sweep", values=[X0 + d for d in deltas])
    out: MetaResult = ex.evaluate("m", "results")
    assert out.values == pytest.approx(list(X0 + deltas))
    assert out.measures == pytest.approx([expected(d) for d in deltas], abs=2e-3)
    assert ex.calls[("src", "out")] == 1  # upstream of the mirror: cached
    assert ex.calls[("mB", "out")] == len(deltas)
    assert ex.graph.nodes["mB"].params["point"] == [20.0, 0, 0]  # restored afterwards


def test_w25_replicates_average_and_optimum() -> None:
    n = 8
    deltas = [k * WAVELENGTH_MM / 2 / n for k in range(n)]  # one period of cos 2k₀δ
    ex = with_meta("meta.replicates", values=[X0 + d for d in deltas])
    out: MetaResult = ex.evaluate("m", "results")
    ny, nx = out.mean.values.shape
    centre = out.mean.values[ny // 2 - 1 : ny // 2 + 2, nx // 2 - 1 : nx // 2 + 2].mean()
    assert centre == pytest.approx(2 * I_ARM, abs=2e-3)
    ex = with_meta(
        "meta.optimise",
        bounds=[X0 + WAVELENGTH_MM / 4, X0 + 3 * WAVELENGTH_MM / 4],
        goal="maximum",
        tolerance=1e-7,
    )
    best: MetaResult = ex.evaluate("m", "results")
    assert best.best == pytest.approx(X0 + WAVELENGTH_MM / 2, abs=1e-3 * WAVELENGTH_MM)


def test_meta_errors() -> None:
    ex = with_meta("meta.sweep", values=[X0], parameter="nope.point[0]")
    with pytest.raises(ExecutionError, match="exec.param"):
        ex.evaluate("m", "results")
    ex = with_meta("meta.sweep", values=[X0], collect="m.results")
    with pytest.raises(ExecutionError, match="exec.cycle"):
        ex.evaluate("m", "results")
