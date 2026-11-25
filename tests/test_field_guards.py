"""Guards of the field propagation (#135, case 1): Fieldtracing 0.2.x's curvature fit that looped
forever on the camera of the FrObSy exposure head with the head shifted, and a grid too large
for the memory refused with exec.field instead of asked for (17.8 GiB)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import field
from rtt_fusion.graph.light import CoherenceKey

DATA = Path(__file__).parent / "data" / "curvature_loop_135.npz"


def test_the_curvature_fit_ends() -> None:
    """The samples (2650 points of ∇ψ) on which the fit of Fieldtracing 0.2.0 solved the same
    linear program forever: with the fix it ends after a few, and its bound holds everywhere."""
    pytest.importorskip("fieldtracing")
    import fieldtracing.selection as selection

    if not hasattr(selection, "_curvature_lp"):  # a Fieldtracing before the subset loop (0.1)
        pytest.skip("this Fieldtracing fits the curvature without the subset loop")

    from rtt_fusion.graph._fieldtracing_fix import install

    install()
    data = np.load(DATA)
    p, g = data["points"], data["gradient"]
    calls = []
    lp = selection._curvature_lp

    def counting(pp: Any, gg: Any) -> Any:
        calls.append(len(pp))
        if len(calls) > 50:
            raise AssertionError("the curvature fit does not end")
        return lp(pp, gg)

    selection._curvature_lp = counting
    try:
        c, bound = selection.fit_curvature(p, g)
    finally:
        selection._curvature_lp = lp
    assert len(calls) < 10, calls
    residual = np.max(np.abs(g - p @ c), axis=1)
    assert bound >= residual.max() - 1e-12
    assert bound == pytest.approx(140.6, rel=1e-3)  # the optimum the subsets reached


def test_a_grid_beyond_the_memory_is_refused() -> None:
    """A Debye plan of 19944 samples per axis (#135) – or any plan whose grid would take more than
    FIELD_GRID_BYTES – raises FieldError (exec.field at the detector) before it runs."""
    ran: list[bool] = []
    plan = SimpleNamespace(
        method="debye",
        parameters={"samples": 4986.0},
        run=lambda f: ran.append(True),
    )
    assert field.grid_bytes(plan) > field.FIELD_GRID_BYTES
    fb = SimpleNamespace(key=CoherenceKey("source", 0, 0), field=None)
    with pytest.raises(field.FieldError, match="4986 samples per axis"):
        field._run(fb, 164.5, plan, [])  # type: ignore[arg-type]
    assert not ran
    small = SimpleNamespace(method="debye", parameters={"samples": 768.0})
    assert field.grid_bytes(small) < field.FIELD_GRID_BYTES  # the largest Debye about its axis
