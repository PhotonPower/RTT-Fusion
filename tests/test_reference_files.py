"""R7: reference files for the C++ port (reference/*.result.json).

Every case of tools/reference_cases.py is recomputed and compared with its file within the stored
tolerance; independent values hold within their physical tolerance; the encoding round-trips.
Cases that need FourierGrating are skipped without it (RTT_FUSION_REQUIRE_MODULES=1: fail).
"""

from __future__ import annotations

import importlib.util
import math
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.results import FIELD_REFERENCE, decode, encode, load_json, to_json

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

from reference_cases import CASES, NEEDS_FOURIERGRATING  # noqa: E402

HAS_FG = importlib.util.find_spec("fouriergrating") is not None


def _skip_without_fg(name: str) -> None:
    if name in NEEDS_FOURIERGRATING and not HAS_FG:
        if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
            raise RuntimeError("fouriergrating required")
        pytest.skip("fouriergrating not installed")


def _close(a: Any, b: Any, rel: float, abs_: float, where: str) -> None:
    if isinstance(a, dict):
        assert isinstance(b, dict) and set(a) == set(b), where
        for key in a:
            _close(a[key], b[key], rel, abs_, f"{where}/{key}")
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), where
        for i, (x, y) in enumerate(zip(a, b, strict=True)):
            _close(x, y, rel, abs_, f"{where}/{i}")
    elif isinstance(a, np.ndarray):
        assert isinstance(b, np.ndarray) and a.shape == b.shape, where
        if a.dtype.kind in "fc":
            assert np.allclose(a, b, rtol=max(rel, 1e-10), atol=abs_, equal_nan=True), where
        else:
            assert np.array_equal(a, b), where
    elif isinstance(a, float | complex):
        assert abs(a - b) <= max(abs_, max(rel, 1e-10) * abs(b)) or (
            isinstance(a, float) and math.isnan(a) and math.isnan(b)
        ), where
    else:
        assert a == b, where


def _load(name: str) -> dict[str, Any]:
    type_name, data = load_json((ROOT / "reference" / f"{name}.result.json").read_text("utf-8"))
    assert type_name == FIELD_REFERENCE
    return data


@pytest.mark.parametrize("name", list(CASES))
def test_file_is_up_to_date(name: str) -> None:
    _skip_without_fg(name)
    stored = _load(name)
    fresh = decode(encode(CASES[name]()))
    tol = dict(fresh["tolerance"])
    inputs_rel, inputs_abs = 1e-12, 0.0
    if name in NEEDS_FOURIERGRATING:
        # FourierGrating's eigen solvers round differently between builds (BLAS, compiler): the
        # regenerated table may differ in the 1e-10 range from the stored one; the port itself is
        # checked against the stored inputs (test_derived_values_from_stored_inputs)
        tol = {"rel": max(tol["rel"], 1e-8), "abs": max(tol["abs"], 1e-9)}
        inputs_rel, inputs_abs = 1e-8, 1e-9
    _close(fresh["expected"], stored["expected"], tol["rel"], tol["abs"], f"{name}/expected")
    _close(fresh["inputs"], stored["inputs"], inputs_rel, inputs_abs, f"{name}/inputs")
    for key in ("case", "description", "source", "parameters", "tolerance"):
        _close(fresh[key], stored[key], 1e-12, 0.0, f"{name}/{key}")


def test_derived_values_from_stored_inputs() -> None:
    # interpolation from the stored table reproduces the stored values exactly (no FourierGrating)
    from rtt_fusion.couplings import jones_table as jt

    d = _load("jones_table_lamellar")
    p, i = d["parameters"], d["inputs"]
    table = jt.JonesTable(
        tuple(p["orders"]),
        p["side"],
        i["wavelengths_um"],
        i["sigma_x"],
        i["sigma_y"],
        p["period_um"],
        p["n_in"],
        p["n_out"],
        i["jones"],
    )
    for q, expected in zip(i["queries"], d["expected"]["interpolated"], strict=True):
        got = np.array([jt.interpolate(table, m, *q) for m in p["orders"]])
        assert np.allclose(got, expected, rtol=1e-12, atol=d["tolerance"]["abs"])


@pytest.mark.parametrize("name", list(CASES))
def test_independent_values(name: str) -> None:
    data = _load(name)
    ref, tol = data["reference"], data["reference_tolerance"]
    exp = data["expected"]
    if name == "ports_grid_to_spectrum":
        assert np.max(np.abs(ref["single_wave_amplitude"] - np.array([0.3, 0.7j]))) < tol["abs"]
    elif name == "hemisphere_quadrature":
        assert abs(np.sum(exp["forward_weights"]) - ref["sum_weights"]) < tol["abs"]
        cos_int = np.sum(exp["forward_weights"] * exp["forward_directions"][:, 2])
        assert abs(cos_int - ref["integral_cos_theta"]) < tol["abs"]
    elif name in ("focused_beam_spectrum", "far_field_spectrum_dipole"):
        assert abs(exp["power"] / ref["power"] - 1) < tol["rel"]
    elif name == "jones_table_lamellar":
        assert np.max(np.abs(exp["interpolated"] - ref["direct"])) < tol["abs"]
    else:
        assert ref == {}


@pytest.mark.parametrize("name", list(CASES))
def test_round_trip(name: str) -> None:
    data = _load(name)
    again = load_json(to_json(FIELD_REFERENCE, data))[1]
    _close(data, again, 0.0, 0.0, name)
