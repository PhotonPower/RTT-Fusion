"""R6: prototype of the tabulated diffraction interaction (Jones matrix per order).

Grating of pilot K1 (lamellar, Λ = 0.7 µm, n = 2.1 / air). Checks of the proposed rules
(docs/adr-entwurf-beugungstabelle.md): PRT matrix (P k_in = k_m, transversality, energy over all
orders = 1 for any polarisation), normal incidence without basis degeneracy, power normalisation
against the efficiencies of K1, phase reference at the grating origin, interpolation accuracy and
detection of Rayleigh anomalies.
"""

from __future__ import annotations

import copy
import os
from typing import Any

import numpy as np
import pytest

try:
    import fouriergrating  # noqa: F401
except ImportError:
    if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
        raise
    pytest.skip("fouriergrating not installed", allow_module_level=True)

from test_pilot_k1 import CONFIG

from rtt_fusion.adapters import fouriergrating as fga
from rtt_fusion.couplings import jones_table as jt
from rtt_fusion.couplings.grating_table import order_efficiencies

LAM = 0.6
PERIOD = 0.7


def _handle(config: Any = CONFIG, lams: tuple[float, ...] = (LAM,)) -> fga.Handle:
    return fga.prepare(fga.Structure(config), list(lams))


def _k(sigma: tuple[float, float], n: float, sign: float = 1.0) -> np.ndarray:  # type: ignore[type-arg]
    return np.array([sigma[0], sigma[1], sign * np.sqrt(n * n - sigma[0] ** 2 - sigma[1] ** 2)]) / n


@pytest.mark.parametrize("sigma", [(0.0, 0.0), (0.07, 0.05), (-0.1, 0.12)])
def test_prt_matrices_conserve_energy_for_any_polarisation(sigma: tuple[float, float]) -> None:
    h = _handle()
    k_in = _k(sigma, 1.0)
    rng = np.random.default_rng(1)
    tables = {side: jt.direct_jones(h, LAM, *sigma, side) for side in ("transmitted", "reflected")}
    for _ in range(3):
        e = rng.normal(size=3) + 1j * rng.normal(size=3)
        e = e - np.dot(k_in, e) * k_in
        e = e / np.linalg.norm(e)
        total = 0.0
        for side, table in tables.items():
            reflected = side == "reflected"
            for m, jones in table.items():
                s_out = (sigma[0] + m * LAM / PERIOD, sigma[1])
                p = jt.prt_matrix(jones, sigma, s_out, 1.0, 1.0, reflected)
                k_m = _k(s_out, 1.0, -1.0 if reflected else 1.0)
                assert np.allclose(p @ k_in, k_m, atol=1e-13)
                e_out = p @ e
                assert abs(np.dot(k_m, e_out)) < 1e-12
                total += float(np.vdot(e_out, e_out).real)
        assert abs(total - 1.0) < 1e-9


def test_normal_incidence_is_anisotropic_and_well_defined() -> None:
    jones = jt.direct_jones(_handle(), LAM, 0.0, 0.0)
    j0 = jones[0]
    # TE (E along the grooves, y) and TM (E along x) differ; no coupling at normal incidence
    assert abs(abs(j0[0, 0]) - abs(j0[1, 1])) > 1e-2
    assert abs(j0[0, 1]) < 1e-12 and abs(j0[1, 0]) < 1e-12


@pytest.mark.parametrize("sigma", [(0.0, 0.0), (0.0805, 0.0)])
def test_power_normalisation_equals_efficiencies(sigma: tuple[float, float]) -> None:
    h = _handle()
    k0 = 2 * np.pi / (LAM * 1e-3)
    jones = jt.direct_jones(h, LAM, *sigma)
    k_in = _k(sigma, 1.0)
    for pol in ("s", "p"):
        eff = order_efficiencies(h, LAM, sigma[0] * k0, sigma[1] * k0, pol)["transmitted"]
        if sigma == (0.0, 0.0):
            e = np.array([0.0, 1.0, 0.0]) if pol == "s" else np.array([1.0, 0.0, 0.0])
        else:  # incidence in the x-z plane: s = k × z ∝ −y, p = k × s
            s_vec = np.cross(k_in, [0.0, 0.0, 1.0])
            s_vec = s_vec / np.linalg.norm(s_vec)
            e = s_vec if pol == "s" else np.cross(k_in, s_vec)
        for m, jm in jones.items():
            s_out = (sigma[0] + m * LAM / PERIOD, sigma[1])
            e_out = jt.prt_matrix(jm, sigma, s_out, 1.0, 1.0) @ e
            assert abs(np.vdot(e_out, e_out).real - eff[m]) < 1e-12


def test_phase_reference_is_the_grating_origin() -> None:
    shifted = copy.deepcopy(CONFIG)
    x0 = 0.1
    shifted["layers"][0]["breakpoints_um"] = [0.0, x0, 0.3 + x0, PERIOD]  # type: ignore[index]
    shifted["layers"][0]["materials"] = ["air", "hi", "air"]  # type: ignore[index]
    a = jt.direct_jones(_handle(), LAM, 0.05, 0.02)
    b = jt.direct_jones(_handle(shifted), LAM, 0.05, 0.02)
    for m in a:
        assert np.max(np.abs(b[m] - a[m] * np.exp(-1j * m * 2 * np.pi * x0 / PERIOD))) < 1e-10


def test_interpolation_and_rayleigh_anomaly() -> None:
    lams = np.linspace(0.58, 0.62, 4)
    h = _handle(lams=(*lams, 0.6012))
    table = jt.build_table(
        h, lams, np.linspace(-0.1, 0.1, 5), np.array([-0.06, -0.02, 0.02, 0.06]), (-1, 0, 1)
    )
    q = (0.6012, 0.031, 0.007)
    direct = jt.direct_jones(h, *q)
    assert not jt.near_rayleigh_anomaly(table, *q)
    for m in (-1, 0, 1):
        assert np.max(np.abs(jt.interpolate(table, m, *q) - direct[m])) < 5e-3


def test_rayleigh_anomaly_detection() -> None:
    # order ±1 grazes at |σ_x ± λ/Λ| = 1, i.e. σ_x ≈ ±0.143 at λ = 0.6; only the axes matter here
    table = jt.JonesTable(
        (0,),
        "transmitted",
        np.linspace(0.58, 0.62, 5),
        np.linspace(0.0, 0.3, 13),
        np.linspace(-0.06, 0.06, 4),
        PERIOD,
        1.0,
        1.0,
        np.zeros((1, 5, 13, 4, 2, 2), dtype=np.complex128),
    )
    assert jt.near_rayleigh_anomaly(table, 0.6, 0.140, 0.0)
    assert jt.near_rayleigh_anomaly(table, 0.6, 0.150, 0.01)
    assert not jt.near_rayleigh_anomaly(table, 0.6, 0.05, 0.0)
    assert not jt.near_rayleigh_anomaly(table, 0.6, 0.24, 0.0)


def test_w62_table_from_apply_many_equals_point_by_point() -> None:
    """L4 O2: the table built with one apply_many per wavelength equals direct_jones point by
    point, bitwise (FourierGrating still solves each incidence on its own)."""
    pytest.importorskip("fouriergrating")
    handle = _handle()
    sx, sy = np.array([-0.2, 0.0, 0.15]), np.array([0.0, 0.1])
    table = jt.build_table(handle, np.array([LAM]), sx, sy, (-1, 0, 1))
    for b, x in enumerate(sx):
        for c, y in enumerate(sy):
            d = jt.direct_jones(handle, LAM, float(x), float(y))
            for o, m in enumerate((-1, 0, 1)):
                if m in d:
                    assert np.array_equal(table.jones[o, 0, b, c], d[m])
                else:
                    assert np.all(np.isnan(table.jones[o, 0, b, c]))
