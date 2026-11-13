"""Grating as an operator on fields (grating ADR, stage F1: GF1, GF2, interpolated transfer).

References: `fouriergrating.solve` for the plane wave at the beam's direction (efficiency per
order, the GUI's and the plugin's numbers), the grating's own exact path for the interpolated one.
Powers are measured spectrally (Σ Re k_z·|E|² per plane wave of each order's field).
"""

from __future__ import annotations

import copy
from typing import Any

import numpy as np
import pytest

from rtt_fusion.couplings import grating_field as gf
from rtt_fusion.graph import Graph, Node, builtin_registry, optics, spectral
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.field import rays_to_field
from rtt_fusion.graph.light import Plane, RayBundle

fg = pytest.importorskip("fouriergrating")
pytest.importorskip("fieldtracing")

LAM = 0.55
LAMELLAR = {
    "dimension": "1D",
    "solver": "rcwa",
    "truncation": {"Nx": 8},
    "period_um": {"x": 0.7},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}},
    "superstrate": "air",
    "substrate": "glass",
    "layers": [
        {
            "type": "lamellar1d",
            "thickness_um": 0.25,
            "breakpoints_um": [0, 0.35, 0.7],
            "materials": ["glass", "air"],
        }
    ],
}


def _field(theta_deg: float, radius: float, rings: int) -> Any:
    """A collimated beam at θ in the x-z plane (polarised in that plane), as a field in the
    grating's plane z = 0 and frame (x = grating vector)."""
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": radius,
                "rings": rings,
                "center": [0, 0, 1.0],
                "direction": [0, 0, -1],
                "wavelengths_um": [LAM],
                "field_angles_deg": [[theta_deg, 0.0]],
            },
        )
    )
    rays: RayBundle = Executor(g, reg, ev).evaluate("src", "out")[0]
    log: list[str] = []
    target = gf.onto_plane(rays_to_field(rays), Plane.make([0, 0, 0], [0, 0, 1]), False, log)
    return gf.in_frame(target.field, gf.frame(target, [1.0, 0.0, 0.0]).q)


def _power(field: Any, n: float) -> float:
    s = spectral.sample(field, 4096)
    sp = spectral.spectrum(s.residual)
    kxx, kyy = s.grid.frequencies()
    k = 2 * np.pi / (LAM * 1e-3) * n
    kx, ky = kxx + s.grid.carrier[0], kyy + s.grid.carrier[1]
    kz = np.sqrt((k**2 - kx**2 - ky**2).astype(complex)).real
    return float(np.sum(kz[..., None] * np.abs(sp) ** 2))


def _plane_wave_efficiencies(field: Any) -> dict[tuple[str, tuple[int, int]], float]:
    """FourierGrating's efficiencies for the plane wave at the field's mean direction (TM)."""
    sin_x = float(np.asarray(field.gradient).mean(axis=0)[0] / (2 * np.pi / (LAM * 1e-3)))
    cfg = copy.deepcopy(LAMELLAR)
    cfg["incidence"] = {
        "lambda_um": LAM,
        "phi_deg": float(np.degrees(np.arcsin(abs(sin_x)))),
        "psi_deg": 0.0 if sin_x >= 0 else 180.0,
        "delta_deg": 90.0,
    }
    res = fg.solve(cfg)
    out = {("front", (int(o["n"]), 0)): float(o["eff"]) for o in res.reflected}
    out |= {("back", (int(o["n"]), 0)): float(o["eff"]) for o in res.transmitted}
    return out


def _efficiencies(field: Any, applied: gf.Applied) -> dict[tuple[str, tuple[int, int]], float]:
    p0 = _power(field, 1.0)
    return {
        key: _power(f, 1.0 if key[0] == "front" else 1.5) / p0 for key, f in applied.fields.items()
    }


def _handle() -> Any:
    from rtt_fusion.adapters import fouriergrating as fga

    return fga.prepare(fga.Structure(LAMELLAR), (LAM,))


def test_gf1_orders_of_a_beam_are_the_plane_wave_efficiencies() -> None:
    """GF1: a collimated beam (radius 0.15 mm, λ/D ≈ 2e−3) on a lamellar grating, exact path: the
    power of every order is FourierGrating's efficiency of the plane wave; energy kept; every
    order a field with the carrier k_c + nK."""
    field = _field(10.0, 0.15, 5)
    applied = gf.apply(_handle(), field)
    eff = _efficiencies(field, applied)
    ref = _plane_wave_efficiencies(field)
    assert set(eff) == set(ref)  # the same propagating orders
    worst = max(abs(eff[k] - ref[k]) for k in ref)
    assert worst < 2e-5, (worst, eff, ref)
    assert sum(eff.values()) == pytest.approx(1.0, abs=3e-5)
    k_step = 2e3 * np.pi / 0.7
    carrier = np.asarray(field.gradient).mean(axis=0)
    for (side, (n, _)), f in applied.fields.items():
        g = np.asarray(f.gradient).mean(axis=0)
        expected = carrier[0] + n * k_step
        assert g[0] == pytest.approx(expected if side == "back" else expected, rel=1e-9)


def test_gf2_a_wider_beam_comes_closer_to_the_plane_wave() -> None:
    """GF2: the deviation from the plane wave's efficiencies shrinks as the beam gets wider
    (0.15 → 0.3 mm radius). Interpolated transfer: the wide beam's 128² spectrum is verified on a
    33² sub-grid; the narrow one's 32² spectrum is too coarse for 1e−8 and takes the exact path."""
    handle = _handle()
    errors = []
    for radius, rings in ((0.15, 5), (0.3, 12)):
        field = _field(10.0, radius, rings)
        log: list[str] = []
        applied = gf.apply(handle, field, transfer="interpolated", log=log)
        assert any(("verified" if radius > 0.2 else "exact") in line for line in log), log
        eff = _efficiencies(field, applied)
        ref = _plane_wave_efficiencies(field)
        errors.append(max(abs(eff[k] - ref[k]) for k in ref))
    assert errors[1] < errors[0] < 2e-5, errors


def test_interpolated_transfer_equals_the_exact_one() -> None:
    """The interpolated transfer, verified to a tolerance of |T|, gives the exact path's fields
    within it (tolerance 1e−6: a 17² sub-grid of the 32² spectrum)."""
    handle = _handle()
    field = _field(10.0, 0.15, 5)
    exact = gf.apply(handle, field)
    log: list[str] = []
    interp = gf.apply(handle, field, transfer="interpolated", tolerance=1e-6, log=log)
    assert any("verified" in line for line in log), log
    assert set(exact.fields) == set(interp.fields)
    scale = max(float(np.max(np.abs(f.amplitude))) for f in exact.fields.values())
    for key, f in exact.fields.items():
        diff = np.abs(np.asarray(interp.fields[key].amplitude) - np.asarray(f.amplitude))
        assert float(diff.max()) < 2e-6 * scale, key


def test_interpolation_falls_back_near_a_rayleigh_anomaly() -> None:
    """At 12.36° the reflected order +1 grazes (sin θ + λ/Λ = 1): it appears inside the beam's
    spectrum, the interpolation is refused and the exact path taken."""
    handle = _handle()
    field = _field(float(np.degrees(np.arcsin(1 - LAM / 0.7))), 0.15, 5)
    log: list[str] = []
    applied = gf.apply(handle, field, transfer="interpolated", log=log)
    assert any("exact" in line for line in log), log
    assert ("back", (0, 0)) in applied.fields


def test_orders_can_be_chosen() -> None:
    field = _field(10.0, 0.15, 5)
    applied = gf.apply(_handle(), field, orders=lambda n, m: n == 0, transfer="interpolated")
    assert set(applied.fields) == {("front", (0, 0)), ("back", (0, 0))}
