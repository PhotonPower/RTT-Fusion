"""Reference cases for the C++ port of rtt-fusion (stage R7; Fusion ADR draft point 9).

Each case stores inputs, the values of the Python prototype (`expected`, same algorithm, tight
tolerance) and, where available, an independent value (`reference`, physical tolerance).
Written by tools/export_reference.py to reference/<case>.result.json, checked by
tests/test_reference_files.py. Solver-free cases (ports, conversions, interpolation) are what
`rtt-fusion` implements; adapter cases document the conventions of the external adapters.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from typing import Any

import numpy as np

from rtt_fusion.convert import far_field_to_spectrum, hemisphere
from rtt_fusion.couplings.focus import richards_wolf
from rtt_fusion.ports import GriddedField, Port, power, spectrum_power, to_plane_waves
from rtt_fusion.results import field_reference


def _tight(*arrays: Any) -> dict[str, float]:
    scale = max(float(np.max(np.abs(np.asarray(a)))) for a in arrays) or 1.0
    return {"rel": 1e-12, "abs": 1e-13 * scale}


def ports_grid_to_spectrum() -> dict[str, Any]:
    port = Port(0.0, 1, 1.5 + 0j, 0.5)
    nx, ny, dx, dy = 12, 8, 0.5e-3 / 8, 0.5e-3 / 6
    rng = np.random.default_rng(11)
    field = rng.normal(size=(ny, nx, 2)) + 1j * rng.normal(size=(ny, nx, 2))
    g = GriddedField(port, dx, dy, field)
    s = to_plane_waves(g)
    # independent: a single grid plane wave comes back with its amplitude
    x = (np.arange(nx) - nx // 2) * dx
    y = (np.arange(ny) - ny // 2) * dy
    kx, ky = 2 * np.pi * 3 / (nx * dx), -2 * np.pi * 2 / (ny * dy)
    phase = np.exp(1j * (kx * x[None, :] + ky * y[:, None]))
    single = to_plane_waves(
        GriddedField(port, dx, dy, np.stack([0.3 * phase, 0.7j * phase], -1)), threshold=1e-9
    )
    return field_reference(
        "ports_grid_to_spectrum",
        "Gridded port field (G) to plane waves (S) on the centred grid, Ez and power.",
        "rtt_fusion.ports (to_plane_waves, power, PlaneWaves.ez); docs/konventionen.md",
        {
            "z_mm": 0.0,
            "direction": 1,
            "index": 1.5,
            "wavelength_um": 0.5,
            "nx": nx,
            "ny": ny,
            "dx_mm": dx,
            "dy_mm": dy,
        },
        {"field": field},
        {
            "kx_per_mm": s.kx,
            "ky_per_mm": s.ky,
            "amplitude": s.field,
            "ez": s.ez(),
            "power": power(g),
        },
        _tight(s.field, s.ez()),
        {
            "single_wave_amplitude": single.field[0],
            "single_wave_kx_per_mm": kx,
            "single_wave_ky_per_mm": ky,
        },
        {"abs": 1e-13},
    )


def hemisphere_quadrature() -> dict[str, Any]:
    dirs, w = hemisphere(6, 8, True)
    back, wb = hemisphere(6, 8, False)
    return field_reference(
        "hemisphere_quadrature",
        "Directions and solid-angle weights: Gauss-Legendre in cos(theta), uniform phi.",
        "rtt_fusion.convert.hemisphere",
        {"n_theta": 6, "n_phi": 8},
        {},
        {
            "forward_directions": dirs,
            "forward_weights": w,
            "backward_directions": back,
            "backward_weights": wb,
        },
        _tight(dirs, w),
        {"sum_weights": 2 * np.pi, "integral_cos_theta": np.pi},
        {"abs": 1e-13},
    )


def focused_beam_spectrum() -> dict[str, Any]:
    port = Port(0.0, 1, 1.0 + 0j, 0.6)
    na = 0.8
    beam = richards_wolf(port, na, (0.0, 0.0, 0.002), (1.0, 0.0), 8, 16)
    k = port.k0
    return field_reference(
        "focused_beam_spectrum",
        "Richards-Wolf beam (NA 0.8, x-polarised, uniform pupil) as weighted plane waves; power.",
        "rtt_fusion.couplings.focus.richards_wolf, rtt_fusion.ports.spectrum_power",
        {
            "z_mm": 0.0,
            "index": 1.0,
            "wavelength_um": 0.6,
            "na": na,
            "focus_mm": [0.0, 0.0, 0.002],
            "n_theta": 8,
            "n_phi": 16,
        },
        {},
        {
            "kx_per_mm": beam.kx,
            "ky_per_mm": beam.ky,
            "amplitude": beam.field,
            "solid_angle": beam.solid_angle,
            "ez": beam.ez(),
            "power": spectrum_power(beam),
            "power_na_0p5": spectrum_power(beam, 0.5),
        },
        _tight(beam.field, beam.ez()),
        # |a|² = cos θ: (2π)²/k ∫ cos θ dΩ over the cap = (2π)²/k · π NA² (exact for Gauss-Legendre)
        {"power": (2 * np.pi) ** 2 / k * np.pi * na**2},
        {"rel": 1e-12},
    )


def far_field_spectrum_dipole() -> dict[str, Any]:
    # far field of an x dipole at c, transverse part of p: E∞ = p − (k̂·p) k̂ (arbitrary scale)
    port = Port(0.004, 1, 1.0 + 0j, 0.6)
    center = np.array([0.0005, -0.0002, 0.002])
    dirs, w = hemisphere(10, 20, True)
    p = np.array([1.0, 0.0, 0.0])
    e_inf = (p[None, :] - (dirs @ p)[:, None] * dirs).astype(np.complex128)
    e_inf += 1e-3 * dirs  # a small radial part, as from a numerical far field (removed)
    s = far_field_to_spectrum(port, center, dirs, w, e_inf)
    k = port.k0
    return field_reference(
        "far_field_spectrum_dipole",
        "Far field of an x dipole beside the axis to plane waves at a port behind it (Weyl).",
        "rtt_fusion.convert.far_field_to_spectrum; adapter CliffordBem",
        {
            "z_mm": 0.004,
            "direction": 1,
            "index": 1.0,
            "wavelength_um": 0.6,
            "center_mm": center.tolist(),
            "n_theta": 10,
            "n_phi": 20,
        },
        {"directions": dirs, "weights": w, "e_inf": e_inf},
        {"kx_per_mm": s.kx, "ky_per_mm": s.ky, "amplitude": s.field, "power": spectrum_power(s)},
        _tight(s.field),
        # k ∫_hemisphere |p⊥|² dΩ = k · 4π/3 for |p| = 1
        {"power": k * 4 * np.pi / 3},
        {"rel": 1e-12},
    )


def _fg_config() -> dict[str, Any]:
    return {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 15, "Nz": 0},
        "period_um": {"x": 0.7, "z": 0.0},
        "materials": {"air": {"n": 1.0}, "hi": {"n": 2.1}},
        "superstrate": "air",
        "substrate": "air",
        "layers": [
            {
                "type": "lamellar1d",
                "thickness_um": 0.25,
                "breakpoints_um": [0.0, 0.3, 0.7],
                "materials": ["hi", "air"],
            }
        ],
    }


def jones_table_lamellar() -> dict[str, Any]:
    from rtt_fusion.adapters import fouriergrating as fga
    from rtt_fusion.couplings import jones_table as jt

    lams = np.linspace(0.58, 0.62, 4)
    queries = np.array([[0.6012, 0.031, 0.007], [0.5893, -0.064, -0.012], [0.6, 0.0, 0.0]])
    handle = fga.prepare(fga.Structure(_fg_config()), [*lams.tolist(), 0.6012, 0.5893, 0.6])
    table = jt.build_table(
        handle, lams, np.linspace(-0.1, 0.1, 5), np.linspace(-0.06, 0.06, 4), (-1, 0, 1)
    )
    interpolated = np.array([[jt.interpolate(table, m, *q) for m in (-1, 0, 1)] for q in queries])
    direct = np.array([[jt.direct_jones(handle, *q)[m] for m in (-1, 0, 1)] for q in queries])
    sigma = (0.031, 0.007)
    prt = np.array(
        [
            jt.prt_matrix(
                jt.direct_jones(handle, 0.6012, *sigma)[m],
                sigma,
                (sigma[0] + m * 0.6012 / 0.7, sigma[1]),
                1.0,
                1.0,
            )
            for m in (-1, 0, 1)
        ]
    )
    return field_reference(
        "jones_table_lamellar",
        "Tabulated diffraction interaction: table, cubic interpolation, PRT matrices (grating K1).",
        "rtt_fusion.couplings.jones_table; docs/adr-entwurf-beugungstabelle.md",
        {
            "period_um": 0.7,
            "n_in": 1.0,
            "n_out": 1.0,
            "orders": [-1, 0, 1],
            "side": "transmitted",
            "fouriergrating_config": json.dumps(_fg_config(), sort_keys=True),
        },
        {
            "wavelengths_um": table.wavelengths_um,
            "sigma_x": table.sigma_x,
            "sigma_y": table.sigma_y,
            "jones": table.jones,
            "queries": queries,
        },
        {
            "interpolated": interpolated,
            "prt_at_query_0": prt,
            "near_anomaly": np.array([jt.near_rayleigh_anomaly(table, *q) for q in queries]),
        },
        _tight(interpolated, prt),
        {"direct": direct},
        {"abs": 5e-3},  # cubic on this coarse table (J5)
    )


def fouriergrating_stack() -> dict[str, Any]:
    from rtt_fusion.adapters import fouriergrating as fga
    from rtt_fusion.ports import PlaneWaves

    cfg = {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 3, "Nz": 0},
        "period_um": {"x": 0.7, "z": 0.0},
        "materials": {
            "air": {"n": 1.0},
            "glass": {"n": 1.5},
            "hi": {"n": 2.1},
            "lossy": {"n": [1.8, 0.05]},
            "mid": {"n": 1.45},
        },
        "superstrate": "air",
        "substrate": "glass",
        "layers": [
            {"type": "homogeneous", "thickness_um": 0.13, "material": "mid"},
            {"type": "homogeneous", "thickness_um": 0.21, "material": "lossy"},
            {"type": "homogeneous", "thickness_um": 0.09, "material": "hi"},
        ],
    }
    handle = fga.prepare(fga.Structure(cfg), [0.6])
    k = 2 * np.pi / 0.6e-3 * np.sin(0.4)
    inc = PlaneWaves(
        fga.front_port(handle, 0.6),
        np.array([k * np.cos(0.3)]),
        np.array([k * np.sin(0.3)]),
        np.array([[0.4 - 0.1j, 0.2 + 0.7j]]),
    )
    out = fga.apply(handle, inc)
    front, back = out["front"], out["back"]
    assert front.orders is not None and back.orders is not None
    return field_reference(
        "fouriergrating_stack",
        "FourierGrating adapter: planar stack, oblique with azimuth, reflected and transmitted.",
        "rtt_fusion.adapters.fouriergrating (test F2); FourierGrating 0.16.0",
        {
            "fouriergrating_config": json.dumps(cfg, sort_keys=True),
            "wavelength_um": 0.6,
            "thickness_mm": handle.thickness_mm,
        },
        {"incident_kx_per_mm": inc.kx, "incident_ky_per_mm": inc.ky, "incident": inc.field},
        {
            "front_orders": front.orders,
            "front_kx_per_mm": front.kx,
            "front_amplitude": front.field,
            "back_orders": back.orders,
            "back_kx_per_mm": back.kx,
            "back_amplitude": back.field,
        },
        {"rel": 1e-9, "abs": 1e-11},
    )


CASES: dict[str, Callable[[], dict[str, Any]]] = {
    f.__name__: f
    for f in (
        ports_grid_to_spectrum,
        hemisphere_quadrature,
        focused_beam_spectrum,
        far_field_spectrum_dipole,
        jones_table_lamellar,
        fouriergrating_stack,
    )
}
NEEDS_FOURIERGRATING = {"jones_table_lamellar", "fouriergrating_stack"}
