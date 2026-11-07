"""Pilot K4 (R10): biosensor with an index profile from molecular dynamics (BioOptoMD, data port N).

Kretschmann surface-plasmon sensor at 633 nm: glass prism (1.515) / gold 50 nm / adlayer n(z) /
water (1.333), light from the prism (BioOptoMD's substrate side). The adlayer profiles are
written in BioOptoMD's export format (CSV of `analysis.export_refractive_index`; BioOptoMD needs
CUDA, so the files are synthetic here). Checks:

1. the layered stack from `profile_layers` gives the same reflectance through the FourierGrating
   adapter and an independent transfer matrix (< 1e-9);
2. de Feijter: two profiles with equal excess optical thickness ∫(n − n_w) dz but different
   shapes shift the resonance equally to first order (here within 3 %), and the shift is linear in
   the amplitude (3 %);
3. coarse-graining the 0.05-nm bins to 1-nm layers changes the shift by < 1 %.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from rtt_fusion.adapters import biooptomd as bmd

LAM = 0.633  # µm
K0 = 2 * np.pi / (LAM * 1e-3)
N_PRISM, N_WATER = 1.515, 1.333
N_GOLD = 0.197 + 3.44j  # exp(-iωt), Im > 0
D_GOLD_NM = 50.0
DZ = 0.05


def _write_profile(path: Path, shape: str, amplitude: float) -> bmd.Profile:
    z = (np.arange(400) + 0.5) * DZ
    if shape == "exp":
        dn = amplitude * np.exp(-z / 3.0)
    else:
        dn = np.where(z < 6.0, amplitude * 0.5, 0.0)
    n = np.vstack([N_WATER + dn, N_WATER + 0.98 * dn])  # two replicas
    header = "BioOptoMD RefractiveIndexProfiler (synthetic test data)\nz_nm,n_replica_0,n_replica_1"
    np.savetxt(path, np.column_stack([z, n.T]), delimiter=",", header=header, comments="# ")
    return bmd.read_profile_csv(path, DZ, N_WATER, LAM)


def _rp(layers: list[tuple[float, complex]], theta: float) -> float:
    """p reflectance of prism | layers | water, transfer matrix (Byrnes, tangential E, H)."""
    stack = [(N_PRISM, 0.0), *((n, d) for d, n in layers), (N_WATER, 0.0)]
    kx = N_PRISM * K0 * np.sin(theta)
    kz = [np.sqrt(complex((n * K0) ** 2 - kx * kx)) for n, _ in stack]
    kz = [k if k.imag >= 0 else -k for k in kz]
    eta = [(n * K0) ** 2 / k for (n, _), k in zip(stack, kz, strict=True)]  # p: ε k0²/k_z
    m = np.eye(2, dtype=complex)
    for j in range(1, len(stack) - 1):
        d = stack[j][1]
        c, s = np.cos(kz[j] * d), np.sin(kz[j] * d)
        m = m @ np.array([[c, -1j * s / eta[j]], [-1j * eta[j] * s, c]])
    b, c = m @ np.array([1.0, eta[-1]])
    r = (eta[0] * b - c) / (eta[0] * b + c)
    return float(abs(r) ** 2)


def _layers(profile: bmd.Profile | None, coarse: int = 1) -> list[tuple[float, complex]]:
    out: list[tuple[float, complex]] = [(D_GOLD_NM * 1e-6, N_GOLD)]
    if profile is None:
        return out
    t, n = bmd.profile_layers(profile)
    if coarse > 1:  # block averages of `coarse` bins
        k = len(n) // coarse * coarse
        n = np.concatenate([n[:k].reshape(-1, coarse).mean(axis=1), n[k:]])
        t = np.concatenate([np.full(k // coarse, coarse * DZ), t[k:]])
    out += [(d, complex(v)) for d, v in bmd.stack_along_z(t, n, "substrate")]
    return out


def _resonance(layers: list[tuple[float, complex]]) -> float:
    res = minimize_scalar(
        lambda th: _rp(layers, th),
        bounds=(np.radians(66), np.radians(76)),
        method="bounded",
        options={"xatol": 1e-10},
    )
    return float(res.x)


def test_profile_reading_and_layers(tmp_path: Path) -> None:
    p = _write_profile(tmp_path / "n.csv", "exp", 0.12)
    assert p.n.shape == (2, 400) and abs(p.dz_nm - DZ) < 1e-15
    t, n = bmd.profile_layers(p, replica=0)
    assert abs(t.sum() - 400 * DZ) < 1e-9  # exponential never reaches n_w exactly
    assert abs(bmd.excess_optical_thickness_nm(p, 0) - 0.12 * 3 * (1 - np.exp(-20 / 3))) < 1e-4
    box = _write_profile(tmp_path / "b.csv", "box", 0.12)
    tb, nb = bmd.profile_layers(box, replica=0)
    assert len(nb) == 1 and abs(tb[0] - 6.0) < 1e-9  # merged, trailing water removed
    assert bmd.stack_along_z(tb, nb, "solution") == bmd.stack_along_z(tb, nb, "substrate")[::-1]


def test_fouriergrating_adapter_equals_transfer_matrix(tmp_path: Path) -> None:
    pytest.importorskip("fouriergrating")
    from rtt_fusion.adapters import fouriergrating as fga
    from rtt_fusion.ports import PlaneWaves

    p = _write_profile(tmp_path / "n.csv", "box", 0.24)
    layers = _layers(p)
    materials = {
        "prism": {"n": N_PRISM},
        "water": {"n": N_WATER},
        "gold": {"n": [N_GOLD.real, N_GOLD.imag]},
    }
    fg_layers = []
    for i, (d, n) in enumerate(reversed(layers)):  # FourierGrating lists from the substrate up
        name = "gold" if n == N_GOLD else f"l{i}"
        if name != "gold":
            materials[name] = {"n": float(n.real)}
        fg_layers.append({"type": "homogeneous", "thickness_um": d * 1e3, "material": name})
    cfg = {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 0, "Nz": 0},
        "period_um": {"x": 1.0, "z": 0.0},
        "materials": materials,
        "superstrate": "prism",
        "substrate": "water",
        "layers": fg_layers,
    }
    h = fga.prepare(fga.Structure(cfg), [LAM])
    for deg in (60.0, 71.0):
        kx = N_PRISM * K0 * np.sin(np.radians(deg))
        out = fga.apply(
            h,
            PlaneWaves(
                fga.front_port(h, LAM),
                np.array([kx]),
                np.array([0.0]),
                np.array([[1.0, 0.0]], dtype=np.complex128),
            ),
        )  # p in x-z
        e_r = out["front"].field[0]
        kz = np.sqrt((N_PRISM * K0) ** 2 - kx * kx)
        r_fg = (abs(e_r[0]) ** 2) * (1 + (kx / kz) ** 2) / (1 + (kx / kz) ** 2)
        assert abs(r_fg - _rp(layers, np.radians(deg))) < 1e-9


def test_de_feijter_and_linearity(tmp_path: Path) -> None:
    bare = _resonance(_layers(None))
    shift: dict[str, float] = {}
    excess: dict[str, float] = {}
    for name, shape, amp in (("exp", "exp", 0.12), ("box", "box", 0.12), ("exp2", "exp", 0.24)):
        p = _write_profile(tmp_path / f"{name}.csv", shape, amp)
        shift[name] = _resonance(_layers(p)) - bare
        excess[name] = bmd.excess_optical_thickness_nm(p)
        shift[name + "_coarse"] = _resonance(_layers(p, coarse=20)) - bare
    assert shift["exp"] > 0  # an adlayer moves the resonance to larger angles
    # shift per excess optical thickness independent of the profile shape (measured 0.07 %)
    assert abs((shift["box"] / excess["box"]) / (shift["exp"] / excess["exp"]) - 1) < 0.005
    assert abs(shift["exp2"] / shift["exp"] - 2) < 0.06
    for name in ("exp", "box", "exp2"):
        assert abs(shift[name + "_coarse"] / shift[name] - 1) < 0.01


def test_dipoles_from_tracker() -> None:
    # two frames: R_F = (1, 2, 3) nm, μ̂ along (0, 0.6, 0.8), d = 4 nm; surface at z = 0.01 mm
    frame = [1.0, 2.0, 3.0, 0.0, 0.6, 0.8, 4.0]
    series = np.array([[frame, frame]])
    up = bmd.dipoles_from_tracker(series, 0.01, "substrate")
    assert np.allclose(up.positions_mm[0], [1e-6, 2e-6, 0.01 + 4e-6])
    assert np.allclose(up.orientations[0], [0.0, 0.6, 0.8])
    down = bmd.dipoles_from_tracker(series, 0.01, "solution")
    assert np.allclose(down.positions_mm[1], [1e-6, 2e-6, 0.01 - 4e-6])
    assert np.allclose(down.orientations[1], [0.0, 0.6, -0.8])
