"""R1: conformance of the FourierGrating adapter (docs/schnittstelle.md, Fusion plan 6.5).

Independent references: vector Fresnel operator in Convention A (built here from s/p bases and
Byrnes' coefficients with angles), transfer matrix of tangential fields for planar stacks, the
grating equation, energy conservation, reciprocity, and the BWPM adapter as a second solver.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

try:
    import fouriergrating  # noqa: F401
except ImportError:
    if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
        raise
    pytest.skip("fouriergrating not installed", allow_module_level=True)

from rtt_fusion.adapters import fouriergrating as fga
from rtt_fusion.ports import PlaneWaves

LAM = 0.6  # µm
K0 = 2 * np.pi / (LAM * 1e-3)


def _config(
    layers: list[dict[str, object]],
    sub: str = "glass",
    sup: str = "air",
    period: float = 0.7,
    nx: int = 3,
) -> dict[str, object]:
    return {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": nx, "Nz": 0},
        "period_um": {"x": period, "z": 0.0},
        "materials": {
            "air": {"n": 1.0},
            "glass": {"n": 1.5},
            "hi": {"n": 2.1},
            "lossy": {"n": [1.8, 0.05]},
            "mid": {"n": 1.45},
        },
        "superstrate": sup,
        "substrate": sub,
        "layers": layers,
    }


def _incident(
    handle: fga.Handle, theta: float, azimuth: float, e: tuple[complex, complex]
) -> PlaneWaves:
    port = fga.front_port(handle, LAM)
    k = handle.front_index.real * K0 * np.sin(theta)
    return PlaneWaves(
        port,
        np.array([k * np.cos(azimuth)]),
        np.array([k * np.sin(azimuth)]),
        np.array([e], dtype=np.complex128),
    )


def _zero_order(waves: PlaneWaves) -> np.ndarray:  # type: ignore[type-arg]
    assert waves.orders is not None
    j = [i for i, o in enumerate(waves.orders) if tuple(o) == (0, 0)]
    assert len(j) == 1
    return np.asarray(waves.field[j[0]])


def _fresnel_vector(e_perp: np.ndarray, kx: float, ky: float, n1: float, n2: float):  # type: ignore[no-untyped-def,type-arg]
    """Transverse (E_x, E_y) of transmitted and reflected wave, Convention A (as BWPM K1)."""
    kp2 = kx * kx + ky * ky
    k1z, k2z = np.sqrt(complex(n1**2 * K0**2 - kp2)), np.sqrt(complex(n2**2 * K0**2 - kp2))
    k_in, k_t, k_r = np.array([kx, ky, k1z]), np.array([kx, ky, k2z]), np.array([kx, ky, -k1z])
    e_in = np.array([e_perp[0], e_perp[1], -(kx * e_perp[0] + ky * e_perp[1]) / k1z])
    if kp2 == 0:
        s = np.array([0.0, 1.0, 0.0])  # any s works at normal incidence
    else:
        s = np.cross(k_in, [0.0, 0.0, 1.0]).real
        s = s / np.linalg.norm(s)
    p_in, p_t, p_r = (np.cross(k / (n * K0), s) for k, n in ((k_in, n1), (k_t, n2), (k_r, n1)))
    c1, c2 = k1z / (n1 * K0), k2z / (n2 * K0)
    r_s = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
    r_p = (n2 * c1 - n1 * c2) / (n2 * c1 + n1 * c2)
    t_s = 2 * n1 * c1 / (n1 * c1 + n2 * c2)
    t_p = 2 * n1 * c1 / (n2 * c1 + n1 * c2)
    a_s, a_p = e_in @ s, e_in @ p_in
    return (t_s * a_s * s + t_p * a_p * p_t)[:2], (r_s * a_s * s + r_p * a_p * p_r)[:2]


@pytest.mark.parametrize(
    ("theta", "azimuth"), [(0.0, 0.0), (0.0, 1.0), (0.5, 0.0), (0.5, 0.8), (1.2, 2.5)]
)
def test_single_interface_equals_vector_fresnel(theta: float, azimuth: float) -> None:
    handle = fga.prepare(fga.Structure(_config([])), [LAM])
    e = (0.6 + 0.2j, -0.3 + 0.5j)
    inc = _incident(handle, theta, azimuth, e)
    out = fga.apply(handle, inc)
    t_ref, r_ref = _fresnel_vector(np.array(e), float(inc.kx[0]), float(inc.ky[0]), 1.0, 1.5)
    assert np.max(np.abs(_zero_order(out["front"]) - r_ref)) < 1e-10
    assert np.max(np.abs(_zero_order(out["back"]) - t_ref)) < 1e-10


def _transfer(
    stack: list[tuple[complex, float]], n_in: complex, n_out: complex, kp2: float, pol: str
) -> tuple[complex, complex]:
    """(r, t) of tangential fields; phases at the first and last interface (FG references)."""
    layers = [(n_in, 0.0), *stack, (n_out, 0.0)]
    kz = [np.sqrt(complex(n * n * K0 * K0 - kp2)) for n, _ in layers]
    eta = [k if pol == "s" else n * n * K0 * K0 / k for (n, _), k in zip(layers, kz, strict=True)]
    a_plus, a_minus = 1.0 + 0j, 0.0 + 0j  # at the last interface, in the exit medium
    for j in range(len(layers) - 1, 0, -1):
        e_t, h_t = a_plus + a_minus, eta[j] * (a_plus - a_minus)
        a_plus, a_minus = 0.5 * (e_t + h_t / eta[j - 1]), 0.5 * (e_t - h_t / eta[j - 1])
        d = layers[j - 1][1] * 1e-3
        if j - 1 > 0:  # through layer j - 1 to its upper interface
            a_plus, a_minus = (
                a_plus * np.exp(-1j * kz[j - 1] * d),
                a_minus * np.exp(1j * kz[j - 1] * d),
            )
    return complex(a_minus / a_plus), complex(1.0 / a_plus)


STACK_CFG = [  # FG lists layers from the substrate upwards
    {"type": "homogeneous", "thickness_um": 0.13, "material": "mid"},
    {"type": "homogeneous", "thickness_um": 0.21, "material": "lossy"},
    {"type": "homogeneous", "thickness_um": 0.09, "material": "hi"},
]
STACK_TOP_DOWN = [(2.1, 0.09), (1.8 + 0.05j, 0.21), (1.45, 0.13)]


@pytest.mark.parametrize(("theta", "azimuth"), [(0.0, 0.0), (0.4, 0.3), (0.9, 2.0)])
def test_planar_stack_against_transfer_matrix(theta: float, azimuth: float) -> None:
    handle = fga.prepare(fga.Structure(_config(STACK_CFG)), [LAM])
    assert abs(handle.thickness_mm - 0.43e-3) < 1e-15
    e = (0.4 - 0.1j, 0.2 + 0.7j)
    inc = _incident(handle, theta, azimuth, e)
    out = fga.apply(handle, inc)
    kx, ky = float(inc.kx[0]), float(inc.ky[0])
    kp2 = kx * kx + ky * ky
    s_dir = np.array([-ky, kx]) / np.sqrt(kp2) if kp2 else np.array([0.0, 1.0])
    p_dir = np.array([kx, ky]) / np.sqrt(kp2) if kp2 else np.array([1.0, 0.0])
    r_s, t_s = _transfer(STACK_TOP_DOWN, 1.0, 1.5, kp2, "s")
    r_p, t_p = _transfer(STACK_TOP_DOWN, 1.0, 1.5, kp2, "p")
    ev = np.array(e)
    expected_r = r_s * (ev @ s_dir) * s_dir + r_p * (ev @ p_dir) * p_dir
    expected_t = t_s * (ev @ s_dir) * s_dir + t_p * (ev @ p_dir) * p_dir
    assert np.max(np.abs(_zero_order(out["front"]) - expected_r)) < 1e-9
    assert np.max(np.abs(_zero_order(out["back"]) - expected_t)) < 1e-9


GRATING = [
    {
        "type": "lamellar1d",
        "thickness_um": 0.25,
        "breakpoints_um": [0.0, 0.3, 0.7],
        "materials": ["hi", "air"],
    }
]


@pytest.mark.parametrize(("theta", "azimuth"), [(0.2, 0.0), (0.35, 0.6)])
def test_grating_equation_and_energy(theta: float, azimuth: float) -> None:
    handle = fga.prepare(fga.Structure(_config(GRATING, nx=12)), [LAM])
    inc = _incident(handle, theta, azimuth, (0.3 + 0.4j, -0.5))
    out = fga.apply(handle, inc)
    period_mm = 0.7e-3
    for side in ("front", "back"):
        waves = out[side]
        assert len(waves.kx) >= 2 and waves.orders is not None
        for j, (n, m) in enumerate(waves.orders):
            assert m == 0
            assert abs(waves.kx[j] - (inc.kx[0] + n * 2 * np.pi / period_mm)) < 1e-6
            assert abs(waves.ky[j] - inc.ky[0]) < 1e-6
    p_in = float(inc.power_density()[0])
    p_out = float(np.sum(out["front"].power_density()) + np.sum(out["back"].power_density()))
    assert abs(p_out / p_in - 1.0) < 1e-9  # lossless grating


def test_reciprocity_of_planar_transmission() -> None:
    # vacuum on both sides, asymmetric absorbing stack: t equal from both sides, r not
    cfg_up = _config(STACK_CFG, sub="air")
    cfg_down = _config(list(reversed(STACK_CFG)), sub="air")
    results = []
    for cfg in (cfg_up, cfg_down):
        handle = fga.prepare(fga.Structure(cfg), [LAM])
        results.append(fga.apply(handle, _incident(handle, 0.5, 0.7, (0.6, 0.3j))))
    t = [_zero_order(r["back"]) for r in results]
    r = [_zero_order(x["front"]) for x in results]
    assert np.max(np.abs(t[0] - t[1])) < 1e-10
    assert np.max(np.abs(r[0] - r[1])) > 1e-3


def test_hash_and_errors() -> None:
    a = fga.prepare(fga.Structure(_config(GRATING)), [LAM])
    b = fga.prepare(fga.Structure(_config(GRATING)), [LAM])
    c = fga.prepare(fga.Structure(_config(STACK_CFG)), [LAM])
    assert a.content_hash == b.content_hash != c.content_hash
    with pytest.raises(ValueError):
        fga.prepare(fga.Structure(_config(GRATING, sub="lossy")), [LAM])
    with pytest.raises(ValueError):
        fga.apply(
            a,
            PlaneWaves(
                fga.front_port(a, LAM, -1), np.zeros(1), np.zeros(1), np.ones((1, 2), complex)
            ),
        )


def test_cross_check_with_bwpm() -> None:
    # same planar stack with both adapters: layers multiple of dz, BWPM's back port lies one slab
    # (substrate) behind the last interface
    bwpm_adapter = pytest.importorskip("bwpm.adapter")
    dz = LAM * 1e-3 / 40
    layers_dz = [(2.1, 6), (1.45, 10)]  # top down: thickness in slabs
    cfg = _config(
        [
            {"type": "homogeneous", "thickness_um": n * dz * 1e3, "material": m}
            for m, n in (("mid", 10), ("hi", 6))
        ]
    )
    fg_handle = fga.prepare(fga.Structure(cfg), [LAM])
    media = [complex(n) for n, count in layers_dz for _ in range(count)] + [1.5 + 0j]
    unique = sorted(set(media), key=lambda c: c.real)
    n_grid = 80  # |k⊥| = k0·sqrt(mx² + my²)·40/n_grid
    labels = np.array([np.full((n_grid, n_grid), unique.index(m)) for m in media])
    structure = bwpm_adapter.Structure(n_grid, n_grid, dz, dz, dz, tuple(unique), labels, 1.0 + 0j)
    bw_handle = bwpm_adapter.prepare(structure, [LAM])
    kx = 2 * np.pi * 1 / (n_grid * dz)  # 0.5 k0
    ky = 2 * np.pi * 1 / (n_grid * dz)
    e = (0.5 + 0.1j, -0.2 + 0.4j)
    bw_out = bwpm_adapter.apply(
        bw_handle,
        bwpm_adapter.plane_wave(bwpm_adapter.front_port(bw_handle, LAM), 1, 1, e),
        iterations=60,
        backend="auto",
    )
    fg_inc = PlaneWaves(
        fga.front_port(fg_handle, LAM),
        np.array([kx]),
        np.array([ky]),
        np.array([e], dtype=np.complex128),
    )
    fg_out = fga.apply(fg_handle, fg_inc)
    centre = (n_grid // 2, n_grid // 2)
    kz_sub = np.sqrt(complex((1.5 * K0) ** 2 - kx * kx - ky * ky))
    assert np.max(np.abs(bw_out["front"].field[centre] - _zero_order(fg_out["front"]))) < 1e-9
    assert (
        np.max(
            np.abs(
                bw_out["back"].field[centre]
                - _zero_order(fg_out["back"]) * np.exp(1j * kz_sub * dz)
            )
        )
        < 1e-9
    )
