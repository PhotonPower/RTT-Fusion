"""R4: conformance of the CliffordBem adapter (sphere, flat elements).

References: the plane-wave solution of CliffordBem itself (phase and polarisation mapping of the
adapter), Mie series written here (tests/mie_reference.py, Bohren & Huffman) as plausibility at
the accuracy of flat elements, energy (σ_sca = σ_ext for a lossless sphere) and translation
covariance of the port amplitudes.
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pytest

try:
    import cliffordbem as cb
except ImportError:
    if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
        raise
    pytest.skip("cliffordbem not installed", allow_module_level=True)

from mie_reference import S12, qext

from rtt_fusion.adapters import cliffordbem as cba
from rtt_fusion.ports import PlaneWaves, Port, spectrum_power

LAM, R, N_P = 0.6, 0.25, 1.5  # µm
K_UM = 2 * np.pi / LAM
CENTER = (0.0, 0.0, 0.002)  # mm


def _mesh() -> Any:
    m0 = cb.make_icosphere(6)
    return cb.TriangleMesh(np.asarray(m0.points) * R, m0.triangles)


@pytest.fixture(scope="module")
def handle() -> Any:
    return cba.prepare(cba.Structure(_mesh(), N_P**2, 1.0, CENTER, 0.0, 0.004), [LAM])


def _plane(h: Any) -> PlaneWaves:
    return PlaneWaves(
        cba.front_port(h, LAM),
        np.array([0.0]),
        np.array([0.0]),
        np.array([[1.0, 0.0]], dtype=np.complex128),
    )


@pytest.fixture(scope="module")
def solution(handle: Any) -> Any:
    return cba.scatter(handle, _plane(handle))


def test_plane_wave_equals_cliffordbem(handle: Any, solution: Any) -> None:
    # the adapter's incident wave has phase 1 at the front port, i.e. e^(ik·2 µm) at the centre
    problem = cba._problem(handle, LAM)
    ref = problem.solve_plane_wave(d=(0, 0, 1), p=(1, 0, 0), options=cb.SolveOptions(tol=1e-10))
    phase = np.exp(1j * K_UM * CENTER[2] * 1000)
    assert np.max(np.abs(solution.h - phase * np.asarray(ref.h))) < 1e-7 * np.max(np.abs(ref.h))


def test_mie_plausibility_and_energy(handle: Any, solution: Any) -> None:
    phase = np.exp(1j * K_UM * CENTER[2] * 1000)
    forward = (
        cb.far_field_E(handle.structure.mesh, solution.h_scat, K_UM, 1.0, np.array([[0, 0, 1.0]]))[
            0
        ]
        / phase
    )
    sigma_ext = 4 * np.pi / K_UM * np.imag(forward[0])
    sigma_mie = qext(N_P, K_UM * R) * np.pi * R**2
    assert abs(sigma_ext / sigma_mie - 1) < 0.04  # flat elements, 720 triangles: 3 %
    theta = np.array([0.3, 1.0])
    _, s2 = S12(N_P, K_UM * R, theta)
    for t, s in zip(theta, s2, strict=True):
        d = np.array([np.sin(t), 0.0, np.cos(t)])
        e = cb.far_field_E(handle.structure.mesh, solution.h_scat, K_UM, 1.0, d[None])[0] / phase
        e_theta = e @ np.array([np.cos(t), 0.0, -np.sin(t)])
        assert abs(e_theta - 1j / K_UM * s) < 0.04 * abs(1j / K_UM * s)
    out = cba.apply(handle, _plane(handle))
    # scattered power from the port spectra (W in the units of power_density per unit intensity)
    p_sca = spectrum_power(out["front"]) + spectrum_power(out["back_scattered"])
    k_mm = K_UM * 1000
    sigma_sca = p_sca / k_mm * 1e6  # intensity k|E|² = k_mm; area mm² → µm²
    assert abs(sigma_sca / sigma_ext - 1) < 0.02  # lossless sphere


def test_translation_covariance(handle: Any) -> None:
    out = cba.apply(handle, _plane(handle), 16, 32)
    shift_mm = np.array([0.3e-3, -0.1e-3, 0.4e-3])
    moved = cba.prepare(
        cba.Structure(_mesh(), N_P**2, 1.0, tuple(np.add(CENTER, shift_mm)), 0.0, 0.004), [LAM]
    )
    out2 = cba.apply(moved, _plane(moved), 16, 32)
    k_in = np.array([0.0, 0.0, K_UM * 1000])
    for side in ("front", "back_scattered"):
        a, b = out[side], out2[side]
        kz = np.real(a.port.kz(a.kx, a.ky)) * a.port.direction
        k_out = np.column_stack([a.kx, a.ky, kz])
        factor = np.exp(1j * ((k_in[None, :] - k_out) @ shift_mm))
        expected = a.field * factor[:, None]
        assert np.max(np.abs(b.field - expected)) < 1e-5 * np.max(np.abs(a.field))


def test_back_port_carries_the_incident_wave(handle: Any) -> None:
    out = cba.apply(handle, _plane(handle), 8, 16)
    back = out["back"]
    assert abs(back.field[0, 0] - np.exp(1j * K_UM * 1000 * 0.004)) < 1e-12
    with pytest.raises(ValueError):
        cba.prepare(cba.Structure(_mesh(), N_P**2, 1.0, CENTER, 0.0019, 0.004), [LAM])


# ---------------------------------------------------------------------------- L4 / L5


def _incident(port: Any, kx: float, ky: float, e: tuple[complex, complex]) -> PlaneWaves:
    return PlaneWaves(port, np.array([kx]), np.array([ky]), np.array([e], dtype=np.complex128))


def _random_incidents(h: Any, n: int, seed: int = 7) -> list[PlaneWaves]:
    rng = np.random.default_rng(seed)
    k_mm = K_UM * 1000
    out = []
    for _ in range(n):
        s = 0.3 * np.sqrt(rng.random())  # sin θ up to 0.3
        phi = 2 * np.pi * rng.random()
        e = rng.normal(size=2) + 1j * rng.normal(size=2)
        out.append(
            _incident(
                cba.front_port(h, LAM), s * k_mm * np.cos(phi), s * k_mm * np.sin(phi), (e[0], e[1])
            )
        )
    return out


def _fresh() -> Any:
    return cba.prepare(cba.Structure(_mesh(), N_P**2, 1.0, CENTER, 0.0, 0.004), [LAM])


def test_l4_apply_many_equals_apply() -> None:
    plain, many = _fresh(), _fresh()
    incidents = _random_incidents(plain, 20)
    one_by_one = [cba.apply(plain, inc, 12, 24) for inc in incidents]
    exact = cba.apply_many(_fresh(), incidents, 12, 24, recycle=False)
    recycled = cba.apply_many(many, incidents, 12, 24)
    assert LAM in many._recycling and LAM not in plain._recycling
    for a, b, c in zip(one_by_one, exact, recycled, strict=True):
        for key in a:
            scale = np.abs(a[key].field).max()
            # the same solves; CliffordBem's H-matrix sums in an unspecified order with several
            # threads (its docs: bitwise only with OMP_NUM_THREADS=1), so to the solver tolerance
            assert np.abs(a[key].field - b[key].field).max() < 1e-6 * scale
            assert np.abs(a[key].field - c[key].field).max() < 1e-6 * scale  # solver tolerance


def test_l5_light_from_the_back_mirrors_the_front() -> None:
    # the sphere lies midway between the ports and its mesh is symmetric under z → −z: light
    # coming in at the back (−z) gives the mirror image of light coming in at the front (+z)
    h = _fresh()
    front = cba.apply(h, _incident(cba.front_port(h, LAM), 300.0, -200.0, (1.0, 0.5j)), 16, 32)
    back_port = Port(0.004, -1, 1.0 + 0j, LAM)
    back = cba.apply(h, _incident(back_port, 300.0, -200.0, (1.0, 0.5j)), 16, 32)
    assert sorted(back) == ["back", "front", "front_scattered"]
    assert back["back"].port.direction == 1 and back["front"].port.direction == -1
    for f_key, b_key in (
        ("front", "back"),
        ("back", "front"),
        ("back_scattered", "front_scattered"),
    ):
        a, b = front[f_key], back[b_key]
        assert np.allclose(a.kx, b.kx) and np.allclose(a.ky, b.ky)
        assert np.abs(a.field - b.field).max() < 1e-4 * np.abs(a.field).max()
    with pytest.raises(ValueError, match="back port with direction -1"):
        cba.apply(h, _incident(Port(0.004, 1, 1.0 + 0j, LAM), 0.0, 0.0, (1.0, 0.0)))
