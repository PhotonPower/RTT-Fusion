"""R8: conformance of the SpeckleBem adapter (closed bodies), the solver with exp(+jωt).

The conventions are the risk here: time factor (conjugation of ε, incident amplitude and far
field), SI units and phase reference. References: Mie series (tests/mie_reference.py) for an
absorbing sphere, the sign of the absorption (a wrong conjugation turns absorption into gain),
the CliffordBem adapter on the same sphere (second BEM, other conventions), and linearity of
several incident waves.
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pytest

try:
    import specklebem as sb
except ImportError:
    if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
        raise
    pytest.skip("specklebem not installed", allow_module_level=True)

from mie_reference import S12, qext

from rtt_fusion.adapters import specklebem as sba
from rtt_fusion.ports import PlaneWaves, Port, spectrum_power

LAM, R = 0.6, 0.25  # µm
EPS = 2.25 + 0.3j  # absorbing, our convention
K_UM = 2 * np.pi / LAM
CENTER = (0.0, 0.0, 0.002)


@pytest.fixture(scope="module")
def handle() -> Any:
    mesh = sb.make_sphere(R * 1e-6, LAM * 1e-6 / 10)
    return sba.prepare(sba.Structure(mesh, EPS, 1.0, CENTER, 0.0, 0.004), [LAM])


def _wave(h: Any, e: tuple[complex, complex], kx: float = 0.0) -> PlaneWaves:
    return PlaneWaves(
        sba.front_port(h, LAM), np.array([kx]), np.array([0.0]), np.array([e], dtype=np.complex128)
    )


def test_describe_and_input_convention() -> None:
    d = sba.describe()
    assert "exp(+j omega t)" in d["conventions_of_solver"][0]
    mesh = sb.make_sphere(R * 1e-6, LAM * 1e-6 / 6)
    with pytest.raises(ValueError):  # ε in the convention of SpeckleBem is refused
        sba.prepare(sba.Structure(mesh, np.conj(EPS), 1.0, CENTER, 0.0, 0.004), [LAM])


def test_far_field_against_mie_absorbing_sphere(handle: Any) -> None:
    pol = np.array([1.0, 0.3j]) / np.linalg.norm([1.0, 0.3j])
    sol = sba.scatter(handle, _wave(handle, (pol[0], pol[1])))
    theta = np.array([0.3, 1.0])
    dirs = np.column_stack([np.sin(theta), 0 * theta, np.cos(theta)])
    phase = np.exp(1j * K_UM * CENTER[2] * 1000)  # incident phase at the centre
    e = sol.far_field(dirs) * 1e3 / phase  # r in µm for the comparison with Mie (1/k in µm)
    s1, s2 = S12(np.sqrt(EPS), K_UM * R, theta)
    for j, t in enumerate(theta):
        e_theta = e[j] @ np.array([np.cos(t), 0.0, -np.sin(t)])
        e_phi = e[j] @ np.array([0.0, 1.0, 0.0])
        assert abs(e_theta - 1j / K_UM * s2[j] * pol[0]) < 0.03 * abs(1j / K_UM * s2[j])
        assert abs(e_phi - 1j / K_UM * s1[j] * pol[1]) < 0.03 * abs(1j / K_UM * s1[j] * pol[1])


def test_absorption_is_positive(handle: Any) -> None:
    out = sba.apply(handle, _wave(handle, (1.0, 0.0)), 24, 48)
    k_mm = K_UM * 1000
    sigma_sca = (spectrum_power(out["front"]) + spectrum_power(out["back_scattered"])) / k_mm * 1e6
    sol = sba.scatter(handle, _wave(handle, (1.0, 0.0)))
    forward = (
        sol.far_field(np.array([[0.0, 0.0, 1.0]]))[0] * 1e3 / np.exp(1j * K_UM * CENTER[2] * 1000)
    )
    sigma_ext = 4 * np.pi / K_UM * np.imag(forward[0])
    mie_ext = qext(np.sqrt(EPS), K_UM * R) * np.pi * R**2
    assert abs(sigma_ext / mie_ext - 1) < 0.03
    # without the conjugation the sphere would amplify: σ_abs = σ_ext − σ_sca < 0
    assert sigma_ext - sigma_sca > 0.2 * sigma_ext


def test_cross_check_with_cliffordbem(handle: Any) -> None:
    cb = pytest.importorskip("cliffordbem")
    from rtt_fusion.adapters import cliffordbem as cba

    m0 = cb.make_icosphere(8)
    mesh = cb.TriangleMesh(np.asarray(m0.points) * R, m0.triangles)
    cb_handle = cba.prepare(cba.Structure(mesh, EPS, 1.0, CENTER, 0.0, 0.004), [LAM])
    kx = 0.3 * K_UM * 1000  # oblique, 17.5°
    e = (0.8, 0.2 + 0.4j)
    a = sba.apply(handle, _wave(handle, e, kx), 16, 32)
    b = cba.apply(
        cb_handle,
        PlaneWaves(
            cba.front_port(cb_handle, LAM),
            np.array([kx]),
            np.array([0.0]),
            np.array([e], dtype=np.complex128),
        ),
        16,
        32,
    )
    for side in ("front", "back_scattered"):
        scale = np.max(np.abs(b[side].field))
        assert np.max(np.abs(a[side].field - b[side].field)) < 0.04 * scale
        assert abs(spectrum_power(a[side]) / spectrum_power(b[side]) - 1) < 0.04


def test_superposition_of_two_waves(handle: Any) -> None:
    e1, e2 = (1.0, 0.0), (0.0, 0.5j)
    kx = 0.2 * K_UM * 1000
    port = sba.front_port(handle, LAM)
    both = PlaneWaves(
        port, np.array([0.0, kx]), np.array([0.0, 0.0]), np.array([e1, e2], dtype=np.complex128)
    )
    dirs = np.array([[0.0, 0.0, 1.0], [0.6, 0.0, -0.8]])
    sum_ff = sba.scatter(handle, both).far_field(dirs)
    one = sba.scatter(handle, _wave(handle, e1)).far_field(dirs)
    two = sba.scatter(
        handle,
        PlaneWaves(port, np.array([kx]), np.array([0.0]), np.array([e2], dtype=np.complex128)),
    ).far_field(dirs)
    assert np.max(np.abs(sum_ff - one - two)) < 1e-6 * np.max(np.abs(sum_ff))


# ---------------------------------------------------------------------------- L4 / L5


def _incident(port: Any, kx: float, ky: float, e: tuple[complex, complex]) -> PlaneWaves:
    return PlaneWaves(port, np.array([kx]), np.array([ky]), np.array([e], dtype=np.complex128))


def test_l4_apply_many_equals_apply(handle: Any) -> None:
    rng = np.random.default_rng(11)
    k_mm = K_UM * 1000
    incidents = []
    for _ in range(5):
        s, phi = 0.3 * np.sqrt(rng.random()), 2 * np.pi * rng.random()
        e = rng.normal(size=2) + 1j * rng.normal(size=2)
        incidents.append(
            _incident(
                sba.front_port(handle, LAM),
                s * k_mm * np.cos(phi),
                s * k_mm * np.sin(phi),
                (e[0], e[1]),
            )
        )
    many = sba.apply_many(handle, incidents, 12, 24)
    for inc, b in zip(incidents, many, strict=True):
        a = sba.apply(handle, inc, 12, 24)
        for key in a:
            assert np.abs(a[key].field - b[key].field).max() < 1e-6 * np.abs(a[key].field).max()


def test_l5_light_from_the_back_mirrors_the_front(handle: Any) -> None:
    # the sphere lies midway between the ports: light coming in at the back (−z) gives the
    # mirror image of light coming in at the front (+z), up to the mesh, which need not be
    # symmetric under z → −z (a wrong direction would differ by order one)
    front = sba.apply(handle, _incident(sba.front_port(handle, LAM), 0.0, 0.0, (1.0, 0.0)), 16, 32)
    back = sba.apply(
        handle, _incident(Port(0.004, -1, 1.0 + 0j, LAM), 0.0, 0.0, (1.0, 0.0)), 16, 32
    )
    assert sorted(back) == ["back", "front", "front_scattered"]
    assert back["back"].port.direction == 1 and back["front"].port.direction == -1
    for f_key, b_key in (("front", "back"), ("back_scattered", "front_scattered")):
        a, b = front[f_key], back[b_key]
        assert np.allclose(a.kx, b.kx) and np.allclose(a.ky, b.ky)
        assert np.abs(a.field - b.field).max() < 2e-2 * np.abs(a.field).max()
        pa, pb = a.power_density().sum(), b.power_density().sum()
        assert abs(pb / pa - 1) < 2e-2
    with pytest.raises(ValueError, match="back port with direction -1"):
        sba.apply(handle, _incident(Port(0.004, 1, 1.0 + 0j, LAM), 0.0, 0.0, (1.0, 0.0)))
