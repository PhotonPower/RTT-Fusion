"""A dipole emitter in a planar multilayer (pilot K5, stage K5a, `couplings.dipole_layers`).

References, independent of the code under test (W94):
- (a) homogeneous medium (an interface between equal indices): P/P₀ = 1, half up, half down;
  the angular spectrum is the free dipole's far field E∞ = k₀²/(4π) (I − k̂k̂) p e^(ikr)/r written
  out, as A = −k E∞ r e^(−ikr)/(2πi) with the phase of the dipole's offset from the interface;
- (b) a near-perfect mirror (n = i·10⁵): the total rate of the image dipole (p_img =
  (−p_x, −p_y, p_z) at −d, its full near field written out; Drexhage, Chance–Prock–Silbey), the
  error falling like 1/|n| (skin depth);
- (c) energy: the near-field rate (Chance–Prock–Silbey, from the stack's field at the dipole)
  equals the power leaving through both half-spaces (far field) for lossless stacks – water/glass
  and an adlayer with the dipole inside it – to 1e-9; with a gold film the difference is the
  absorption (> 0);
- (d) the power of the angular spectrum in the glass, integrated over the solid angle with its
  own quadrature, equals the flux into the glass (also the share within an NA);
- (e) the κ quadrature has converged (doubling changes < 1e-6 relative, gold film included).
"""

from __future__ import annotations

import numpy as np
import pytest

from rtt_fusion._types import FloatArray
from rtt_fusion.couplings import dipole_layers as dl

LAM = 0.6
K0 = 2 * np.pi / (LAM * 1e-3)
PX, PZ = [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]


def _directions(theta: FloatArray, phi: FloatArray) -> FloatArray:
    return np.stack([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), -np.cos(theta)], -1)


@pytest.mark.parametrize("n", [1.0, 1.33])
def test_w94_homogeneous_medium_is_the_free_dipole(n: float) -> None:
    st = dl.Stack(n, (), n, 0.0)
    for p in (PX, PZ, [0.3, -0.5, 0.8]):
        r = dl.rates(st, LAM, 2e-4, p)
        assert r.total == pytest.approx(1.0, abs=1e-9)
        assert r.up == pytest.approx(0.5, abs=1e-9) and r.down == pytest.approx(0.5, abs=1e-9)
    pv = np.array([0.3, -0.5, 0.8])
    z0, k = 3e-4, n * K0
    theta, phi = np.array([0.1, 0.7, 1.2, 1.5]), np.array([0.3, 2.0, 4.0, 5.5])
    a = dl.angular_spectrum(dl.Stack(n, (), n, 0.0), LAM, z0, pv, theta, phi, "bottom")
    kh = _directions(theta, phi)
    e_inf = K0**2 / (4 * np.pi) * (pv[None, :] - kh * (kh @ pv)[:, None])  # × e^(ikr)/r
    want = -k * e_inf / (2j * np.pi) * np.exp(1j * k * np.cos(theta) * z0)[:, None]
    assert np.abs(a - want).max() < 1e-12 * np.abs(want).max()


def _image_rate(p: list[float], d: float) -> float:
    """1 + 6π/(k³|p|²) Im(p*·E_img(r0)), E_img the full field of the image dipole."""
    pv = np.asarray(p)
    img = pv * np.array([-1.0, -1.0, 1.0])
    r, nh = 2 * d, np.array([0.0, 0.0, 1.0])
    e = (
        (
            K0**2 * np.cross(np.cross(nh, img), nh) / r
            + (3 * nh * (nh @ img) - img) * (1 / r**3 - 1j * K0 / r**2)
        )
        * np.exp(1j * K0 * r)
        / (4 * np.pi)
    )
    return float(1 + 6 * np.pi / (K0**3 * (pv @ pv)) * np.imag(pv @ e))


@pytest.mark.parametrize("p", [PX, PZ])
def test_w94_mirror_is_the_image_dipole(p: list[float]) -> None:
    for d in (5e-5, 2e-4, 1e-3):
        want = _image_rate(p, d)
        errs = [dl.rates(dl.Stack(n, (), 1.0, 0.0), LAM, d, p).total - want for n in (3e4j, 1e5j)]
        assert abs(errs[1]) < 1e-4  # measured ≤ 8.3e-5 (vertical), 9.4e-6 (horizontal)
        assert abs(errs[1]) < 0.4 * abs(errs[0])  # ~1/|n|: the skin depth, not the method


STACKS = {
    "water/glass, 10 nm above": (dl.Stack(1.52, (), 1.33, 0.0), 1e-5),
    "water/glass, 1 µm above": (dl.Stack(1.52, (), 1.33, 0.0), 1e-3),
    "in an adlayer": (dl.Stack(1.52, ((1.45, 6e-6), (1.40, 4e-6)), 1.33, 0.0), 3e-6),
    "below a film, in glass": (dl.Stack(1.52, ((2.0, 1e-4),), 1.0, 0.0), -5e-5),
}


@pytest.mark.parametrize("case", list(STACKS))
def test_w94_near_field_rate_equals_radiated_power(case: str) -> None:
    st, z0 = STACKS[case]
    for p in (PX, PZ):
        r = dl.rates(st, LAM, z0, p)
        assert abs(r.absorbed) < 1e-9 * r.total, (case, p, r)


def test_w94_gold_film_absorbs_and_converges() -> None:
    st = dl.Stack(1.515, ((0.197 + 3.44j, 50e-6),), 1.333, -50e-6)
    for p in (PX, PZ):
        a = dl.rates(st, 0.633, 1e-5, p, na=1.49)
        b = dl.rates(st, 0.633, 1e-5, p, na=1.49, n_per=2 * dl.N_KAPPA)
        assert a.total == pytest.approx(b.total, rel=1e-6)
        assert a.down == pytest.approx(b.down, rel=1e-6)
        assert a.absorbed > 0.5 * a.total  # 10 nm above gold: quenched
        assert 0 < a.collected <= a.down


def _spectrum_power(st: dl.Stack, z0: float, p: list[float], s_lo: float, s_hi: float) -> float:
    """(n/2)(2π/k)² ∫ |A|² dΩ over the glass for sin θ ∈ [s_lo, s_hi] (Gauss–Legendre in θ,
    uniform φ), / P₀."""
    n = float(np.real(st.n_bottom))
    k = n * K0
    x, w = np.polynomial.legendre.leggauss(200)
    t_lo, t_hi = np.arcsin(s_lo), np.arcsin(s_hi)
    t = t_lo + (t_hi - t_lo) * (x + 1) / 2
    wt = w * (t_hi - t_lo) / 2 * np.sin(t)
    nphi = 8
    phi = np.arange(nphi) * 2 * np.pi / nphi
    tt, pp = np.meshgrid(t, phi, indexing="ij")
    a = dl.angular_spectrum(st, LAM, z0, p, tt.ravel(), pp.ravel(), "bottom")
    dens = np.sum(np.abs(a) ** 2, -1).reshape(tt.shape)
    power = n / 2 * (2 * np.pi / k) ** 2 * np.sum(dens * wt[:, None]) * 2 * np.pi / nphi
    n_src = float(np.real(st.indices[st.medium_at(z0)]))
    return float(power / (K0**4 * n_src * float(np.dot(p, p)) / (12 * np.pi)))


def test_w94_angular_spectrum_carries_the_flux() -> None:
    """Sub- and supercritical light (split at the critical angle, where |A| has a kink)."""
    st, z0 = dl.Stack(1.52, (), 1.33, 0.0), 2e-5
    crit = 1.33 / 1.52
    for p in (PX, PZ):
        r = dl.rates(st, LAM, z0, p, na=1.3)
        sub = _spectrum_power(st, z0, p, 0.0, crit)
        sup = _spectrum_power(st, z0, p, crit, 1.0)
        assert sub + sup == pytest.approx(r.down, rel=1e-6)
        within = sub + _spectrum_power(st, z0, p, crit, 1.3 / 1.52)
        assert within == pytest.approx(r.collected, rel=1e-6)
