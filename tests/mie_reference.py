"""Mie series for a homogeneous sphere (Bohren & Huffman, Absorption and Scattering of Light by
Small Particles, Wiley 1983, ch. 4): coefficients a_n, b_n, amplitudes S1, S2, Q_ext.

Time factor e^(−iωt); for x-polarised incidence along +z the far field is
E_θ = (i/k) S2 cos φ, E_φ = −(i/k) S1 sin φ (E_s ≈ E∞ e^(ikr)/r). Independent of every module.
"""

from __future__ import annotations

import numpy as np
from scipy.special import spherical_jn, spherical_yn


def mie_ab(m: complex, x: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    nmax = int(x + 4 * x ** (1 / 3) + 2)
    n = np.arange(1, nmax + 1)
    mx = m * x
    jn, yn = spherical_jn(n, x), spherical_yn(n, x)
    jnp, ynp = spherical_jn(n, x, derivative=True), spherical_yn(n, x, derivative=True)
    hn, hnp = jn + 1j * yn, jnp + 1j * ynp
    jm, jmp = spherical_jn(n, mx), spherical_jn(n, mx, derivative=True)
    psi_x, dpsi_x = x * jn, jn + x * jnp
    xi_x, dxi_x = x * hn, hn + x * hnp
    psi_m, dpsi_m = mx * jm, jm + mx * jmp
    a = (m * psi_m * dpsi_x - psi_x * dpsi_m) / (m * psi_m * dxi_x - xi_x * dpsi_m)
    b = (psi_m * dpsi_x - m * psi_x * dpsi_m) / (psi_m * dxi_x - m * xi_x * dpsi_m)
    return n, a, b


def S12(m: complex, x: float, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    n, a, b = mie_ab(m, x)
    mu = np.cos(theta)
    pi = np.zeros((len(n), np.size(theta)))
    tau = np.zeros_like(pi)
    pi[0], tau[0] = 1.0, mu
    if len(n) > 1:
        pi[1] = 3 * mu
        tau[1] = 2 * mu * pi[1] - 3 * pi[0]
    for k in range(2, len(n)):
        nn = k + 1
        pi[k] = ((2 * nn - 1) * mu * pi[k - 1] - nn * pi[k - 2]) / (nn - 1)
        tau[k] = nn * mu * pi[k] - (nn + 1) * pi[k - 1]
    f = ((2 * n + 1) / (n * (n + 1)))[:, None]
    s1 = np.sum(f * (a[:, None] * pi + b[:, None] * tau), axis=0)
    s2 = np.sum(f * (a[:, None] * tau + b[:, None] * pi), axis=0)
    return s1, s2


def qext(m: complex, x: float) -> float:
    n, a, b = mie_ab(m, x)
    return float(2 / x**2 * np.sum((2 * n + 1) * np.real(a + b)))
