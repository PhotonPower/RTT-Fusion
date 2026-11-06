"""Conversions between port representations that do not involve a solver (for the C++ port).

`hemisphere`: quadrature of directions; `far_field_to_spectrum`: radiated far field of a body →
plane waves at a port (Weyl expansion, docs/konventionen.md and adapter CliffordBem).
"""

from __future__ import annotations

import numpy as np

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.ports import PlaneWaves, Port


def hemisphere(n_theta: int, n_phi: int, forward: bool) -> tuple[FloatArray, FloatArray]:
    """Unit directions (N, 3) and weights dΩ (N,): Gauss–Legendre in cos θ on (0, 1], uniform φ
    with nodes at (j + ½)·2π/n_phi; θ measured from +z (forward) or −z (backward)."""
    x, w = np.polynomial.legendre.leggauss(n_theta)
    mu, wmu = 0.5 * (x + 1.0), 0.5 * w
    phi = 2 * np.pi * (np.arange(n_phi) + 0.5) / n_phi
    mm, pp = np.meshgrid(mu, phi, indexing="ij")
    st = np.sqrt(1 - mm**2)
    sign = 1.0 if forward else -1.0
    dirs = np.column_stack(
        [(st * np.cos(pp)).ravel(), (st * np.sin(pp)).ravel(), sign * mm.ravel()]
    )
    weights = (wmu[:, None] * np.full(n_phi, 2 * np.pi / n_phi)[None, :]).ravel()
    return dirs, weights


def far_field_to_spectrum(
    port: Port,
    center_mm: FloatArray,
    directions: FloatArray,
    weights: FloatArray,
    e_inf: ComplexArray,
) -> PlaneWaves:
    """Field radiated by a body at `center_mm`, given by its far field E∞ (E_s ≈ E∞ e^(ikr)/r, r in
    mm) on quadrature directions of the hemisphere towards the port, as plane waves at the port.

    E_s(r) = (ik/2π) ∫ E∞(k̂) e^(ik k̂·(r − c)) dΩ (Weyl; evanescent waves left out, exact for what
    leaves through NA < n). Wave j: amplitude (ik/2π)·w_j·E∞⊥(k̂_j)·e^(ik k̂_j·(P₀ − c)) with
    P₀ = (0, 0, z_port) and E∞⊥ the part transverse to k̂_j (removes numerical radial parts that
    E_z = −(k_x E_x + k_y E_y)/k_z would amplify near grazing). k = n k₀ of the port medium.
    """
    dirs = np.asarray(directions, dtype=np.float64)
    if np.any(np.sign(dirs[:, 2]) != port.direction):
        raise ValueError("directions must point through the port (sign of k_z = direction)")
    e = np.asarray(e_inf, dtype=np.complex128)
    e = e - np.sum(e * dirs, axis=1, keepdims=True) * dirs
    k = float(np.real(port.index)) * port.k0
    p0 = np.array([0.0, 0.0, port.z]) - np.asarray(center_mm, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    amp = (1j * k / (2 * np.pi)) * w[:, None] * e * np.exp(1j * k * (dirs @ p0))[:, None]
    return PlaneWaves(
        port, dirs[:, 0] * k, dirs[:, 1] * k, amp[:, :2].astype(np.complex128), solid_angle=w
    )
