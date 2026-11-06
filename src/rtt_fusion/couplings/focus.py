"""Focused vector beam as plane waves (Richards–Wolf), for pilots K3 and the objective side.

Aplanatic objective in a medium of index n, numerical aperture NA, pupil polarisation e_in
(x, y, 0), uniform pupil illumination. Angular spectrum per solid angle (CliffordBem
docs/results_beams.md, Richards–Wolf; same form as Fieldtracing's Debye integral):

    a(θ, φ) = √cos θ [(e_in·φ̂) φ̂ + (e_in·ρ̂) θ̂],  θ ≤ arcsin(NA/n),

E(r) = ∫ a e^(ik k̂·(r − r_f)) dΩ. Discretised with Gauss–Legendre in cos θ on
[cos θ_max, 1] and a uniform φ grid; the waves are referred to the port origin (0, 0, z_port), so
the beam focuses at r_f.
"""

from __future__ import annotations

import numpy as np

from rtt_fusion.ports import PlaneWaves, Port


def richards_wolf(
    port: Port,
    na: float,
    focus_mm: tuple[float, float, float],
    pupil_polarisation: tuple[complex, complex] = (1.0, 0.0),
    n_theta: int = 48,
    n_phi: int = 96,
) -> PlaneWaves:
    """Plane waves of the focused beam on the port (direction +1, medium of the port)."""
    n = float(np.real(port.index))
    if not 0 < na < n:
        raise ValueError("0 < NA < n required")
    cmax = np.sqrt(1 - (na / n) ** 2)
    x, w = np.polynomial.legendre.leggauss(n_theta)
    mu = cmax + (1 - cmax) * 0.5 * (x + 1)
    wmu = (1 - cmax) * 0.5 * w
    phi = 2 * np.pi * (np.arange(n_phi) + 0.5) / n_phi
    mm, pp = np.meshgrid(mu, phi, indexing="ij")
    ct, st = mm.ravel(), np.sqrt(1 - mm.ravel() ** 2)
    cp, sp = np.cos(pp.ravel()), np.sin(pp.ravel())
    weights = (wmu[:, None] * np.full(n_phi, 2 * np.pi / n_phi)[None, :]).ravel()
    e_in = np.array([pupil_polarisation[0], pupil_polarisation[1], 0.0], dtype=complex)
    phi_hat = np.column_stack([-sp, cp, np.zeros_like(cp)])
    rho_hat = np.column_stack([cp, sp, np.zeros_like(cp)])
    theta_hat = np.column_stack([ct * cp, ct * sp, -st])
    a = np.sqrt(ct)[:, None] * (
        (phi_hat @ e_in)[:, None] * phi_hat + (rho_hat @ e_in)[:, None] * theta_hat
    )
    k = n * port.k0
    dirs = np.column_stack([st * cp, st * sp, ct])
    shift = np.asarray(focus_mm) - np.array([0.0, 0.0, port.z])
    amp = a * weights[:, None] * np.exp(-1j * k * (dirs @ shift))[:, None]
    return PlaneWaves(
        port, dirs[:, 0] * k, dirs[:, 1] * k, amp[:, :2].astype(np.complex128), solid_angle=weights
    )
