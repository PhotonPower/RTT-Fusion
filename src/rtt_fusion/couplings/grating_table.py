"""Coupling K1: diffraction efficiencies from FourierGrating as a table in Raytatouille (stage R3).

Raytatouille models a grating as a thin phase layer (ADR 0025): the order m leaves along the
local grating equation, and an optional `diffraction_efficiency` list gives the power fraction η
of each order. η is a plain number per surface: independent of wavelength, angle of incidence and
polarisation (ADR 0025, point 5), and it enters only the ray weight, after the surface's own
interaction (Fresnel, coating, ideal) of order 0.

This module computes η per transmitted order with the FourierGrating adapter for one wavelength,
one direction of incidence and one polarisation (or the mean of s and p for unpolarised light)
and writes it into a Raytatouille system with a JSON patch. Order m of Raytatouille is order n = m
of FourierGrating when the grating vector of the surface is the local x axis (orientation 0) and
FourierGrating's x axis is the same x axis (frame of the adapter: x′ = x); test K1 checks this.

Use it where the surface's own interaction transmits fully (grating between equal media, as in
the pilot); otherwise η must be divided by the order-0 power factor of that interaction.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from rtt_fusion.adapters import fouriergrating as fga
from rtt_fusion.ports import PlaneWaves


def _unit_polarisation(kx: float, ky: float, kz: float, pol: str) -> tuple[complex, complex]:
    """Transverse (E_x, E_y) of a unit s or p wave (Convention A, N = +z)."""
    kp = np.hypot(kx, ky)
    if kp == 0:  # normal incidence: s along y, p along x (any orthonormal pair)
        return (0.0, 1.0) if pol == "s" else (1.0, 0.0)
    s = np.array([ky, -kx]) / kp  # k × N with N = +z
    if pol == "s":
        return complex(s[0]), complex(s[1])
    k = np.sqrt(kx * kx + ky * ky + kz * kz)
    # p = k̂ × s; its transverse part is (k_z/|k|)·(k_x, k_y)/k⊥ up to orientation
    p = np.cross(np.array([kx, ky, kz]) / k, np.array([s[0], s[1], 0.0]))
    return complex(p[0]), complex(p[1])


def order_efficiencies(
    handle: fga.Handle,
    wavelength_um: float,
    kx: float,
    ky: float,
    polarisation: str = "unpolarised",
) -> dict[str, dict[int, float]]:
    """Power fractions per order, {"transmitted": {n: η}, "reflected": {n: η}} for 1D gratings.

    Args:
        handle: prepared FourierGrating structure.
        wavelength_um: vacuum wavelength.
        kx, ky: transverse wave vector of the incident wave in the superstrate (1/mm).
        polarisation: "s", "p" or "unpolarised" (mean of s and p).
    """
    if polarisation == "unpolarised":
        s = order_efficiencies(handle, wavelength_um, kx, ky, "s")
        p = order_efficiencies(handle, wavelength_um, kx, ky, "p")
        return {
            side: {
                n: 0.5 * (s[side].get(n, 0.0) + p[side].get(n, 0.0))
                for n in set(s[side]) | set(p[side])
            }
            for side in s
        }
    port = fga.front_port(handle, wavelength_um)
    kz = float(np.real(port.kz(np.array(kx), np.array(ky))))
    e = _unit_polarisation(kx, ky, kz, polarisation)
    incident = PlaneWaves(port, np.array([kx]), np.array([ky]), np.array([e], dtype=np.complex128))
    p_in = float(incident.power_density()[0])
    out = fga.apply(handle, incident)
    result: dict[str, dict[int, float]] = {}
    for side, name in (("back", "transmitted"), ("front", "reflected")):
        waves = out[side]
        assert waves.orders is not None
        density = waves.power_density()
        result[name] = {int(o[0]): float(density[j]) / p_in for j, o in enumerate(waves.orders)}
    return result


def with_efficiencies(system: Any, surface_id: str, table: dict[int, float]) -> Any:
    """Raytatouille system with `diffraction_efficiency` of the surface set from the table."""
    import raytatouille as rt

    pointer = system.locate_surface(surface_id)
    if pointer is None:
        raise ValueError(f"no surface {surface_id!r}")
    value = [
        {"order": int(m), "efficiency": float(min(max(eta, 0.0), 1.0))}
        for m, eta in sorted(table.items())
    ]
    return rt.apply_patch(
        system, [{"op": "add", "path": pointer + "/diffraction_efficiency", "value": value}]
    )
