"""Adapter for hp-FEM gratings (conical Nédélec/H1 solver, `hpfem.grating`), stage R9.

Interface: docs/schnittstelle.md; representation S (diffraction orders), like FourierGrating.

Conventions of hp-FEM (hp-FEM docs/python.md, `hpfem.grating.jones`): time factor e^(−iωt) (as
ours), SI metres; solver frame x = period, y = stack normal (light from +y), z = along the lines,
so the RTT frame is x′ = x, y′ = z, z′ = −y (as for FourierGrating); incidence (θ, φ) with
k = k(sin θ cos φ, −cos θ, sin θ sin φ). `grating.jones` gives power-normalised Jones matrices in
the basis ŝ = k̂ × ŷ normalised (ẑ at normal incidence), p̂ = ŝ × k̂ – both differ from
Convention A in sign –, referred to the left cell edge x_min and to the top (reflected,
incident) and bottom (transmitted) interface of the background stack. The adapter rebuilds
Cartesian field matrices with hp-FEM's own basis (`hpfem.conical_polarisation`), removes the
power normalisation and moves the references to x = 0 and to the port planes y_front, y_back by
the plane-wave phases.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from rtt_fusion import __version__
from rtt_fusion.ports import PlaneWaves, Port

ADAPTER_VERSION = "0.1.0"
MM_PER_M = 1e3


def _hp() -> Any:
    import hpfem

    return hpfem


@dataclass(frozen=True)
class Structure:
    """Grating unit cell for `hpfem.grating.jones`.

    Attributes:
        mesh, materials, stack: as `hpfem.grating.solve` (SI, solver frame); the x origin of the
            mesh is x′ = 0 of the ports.
        n_cover, n_substrate: real indices of the incidence and exit medium.
        y_front, y_back: port planes in the solver frame (m), above and below the structure.
        order: polynomial order; options: further keyword arguments of `grating.solve`.
    """

    mesh: Any
    materials: dict[int, Any]
    stack: Any
    n_cover: float
    n_substrate: float
    y_front: float
    y_back: float
    order: int = 5
    options: tuple[tuple[str, Any], ...] = ()


@dataclass
class Handle:
    structure: Structure
    wavelengths_um: tuple[float, ...]
    content_hash: str
    _cache: dict[tuple[float, float, float], Any] = field(default_factory=dict, repr=False)


def describe() -> dict[str, Any]:
    return {
        "name": "hp-FEM",
        "version": str(getattr(_hp(), "__version__", "0.4.0")),
        "adapter_version": ADAPTER_VERSION,
        "method": "hp finite elements (Nédélec + H1), conical 2.5D grating cell, PML, Bloch",
        "representation": "S",
        "ports": [
            {"name": "front", "plane": "y_front -> z' = 0", "in": "+z", "out": "-z"},
            {"name": "back", "plane": "y_back -> z' = y_front - y_back", "out": "+z"},
        ],
        "linear": True,
        "reciprocal": True,
        "operator": False,
        "deterministic": {"cpu": "direct solver: bitwise"},
        "conventions_of_solver": [
            "SI metres; frame x' = x, y' = z, z' = -y",
            "Jones in s = k x y (z at normal incidence), p = s x k: rebuilt as Cartesian fields",
            "power normalisation removed; references moved from x_min and the stack interfaces",
        ],
        "validity": [
            "convergence in mesh size and order is the user's responsibility (reference K7)",
            "isotropic cover and substrate",
        ],
    }


def prepare(structure: Structure, wavelengths_um: list[float] | tuple[float, ...]) -> Handle:
    if structure.y_front < structure.y_back:
        raise ValueError("y_front must lie above y_back")
    if not wavelengths_um:
        raise ValueError("at least one wavelength is required")
    h = hashlib.sha256()
    h.update(f"rtt_fusion {__version__} hpfem-adapter {ADAPTER_VERSION}".encode())
    h.update(np.asarray(structure.mesh.vertices, dtype=np.float64).tobytes())
    h.update(
        repr(
            (
                structure.n_cover,
                structure.n_substrate,
                structure.y_front,
                structure.y_back,
                structure.order,
                structure.options,
                sorted(structure.materials),
            )
        ).encode()
    )
    h.update(np.asarray(wavelengths_um, dtype=np.float64).tobytes())
    return Handle(structure, tuple(float(w) for w in wavelengths_um), h.hexdigest())


def front_port(handle: Handle, wavelength_um: float, direction: int = 1) -> Port:
    return Port(0.0, direction, complex(handle.structure.n_cover), wavelength_um)


def back_port(handle: Handle, wavelength_um: float) -> Port:
    s = handle.structure
    return Port((s.y_front - s.y_back) * MM_PER_M, 1, complex(s.n_substrate), wavelength_um)


def _basis(k: np.ndarray, pol: Any) -> np.ndarray:  # type: ignore[type-arg]
    hp = _hp()
    return np.asarray(hp.conical_polarisation(k, np.array([0.0, 1.0, 0.0]), pol), dtype=complex)


def _field_matrices(handle: Handle, wavelength_um: float, theta: float, phi: float) -> Any:
    """Per side and order m: (k_t, β, k_n, 3×3 Cartesian matrix in the solver frame) such that
    E_out(x = 0, port plane) = M · E_in(x = 0, front plane)."""
    key = (wavelength_um, theta, phi)
    if key in handle._cache:
        return handle._cache[key]
    hp = _hp()
    s = handle.structure
    omega = hp.units.angular_frequency(wavelength=wavelength_um * 1e-6)
    res = hp.grating.jones(
        s.mesh, s.materials, s.stack, theta, phi, omega, s.order, **dict(s.options)
    )
    wave = res.s.wave
    x_min = res.s.x_min
    geometry = res.s.problem.setup.background
    top, bottom = geometry.top, geometry.bottom
    kx, ky, beta = float(wave.kx), float(np.real(wave.ky)), float(wave.beta)
    k_in = np.array([kx, -ky, beta])
    b_in = np.column_stack([_basis(k_in, hp.Polarisation.S), _basis(k_in, hp.Polarisation.P)])
    out: dict[str, dict[int, Any]] = {"R": {}, "T": {}}
    for side, orders in (("R", res.R_orders), ("T", res.T_orders)):
        gratings = res.s.R_orders if side == "R" else res.s.T_orders
        kts = {o.m: (float(o.kt), complex(o.kn)) for o in gratings}
        for jo in orders:
            kt, kn = kts[jo.m]
            kn_r = float(np.real(kn))
            sign = 1.0 if side == "R" else -1.0
            k_out = np.array([kt, sign * kn_r, beta])
            b_out = np.column_stack(
                [_basis(k_out, hp.Polarisation.S), _basis(k_out, hp.Polarisation.P)]
            )
            field_jones = np.asarray(jo.jones) * np.sqrt(ky / kn_r)  # undo power normalisation
            m3 = b_out @ field_jones @ b_in.conj().T
            # references: x_min → 0; incident top → y_front; outgoing top/bottom → port plane
            phase = np.exp(1j * (kx - kt) * x_min)
            phase *= np.exp(
                -1j * ky * (top - s.y_front)
            )  # E_in(top) = E_in(y_f) e^(−i k_y (top − y_f))
            if side == "R":
                phase *= np.exp(1j * kn_r * (s.y_front - top))
            else:
                phase *= np.exp(-1j * kn_r * (s.y_back - bottom))
            out[side][jo.m] = (kt, beta, kn_r, m3 * phase)
    handle._cache[key] = out
    return out


def apply(handle: Handle, incident: PlaneWaves) -> dict[str, PlaneWaves]:
    """Plane waves on the front port → reflected (front, −z) and transmitted (back, +z) orders."""
    p = incident.port
    s = handle.structure
    if p.direction != 1 or p.z != 0.0 or complex(p.index) != complex(s.n_cover):
        raise ValueError("incident waves must be on the front port with direction +1")
    k_mm = s.n_cover * p.k0
    ez = incident.ez()
    sides: dict[str, dict[tuple[float, float], Any]] = {"front": {}, "back": {}}
    for j in range(len(incident.kx)):
        kx_p, ky_p = float(incident.kx[j]), float(incident.ky[j])
        sin_t = np.hypot(kx_p, ky_p) / k_mm
        theta = float(np.arcsin(min(sin_t, 1.0)))
        phi = float(np.arctan2(ky_p, kx_p)) if sin_t > 0 else 0.0
        e_rtt = np.array([incident.field[j, 0], incident.field[j, 1], ez[j]])
        e_hp = np.array([e_rtt[0], -e_rtt[2], e_rtt[1]])  # (x, y, z) = (x′, −z′, y′)
        mats = _field_matrices(handle, p.wavelength_um, theta, phi)
        for side, key in (("front", "R"), ("back", "T")):
            for m, (kt, beta, _, m3) in mats[key].items():
                e = m3 @ e_hp
                e_out = np.array([e[0], e[2], -e[1]])  # → (x′, y′, z′)
                kxy = (round(kt / MM_PER_M, 6), round(beta / MM_PER_M, 6))
                entry = sides[side].setdefault(kxy, [np.zeros(3, dtype=complex), m])
                entry[0] = entry[0] + e_out
    result = {}
    for side, port in (
        ("front", front_port(handle, p.wavelength_um, -1)),
        ("back", back_port(handle, p.wavelength_um)),
    ):
        keys = sorted(sides[side])
        result[side] = PlaneWaves(
            port,
            np.array([k[0] for k in keys]),
            np.array([k[1] for k in keys]),
            np.array([sides[side][k][0][:2] for k in keys], dtype=np.complex128).reshape(-1, 2),
            np.array([[sides[side][k][1], 0] for k in keys], dtype=np.int64).reshape(-1, 2),
        )
    return result


def estimate(handle: Handle) -> dict[str, float]:
    n = int(getattr(handle.structure.mesh, "num_cells", 1000))
    return {
        "seconds": 2e-3 * n * handle.structure.order**2,
        "memory_bytes": 4e3 * n * handle.structure.order**3,
    }
