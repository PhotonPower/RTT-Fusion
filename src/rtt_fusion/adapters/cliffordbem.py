"""Adapter for CliffordBem (boundary elements for particles), stage R4.

Interface: docs/schnittstelle.md. Representation S (plane waves) at two port planes in front of
(z_front) and behind (z_back) the particle; the particle is a closed body with its centre at
`center_mm` in the RTT frame.

Conversions (only here):

- units: CliffordBem lengths in µm (L = 1 µm), ω = k₀ = 2π/λ[µm]; ports in mm;
- incident plane waves (phase 1 at the front-port origin (0, 0, z_front)) → one
  `AngularSpectrumField` with focus at that origin in body coordinates (CliffordBem docs/python.md:
  E = Σ a_j exp(i k k̂_j·(x − x_f)), time factor e^(−iωt) as in Raytatouille);
- scattered field → plane waves: the far field E∞ (E_s ≈ E∞ e^(ikr)/r, CliffordBem `far_field_E`,
  consistent with σ_ext = 4π/k·Im(p*·E∞(d))) gives the exact angular spectrum of the radiated
  field in a half space beyond the body (Weyl): E_s(r) = (ik/2π) ∫ E∞(k̂) e^(ik k̂·(r − c)) dΩ over
  the hemisphere of propagation; discretised with Gauss-Legendre in cos θ and a uniform φ grid.
  Each output wave carries amplitude (ik/2π)·w_j·E∞(k̂_j)·e^(ik k̂_j·(P₀ − c)) with P₀ the port
  origin; evanescent waves are left out (ports ≥ about a wavelength from the body).

Outputs: "front" = backscattered waves (direction −1); "back" = incident waves propagated to the
back port plus forward-scattered waves (direction +1).
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from rtt_fusion import __version__
from rtt_fusion._types import ComplexArray
from rtt_fusion.convert import far_field_to_spectrum, hemisphere
from rtt_fusion.ports import PlaneWaves, Port

ADAPTER_VERSION = "0.1.0"
UM_PER_MM = 1000.0


def _cb() -> Any:
    import cliffordbem

    return cliffordbem


@dataclass(frozen=True)
class Structure:
    """Particle and ports.

    Attributes:
        mesh: CliffordBem `TriangleMesh` in µm, centred on the particle centre.
        eps_particle: relative permittivity of the particle.
        n_background: real refractive index around the particle (and at the ports).
        center_mm: particle centre in the RTT frame (mm).
        z_front_mm, z_back_mm: port planes (mm), in front of and behind the particle.
        tol: GMRES tolerance of the BEM solve.
    """

    mesh: Any
    eps_particle: complex
    n_background: float
    center_mm: tuple[float, float, float]
    z_front_mm: float
    z_back_mm: float
    tol: float = 1e-8


@dataclass
class Handle:
    structure: Structure
    wavelengths_um: tuple[float, ...]
    content_hash: str
    _problems: dict[float, Any] = field(default_factory=dict, repr=False)


def describe() -> dict[str, Any]:
    cb = _cb()
    return {
        "name": "CliffordBem",
        "version": str(cb.__version__),
        "adapter_version": ADAPTER_VERSION,
        "method": "boundary elements in Cl3(C), H-matrices, GMRES",
        "representation": "S",
        "ports": [
            {"name": "front", "plane": "z_front", "in": "+z", "out": "-z (backscattered)"},
            {"name": "back", "plane": "z_back", "out": "+z (incident + forward scattered)"},
        ],
        "linear": True,
        "reciprocal": True,
        "operator": False,
        "deterministic": {"cpu": "OpenMP reductions: in tolerance"},
        "conventions_of_solver": [
            "lengths in um (L = 1 um), omega = k0 = 2 pi / lambda[um]",
            "incident field: AngularSpectrumField, focus = front-port origin",
            "far field E_s ~ E_inf exp(ikr)/r -> plane waves via the Weyl expansion",
        ],
        "validity": [
            "ports outside the smallest slab containing the body; evanescent waves left out",
            "flat triangles: about 1-3 % in cross sections at 720-2000 triangles (Mie)",
        ],
    }


def prepare(structure: Structure, wavelengths_um: list[float] | tuple[float, ...]) -> Handle:
    s = structure
    points = np.asarray(s.mesh.points, dtype=np.float64)
    zmin = s.center_mm[2] + points[:, 2].min() / UM_PER_MM
    zmax = s.center_mm[2] + points[:, 2].max() / UM_PER_MM
    if not (s.z_front_mm < zmin and s.z_back_mm > zmax):
        raise ValueError("ports must lie in front of and behind the body")
    if not wavelengths_um:
        raise ValueError("at least one wavelength is required")
    h = hashlib.sha256()
    h.update(
        f"rtt_fusion {__version__} cbem-adapter {ADAPTER_VERSION} {_cb().__version__}".encode()
    )
    h.update(points.tobytes())
    h.update(np.asarray(s.mesh.triangles, dtype=np.int64).tobytes())
    h.update(
        repr(
            (complex(s.eps_particle), s.n_background, s.center_mm, s.z_front_mm, s.z_back_mm, s.tol)
        ).encode()
    )
    h.update(np.asarray(wavelengths_um, dtype=np.float64).tobytes())
    return Handle(s, tuple(float(w) for w in wavelengths_um), h.hexdigest())


def front_port(handle: Handle, wavelength_um: float, direction: int = 1) -> Port:
    s = handle.structure
    return Port(s.z_front_mm, direction, complex(s.n_background), wavelength_um)


def back_port(handle: Handle, wavelength_um: float) -> Port:
    s = handle.structure
    return Port(s.z_back_mm, 1, complex(s.n_background), wavelength_um)


def _problem(handle: Handle, wavelength_um: float) -> Any:
    if wavelength_um not in handle.wavelengths_um:
        raise ValueError("wavelength was not prepared")
    if wavelength_um not in handle._problems:
        cb = _cb()
        s = handle.structure
        omega = 2 * np.pi / wavelength_um
        handle._problems[wavelength_um] = cb.ScatteringProblem(
            s.mesh,
            cb.Medium(eps=complex(s.eps_particle)),
            omega,
            outer=cb.Medium(eps=s.n_background**2),
        )
    return handle._problems[wavelength_um]


@dataclass(frozen=True)
class Scattering:
    """BEM solution of one incidence (body coordinates in µm).

    Attributes:
        field: the CliffordBem incident field (AngularSpectrumField).
        rhs: incident trace b; h: total trace; h_scat = h − b.
        k_um: wavenumber in the background (1/µm).
    """

    field: Any
    rhs: ComplexArray
    h: ComplexArray
    k_um: float

    @property
    def h_scat(self) -> ComplexArray:
        return np.asarray(self.h - self.rhs)


def scatter(handle: Handle, incident: PlaneWaves) -> Scattering:
    """Solve the BEM for plane waves incident on the front port."""
    cb = _cb()
    s = handle.structure
    p = incident.port
    if p.direction != 1 or p.z != s.z_front_mm or complex(p.index) != complex(s.n_background):
        raise ValueError("incident waves must be on the front port with direction +1")
    problem = _problem(handle, p.wavelength_um)
    k_um = s.n_background * 2 * np.pi / p.wavelength_um
    kz = p.kz(incident.kx, incident.ky)
    if np.any(np.abs(np.imag(kz)) > 0):
        raise ValueError("incident waves must propagate")
    k_mm = k_um * UM_PER_MM
    dirs = np.column_stack([incident.kx, incident.ky, np.real(kz)]) / k_mm
    amps = np.column_stack([incident.field[:, 0], incident.field[:, 1], incident.ez()])
    origin = (np.array([0.0, 0.0, s.z_front_mm]) - np.asarray(s.center_mm)) * UM_PER_MM
    medium = cb.Medium(eps=s.n_background**2)
    fld = cb.AngularSpectrumField(medium, 2 * np.pi / p.wavelength_um, dirs, amps, focus=origin)
    b = np.asarray(fld.project(s.mesh, medium))
    r = problem.solve_rhs(b, cb.SolveOptions(tol=s.tol))
    return Scattering(fld, b, np.asarray(r.h), k_um)


def scattered_plane_waves(
    handle: Handle, sol: Scattering, port: Port, n_theta: int = 48, n_phi: int = 96
) -> PlaneWaves:
    """Scattered field at a port as plane waves (`rtt_fusion.convert.far_field_to_spectrum`)."""
    cb = _cb()
    s = handle.structure
    dirs, w = hemisphere(n_theta, n_phi, port.direction == 1)
    e_um = np.asarray(
        cb.far_field_E(s.mesh, sol.h_scat, sol.k_um, s.n_background**2, dirs), dtype=np.complex128
    )
    # E∞ with r in µm → r in mm: E∞/r is invariant, so E∞[mm] = E∞[µm]·1e-3
    return far_field_to_spectrum(port, np.asarray(s.center_mm), dirs, w, e_um / UM_PER_MM)


def apply(
    handle: Handle, incident: PlaneWaves, n_theta: int = 48, n_phi: int = 96
) -> dict[str, PlaneWaves]:
    """Plane waves on the front port → backscattered (front) and transmitted + scattered (back)."""
    sol = scatter(handle, incident)
    p = incident.port
    front = scattered_plane_waves(
        handle, sol, front_port(handle, p.wavelength_um, -1), n_theta, n_phi
    )
    back_s = scattered_plane_waves(handle, sol, back_port(handle, p.wavelength_um), n_theta, n_phi)
    dz = handle.structure.z_back_mm - handle.structure.z_front_mm
    inc_back = (
        np.asarray(incident.field) * np.exp(1j * p.kz(incident.kx, incident.ky) * dz)[:, None]
    )
    back = PlaneWaves(
        back_s.port,
        np.concatenate([incident.kx, back_s.kx]),
        np.concatenate([incident.ky, back_s.ky]),
        np.concatenate([inc_back, back_s.field]).astype(np.complex128),
    )
    # "back" mixes the incident waves (discrete or with their own weights) and the scattered
    # spectrum, so it carries no common solid-angle weights; "back_scattered" carries them.
    return {"front": front, "back": back, "back_scattered": back_s}


def estimate(handle: Handle) -> dict[str, float]:
    n = len(handle.structure.mesh.triangles)
    return {"seconds": 4e-3 * n * np.log(max(n, 2)), "memory_bytes": 2e4 * n * np.log(max(n, 2))}
