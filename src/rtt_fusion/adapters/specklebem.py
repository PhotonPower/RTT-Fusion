"""Adapter for SpeckleBem (surface integral equations, MLFMM), stage R8: closed bodies.

Interface: docs/schnittstelle.md; representation S as for CliffordBem (ports in front of and
behind the body, far field → port spectra with `rtt_fusion.convert.far_field_to_spectrum`).

SpeckleBem uses the conventions of RF engineering (SpeckleBem docs/06_conventions.md), the
largest difference among the solvers of the Fusion:

- time factor exp(+jωt): a physical field Re{E₊ e^(jωt)} is Re{E₋ e^(−iωt)} with E₋ = conj(E₊).
  So every complex quantity crosses the boundary **complex conjugated**: the permittivity
  (passive media have Im ε ≤ 0 there, Im ε ≥ 0 here), the incident amplitude (a plane wave
  e0 exp(−j k k̂·r) carries e0 = conj(a) for our a exp(i k k̂·r)) and the far field
  (E_s ≈ F e^(−jkr)/r → E∞ = conj(F));
- SI units: lengths in m, so F[mm-convention] = F[m]·1e3 (E∞/r is invariant);
- phase reference of the plane wave at the origin of the mesh (body centre).

Each incident plane wave is one solve; further waves reuse the assembled operator
(`Simulation.set_excitation`). Rough surfaces (open meshes with PML) are not covered: their
illumination must be a Gaussian beam in SpeckleBem (no arbitrary spectra yet, Fusion plan D12).

Light may also come in at the back port (direction −1, L5): then "back" holds the backscattered
waves (direction +1) and "front" the transmitted and forward-scattered ones (direction −1), with
"front_scattered" alone. `apply_many` (L4) applies many incidences with the operator of each
wavelength assembled once.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from rtt_fusion import __version__
from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.convert import far_field_to_spectrum, hemisphere
from rtt_fusion.ports import PlaneWaves, Port

ADAPTER_VERSION = "0.1.0"
M_PER_MM = 1e-3


def _sb() -> Any:
    import specklebem

    return specklebem


@dataclass(frozen=True)
class Structure:
    """Closed body and ports.

    Attributes:
        mesh: SpeckleBem `TriangleMesh` in metres, centred on the body centre.
        eps_particle: relative permittivity in **our** convention (Im ε ≥ 0 for absorption).
        n_background: real refractive index around the body.
        center_mm, z_front_mm, z_back_mm: body centre and port planes in the RTT frame (mm).
        tol: GMRES tolerance; compression: "dense" or "mlfmm".
    """

    mesh: Any
    eps_particle: complex
    n_background: float
    center_mm: tuple[float, float, float]
    z_front_mm: float
    z_back_mm: float
    tol: float = 1e-8
    compression: str = "dense"


@dataclass
class Handle:
    structure: Structure
    wavelengths_um: tuple[float, ...]
    content_hash: str
    _sims: dict[float, Any] = field(default_factory=dict, repr=False)


def describe() -> dict[str, Any]:
    return {
        "name": "SpeckleBem",
        "version": "0.1.0",
        "adapter_version": ADAPTER_VERSION,
        "method": "surface integral equations (PMCHWT/ICTF/MCTF), RWG, dense or MLFMM",
        "representation": "S",
        "ports": [
            {"name": "front", "plane": "z_front", "in": "+z", "out": "-z (backscattered)"},
            {"name": "back", "plane": "z_back", "out": "+z (incident + forward scattered)"},
        ],
        "linear": True,
        "reciprocal": True,
        "operator": False,
        "deterministic": {"cpu": "OpenMP: in tolerance"},
        "conventions_of_solver": [
            "time factor exp(+j omega t): permittivity, incident amplitude and far field "
            "complex conjugated",
            "SI metres: far field F[mm] = F[m] * 1e3",
            "plane-wave phase at the mesh origin (body centre)",
        ],
        "validity": [
            "closed meshes only (box route); rough surfaces need a Gaussian beam (D12)",
            "flat RWG elements: about 1-3 % in far fields at lambda/12 (Mie, absorbing sphere)",
        ],
    }


def prepare(structure: Structure, wavelengths_um: list[float] | tuple[float, ...]) -> Handle:
    s = structure
    vertices = np.asarray(s.mesh.vertices, dtype=np.float64)
    zmin = s.center_mm[2] + vertices[:, 2].min() / M_PER_MM
    zmax = s.center_mm[2] + vertices[:, 2].max() / M_PER_MM
    if not (s.z_front_mm < zmin and s.z_back_mm > zmax):
        raise ValueError("ports must lie in front of and behind the body")
    if complex(s.eps_particle).imag < 0:
        raise ValueError("eps_particle in the convention exp(-i omega t): Im eps >= 0 if passive")
    if not wavelengths_um:
        raise ValueError("at least one wavelength is required")
    h = hashlib.sha256()
    h.update(f"rtt_fusion {__version__} speckle-adapter {ADAPTER_VERSION}".encode())
    h.update(vertices.tobytes())
    h.update(np.asarray(s.mesh.triangles, dtype=np.int64).tobytes())
    h.update(
        repr(
            (
                complex(s.eps_particle),
                s.n_background,
                s.center_mm,
                s.z_front_mm,
                s.z_back_mm,
                s.tol,
                s.compression,
            )
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


def _plane_wave(
    handle: Handle, wavelength_um: float, direction: FloatArray, e_ours: ComplexArray
) -> Any:
    sb = _sb()
    background = sb.Material(eps_r=complex(handle.structure.n_background**2))
    return sb.PlaneWave(
        wavelength_um * 1e-6,
        np.asarray(direction, dtype=np.float64),
        np.conj(np.asarray(e_ours, dtype=np.complex128)),
        background,
    )


@dataclass
class Scattering:
    """Incident waves of one call; far fields by superposition of one solve per wave."""

    handle: Handle
    sim: Any
    excitations: list[Any]
    solved: int = -1  # index of the excitation whose currents the Simulation holds

    def far_field(self, directions: FloatArray) -> ComplexArray:
        """E∞ in our convention (E_s ≈ E∞ e^(ikr)/r, r in mm) for unit directions (N, 3)."""
        d = np.asarray(directions, dtype=np.float64)
        total = np.zeros((len(d), 3), dtype=np.complex128)
        for j, excitation in enumerate(self.excitations):
            if self.solved != j:
                self.sim.set_excitation(excitation)
                self.sim.solve(tol=self.handle.structure.tol)
                self.solved = j
            total += np.conj(np.asarray(self.sim.far_field(d))) / M_PER_MM
        return total


def _incident_side(handle: Handle, p: Port) -> str:
    """ "front" (front port, direction +1) or "back" (back port, direction −1)."""
    s = handle.structure
    if complex(p.index) == complex(s.n_background):
        if p.direction == 1 and p.z == s.z_front_mm:
            return "front"
        if p.direction == -1 and p.z == s.z_back_mm:
            return "back"
    raise ValueError(
        "incident waves must be on the front port with direction +1 or on the back port with "
        "direction -1"
    )


def scatter(handle: Handle, incident: PlaneWaves) -> Scattering:
    """Solve for plane waves incident on the front port (+z) or the back port (−z), one solve
    per wave."""
    sb = _sb()
    s = handle.structure
    p = incident.port
    _incident_side(handle, p)
    if p.wavelength_um not in handle.wavelengths_um:
        raise ValueError("wavelength was not prepared")
    kz = p.kz(incident.kx, incident.ky)
    if np.any(np.abs(np.imag(kz)) > 0):
        raise ValueError("incident waves must propagate")
    k = s.n_background * p.k0
    dirs = np.column_stack([incident.kx, incident.ky, p.direction * np.real(kz)]) / k
    amps = np.column_stack([incident.field[:, 0], incident.field[:, 1], incident.ez()])
    shift = np.asarray(s.center_mm) - np.array([0.0, 0.0, p.z])
    excitations: list[Any] = []
    sim = handle._sims.get(p.wavelength_um)
    for j in range(len(dirs)):
        a = amps[j] * np.exp(1j * k * (dirs[j] @ shift))  # amplitude at the body centre
        exc = _plane_wave(handle, p.wavelength_um, dirs[j], a)
        if sim is None:
            sim = sb.Simulation(
                s.mesh,
                exc,
                object=sb.Material(eps_r=np.conj(complex(s.eps_particle))),
                compression=s.compression,
            )
            sim.assemble()
            handle._sims[p.wavelength_um] = sim
        excitations.append(exc)
    return Scattering(handle, sim, excitations)


def scattered_plane_waves(
    handle: Handle, sol: Scattering, port: Port, n_theta: int = 48, n_phi: int = 96
) -> PlaneWaves:
    dirs, w = hemisphere(n_theta, n_phi, port.direction == 1)
    return far_field_to_spectrum(
        port, np.asarray(handle.structure.center_mm), dirs, w, sol.far_field(dirs)
    )


def apply(
    handle: Handle, incident: PlaneWaves, n_theta: int = 48, n_phi: int = 96
) -> dict[str, PlaneWaves]:
    """Plane waves on the front port → backscattered (front), incident + forward scattered (back)
    and the forward-scattered spectrum alone (back_scattered, with solid-angle weights); on the
    back port (direction −1) the same mirrored: back, front and front_scattered."""
    side = _incident_side(handle, incident.port)
    sol = scatter(handle, incident)
    p = incident.port
    # both hemispheres in one far-field call: one solve per incident wave
    d_front, w_front = hemisphere(n_theta, n_phi, False)
    d_back, w_back = hemisphere(n_theta, n_phi, True)
    e_inf = sol.far_field(np.vstack([d_front, d_back]))
    c = np.asarray(handle.structure.center_mm)
    front = far_field_to_spectrum(
        front_port(handle, p.wavelength_um, -1), c, d_front, w_front, e_inf[: len(d_front)]
    )
    back_s = far_field_to_spectrum(
        back_port(handle, p.wavelength_um), c, d_back, w_back, e_inf[len(d_front) :]
    )
    reflected, far_s = (front, back_s) if side == "front" else (back_s, front)
    dz = handle.structure.z_back_mm - handle.structure.z_front_mm  # path between the ports
    inc_far = np.asarray(incident.field) * np.exp(1j * p.kz(incident.kx, incident.ky) * dz)[:, None]
    through = PlaneWaves(
        far_s.port,
        np.concatenate([incident.kx, far_s.kx]),
        np.concatenate([incident.ky, far_s.ky]),
        np.concatenate([inc_far, far_s.field]).astype(np.complex128),
    )
    other = "back" if side == "front" else "front"
    return {side: reflected, other: through, f"{other}_scattered": far_s}


def apply_many(
    handle: Handle, incidents: list[PlaneWaves], n_theta: int = 48, n_phi: int = 96
) -> list[dict[str, PlaneWaves]]:
    """`apply` for many incidences (L4): the operator of each wavelength is assembled once and
    kept in the handle (`set_excitation` per wave); results equal `apply`."""
    return [apply(handle, inc, n_theta, n_phi) for inc in incidents]


def estimate(handle: Handle) -> dict[str, float]:
    n = 1.5 * len(handle.structure.mesh.triangles)  # RWG edges
    return {"seconds": 6e-10 * (2 * n) ** 2 * 100, "memory_bytes": 16.0 * (2 * n) ** 2}
