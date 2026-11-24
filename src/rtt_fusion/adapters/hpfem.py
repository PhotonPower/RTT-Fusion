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
import importlib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from rtt_fusion import __version__
from rtt_fusion.ports import PlaneWaves, Port

ADAPTER_VERSION = "0.2.0"
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

    # what the Jones table (`couplings.jones_table`) asks of a prepared grating, as of
    # FourierGrating's handle
    @property
    def front_index(self) -> complex:
        """Index of the incidence medium (front port)."""
        return complex(self.structure.n_cover)

    @property
    def back_index(self) -> complex:
        """Index of the exit medium (back port)."""
        return complex(self.structure.n_substrate)

    @property
    def period_um(self) -> float:
        """Grating period Λ (µm): the width of the cell's mesh."""
        x = np.asarray(self.structure.mesh.vertices, dtype=np.float64)[:, 0]
        return float(np.ptp(x)) * 1e6


#: solver options of a job (`solver` of hp-FEM's job schema v1) that `grating.jones` takes
_OPTIONS = ("bottom", "orders_max", "snap_tolerance", "pml_target", "check")


def _real_index(material: Any, where: str, lossless: bool) -> float:
    eps = complex(np.asarray(material.eps_tensor)[0, 0])
    if lossless and abs(eps.imag) > 1e-12:
        raise ValueError(f"{where} must be lossless (eps = {eps:.6g})")
    return float(np.sqrt(eps).real)


def structure_from_job(
    job: Mapping[str, Any], wavelength_um: float, base: str | Path | None = None
) -> Structure:
    """The grating cell of an hp-FEM job document (schema v1, problem "grating", as hp-FEM's
    `python -m hpfem.run` and the hp-FEM plugin of Raytatouille Fusion write it) at one
    wavelength, built with hp-FEM's own job reader (`hpfem.run`, so a job computes here as it
    does there): unit cell and mesh (structured, gmsh or a mesh file relative to `base`), the
    materials by tag, the background stack and the solver options (order, PML, bottom, orders).
    The job's incidence, sweep and maps are not used – the incident waves set them.

    Ports (Convention of the R9 adapter): the front plane on the highest, the back plane on the
    lowest point of the structure and the stack – the scatterers are the cells whose material
    differs from the stack's at their centroid, as hp-FEM finds them –, so that a block's surface
    point is the top of the grating. The stack's media are taken at `wavelength_um` (a Jones
    table of a few nm around it neglects their dispersion there; the cells' materials follow the
    wavelength of every solve). A lossy substrate has no transmitted orders.

    Raises ValueError for a job that is not a grating, a lossy incidence medium or an error of
    hp-FEM's reader.
    """
    hp = _hp()
    run = importlib.import_module("hpfem.run")
    if str(job.get("problem", "grating")) != "grating":
        raise ValueError(f"problem {job.get('problem')!r}: the grating block needs 'grating'")
    unit = float(job.get("unit", 1.0))
    omega = float(hp.units.angular_frequency(wavelength=wavelength_um * 1e-6))
    try:
        cell = run._unit_cell(job["model"], unit)
        materials = {
            int(t): run._material(m, f"materials[{t}]") for t, m in dict(job["materials"]).items()
        }
        stack = run._stack(job["stack"], omega, unit)
        mesh = run._build_mesh(job, cell, materials, omega, unit, Path(base or Path.cwd()))
    except KeyError as e:
        raise ValueError(f"hp-FEM job: missing {e}") from e
    n_cover = _real_index(stack.incidence_medium, "the incidence medium", True)
    n_substrate = _real_index(stack.substrate, "the substrate", False)
    lo, hi = _structure_extent(mesh, materials, stack, omega)
    solver = dict(job.get("solver", {}))
    options: dict[str, Any] = {k: solver[k] for k in _OPTIONS if k in solver}
    pml = solver.get("pml")
    if isinstance(pml, Mapping):
        options["pml"] = {str(k): float(v) * unit for k, v in pml.items()}
    return Structure(
        mesh,
        materials,
        stack,
        n_cover,
        n_substrate,
        max(hi, float(stack.top)),
        min(lo, float(stack.bottom)),
        int(solver.get("order", 4)),
        tuple(sorted(options.items())),
    )


def _structure_extent(mesh: Any, materials: dict[int, Any], stack: Any, omega: float) -> Any:
    """Lowest and highest vertex y (m) of the cells whose material differs from the stack's at
    their centroid (hp-FEM's rule, `grating._structure_extent`); (+inf, −inf) for none."""
    at = {t: (m.at(omega) if hasattr(m, "at") else m) for t, m in materials.items()}
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    centroids = np.asarray(mesh.cell_centroids, dtype=np.float64)
    tags = np.asarray(mesh.cell_tags)
    lo, hi = np.inf, -np.inf
    for c in range(mesh.num_cells):
        material = at[int(tags[c])]
        background = stack.material_at([float(centroids[c, 0]), float(centroids[c, 1])])
        if material.mu_r != background.mu_r or not np.array_equal(
            np.asarray(material.eps_tensor), np.asarray(background.eps_tensor)
        ):
            ys = vertices[np.asarray(mesh.cell_vertices(c), dtype=int), 1]
            lo, hi = min(lo, float(ys.min())), max(hi, float(ys.max()))
    return lo, hi


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


def apply_many(handle: Handle, incidents: list[PlaneWaves]) -> list[dict[str, PlaneWaves]]:
    """`apply` for many incidences (L4): the Jones matrices of each (λ, θ, φ) are computed once
    and kept in the handle; a conical operator in hp-FEM (plan L4, stage O3) will share the
    factorisation between s and p."""
    return [apply(handle, inc) for inc in incidents]


def estimate(handle: Handle) -> dict[str, float]:
    n = int(getattr(handle.structure.mesh, "num_cells", 1000))
    return {
        "seconds": 2e-3 * n * handle.structure.order**2,
        "memory_bytes": 4e3 * n * handle.structure.order**3,
    }
