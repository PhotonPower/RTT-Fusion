"""Particle as the sample (stage R19c, concept section 4): CliffordBem in the working graph.

`fusion.particle` has the sockets of a mirror ("in", "out"), so a sample mirror can be replaced
by it with the links kept. The incident rays converge onto the particle (e.g. the focus of an
objective); the block turns them into an angular spectrum (representation S), lets CliffordBem
scatter it (`adapters.cliffordbem`, stage R4) and sends the backscattered spectrum on as rays
leaving the particle – back into the objective.

R -> S (`rays_to_spectrum`): every ray is one plane wave of the focused field, the Debye /
Richards–Wolf picture of a converging beam (fusion plan 6.3, "R -> H -> S" collapsed for a
focus): direction k̂_j, solid angle ΔΩ_j from the triangulated direction cosines (dΩ = dσ/cos θ,
a third of the adjacent triangles per ray), power P_j of the ray from its tube at the source,
    A_j = −i (k/2π) sqrt(P_j/ΔΩ_j) · ê_j · exp(i k₀ OPL_j(F)),   field_j = A_j ΔΩ_j,
with ê_j the ray's polarisation and OPL_j(F) its optical path to the closest approach to the
particle centre F. Then E(r) = Σ field_j exp(ik k̂_j·(r − F)) and `ports.spectrum_power` equals
k·Σ P_j (the flux convention of `ports`).

S -> R (`spectrum_to_rays`): the inverse – a ray per plane wave starting at F along k̂_j with the
power (2π/k)² |field_j|²/ΔΩ_j and the phase of the far field E∞ = −i (2π/k) A (CliffordBem's
Weyl expansion E = (ik/2π) ∫ E∞ e^{ik k̂·(r−c)} dΩ); the source coordinates of the new rays are
their direction cosines, so that downstream ray tubes carry the right power.

Collimated light (all rays parallel, e.g. Köhler illumination of the sample) is one plane wave
with the field of the ray through the particle centre.

Plane-wave spectra (S, stage R30d, docs/plan-s-spektren.md): a `SpectrumBundle` on the input is
the incident angular spectrum as it is – no internal R -> S (`incident_waves`: the spectrum moved
onto CliffordBem's front port in the spectrum, `spectrum.propagate`, an exact phase per wave; the
particle's axis is the spectrum's normal). The backscattered far field leaves as a spectrum too
(`scattered_spectrum`): the waves of CliffordBem's front port referred to the particle centre, in
the mirrored frame (u, −v, −n) as behind a grating, with the sampling of a point source at the
centre (spectral phase ψ̃ = 0, ρ = 0, Hessian 0: the far field E∞ is the residual) moved in the
spectrum onto a plane `OUT_WAVELENGTHS` wavelengths in front of the particle, so that S -> H there
is the inverse PFT of the spherical wave E∞ e^(ikr)/r (exact for it) and S -> R gives rays from
that plane, away from the centre. Rays in, rays out (the Debye way above); spectra in, spectra
out.

Each coherence group is scattered on its own (the BEM problem is linear; groups add as
intensities). The particle: a sphere of radius `radius_um` and index `index` (complex allowed)
in a medium of index `n_medium`, meshed as an icosphere of level `refinement` (CliffordBem
`make_icosphere`, as in pilot K3).
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np
from scipy.spatial import Delaunay, QhullError

from rtt_fusion._types import FloatArray
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.field import FieldError
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle
from rtt_fusion.graph.model import (
    Category,
    Node,
    NodeType,
    ParamSpec,
    PortSpec,
    light_in,
    light_out,
)
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.graph.spectrum import CUT_WARNING, SpectrumBundle, propagate
from rtt_fusion.ports import PlaneWaves, Port

UM_PER_MM = 1000.0
#: the scattered spectrum leaves on a plane this many wavelengths (in the medium) in front of the
#: particle: S -> H there is the inverse PFT of the spherical wave E∞ e^(ikr)/r (exact for it; the
#: near field, 1/(kD) of it, and the evanescent waves are not in the far field anyway), and the
#: blocks behind the particle – an objective, a detector – lie beyond it
OUT_WAVELENGTHS = 100.0


def _vertex_areas(points: FloatArray) -> FloatArray:
    """A third of the areas of the adjacent Delaunay triangles per point."""
    tri = Delaunay(points).simplices
    a, b, c = points[tri[:, 0]], points[tri[:, 1]], points[tri[:, 2]]
    area = 0.5 * np.abs(
        (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    )
    out = np.zeros(len(points))
    for k in range(3):
        np.add.at(out, tri[:, k], area / 3.0)
    return out


def _frame(axis: FloatArray) -> FloatArray:
    """Rows u, v, n with n = axis: global -> local is x @ R.T."""
    p = Plane.make(np.zeros(3), axis)
    return np.asarray(np.stack([p.u, p.v, p.normal]))


def ray_powers(bundle: RayBundle, alive: Any) -> FloatArray:
    """Power of every living ray: irradiance·|amplitude·prt·jones|²·(source area of the ray)."""
    src = np.asarray(bundle.source_xy)[alive]
    areas = _vertex_areas(src)
    field = np.asarray(bundle.amplitude)[alive][:, None] * bundle.field_direction()[alive]
    cos_in = np.broadcast_to(np.asarray(bundle.source_cos, dtype=float), (bundle.size,))[alive]
    return np.asarray(
        bundle.irradiance * np.sum(np.abs(field) ** 2, axis=1) * areas * cos_in, dtype=float
    )


def rays_to_spectrum(
    bundle: RayBundle, centre: FloatArray, axis: FloatArray, port_z_mm: float, n: float = 1.0
) -> tuple[PlaneWaves, FloatArray]:
    """The converging rays as plane waves on a port (local frame: z = axis, origin = centre).

    Returns the plane waves (field_j = A_j ΔΩ_j, with `solid_angle`) and the rotation global ->
    local (rows u, v, n)."""
    rot = _frame(axis)
    alive = np.flatnonzero(np.asarray(bundle.alive))
    if alive.size < 6:
        raise ExecutionError("exec.particle", "too few rays reach the particle")
    d = np.asarray(bundle.direction)[alive] @ rot.T
    if np.any(d[:, 2] <= 0):
        raise ExecutionError("exec.particle", "rays must run towards +axis onto the particle")
    port = Port(port_z_mm, 1, complex(n), bundle.wavelength_um)
    k0 = bundle.k0
    k = n * k0
    if float(np.ptp(d, axis=0).max()) < 1e-9:
        # collimated light (Köhler illumination of the sample): one plane wave with the field of
        # the ray that passes closest to the particle centre
        pos = np.asarray(bundle.position)[alive]
        d_glob = np.asarray(bundle.direction)[alive]
        t = np.sum((centre - pos) * d_glob, axis=1)
        miss = np.linalg.norm(pos + t[:, None] * d_glob - centre, axis=1)
        j = int(np.argmin(miss))
        e0 = bundle.field_direction()[alive][j] * np.asarray(bundle.amplitude)[alive][j]
        e0 = (
            e0
            * np.sqrt(bundle.irradiance)
            * np.exp(1j * k0 * (np.asarray(bundle.opl)[alive][j] + n * t[j]))
        )
        e0_loc = e0 @ rot.T
        shift = np.array([0.0, 0.0, -port_z_mm])
        field = e0_loc[None, :] * np.exp(-1j * k * (d[:1] @ shift))[:, None]
        single = PlaneWaves(port, d[:1, 0] * k, d[:1, 1] * k, field[:, :2].astype(np.complex128))
        return single, rot
    sigma = d[:, :2]
    try:
        cos_t = d[:, 2]
        d_omega = _vertex_areas(sigma) / cos_t
    except QhullError as err:
        raise ExecutionError("exec.particle", f"directions do not span a cone: {err}") from err
    power = ray_powers(bundle, alive)
    pos = np.asarray(bundle.position)[alive]
    d_glob = np.asarray(bundle.direction)[alive]
    t = np.sum((centre - pos) * d_glob, axis=1)  # to the closest approach to F
    opl_f = np.asarray(bundle.opl)[alive] + n * t
    e = bundle.field_direction()[alive]
    e = e / np.maximum(np.linalg.norm(e, axis=1), 1e-300)[:, None]
    e_loc = e @ rot.T
    amp = (
        -1j
        * (k / (2 * np.pi))
        * np.sqrt(power / np.maximum(d_omega, 1e-300))[:, None]
        * e_loc
        * np.exp(1j * k0 * opl_f)[:, None]
    )
    # waves referred to the port origin (0, 0, z_port), the beam focuses at the local origin
    shift = np.array([0.0, 0.0, -port_z_mm])
    field = amp * d_omega[:, None] * np.exp(-1j * k * (d @ shift))[:, None]
    waves = PlaneWaves(
        port, d[:, 0] * k, d[:, 1] * k, field[:, :2].astype(np.complex128), solid_angle=d_omega
    )
    return waves, rot


def spectrum_to_rays(
    waves: PlaneWaves,
    centre: FloatArray,
    rot: FloatArray,
    template: RayBundle,
    path: tuple[str, ...],
) -> RayBundle:
    """Plane waves of a continuous spectrum (with `solid_angle`) as rays leaving the centre."""
    if waves.solid_angle is None:
        raise ExecutionError("exec.particle", "scattered waves without solid angles")
    port = waves.port
    n = float(np.real(port.index))
    k = n * port.k0
    kz = np.real(port.kz(waves.kx, waves.ky)) * port.direction
    d = np.column_stack([waves.kx, waves.ky, kz]) / k
    w = np.asarray(waves.solid_angle)
    field = np.column_stack([waves.field[:, 0], waves.field[:, 1], waves.ez()])
    # back to the particle centre: the waves are referred to the port origin (0, 0, z_port)
    shift = np.array([0.0, 0.0, port.z])
    a = field / w[:, None] * np.exp(-1j * k * (d @ shift))[:, None]  # A_j at the centre
    e_inf = -1j * (2 * np.pi / k) * a  # far field of the scattered wave per direction
    power = (2 * np.pi / k) ** 2 * np.sum(np.abs(a) ** 2, axis=1) * w
    d_glob = d @ rot
    e_glob = e_inf @ rot
    big = np.argmax(np.abs(e_glob), axis=1)
    alpha = np.angle(e_glob[np.arange(len(d)), big])
    e_ref = np.real(e_glob * np.exp(-1j * alpha)[:, None])
    norm = np.linalg.norm(e_ref, axis=1)
    e_ref = e_ref / np.where(norm > 0, norm, 1.0)[:, None]
    sigma = d_glob @ rot[:2].T  # direction cosines in the local u, v axes
    cos_t = np.abs(d[:, 2])
    area_sigma = w * cos_t  # dσ = cos θ dΩ
    scale = np.sqrt(power / np.maximum(area_sigma, 1e-300))
    jones = np.array([1.0, 0.0, 0.0], dtype=complex)
    prt = np.einsum("ni,j->nij", e_ref, jones).astype(complex)
    phase = np.exp(1j * alpha)
    return RayBundle(
        key=template.key,
        wavelength_um=template.wavelength_um,
        position=np.tile(centre, (len(d), 1)),
        direction=d_glob,
        opl=np.zeros(len(d)),
        amplitude=(scale * phase).astype(complex),
        source_xy=sigma,
        alive=power > 0,
        irradiance=1.0,
        path=path,
        jones=jones,
        prt=prt,
        source_cos=1.0,
        source_frame=None,
        index=n,
        source_index=n,
    )


def incident_waves(
    sb: SpectrumBundle, centre: FloatArray, port_z_mm: float, n: float, notes: list[str]
) -> tuple[PlaneWaves, FloatArray]:
    """S at the particle (R30d): the incident spectrum on CliffordBem's front port, local frame
    (u, v of the spectrum's plane, z = its normal, origin = `centre`), port at z = `port_z_mm`.

    The spectrum is moved onto the port plane in the spectrum (`spectrum.propagate`, plan S4:
    e^(i k·Δr) per wave, exact; its sampling is not needed there and left behind); waves that
    are evanescent there (decayed on the way) are left out with `exec.spectrum_cut`. Discrete
    waves stay discrete, a continuous spectrum keeps its dΩ. Returns the plane waves and the
    rotation global -> local (rows u, v, n)."""
    if abs(sb.refractive_index - n) > 1e-9:
        raise ExecutionError(
            "exec.particle",
            f"the spectrum runs in a medium of index {sb.refractive_index:g}, the particle sits "
            f"in n_medium = {n:g}",
        )
    plane = sb.plane
    axis = np.asarray(plane.normal, dtype=float)
    port_plane = Plane(np.asarray(centre) + port_z_mm * axis, axis, plane.u, plane.v)
    moved = propagate(replace(sb, support=None), port_plane, notes)
    w = moved.waves
    kz = w.port.kz(w.kx, w.ky)
    keep = (np.imag(kz) == 0.0) & (np.real(kz) > 0.0)
    if not keep.all():
        e = np.column_stack([np.asarray(w.field), w.ez()])
        weight = np.sum(np.abs(e) ** 2, axis=1)
        if w.solid_angle is not None:
            weight = weight / np.asarray(w.solid_angle)
        share = float(weight[~keep].sum() / (weight.sum() or 1.0))
        if share > CUT_WARNING:
            notes.append(
                f"exec.spectrum_cut: {share:.2e} of the power left out of the spectrum "
                "(evanescent at the particle)"
            )
    if not keep.any():
        raise ExecutionError("exec.particle", "no propagating wave reaches the particle")
    waves = PlaneWaves(
        Port(port_z_mm, 1, complex(n), sb.wavelength_um),
        np.ascontiguousarray(np.asarray(w.kx)[keep]),
        np.ascontiguousarray(np.asarray(w.ky)[keep]),
        np.ascontiguousarray(np.asarray(w.field)[keep], dtype=np.complex128),
        orders=None if w.orders is None else np.asarray(w.orders)[keep],
        solid_angle=None if w.solid_angle is None else np.asarray(w.solid_angle)[keep],
    )
    return waves, np.asarray(np.stack([plane.u, plane.v, axis]))


def out_distance(sb: SpectrumBundle) -> float:
    """Distance (mm) of the scattered spectrum's plane in front of the particle:
    `OUT_WAVELENGTHS` wavelengths in the medium."""
    return OUT_WAVELENGTHS * sb.wavelength_um / UM_PER_MM / sb.refractive_index


def scattered_spectrum(
    waves: PlaneWaves,
    centre: FloatArray,
    rot: FloatArray,
    distance: float,
    template: SpectrumBundle,
    path: tuple[str, ...],
    notes: list[str] | None = None,
) -> SpectrumBundle:
    """The backscattered far field as a spectrum (R30d): `waves` on CliffordBem's front port
    (direction −1, with `solid_angle`, local frame `rot` about `centre`) referred to the centre
    (factor e^(−i k_z z_port)), in the mirrored frame (u, −v, −n) – the light now runs along −n
    – with the sampling of a point source at the centre: ψ̃ = 0, −ρ = 0 and Hessian 0, the
    residual field_j/Δ²k_j (∝ the far field E∞, slowly varying for a small body). Then moved in
    the spectrum (`spectrum.propagate`, S4) onto the plane `distance` in front of the centre:
    ψ̃ = k·Δr, −ρ = the ray shift −(k∥/k_n)·distance, Hessian −distance·∂²k_n/∂k∥², regular for
    distance > 0, so that S -> H there is the inverse PFT of the spherical wave."""
    from fieldtracing.field import GridlessField

    if waves.solid_angle is None:
        raise ExecutionError("exec.particle", "scattered waves without solid angles")
    port = waves.port
    n = float(np.real(port.index))
    k = n * port.k0
    kz = np.real(port.kz(waves.kx, waves.ky))  # |k_z|; the waves run along −n
    field = np.column_stack([np.asarray(waves.field), waves.ez()])
    field = field * np.exp(-1j * port.direction * kz * port.z)[:, None]  # e^(−i k_z z), k_z < 0
    mirror = np.array([1.0, -1.0, -1.0])
    e = field * mirror
    kx = np.ascontiguousarray(np.asarray(waves.kx, dtype=float))
    ky = np.ascontiguousarray(-np.asarray(waves.ky, dtype=float))
    w = np.asarray(waves.solid_angle, dtype=float)
    u, v, axis = (np.asarray(r, dtype=float) for r in rot)
    at_centre = Plane(np.asarray(centre, dtype=float), -axis, u, -v)
    out_port = Port(0.0, 1, complex(n), port.wavelength_um)
    count = len(kx)
    support = GridlessField(
        np.column_stack([kx, ky]),
        np.zeros(count),
        np.zeros((count, 2)),
        np.zeros((count, 2, 2)),
        np.asarray(e / (w * k * kz)[:, None], dtype=np.complex128),
        float(port.wavelength_um),
        n,
        "frequency",
    )
    sb = SpectrumBundle(
        template.key,
        at_centre,
        PlaneWaves(out_port, kx, ky, np.ascontiguousarray(e[:, :2]), solid_angle=w),
        path,
        template.branch,
        support,
    )
    target = Plane(np.asarray(centre) - distance * axis, -axis, u, -v)
    return propagate(sb, target, notes)


def _structure(node: Node, centre: FloatArray, wavelength_um: float) -> Any:
    from rtt_fusion.adapters import cliffordbem as cba

    try:
        import cliffordbem as cb
    except ImportError as e:
        raise ExecutionError("exec.solver", f"CliffordBem is not installed: {e}", node.id) from e
    radius_um = float(node.params.get("radius_um", 0.25))
    index = complex(node.params.get("index", 1.5))
    n_medium = float(node.params.get("n_medium", 1.0))
    m0 = cb.make_icosphere(int(node.params.get("refinement", 4)))
    mesh = cb.TriangleMesh(np.asarray(m0.points) * radius_um, m0.triangles)
    gap = 2.0 * radius_um / UM_PER_MM + 2.0 * wavelength_um / UM_PER_MM
    return cba.Structure(
        mesh,
        index**2,
        n_medium,
        (0.0, 0.0, 0.0),
        -gap,
        gap,
        float(node.params.get("tol", 1e-8)),
    )


def particle(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
    """Backscatter of every incident bundle: rays as rays (the Debye way, module docstring),
    spectra as spectra (R30d)."""
    from rtt_fusion.adapters import cliffordbem as cba

    centre = np.asarray(node.params.get("point", [0.0, 0.0, 0.0]), dtype=float)
    n_theta = int(node.params.get("n_theta", 40))
    n_phi = int(node.params.get("n_phi", 80))
    # the BEM operator lives in the particle's frame (structure at the origin): its key leaves
    # the position out, a moved particle keeps it (operator store, R26)
    local = {
        k: v
        for k, v in node.params.items()
        if k not in ("point", "pose", "pose_offset", "n_theta", "n_phi")
        and not k.startswith("pose_by_")
    }
    store = ctx.operators if ctx.operators is not None else ctx.shared.setdefault("operators", {})
    # the bundles grouped by operator, each group in one apply_many (L4)
    groups: dict[str, list[tuple[int, Any, Any, Any]]] = {}
    bundles = [b for b in inputs.get("in", []) if isinstance(b, RayBundle | SpectrumBundle)]
    notes: list[str] = []
    try:
        for index, b in enumerate(bundles):
            key = f"particle|{b.wavelength_um!r}|{sorted(local.items(), key=str)!r}"
            if key not in store:
                structure = _structure(node, centre, b.wavelength_um)
                try:
                    store[key] = cba.prepare(structure, (b.wavelength_um,))
                except ValueError as e:
                    raise ExecutionError("exec.solver", f"CliffordBem: {e}", node.id) from e
                ctx.log.append(f"{node.id}: CliffordBem operator prepared")
            handle = store[key]
            z_front, n_bg = handle.structure.z_front_mm, handle.structure.n_background
            if isinstance(b, SpectrumBundle):
                waves, rot = incident_waves(b, centre, z_front, n_bg, notes)
            else:
                alive = np.asarray(b.alive)
                axis = np.asarray(b.direction)[alive].mean(axis=0)
                waves, rot = rays_to_spectrum(b, centre, axis / np.linalg.norm(axis), z_front, n_bg)
            groups.setdefault(key, []).append((index, b, waves, rot))
        results: dict[int, Any] = {}
        for key, group in groups.items():
            handle = store[key]
            try:  # recycling only pays across several solves; one bundle stays as with apply
                applied = cba.apply_many(
                    handle, [g[2] for g in group], n_theta, n_phi, recycle=len(group) > 1
                )
            except (ValueError, RuntimeError) as e:
                raise ExecutionError("exec.solver", f"CliffordBem: {e}", node.id) from e
            for (index, b, waves, rot), result in zip(group, applied, strict=True):
                back = result["front"]  # backscattered, direction −1 on the front port
                path = (*b.path, f"{node.id}[back]")
                if isinstance(b, SpectrumBundle):
                    results[index] = scattered_spectrum(
                        back, centre, rot, out_distance(b), b, path, notes
                    )
                    kind = "spectrum"
                else:
                    results[index] = spectrum_to_rays(back, centre, rot, b, path)
                    kind = "rays"
                ctx.log.append(
                    f"{node.id}: CliffordBem, {len(waves.kx)} incident waves ({kind}), key {b.key}"
                )
            ctx.log.append(f"{node.id}: CliffordBem apply_many, {len(group)} incidence(s)")
    except FieldError as e:
        raise ExecutionError("exec.field", str(e), node.id) from e
    finally:
        ctx.log.extend(f"{node.id}: {line}" for line in notes)
    return [results[i] for i in sorted(results)]


particle.representations = frozenset("RS")  # type: ignore[attr-defined]


TYPE = NodeType(
    "fusion.particle",
    Category.DOMAIN,
    (light_in("in", "RHS"), light_out("out", "RHS", label="backscattered")),
    (
        ParamSpec("point", "point", [0.0, 0.0, 0.0], unit="mm", description="particle centre"),
        ParamSpec("radius_um", "number", 0.25, unit="µm", minimum=0.0),
        ParamSpec("index", "number", 1.5, description="refractive index of the particle"),
        ParamSpec("n_medium", "number", 1.0, description="index of the medium around it"),
        ParamSpec("refinement", "integer", 4, minimum=1, maximum=7, description="icosphere level"),
        ParamSpec("n_theta", "integer", 40, minimum=4, description="far-field quadrature"),
        ParamSpec("n_phi", "integer", 80, minimum=8),
    ),
    (("out", ("in",)),),
    ports=(PortSpec("front", "in", "out", "front"),),
    scatter=(("front", "front", "full", ()),),
    description="Particle (CliffordBem BEM): backscattered light of a focused beam",
)

EVALUATORS: dict[str, Evaluator] = {"fusion.particle": particle}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    registry.add(TYPE)
    return dict(EVALUATORS)


__all__ = [
    "OUT_WAVELENGTHS",
    "CoherenceKey",
    "incident_waves",
    "out_distance",
    "ray_powers",
    "rays_to_spectrum",
    "register",
    "scattered_spectrum",
    "spectrum_to_rays",
]
