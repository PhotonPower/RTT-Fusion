"""Representation H on links: hybrid fields of Fieldtracing (stage R18, fusion plan 6.2/6.3).

A `FieldBundle` is one coherence group as a Fieldtracing `GridlessField` (or, after a propagation
into a focus, a `GriddedField`) on a port plane: coordinates (u, v) of the plane in mm, field
components in the plane's frame (u, v, n), n the propagation axis. The phase ψ is k₀·OPL with the
OPL of the rays, so fields of different paths share one phase reference and interfere.

Operators (table 6.3):

- `rays_to_field` (R → H): the rays are extended in a straight line to the port plane; amplitude
  by flux conservation along the ray tube with the polarisation of the rays (Fieldtracing ADR 0004,
  `fieldtracing.rtt_bridge.field_operator` written out for ray bundles),
      U = a·sqrt(I₀·cos γ_in/(cos γ_out·|J|))·exp(−iπμ/2)·(I − k kᵀ)·prt·jones,
  J and μ from the mesh of the start coordinates (`fieldtracing.mesh.mesh_jacobian_cloud`), the
  Hessian of ψ fitted on that mesh as in `fieldtracing.rtt_bridge.field_from_rays`.
- `field_to_rays` (H → R): one ray per sample along the local wave vector, OPL = ψ/k₀, the field
  as scalar amplitude times a real polarisation direction (exact for linear polarisation).
- `field_image`: the field propagated to a detector plane with Fieldtracing's planner
  (`fieldtracing.chain.plan_propagation`: geometric, Debye, SFT, … by its selection rules) and
  summed per coherence group. Detector and field plane must be parallel (tilted planes later).

Validity: the conditions of the operators in table 6.3 – a dense ray net without caustic at the
port plane (R → H), local Fresnel numbers ≫ 1 (H → R).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import Delaunay, QhullError

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle

if TYPE_CHECKING:
    from fieldtracing.field import GriddedField, GridlessField

# fieldtracing (optional dependency "modules") is imported where fields are made or propagated,
# so that the graph package works without it as long as no link carries H.


class FieldError(ValueError):
    """A conversion or propagation that is not admissible (caustic, tilted planes, …)."""


@dataclass(frozen=True)
class FieldBundle:
    """One coherence group as a hybrid field on a port plane."""

    key: CoherenceKey
    plane: Plane
    field: GridlessField | GriddedField
    path: tuple[str, ...] = ()

    @property
    def wavelength_um(self) -> float:
        return float(self.field.wavelength_um)

    @property
    def k0(self) -> float:
        return float(2.0 * np.pi / (self.wavelength_um * 1e-3))


def _frame(plane: Plane) -> FloatArray:
    """Rows u, v, n: global -> local is x @ R.T, local -> global is x @ R."""
    return np.asarray(np.stack([plane.u, plane.v, plane.normal]))


def port_plane(bundle: RayBundle) -> Plane:
    """The plane through the centroid of the living rays, perpendicular to their mean direction."""
    alive = np.asarray(bundle.alive)
    d = np.asarray(bundle.direction)[alive].mean(axis=0)
    return Plane.make(np.asarray(bundle.position)[alive].mean(axis=0), d / np.linalg.norm(d))


def rays_to_field(bundle: RayBundle, plane: Plane | None = None) -> FieldBundle:
    """R -> H on `plane` (default `port_plane`)."""
    alive = np.asarray(bundle.alive)
    from fieldtracing.field import GridlessField
    from fieldtracing.mesh import fit_derivatives, mesh_jacobian_cloud

    if alive.sum() < 6:
        raise FieldError("too few rays for a field (need a mesh)")
    plane = plane or port_plane(bundle)
    rot = _frame(plane)
    pos = (np.asarray(bundle.position)[alive] - plane.origin) @ rot.T
    d_glob = np.asarray(bundle.direction)[alive]
    d = d_glob @ rot.T
    if np.any(d[:, 2] <= 1e-9):
        raise FieldError("rays do not all cross the port plane forward")
    t = -pos[:, 2] / d[:, 2]
    landing = pos[:, :2] + t[:, None] * d[:, :2]
    n = float(bundle.index)
    opl = np.asarray(bundle.opl)[alive] + n * t  # straight on in the present medium
    start = np.asarray(bundle.source_xy)[alive]
    if bundle.source_frame is not None:
        # express the start coordinates in the port axes along the carried source axes: the
        # rotation or mirroring of the ray map by the path's geometry drops out, and the sign of
        # J and the Maslov index count caustics only
        m = rot[:2] @ np.asarray(bundle.source_frame, dtype=float).T  # (2, 2)
        u_svd, _, vt_svd = np.linalg.svd(m)
        start = start @ (u_svd @ vt_svd).T
    cloud = mesh_jacobian_cloud(start, landing)
    cos_in = np.broadcast_to(np.asarray(bundle.source_cos, dtype=float), (len(alive),))[alive]
    ratio = (
        bundle.irradiance * bundle.source_index * cos_in / (n * d[:, 2] * np.abs(cloud.determinant))
    )
    pol = bundle.field_direction()[alive]
    pol = pol - np.sum(pol * d_glob, axis=1, keepdims=True) * d_glob  # (I − k kᵀ)·prt·jones
    scale = (
        np.asarray(bundle.amplitude)[alive]
        * np.sqrt(np.where(np.isfinite(ratio), ratio, 0.0))
        * np.exp(-0.5j * np.pi * cloud.maslov)
    )
    u_local = (scale[:, None] * pol) @ rot.T
    k0 = bundle.k0
    gradient = k0 * n * d[:, :2]
    d_grad, _ = fit_derivatives(start, gradient, cloud.triangles)
    hessian = np.full((len(start), 2, 2), np.nan)
    ok = cloud.valid & ~cloud.caustic & np.isfinite(ratio)
    ok &= np.all(np.isfinite(d_grad.reshape(len(start), -1)), axis=1)
    hessian[ok] = np.einsum("nij,njk->nik", d_grad[ok], np.linalg.inv(cloud.matrix[ok]))
    ok &= np.all(np.isfinite(hessian.reshape(len(start), -1)), axis=1)
    if ok.sum() < 6:
        raise FieldError("the ray net has a caustic at the port plane (no geometric field there)")
    f = GridlessField(
        landing[ok],
        k0 * opl[ok],
        gradient[ok],
        hessian[ok],
        np.asarray(u_local[ok], dtype=np.complex128),
        bundle.wavelength_um,
        n,
    )
    return FieldBundle(bundle.key, plane, f, bundle.path)


def field_to_rays(fb: FieldBundle) -> RayBundle:
    """H -> R at the field's plane."""
    from fieldtracing.field import GridlessField

    f = fb.field
    if not isinstance(f, GridlessField) or f.domain != "space":
        raise FieldError("rays from a gridded or frequency-domain field are not supported yet")
    rot = _frame(fb.plane)
    pts = np.asarray(f.points)
    position = fb.plane.origin + pts[:, :1] * fb.plane.u + pts[:, 1:] * fb.plane.v
    k_local = f.local_wave_vector() / f.k
    direction = k_local @ rot
    u_glob = np.asarray(f.amplitude) @ rot
    big = np.argmax(np.abs(u_glob), axis=1)
    alpha = np.angle(u_glob[np.arange(len(pts)), big])
    e = np.real(u_glob * np.exp(-1j * alpha)[:, None])
    norm = np.linalg.norm(e, axis=1)
    e = e / np.where(norm > 0, norm, 1.0)[:, None]
    amplitude = np.sum(u_glob * e, axis=1)
    # polarisation matrix: the rotation (u, v, n) of the field plane -> (e, k × e, k), so that
    # prt·jones = e with jones = u and the orientation of the plane's axes is carried on
    jones = fb.plane.u.astype(np.complex128)
    b = np.cross(direction, e)
    prt = (
        np.einsum("ni,j->nij", e, fb.plane.u)
        + np.einsum("ni,j->nij", b, fb.plane.v)
        + np.einsum("ni,j->nij", direction, fb.plane.normal)
    ).astype(np.complex128)
    return RayBundle(
        key=fb.key,
        wavelength_um=fb.wavelength_um,
        position=position,
        direction=direction,
        opl=np.asarray(f.phase) / fb.k0,
        amplitude=np.asarray(amplitude, dtype=np.complex128),
        source_xy=pts.copy(),
        alive=np.ones(len(pts), dtype=bool),
        irradiance=1.0,
        path=fb.path,
        jones=jones,
        prt=prt,
        source_cos=np.asarray(k_local[:, 2]),
        source_frame=np.stack([fb.plane.u, fb.plane.v]),
        index=float(f.refractive_index),
        source_index=float(f.refractive_index),
    )


def field_power(fb: FieldBundle) -> float:
    """∫ n|U|² cos γ dA over the samples of a gridless field (W), by triangles."""
    from fieldtracing.field import GridlessField

    f = fb.field
    if not isinstance(f, GridlessField):
        g = f
        return float(np.sum(np.abs(g.values) ** 2) * g.step**2)
    pts = np.asarray(f.points)
    try:
        tri = Delaunay(pts).simplices
    except QhullError:
        return 0.0
    dens = (
        f.refractive_index
        * np.sum(np.abs(f.amplitude) ** 2, axis=1)
        * (f.local_wave_vector()[:, 2] / f.k)
    )
    a, b, c = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
    area = 0.5 * np.abs(
        (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    )
    return float(np.sum(area * dens[tri].mean(axis=1)))


def _sample(f: GridlessField | GriddedField, q: FloatArray) -> ComplexArray:
    """Physical field (local components) at local points q (N, 2); 0 outside."""
    from fieldtracing.field import GridlessField

    if isinstance(f, GridlessField):
        return np.asarray(f.residual_at(q) * np.exp(1j * f.phase_at(q))[:, None])
    x = np.asarray(f.x)
    values = np.asarray(f.physical())
    out = np.zeros((len(q), 3), dtype=np.complex128)
    for c in range(3):
        for part, unit in ((np.real, 1.0), (np.imag, 1j)):
            interp = RegularGridInterpolator(
                (x, x), part(values[:, :, c]), method="cubic", bounds_error=False, fill_value=0.0
            )
            out[:, c] += unit * interp(np.stack([q[:, 1], q[:, 0]], axis=-1))
    return out


def propagate_to(
    fb: FieldBundle, plane: Plane, window: float | None = None
) -> tuple[FieldBundle, list[str]]:
    """The field on `plane` (Fieldtracing's planner) and the log of the steps. The returned
    bundle keeps the coordinates of a plane parallel to `plane` (its `plane` attribute is the
    frame of the field coordinates, shifted along the normal onto the target)."""
    from fieldtracing.chain import plan_propagation

    n = fb.plane.normal
    notes: list[str] = []
    if abs(float(n @ plane.normal)) < 1.0 - 1e-12:
        # tilted planes: resample geometrically onto a plane parallel to the target (H -> R -> H)
        n_t = plane.normal if float(n @ plane.normal) > 0 else -plane.normal
        fb = rays_to_field(field_to_rays(fb), Plane.make(fb.plane.origin, n_t))
        n = fb.plane.normal
        notes.append(f"{fb.key.source}: field resampled onto a plane parallel to the detector")
    distance = float((plane.origin - fb.plane.origin) @ n)
    if abs(distance) < 1e-9:
        return fb, notes
    plan = plan_propagation(fb.field, distance, focus_window=window)
    result, log = plan.run(fb.field)
    notes.append(f"{fb.key.source}: H propagation {distance:g} mm by {plan.method} ({plan.reason})")
    notes += [f"  {e.operator}: {w}" for e in log for w in e.warnings]
    return FieldBundle(fb.key, fb.plane, result, fb.path), notes


def field_image(
    fields: list[FieldBundle],
    plane: Plane,
    half_u: float,
    half_v: float,
    nx: int,
    ny: int,
    log: list[str] | None = None,
) -> Image:
    """Intensity on the pixel grid of `plane` from fields (amplitudes per coherence group)."""
    u = (np.arange(nx) + 0.5) / nx * 2 * half_u - half_u
    v = (np.arange(ny) + 0.5) / ny * 2 * half_v - half_v
    uu, vv = np.meshgrid(u, v)
    pixels = plane.origin + uu.reshape(-1, 1) * plane.u + vv.reshape(-1, 1) * plane.v
    sums: dict[CoherenceKey, ComplexArray] = {}
    power = 0.0
    for fb in fields:
        power += field_power(fb)
        target, notes = propagate_to(fb, plane, window=2.0 * max(half_u, half_v))
        if log is not None:
            log.extend(notes)
        rot = _frame(target.plane)
        rel = pixels - target.plane.origin
        q = np.stack([rel @ target.plane.u, rel @ target.plane.v], axis=-1)
        e = _sample(target.field, q) @ rot
        sums[fb.key] = np.asarray(sums[fb.key] + e if fb.key in sums else e, dtype=np.complex128)
    values = np.zeros(nx * ny)
    for total in sums.values():
        values = values + np.sum(np.abs(total) ** 2, axis=1)
    return Image(values.reshape(ny, nx), (-half_u, half_u, -half_v, half_v), power, list(sums))
