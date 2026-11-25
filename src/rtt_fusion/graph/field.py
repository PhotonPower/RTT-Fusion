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
  as scalar amplitude times a complex unit polarisation (`light.split_field`, R31: elliptical
  polarisation is kept; before R31 a real direction, exact only for linear polarisation).
- `field_image`: the field propagated to a detector plane with Fieldtracing's planner
  (`fieldtracing.chain.plan_propagation`: geometric, Debye, SFT, … by its selection rules) and
  summed per coherence group. Detector and field plane must be parallel (tilted planes later).

Validity: the conditions of the operators in table 6.3 – a dense ray net without caustic at the
port plane (R → H), local Fresnel numbers ≫ 1 (H → R).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import Delaunay, QhullError

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle, split_field

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
    #: the branch of the light and the position of its last Raytatouille block in it, as for rays
    #: (`RayBundle.branch`, node-wise execution): kept from rays to field and back
    branch: tuple[str, int] | None = None

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
    return FieldBundle(bundle.key, plane, f, bundle.path, bundle.branch)


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
    # the field as a complex amplitude times a complex unit polarisation (R31): elliptical
    # polarisation is kept, linear polarisation gives a real e as before
    amplitude, e = split_field(u_glob)
    # polarisation matrix: the unitary map (u, v, n) of the field plane -> (e, k × e*, k), so
    # that prt·jones = e with jones = u and the orientation of the plane's axes is carried on
    jones = fb.plane.u.astype(np.complex128)
    b = np.cross(direction, np.conj(e))
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
        branch=fb.branch,
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
    x = np.asarray(f.x, dtype=np.float64)
    values = np.asarray(f.physical())
    out = np.zeros((len(q), 3), dtype=np.complex128)
    inside = np.all((q >= x[0]) & (q <= x[-1]), axis=1)
    if not inside.any():
        return out
    # the cubic spline only over the part of the grid the points need (with a margin, where the
    # spline of the whole grid differs by < 1e-9): the grid of a focus is large (2016² for a
    # pupil image), the detector small
    margin = 8
    lo = max(int(np.searchsorted(x, q[inside].min())) - margin, 0)
    hi = min(int(np.searchsorted(x, q[inside].max())) + margin, len(x))
    if hi - lo < 4:
        lo, hi = max(min(lo, len(x) - 4), 0), min(max(hi, 4), len(x))
    xc = x[lo:hi]
    qi = q[inside]
    for c in range(3):
        block = values[lo:hi, lo:hi, c]  # rows: v (q[:, 1]), columns: u (q[:, 0])
        for part, unit in ((np.real, 1.0), (np.imag, 1j)):
            spline = RectBivariateSpline(xc, xc, part(block), kx=3, ky=3, s=0)
            out[inside, c] += unit * spline.ev(qi[:, 1], qi[:, 0])
    return out


def _onto(fb: FieldBundle, plane: Plane, notes: list[str]) -> tuple[FieldBundle, float]:
    """The field on a plane parallel to `plane` (resampled geometrically if tilted) and the
    distance along its normal to `plane`."""
    n = fb.plane.normal
    if abs(float(n @ plane.normal)) < 1.0 - 1e-12:
        # tilted planes: resample geometrically onto a plane parallel to the target (H -> R -> H)
        n_t = plane.normal if float(n @ plane.normal) > 0 else -plane.normal
        rays = field_to_rays(fb)
        # samples whose rays run along the plane or away from it (a wide field, e.g. the far
        # field of a scatterer with nearly grazing directions) do not reach it: left out (R30d)
        forward = np.asarray(rays.direction) @ n_t > 1e-9
        if not forward.all():
            notes.append(
                f"{fb.key.source}: {int((~forward).sum())} of {len(forward)} samples left out "
                "(their rays do not cross the tilted plane forward)"
            )
            rays = rays.evolve(alive=np.asarray(rays.alive) & forward)
        fb = rays_to_field(rays, Plane.make(fb.plane.origin, n_t))
        n = fb.plane.normal
        notes.append(f"{fb.key.source}: field resampled onto a plane parallel to the detector")
    return fb, float((plane.origin - fb.plane.origin) @ n)


def landing(fb: FieldBundle, distance: float) -> FloatArray | None:
    """Where the samples' rays cross the plane `distance` mm on, in the field's coordinates
    (sample + distance · k⊥/k_z); None for a field that is not a cloud in space or has samples
    running backwards."""
    from fieldtracing.field import GridlessField

    f = fb.field
    if not isinstance(f, GridlessField) or f.domain != "space" or len(f.points) == 0:
        return None
    k = f.local_wave_vector()
    if np.any(k[:, 2] <= 0.0):
        return None
    return np.asarray(np.asarray(f.points) + distance * k[:, :2] / k[:, 2:3], dtype=np.float64)


class _OutsideWindow(Exception):
    """A field whose rays all land far outside the detector window (propagate_to, field_sums)."""


def recentred(fb: FieldBundle, distance: float) -> FieldBundle:
    """The field with its coordinates centred where its rays land after `distance` mm (#88).

    Debye puts the focus, Rayleigh–Sommerfeld its output grid on the axis of the field's
    coordinates; RTT-Fusion's port planes put that axis through the bundle's centroid. A bundle
    that converges off that axis (an off-axis field point, one diffraction order of several)
    then focused onto a periodic copy or outside the computed grid. The new origin is the mean
    of the landing points (sample + distance · k⊥/k_z), weighted with |U|²; the samples keep
    their phase, gradient, Hessian and amplitude, the plane its axes (only shifted in itself)."""
    from fieldtracing.field import GridlessField

    f = fb.field
    if not isinstance(f, GridlessField) or f.domain != "space" or len(f.points) == 0:
        return fb
    land = landing(fb, distance)
    w = np.sum(np.abs(np.asarray(f.amplitude)) ** 2, axis=1)
    if land is None or not float(w.sum()) > 0.0:
        return fb
    c = (land * w[:, None]).sum(axis=0) / float(w.sum())
    scale = float(np.max(np.abs(np.asarray(f.points)))) if len(f.points) else 0.0
    if float(np.hypot(c[0], c[1])) <= 1e-12 * max(scale, 1.0):
        return fb
    p = fb.plane
    moved = Plane(
        np.asarray(p.origin) + c[0] * np.asarray(p.u) + c[1] * np.asarray(p.v),
        p.normal,
        p.u,
        p.v,
    )
    shifted = replace(f, points=np.asarray(f.points) - c)
    return FieldBundle(fb.key, moved, shifted, fb.path, fb.branch)


#: output samples per direction of Fieldtracing's Rayleigh–Sommerfeld step (its `n_out` default)
RS_SAMPLES = 64


def plan_onto(
    fb: FieldBundle,
    distance: float,
    window: float | None,
    tolerance: float | None,
    notes: list[str],
) -> Any:
    """Fieldtracing's plan for `distance` mm. A Rayleigh–Sommerfeld step computes a grid of
    RS_SAMPLES² points λ/(8 sin θ) apart around the axis; if that does not cover the image
    window, the plan is Debye (which takes the window) and its estimated error goes to `notes`
    – otherwise the image would be cut to a few micrometres around the axis."""
    from fieldtracing.chain import plan_propagation

    from rtt_fusion.graph._fieldtracing_fix import install

    install()  # Fieldtracing 0.2.x: a curvature fit that can loop forever (#135)

    extra: dict[str, Any] = {} if tolerance is None else {"tolerance": tolerance}
    plan = plan_propagation(fb.field, distance, focus_window=window, **extra)
    if not plan.operators:  # a method without an operator chain yet (Fieldtracing: "mixed")
        method, why = _executable_method(plan.situation, tolerance)
        notes.append(
            f"{fb.key.source}: the planner chose {plan.method} ({plan.reason}), which "
            f"Fieldtracing cannot compute yet – {method} instead ({why})"
        )
        plan = plan_propagation(fb.field, distance, focus_window=window, method=method, **extra)
        if not plan.operators:
            raise FieldError(f"no executable propagation method ({plan.method}) for this field")
    if plan.method != "rayleigh_sommerfeld" or window is None:
        return plan
    covered = RS_SAMPLES * float(plan.parameters.get("step_out_mm", 0.0))
    if covered >= window:
        return plan
    debye = plan_propagation(fb.field, distance, focus_window=window, method="debye", **extra)
    error = float(debye.budget.get("debye", float("nan")))
    eps = float(debye.parameters.get("tolerance", tolerance or 0.02))
    notes.append(
        f"{fb.key.source}: Rayleigh–Sommerfeld would cover {covered * 1e3:.3g} µm of the "
        f"{window:g} mm window – Debye instead, estimated error {error:.2g} (ε = {eps:g})"
    )
    return debye


def _executable_method(situation: Any, tolerance: float | None) -> tuple[str, str]:
    """The method the selection rules after "mixed" give (Fieldtracing's rules 4 and 6): near a
    caustic Debye (if its error 1.5/N_F ≤ ε) or Rayleigh–Sommerfeld, else SFT or the direct FFT,
    whichever needs fewer samples."""
    from fieldtracing.selection import debye_error

    eps = 0.02 if tolerance is None else tolerance
    if situation.target_fresnel * eps < 1.0:
        err = debye_error(situation.fresnel_number)
        if err <= eps:
            return "debye", f"target near a caustic; Debye error {err:.3g} ≤ ε"
        return "rayleigh_sommerfeld", f"target near a caustic; Debye error {err:.3g} > ε"
    if situation.sft_samples < situation.fft_samples:
        return "sft", "fewer samples than the direct FFT"
    return "fft", "fewer samples than the SFT"


def near_focus(
    fb: FieldBundle, plane: Plane, window: float | None = None, tolerance: float | None = None
) -> str | None:
    """Why `plane` needs wave optics, or None: Fieldtracing's selection rule for a target near a
    caustic (rays cross or focus there, 1/F_target > ε), the reason in its words."""
    from fieldtracing.chain import plan_propagation

    from rtt_fusion.graph._fieldtracing_fix import install

    install()  # Fieldtracing 0.2.x: a curvature fit that can loop forever (#135)

    fb, distance = _onto(fb, plane, [])
    if abs(distance) < 1e-9:
        return None
    extra: dict[str, Any] = {} if tolerance is None else {"tolerance": tolerance}
    plan = plan_propagation(fb.field, distance, focus_window=window, **extra)
    if plan.method in ("debye", "rayleigh_sommerfeld"):
        return str(plan.reason)
    if not plan.operators:  # "mixed": the method propagate_to takes instead decides
        method, why = _executable_method(plan.situation, tolerance)
        if method in ("debye", "rayleigh_sommerfeld"):
            return f"{plan.reason}; {why}"
    return None


def propagate_to(
    fb: FieldBundle,
    plane: Plane,
    window: float | None = None,
    tolerance: float | None = None,
    skip_outside: bool = False,
) -> tuple[FieldBundle, list[str]]:
    """The field on `plane` (Fieldtracing's planner) and the log of the steps. The returned
    bundle keeps the coordinates of a plane parallel to `plane` (its `plane` attribute is the
    frame of the field coordinates, shifted along the normal onto the target). `tolerance`: ε of
    the planner's selection rules (None: Fieldtracing's default, 0.02); a larger ε lets a cheaper
    method through, e.g. Debye instead of Rayleigh–Sommerfeld near a focus. With a window, a spot
    beside the field's axis is computed where it lies (#88): Rayleigh–Sommerfeld about where the
    rays land (`recentred`), Debye with a window holding the detector's window and the spot;
    `skip_outside`: a field landing entirely more than a window width beside the window raises
    `_OutsideWindow`."""
    notes: list[str] = []
    fb, distance = _onto(fb, plane, notes)
    if abs(distance) < 1e-9:
        return fb, notes
    if window is None:
        return _run(fb, distance, plan_onto(fb, distance, None, tolerance, notes), notes)
    land = landing(fb, distance)
    du, dv = _in_plane(fb, plane.origin)
    if land is not None and skip_outside:
        # every ray lands more than a window width beside the detector window: its light does
        # not reach the window (and a Debye or RS grid out to it would be huge)
        off = land - np.array([du, dv])
        far = 1.5 * window
        if (
            off[:, 0].min() > far
            or off[:, 0].max() < -far
            or off[:, 1].min() > far
            or off[:, 1].max() < -far
        ):
            raise _OutsideWindow(
                f"{fb.key.source}: the field lands {np.hypot(*off.mean(axis=0)):.3g} mm from "
                f"the window's centre, outside the {window:g} mm window – left out"
            )
    centred = recentred(fb, distance)
    if centred is fb:  # the light lands on the field's axis
        return _run(fb, distance, plan_onto(fb, distance, window, tolerance, notes), notes)
    # (#88) Rayleigh–Sommerfeld and Debye put their output on the axis of the field's
    # coordinates. Rayleigh–Sommerfeld integrates the field itself: centred on where the light
    # lands, its grid holds the spot. Debye works best with the pupil centred on that axis
    # (about the focus, a pupil beside the axis gains a few per cent of power and peak in
    # Fieldtracing): it stays on the pupil's axis with a window that holds the detector's window
    # and the spot, while the spot lies within half the pupil radius of it and that needs no more
    # than DEBYE_SAMPLES samples – otherwise about the focus.
    cu, cv = _in_plane(centred, plane.origin)
    near: list[str] = []
    plan_c = plan_onto(centred, distance, window + 2.0 * max(abs(cu), abs(cv)), tolerance, near)
    if plan_c.method == "rayleigh_sommerfeld":
        return _run(centred, distance, plan_c, notes + near)
    su, sv = _in_plane(fb, centred.plane.origin)  # the spot, from the field's axis
    own = max(window + 2.0 * max(abs(du), abs(dv)), window + 2.0 * max(abs(su), abs(sv)))
    plan = plan_onto(fb, distance, own, tolerance, notes)
    if plan.method == "debye":
        offset = float(np.hypot(su, sv))
        pupil = float(plan.parameters.get("length_u_mm", 0.0))
        samples = float(plan.parameters.get("samples", 0.0))
        # about its own axis Debye takes the directions within the pupil's cone: a spot further
        # off than half the pupil radius (a steep sub-beam) lies outside it
        if offset > 0.5 * pupil or samples > DEBYE_SAMPLES:
            notes.append(
                f"{fb.key.source}: Debye about the focus, {offset:.3g} mm beside the field's axis "
                f"(pupil radius {pupil:.3g} mm"
                + (f", {samples:.0f} samples" if samples > DEBYE_SAMPLES else "")
                + ")"
            )
            return _run(centred, distance, plan_c, notes + near)
    return _run(fb, distance, plan, notes)


def _in_plane(fb: FieldBundle, point: Any) -> tuple[float, float]:
    """A point's coordinates in the field's plane (u, v from its origin)."""
    rel = np.asarray(point, dtype=float) - np.asarray(fb.plane.origin)
    return float(rel @ fb.plane.u), float(rel @ fb.plane.v)


#: most samples across the k-disc for Debye about the field's own axis (#88); beyond, about the
#: focus (a padded grid of 4·n per axis: 768 → 3072², about 450 MB for three components)
DEBYE_SAMPLES = 768


#: most memory a field propagation's grid may take (#135: a Debye grid of 19944² for an order
#: landing far beside the window asked for 17.8 GiB); beyond, the propagation is refused
FIELD_GRID_BYTES = 2 * 2**30


def grid_bytes(plan: Any) -> float:
    """Estimated memory of a plan's largest grid: three complex components on the padded grid
    (Debye: 4·n per axis, the others 2·n), n the plan's samples per axis."""
    n = float(plan.parameters.get("samples", 0.0) or 0.0)
    pad = 4.0 if plan.method == "debye" else 2.0
    return 3 * 16 * (pad * n) ** 2


def _run(
    fb: FieldBundle, distance: float, plan: Any, notes: list[str]
) -> tuple[FieldBundle, list[str]]:
    need = grid_bytes(plan)
    if need > FIELD_GRID_BYTES:
        n = float(plan.parameters.get("samples", 0.0))
        raise FieldError(
            f"{fb.key.source}: the {plan.method} propagation over {distance:g} mm would need a "
            f"grid of {n:.0f} samples per axis ({need / 2**30:.1f} GiB, more than "
            f"{FIELD_GRID_BYTES / 2**30:g} GiB) – the light lands far from the field's axis or "
            "the window is large; a smaller image window, fewer orders or the link on R help"
        )
    result, log = plan.run(fb.field)
    notes.append(f"{fb.key.source}: H propagation {distance:g} mm by {plan.method} ({plan.reason})")
    notes += [f"  {e.operator}: {w}" for e in log for w in e.warnings]
    return FieldBundle(fb.key, fb.plane, result, fb.path, fb.branch), notes


def pixel_centres(plane: Plane, half_u: float, half_v: float, nx: int, ny: int) -> FloatArray:
    """Global positions (ny·nx, 3) of the pixel centres of a detector on `plane`, row by row."""
    u = (np.arange(nx) + 0.5) / nx * 2 * half_u - half_u
    v = (np.arange(ny) + 0.5) / ny * 2 * half_v - half_v
    uu, vv = np.meshgrid(u, v)
    return np.asarray(plane.origin + uu.reshape(-1, 1) * plane.u + vv.reshape(-1, 1) * plane.v)


def field_sums(
    fields: list[FieldBundle],
    plane: Plane,
    half_u: float,
    half_v: float,
    nx: int,
    ny: int,
    log: list[str] | None = None,
    progress: Callable[[float, str], None] | None = None,
    tolerance: float | None = None,
    sums: dict[CoherenceKey, ComplexArray] | None = None,
) -> tuple[dict[CoherenceKey, ComplexArray], float]:
    """The global field (ny·nx, 3) on the pixel centres per coherence group, added to `sums`,
    and the power of the fields (`field_power`)."""
    pixels = pixel_centres(plane, half_u, half_v, nx, ny)
    sums = {} if sums is None else sums
    power = 0.0
    for i, fb in enumerate(fields):
        if progress is not None:
            progress(i / max(len(fields), 1), f"field propagation {i + 1}/{len(fields)}")
        power += field_power(fb)
        try:
            target, notes = propagate_to(
                fb, plane, window=2.0 * max(half_u, half_v), tolerance=tolerance, skip_outside=True
            )
        except _OutsideWindow as outside:
            if log is not None:
                log.append(str(outside))
            continue
        if log is not None:
            log.extend(notes)
        rot = _frame(target.plane)
        rel = pixels - target.plane.origin
        q = np.stack([rel @ target.plane.u, rel @ target.plane.v], axis=-1)
        e = _sample(target.field, q) @ rot
        sums[fb.key] = np.asarray(sums[fb.key] + e if fb.key in sums else e, dtype=np.complex128)
    return sums, power


def image_of(
    sums: dict[CoherenceKey, ComplexArray],
    half_u: float,
    half_v: float,
    nx: int,
    ny: int,
    power: float,
) -> Image:
    """Intensity Σ_groups |E|² of the pixel fields per coherence group (`field_sums`)."""
    values = np.zeros(nx * ny)
    for total in sums.values():
        values = values + np.sum(np.abs(total) ** 2, axis=1)
    return Image(values.reshape(ny, nx), (-half_u, half_u, -half_v, half_v), power, list(sums))


def field_image(
    fields: list[FieldBundle],
    plane: Plane,
    half_u: float,
    half_v: float,
    nx: int,
    ny: int,
    log: list[str] | None = None,
    progress: Callable[[float, str], None] | None = None,
    tolerance: float | None = None,
) -> Image:
    """Intensity on the pixel grid of `plane` from fields (amplitudes per coherence group);
    `tolerance`: ε of the propagation planner (`propagate_to`)."""
    sums, power = field_sums(fields, plane, half_u, half_v, nx, ny, log, progress, tolerance)
    return image_of(sums, half_u, half_v, nx, ny, power)


def field_tolerance(shared: dict[str, Any]) -> float | None:
    """ε of the field propagation from the run's settings (`shared["field_tolerance"]`, e.g. the
    Workflow panel); None or 0: Fieldtracing's default."""
    value = shared.get("field_tolerance")
    return float(value) if value else None
