"""Execution of graphs imported from Raytatouille (stage R15).

R15 runs the Raytatouille blocks of an imported graph as whole paths ("chain fusion"): a detector
block traces every path of the system that ends at it (`graph.meta["paths"]`, written by the
importer), for each field point and wavelength, and images the bundles coherently
(`imaging.coherent_image`). The other Raytatouille blocks only carry the structure; per-block
execution follows when solver blocks sit between them (R19). A muted block takes its paths out of
the computation (bypassing an element needs Raytatouille #380).

Rays: `rt.trace.make_rays` with a hexapolar pupil; for paths that it cannot aim (no stop, folded
and tilted paths) start rays on the same grid over the stop opening – or, without a stop, over the
entrance pupil diameter (`aperture.type = "epd"`) in front of the first surface – parallel along
the field angles. Amplitude √weight (Raytatouille's power-normalised weight), phase k₀·OPL with
Raytatouille's OPL (vacuum-equivalent).

Fields (R18b): if the plan gives the link into the detector the representation H, the rays are
taken right after the last surface before the detector (recorded path, with Raytatouille's
power-normalised polarisation matrix P), turned into hybrid fields (`field.rays_to_field`,
Fieldtracing ADR 0004: U = sqrt(s/2)·sqrt(n_in cos γ_in dA_in/(n_out cos γ_out dA_out))·
exp(−iπμ/2)·P_T·jones with s = weight/(½‖P_T‖²_F)) and propagated to the detector by
Fieldtracing's planner (`field.field_image`). Raytatouille's source is unpolarised: two orthogonal
input polarisations are two coherence groups with half the power each. The refractive indices
come from the recorded OPL (ΔOPL/length of the first and the last segment), the parity of the ray
map from the reflections (direction against the surface axis before and after each event).
Executors of this backend may run R and H (`execute.WITH_FIELDS`); fieldtracing is needed for H.

`shared` of the executor: "system" (raytatouille.System), optional "materials"
(MaterialLibrary), "rings" (pupil rings, default 30), "pixels" ([nx, ny], default [256, 256]),
"gaussian_waist" (Gaussian illumination exp(−ρ²/w²) with w relative to the pupil radius, e.g. a
laser; default uniform – a hard-edged beam, whose Fresnel diffraction H shows and R does not).
The detector block's `half_size` ([hu, hv] mm) overrides the detector aperture as image window.
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.field import FieldBundle, FieldError, field_image, rays_to_field
from rtt_fusion.graph.imaging import coherent_image
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle
from rtt_fusion.graph.model import Node
from rtt_fusion.graph.optics import hexapolar
from rtt_fusion.graph.planner import link_key


def _structural(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[RayBundle]:
    return []  # the light of Raytatouille blocks is computed along whole paths at the detector


def _compiled(ctx: Context) -> Any:
    import raytatouille as rt

    if "compiled" not in ctx.shared:
        system = ctx.shared.get("system")
        if system is None:
            raise ExecutionError(
                "exec.no_system", "the executor needs shared['system'] (rt.System)"
            )
        ctx.shared["compiled"] = rt.compile(
            system, ctx.shared.get("materials", rt.MaterialLibrary())
        )
    return ctx.shared["compiled"]


def _start_rays(
    ctx: Context, compiled: Any, path: str, rings: int, field: int, wavelength: int
) -> Any:
    """Rays of one field point and wavelength: make_rays, or start rays over the stop / EPD."""
    import raytatouille as rt

    try:
        return rt.trace.make_rays(
            compiled,
            rt.trace.HexapolarPupil(rings=rings),
            path=path,
            fields=[field],
            wavelength=wavelength,
        )
    except rt.RaytatouilleError:
        pass  # not aimable (no stop, tilted or folded path): build start rays below
    d = json.loads(ctx.shared["system"].to_json())
    layouts = rt.layout.surfaces(compiled)
    by_id = {L.id: L for L in layouts}
    events = next((p.get("events") for p in d.get("paths", []) if p.get("name") == path), "auto")
    first = by_id[events[0]["surface"]] if isinstance(events, list) and events else layouts[0]
    if first.kind == "stop" and first.aperture is not None and first.aperture[0] == "circular":
        radius, centre, rot = (
            first.aperture[1]["radius"],
            np.asarray(first.translation),
            np.asarray(first.rotation),
        )
    elif d.get("aperture", {}).get("type") == "epd":
        value = d["aperture"]["value"]
        radius = 0.5 * float(value["value"] if isinstance(value, dict) else value)
        centre, rot = np.array([0.0, 0.0, float(first.translation[2]) - 1.0]), np.eye(3)
    else:
        raise ExecutionError(
            "exec.aiming", f"path {path!r} cannot be aimed and has no stop first / EPD"
        )
    if (
        not d.get("object", {}).get("at_infinity", True)
        or d.get("fields", {}).get("type", "angle_deg") != "angle_deg"
    ):
        raise ExecutionError("exec.aiming", "start rays for field angles at infinity only (R15)")
    f = d.get("fields", {}).get("points", [{}])[field]
    direction = np.array(
        [
            np.tan(np.radians(float(f.get("x", 0.0)))),
            np.tan(np.radians(float(f.get("y", 0.0)))),
            1.0,
        ]
    )
    direction /= np.linalg.norm(direction)
    px, py = hexapolar(rings)
    on_plane = np.stack([px * radius, py * radius, np.zeros_like(px)], -1) @ rot.T + centre
    start = on_plane - direction * (1.0 if first.kind != "stop" else 1.0 + radius)
    rays = rt.trace.RayBatch(len(px))
    for name, values in (
        ("pos_x", start[:, 0]),
        ("pos_y", start[:, 1]),
        ("pos_z", start[:, 2]),
        ("dir_x", np.full(len(px), direction[0])),
        ("dir_y", np.full(len(px), direction[1])),
        ("dir_z", np.full(len(px), direction[2])),
        ("pupil_x", px),
        ("pupil_y", py),
    ):
        np.asarray(getattr(rays, name))[:] = values
    np.asarray(rays.field)[:] = field
    np.asarray(rays.wl)[:] = wavelength
    return rays


def detector(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
    import raytatouille as rt

    compiled = _compiled(ctx)
    paths = [
        name
        for name, blocks in ctx.graph.meta.get("paths", {}).items()
        if blocks and blocks[-1] == node.id
    ]
    muted = {n.id for n in ctx.graph.nodes.values() if n.muted}
    rings = int(ctx.shared.get("rings", 30))
    fields_wanted = _input_representation(ctx, node) == "H"
    layouts = list(rt.layout.surfaces(compiled))  # in the order of the surface indices
    waist = float(ctx.shared.get("gaussian_waist", 0.0) or 0.0)
    bundles: list[RayBundle] = []
    fields: list[FieldBundle] = []
    for path in paths:
        if muted & set(ctx.graph.meta["paths"][path]):
            ctx.log.append(f"{node.id}: path {path} left out (muted block on it)")
            continue
        for field in range(compiled.field_count):
            for wl in range(len(compiled.wavelengths_um)):
                rays = _start_rays(ctx, compiled, path, rings, field, wl)
                pupil = np.stack([np.asarray(rays.pupil_x), np.asarray(rays.pupil_y)], -1).copy()
                apod = (
                    np.exp(-np.sum(pupil * pupil, axis=1) / waist**2)
                    if waist > 0
                    else np.ones(len(pupil))
                )
                _, rec = rt.trace.trace(
                    compiled,
                    rays,
                    path=path,
                    record_path=True,
                    max_recorded_rays=max(rays.size, 1),
                    polarization=fields_wanted,
                )
                if fields_wanted:
                    try:
                        fields += _fields(rec, pupil, compiled, layouts, path, field, wl, apod)
                    except FieldError as e:
                        raise ExecutionError("exec.field", f"path {path}: {e}", node.id) from e
                    continue
                count = np.asarray(rec.count)
                idx = np.asarray(rec.ray_indices)
                last = count - 1
                reached = count == rec.slots
                rows = np.arange(len(idx))
                weight = np.asarray(rec.weight)[rows, last]
                bundles.append(
                    RayBundle(
                        key=CoherenceKey("source", field, wl, 0),
                        wavelength_um=float(compiled.wavelengths_um[wl]),
                        position=np.asarray(rec.position)[rows, last],
                        direction=np.asarray(rec.direction)[rows, last],
                        opl=np.asarray(rec.opl)[rows, last],
                        amplitude=(np.sqrt(np.clip(weight, 0.0, None)) * apod[idx]).astype(complex),
                        source_xy=pupil[idx],
                        alive=reached,
                        path=(path,),
                    )
                )
    surface = next((s for s in node.params.get("surfaces", [])), None)
    L = next((x for x in rt.layout.surfaces(compiled) if x.id == surface), None)
    if L is None:
        raise ExecutionError(
            "exec.no_surface", f"detector surface {surface!r} not in the system", node.id
        )
    rot = np.asarray(L.rotation)
    plane = Plane(np.asarray(L.translation, dtype=float), rot[:, 2], rot[:, 0], rot[:, 1])
    if L.aperture is not None and L.aperture[0] == "rectangular":
        hu, hv = L.aperture[1]["half_width_x"], L.aperture[1]["half_width_y"]
    elif L.aperture is not None:
        hu = hv = float(L.aperture[1].get("radius", L.aperture[1].get("semi_axis_x", 1.0)))
    else:
        hu = hv = 1.0
    if node.params.get("half_size"):
        hu, hv = (float(v) for v in node.params["half_size"])
    nx, ny = (int(v) for v in ctx.shared.get("pixels", [256, 256]))
    if fields_wanted:
        try:
            return field_image(fields, plane, float(hu), float(hv), nx, ny, ctx.log)
        except FieldError as e:
            raise ExecutionError("exec.field", str(e), node.id) from e
    return coherent_image(bundles, plane, float(hu), float(hv), nx, ny)


def _input_representation(ctx: Context, node: Node) -> str:
    """The planned representation of the link into the detector (R without a plan)."""
    links = ctx.graph.links_to(node.id)
    if ctx.plan is None or not links:
        return "R"
    return ctx.plan.representations.get(link_key(links[0]), "R")


def _segment_index(rec: Any, rows: Any, a: Any, b: Any) -> float:
    """Refractive index of the segment between slots a and b: ΔOPL / length (median)."""
    pos = np.asarray(rec.position)
    opl = np.asarray(rec.opl)
    length = np.linalg.norm(pos[rows, b] - pos[rows, a], axis=1)
    d_opl = opl[rows, b] - opl[rows, a]
    ok = np.isfinite(length) & np.isfinite(d_opl) & (length > 1e-9)
    return float(np.median(d_opl[ok] / length[ok])) if np.any(ok) else 1.0


def _carried_frame(
    frame: Any,
    direction: Any,
    rows: Any,
    pupil: Any,
    idx: Any,
    port: int,
    rec: Any,
    layouts: list[Any],
) -> Any:
    """The source axes carried along the central ray up to the port slot: reflected at an event
    whose direction changes side of the surface axis, else turned by the smallest rotation."""
    centre = rows[int(np.argmin(np.linalg.norm(np.asarray(pupil)[idx[rows]], axis=1)))]
    f = np.asarray(frame, dtype=float)
    surfaces = list(rec.event_surfaces)
    for j in range(port):
        a, b = direction[centre, j], direction[centre, j + 1]
        if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
            continue
        sid = int(surfaces[j]) if j < len(surfaces) else -1
        axis = np.asarray(layouts[sid].rotation)[:, 2] if 0 <= sid < len(layouts) else None
        if axis is not None and float(a @ axis) * float(b @ axis) < 0:
            n = b - a
            n = n / np.linalg.norm(n)
            f = f - 2.0 * np.outer(f @ n, n)  # reflection
        else:
            f = f @ _rotation(a, b).T
    return f


def _rotation(a: Any, b: Any) -> Any:
    """Smallest rotation taking the unit vector a to b (Rodrigues)."""
    v = np.cross(a, b)
    c = float(a @ b)
    s = float(np.linalg.norm(v))
    if s < 1e-15:
        return np.eye(3) if c > 0 else -np.eye(3)
    k = v / s
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + s * kx + (1 - c) * (kx @ kx)


def _fields(
    rec: Any,
    pupil: Any,
    compiled: Any,
    layouts: list[Any],
    path: str,
    field: int,
    wl: int,
    apod: Any,
) -> list[FieldBundle]:
    """Fields right after the last surface before the detector, one per input polarisation."""
    count = np.asarray(rec.count)
    idx = np.asarray(rec.ray_indices)
    reached = count == rec.slots
    if reached.sum() < 6 or rec.slots < 3:
        return []
    rows = np.flatnonzero(reached)
    port = rec.slots - 2  # state after the event before the detector
    direction = np.asarray(rec.direction)
    prt = np.asarray(rec.prt)[rows, port]
    k_out = direction[rows, port]
    k_in = direction[rows, 0]
    p_t = prt - np.einsum("ni,nj,njk->nik", k_out, k_out, prt)
    norm = 0.5 * np.sum(np.abs(p_t) ** 2, axis=(1, 2))
    weight = np.asarray(rec.weight)[rows, port]
    s = np.where(norm > 0, weight / np.where(norm > 0, norm, 1.0), 0.0)
    # two orthogonal input polarisations transverse to the start direction (unpolarised source)
    k0 = k_in.mean(axis=0)
    k0 = k0 / np.linalg.norm(k0)
    seed = np.eye(3)[int(np.argmin(np.abs(k0)))]  # the axis least parallel to the beam
    x = seed - k0 * float(k0 @ seed)
    x = x / np.linalg.norm(x)
    y = np.cross(k0, x)
    base = RayBundle(
        key=CoherenceKey("source", field, wl, 0),
        wavelength_um=float(compiled.wavelengths_um[wl]),
        position=np.asarray(rec.position)[rows, port],
        direction=k_out,
        opl=np.asarray(rec.opl)[rows, port],
        amplitude=(np.sqrt(np.clip(s / 2.0, 0.0, None)) * np.asarray(apod)[idx[rows]]).astype(
            complex
        ),
        source_xy=np.asarray(pupil)[idx[rows]],
        alive=np.ones(len(rows), dtype=bool),
        path=(path,),
        prt=prt.astype(complex),
        source_frame=_carried_frame(
            np.stack([x, y]), direction, rows, pupil, idx, port, rec, layouts
        ),
        index=_segment_index(rec, rows, port, port + 1),
        source_index=_segment_index(rec, rows, 0, 1),
    )
    out = []
    for p, jones in enumerate((x, y)):
        bundle = base.evolve(key=CoherenceKey("source", field, wl, p), jones=jones.astype(complex))
        out.append(rays_to_field(bundle))
    return out


EVALUATORS: dict[str, Evaluator] = {
    "rtt.source": _structural,
    "rtt.section": _structural,
    "rtt.splitter": _structural,
    "rtt.mirror": _structural,
    "rtt.detector": detector,
}
