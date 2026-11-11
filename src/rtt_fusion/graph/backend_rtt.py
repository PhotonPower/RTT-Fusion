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

`shared` of the executor: "system" (raytatouille.System), optional "materials"
(MaterialLibrary), "rings" (pupil rings, default 30), "pixels" ([nx, ny], default [256, 256]).
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.imaging import coherent_image
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle
from rtt_fusion.graph.model import Node
from rtt_fusion.graph.optics import hexapolar


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
    bundles: list[RayBundle] = []
    for path in paths:
        if muted & set(ctx.graph.meta["paths"][path]):
            ctx.log.append(f"{node.id}: path {path} left out (muted block on it)")
            continue
        for field in range(compiled.field_count):
            for wl in range(len(compiled.wavelengths_um)):
                rays = _start_rays(ctx, compiled, path, rings, field, wl)
                pupil = np.stack([np.asarray(rays.pupil_x), np.asarray(rays.pupil_y)], -1).copy()
                _, rec = rt.trace.trace(
                    compiled, rays, path=path, record_path=True, max_recorded_rays=max(rays.size, 1)
                )
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
                        amplitude=np.sqrt(np.clip(weight, 0.0, None)).astype(complex),
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
    nx, ny = (int(v) for v in ctx.shared.get("pixels", [256, 256]))
    return coherent_image(bundles, plane, float(hu), float(hv), nx, ny)


EVALUATORS: dict[str, Evaluator] = {
    "rtt.source": _structural,
    "rtt.section": _structural,
    "rtt.splitter": _structural,
    "rtt.mirror": _structural,
    "rtt.detector": detector,
}
