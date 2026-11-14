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

Node-wise (R25, ADR ports, poses, branches): every Raytatouille block passes the light of each
branch along its own events (a *run*: the events of the block in a path of the system), the
source starts two input polarisations per field point and wavelength and per branch (aimed along
a system path that starts in the same block), the detector images what arrives. Bundles carry
their branch and position in it (`RayBundle.branch`), so a block passed twice takes the right run
and a block that replaced a Raytatouille block (R19b) acts in between, placed where the replaced
surface was. Rays lost on the way stay in the bundle, marked dead (the image triangulates over the
start grid); a ray that starts on or behind a run's first surface misses it (Raytatouille), as
along a whole path.
On an H link the bundle gets the medium behind its last surface. Chosen with
shared["rtt_execution"] = "nodes", or automatically when a branch of the model is no path of the
system; "paths" (and "auto" otherwise) keeps whole paths. Node-wise equals whole paths to rounding
(W43).

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
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle, rotation_between
from rtt_fusion.graph.model import Node
from rtt_fusion.graph.optics import hexapolar
from rtt_fusion.graph.planner import link_key
from rtt_fusion.graph.pose import angles_xyz, offset_of, rotation_xyz
from rtt_fusion.graph.validate import node_sockets


def _passed(rec: Any) -> Any:
    """Rays that passed every event of a recorded path. Raytatouille gives a ray lost at event j
    count = j + 2 slots and lost_at = j – for the last event that is count = slots too, so the
    loss index decides."""
    return (np.asarray(rec.count) == rec.slots) & (np.asarray(rec.lost_at) == -1)


def _offsets(ctx: Context) -> dict[str, list[float]]:
    """Shifts of Raytatouille blocks per surface id (mm, frame of the parent assembly): the
    parameter `offset`, plus `drive` (set by a phase shifter linked to the block) along the
    block's surface normal (its local +z)."""
    out: dict[str, list[float]] = {}
    for n in ctx.graph.nodes.values():
        if not n.type_id.startswith("rtt.") or n.muted:
            continue
        offset = np.asarray(
            ([float(v) for v in (n.params.get("offset") or [])] + [0.0, 0.0, 0.0])[:3]
        )
        drive = float(n.params.get("drive", 0.0) or 0.0)
        if drive != 0.0:
            offset = offset + drive * _normal_in_parent(ctx, n)
        if np.any(np.abs(offset) > 0):
            for s in n.params.get("surfaces", []):
                out[str(s)] = [float(v) for v in offset]
    return out


def _apply_offsets(d: dict[str, Any], offsets: dict[str, list[float]]) -> dict[str, Any]:
    """Moves every element holding a listed surface by its offset (frame of its parent)."""

    def visit(node: dict[str, Any]) -> None:
        ids = {str(s.get("id")) for s in node.get("surfaces", []) if isinstance(s, dict)}
        hit = next((offsets[i] for i in ids if i in offsets), None)
        if hit is not None:
            pose = node.setdefault("pose", {})
            position = pose.setdefault("position", [0.0, 0.0, 0.0])
            for k in range(3):
                v = position[k]
                if isinstance(v, dict):
                    v["value"] = float(v.get("value", 0.0)) + hit[k]
                else:
                    position[k] = float(v) + hit[k]
        for child in node.get("children", []):
            if isinstance(child, dict):
                visit(child)

    visit(d["root"])
    return d


def _base_json(ctx: Context) -> dict[str, Any]:
    """The system to compute with: generated from the data of the graph's blocks (S2, ADR
    Nachtrag 1: the workflow is the truth), or – for a graph without them – shared["system"]."""
    generated = _generated(ctx)
    if generated is not None:
        return json.loads(generated)  # type: ignore[no-any-return]
    system = ctx.shared.get("system")
    if system is None:
        raise ExecutionError("exec.no_system", "the executor needs shared['system'] (rt.System)")
    d: dict[str, Any] = json.loads(system.to_json())
    return d


def _generated(ctx: Context) -> str | None:
    """The generated system as JSON text, made once per graph state (the blocks' data and the
    branches) and kept in shared."""
    from rtt_fusion.graph.importer import system_from_graph

    fingerprint = json.dumps(
        [
            [
                n.id,
                n.params.get("rtt"),
                n.params.get("runs"),
                n.params.get("rtt_system"),
                n.params.get("rtt_tree"),
                n.params.get("rtt_paths"),
            ]
            for n in ctx.graph.nodes.values()
            if n.type_id.startswith("rtt.") or "rtt" in n.params
        ]
        + [ctx.graph.meta.get("branches")],
        sort_keys=True,
    )
    cached = ctx.shared.get("generated_system")
    if cached is not None and cached[0] == fingerprint:
        return cached[1]  # type: ignore[no-any-return]
    try:
        d = system_from_graph(ctx.graph)
    except ValueError as e:
        raise ExecutionError("exec.system", str(e)) from e
    text = None if d is None else json.dumps(d)
    ctx.shared["generated_system"] = (fingerprint, text)
    return text


def _normal_in_parent(ctx: Context, node: Node) -> Any:
    """The surface normal of a block in the frame of its element's parent assembly: the global
    normal from the layout, the parent frame from three probe translations (cached)."""
    import raytatouille as rt

    cache = ctx.shared.setdefault("drive_normals", {})
    surfaces = [str(s) for s in node.params.get("surfaces", [])]
    if not surfaces:
        raise ExecutionError("exec.param", "drive: the block has no surface", node.id)
    key = surfaces[0]
    if key not in cache:
        materials = ctx.shared.get("materials", rt.MaterialLibrary())

        def surface_pose(d: dict[str, Any]) -> tuple[Any, Any]:
            compiled = rt.compile(rt.System.from_json(json.dumps(d)), materials)
            layout = next((x for x in rt.layout.surfaces(compiled) if x.id == key), None)
            if layout is None:
                raise ExecutionError("exec.param", f"drive: no surface {key!r}", node.id)
            return np.asarray(layout.translation, dtype=float), np.asarray(layout.rotation)

        t0, rot = surface_pose(_base_json(ctx))
        normal = rot[:, 2]
        columns = []
        for k in range(3):
            probe = [0.0, 0.0, 0.0]
            probe[k] = 1.0
            tk, _ = surface_pose(_apply_offsets(_base_json(ctx), {s: probe for s in surfaces}))
            columns.append(tk - t0)
        parent = np.column_stack(columns)  # parent unit vectors in global coordinates
        cache[key] = np.linalg.solve(parent, normal)
    return cache[key]


def _system_json(ctx: Context) -> dict[str, Any]:
    """The system as JSON with the block offsets applied: an element that holds a shifted
    surface moves by the offset (in the frame of its parent assembly); a block with a pose
    offset (R26) moves its tree node rigidly."""
    d = _base_json(ctx)
    offsets = _offsets(ctx)
    if offsets:
        d = _apply_offsets(d, offsets)
    for location, offset in _rigid(ctx).items():
        _apply_rigid(d, location, offset, ctx)
    return d


def _rigid(ctx: Context) -> dict[str, list[float]]:
    """Pose offsets of Raytatouille blocks (six values in the block's local frame) by the
    location of the block's node in the model tree."""
    out: dict[str, list[float]] = {}
    for n in ctx.graph.nodes.values():
        if not n.type_id.startswith("rtt.") or n.muted or "@" in n.id:
            continue
        o = offset_of(n)
        if np.any(o):
            out[str(n.params.get("location", ""))] = [float(v) for v in o]
    return out


def _variation(ctx: Context) -> str:
    """Key of the present geometry changes (offsets, drives, pose offsets) for the caches."""
    return json.dumps({"offsets": _offsets(ctx), "rigid": _rigid(ctx)}, sort_keys=True)


def _apply_rigid(d: dict[str, Any], location: str, offset: list[float], ctx: Context) -> None:
    """Moves the tree node at `location` by the offset in its own frame, rotation about its
    pivot (ADR section 4): with p -> position + pivot + R·(p − pivot), the composed pose keeps
    the pivot and has position + R·d and R·R_offset."""
    node: Any = d
    for part in [p for p in location.split("/") if p]:
        node = node[int(part)] if isinstance(node, list) else node.get(part)
        if node is None:
            raise ExecutionError("exec.pose", f"no node {location!r} in the system")
    owner = next(
        (n for n in ctx.graph.nodes.values() if n.params.get("location") == location), None
    )
    held: set[str] = set()

    def collect(t: dict[str, Any]) -> None:
        held.update(str(s.get("id")) for s in t.get("surfaces", []) if isinstance(s, dict))
        for c in t.get("children", []):
            if isinstance(c, dict):
                collect(c)

    collect(node)
    if owner is not None and held != {str(s) for s in owner.params.get("surfaces", [])}:
        raise ExecutionError(
            "exec.pose",
            "the block is part of a tree node; a pose offset moves whole tree nodes",
            owner.id,
        )
    pose = node.setdefault("pose", {})
    if (
        pose.get("order", "translate_first") != "translate_first"
        or pose.get("reference", "absolute") != "absolute"
    ):
        raise ExecutionError(
            "exec.pose", f"{location}: pose offsets need absolute, " "translate-first poses"
        )

    def values(key: str) -> list[float]:
        raw = pose.get(key, [0.0, 0.0, 0.0])
        return [float(v["value"]) if isinstance(v, dict) else float(v) for v in raw]

    def store(key: str, new: Any) -> None:
        raw = pose.setdefault(key, [0.0, 0.0, 0.0])
        for k in range(3):
            if isinstance(raw[k], dict):
                raw[k]["value"] = float(new[k])
            else:
                raw[k] = float(new[k])

    r = rotation_xyz(np.radians(values("rotation_deg")))
    position = np.asarray(values("position")) + r @ np.asarray(offset[:3])
    store("position", position)
    store("rotation_deg", np.degrees(angles_xyz(r @ rotation_xyz(offset[3:]))))


def _compiled(ctx: Context) -> Any:
    """The compiled system for the present block offsets (cached per offsets)."""
    import raytatouille as rt

    offsets = _offsets(ctx)
    generated = _generated(ctx)
    if generated is not None:  # S2: compiled from the system the workflow makes, per state
        key = f"{hash(generated)}|{_variation(ctx)}"
        cache = ctx.shared.setdefault("compiled_generated", {})
        if key not in cache:
            system = rt.System.from_json(json.dumps(_system_json(ctx)))
            cache[key] = rt.compile(system, ctx.shared.get("materials", rt.MaterialLibrary()))
        return cache[key]
    if not offsets and not _rigid(ctx):
        if "compiled" not in ctx.shared:
            if ctx.shared.get("system") is None:
                raise ExecutionError(
                    "exec.no_system", "the executor needs shared['system'] (rt.System)"
                )
            ctx.shared["compiled"] = rt.compile(
                ctx.shared["system"], ctx.shared.get("materials", rt.MaterialLibrary())
            )
        return ctx.shared["compiled"]
    key = _variation(ctx)
    cache = ctx.shared.setdefault("compiled_by_offsets", {})
    if key not in cache:
        system = rt.System.from_json(json.dumps(_system_json(ctx)))
        cache[key] = rt.compile(system, ctx.shared.get("materials", rt.MaterialLibrary()))
    return cache[key]


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
    d = _system_json(ctx)
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

    if _nodewise(ctx):
        return _detector_nodes(ctx, node, inputs)
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
    steps = max(len(paths) * compiled.field_count * len(compiled.wavelengths_um), 1)
    trace_part = 0.3 if fields_wanted else 1.0  # with fields the propagation takes the rest
    step = 0

    def tick(path: str, field: int, wl: int) -> None:
        nonlocal step
        ctx.report(
            trace_part * step / steps,
            f"{node.name or node.id}: {path}, field {field + 1}/{compiled.field_count}, "
            f"λ {wl + 1}/{len(compiled.wavelengths_um)}",
        )
        step += 1

    for path in paths:
        if muted & set(ctx.graph.meta["paths"][path]):
            ctx.log.append(f"{node.id}: path {path} left out (muted block on it)")
            continue
        foreign = [
            b
            for b in ctx.graph.meta["paths"][path]
            if b in ctx.graph.nodes and not ctx.graph.nodes[b].type_id.startswith("rtt.")
        ]
        if foreign:
            for field in range(compiled.field_count):
                for wl in range(len(compiled.wavelengths_um)):
                    tick(path, field, wl)
                    out = _split_path(ctx, node, path, foreign, field, wl, fields_wanted)
                    if fields_wanted:
                        fields += out
                    else:
                        bundles += out
            continue
        for field in range(compiled.field_count):
            for wl in range(len(compiled.wavelengths_um)):
                tick(path, field, wl)
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
                reached = _passed(rec)
                rows = np.arange(len(idx))
                weight = np.asarray(rec.weight)[rows, last]
                bundles.append(
                    RayBundle(
                        key=CoherenceKey("source", field, wl, 0),
                        wavelength_um=float(compiled.wavelengths_um[wl]),
                        position=np.asarray(rec.position)[rows, last],
                        direction=np.asarray(rec.direction)[rows, last],
                        opl=np.asarray(rec.opl)[rows, last]
                        + _start_opl(rec, rows, _object_index(compiled, wl)),
                        amplitude=(np.sqrt(np.clip(weight, 0.0, None)) * apod[idx]).astype(complex),
                        source_xy=pupil[idx],
                        alive=reached,
                        path=(path,),
                    )
                )
    return _image(ctx, node, compiled, bundles, fields, fields_wanted, trace_part)


def _image(
    ctx: Context,
    node: Node,
    compiled: Any,
    bundles: list[RayBundle],
    fields: list[FieldBundle],
    fields_wanted: bool,
    trace_part: float,
) -> Any:
    """The image on the detector surface: rays coherently (R) or fields propagated (H)."""
    import raytatouille as rt

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
            return field_image(
                fields,
                plane,
                float(hu),
                float(hv),
                nx,
                ny,
                ctx.log,
                lambda part, text: ctx.report(
                    trace_part + (1 - trace_part) * part, f"{node.name or node.id}: {text}"
                ),
            )
        except FieldError as e:
            raise ExecutionError("exec.field", str(e), node.id) from e
    return coherent_image(bundles, plane, float(hu), float(hv), nx, ny)


# ---------------------------------------------------------------------- split paths (R19b)


def _events_of_blocks(ctx: Context, path: str, system_json: dict[str, Any]) -> list[list[Any]]:
    """The events of `path` grouped by the blocks of `graph.meta["paths"][path]`: a Raytatouille
    block takes the consecutive events on its surfaces, any other block (a replaced sample) the
    events up to the next block's surfaces."""
    events = next(
        (p.get("events") for p in system_json.get("paths", []) if p.get("name") == path), None
    )
    if not isinstance(events, list):
        raise ExecutionError(
            "exec.split", f"path {path!r} needs explicit events to be split at a block"
        )
    blocks = ctx.graph.meta["paths"][path]

    def surfaces(block_id: str) -> set[str]:
        n = ctx.graph.nodes.get(block_id)
        if n is None or not n.type_id.startswith("rtt."):
            return set()
        return {str(s) for s in n.params.get("surfaces", [])}

    groups: list[list[Any]] = []
    k = 0
    for i, b in enumerate(blocks):
        own = surfaces(b)
        group: list[Any] = []
        if own:
            while k < len(events) and str(events[k]["surface"]) in own:
                group.append(events[k])
                k += 1
        else:
            nxt = surfaces(blocks[i + 1]) if i + 1 < len(blocks) else set()
            while k < len(events) and (not nxt or str(events[k]["surface"]) not in nxt):
                group.append(events[k])
                k += 1
        groups.append(group)
    if k != len(events):
        raise ExecutionError("exec.split", f"the events of path {path!r} do not follow its blocks")
    return groups


def _split_path(
    ctx: Context,
    detector_node: Node,
    path: str,
    foreign: list[str],
    field: int,
    wl: int,
    fields_wanted: bool,
) -> list[Any]:
    """A path with a non-Raytatouille block (e.g. a grating as the sample): Raytatouille traces up
    to it, the block's evaluator acts on the rays, Raytatouille traces each of its outputs on.
    Two orthogonal input polarisations (unpolarised source). Returns ray bundles at the detector
    (R) or fields after the last surface before it (H)."""
    import raytatouille as rt

    if len(set(foreign)) != 1:
        raise ExecutionError("exec.split", f"path {path!r}: one replaced block per path (R19b)")
    block = ctx.graph.nodes[foreign[0]]
    d = _system_json(ctx)
    groups = _events_of_blocks(ctx, path, d)
    blocks = ctx.graph.meta["paths"][path]
    i = blocks.index(block.id)
    before = [e for g in groups[:i] for e in g]
    own = groups[i]
    after = [e for g in groups[i + 1 :] for e in g]
    if not before or not after or not own:
        raise ExecutionError(
            "exec.split", f"path {path!r}: the block {block.id} needs light before and after it"
        )
    cache = ctx.shared.setdefault("split_systems", {})
    split_key = f"{path}|{_variation(ctx)}"
    if split_key not in cache:
        d["paths"] = [
            *d.get("paths", []),
            {"name": f"{path}#A", "events": before},
            {"name": f"{path}#B", "events": after},
        ]
        system = rt.System.from_json(json.dumps(d))
        cache[split_key] = rt.compile(system, ctx.shared.get("materials", rt.MaterialLibrary()))
    compiled = cache[split_key]
    layouts = list(rt.layout.surfaces(compiled))
    by_id = {L.id: L for L in layouts}
    # the block sits where the replaced surface was (pose of its first event's surface)
    L = by_id.get(str(own[0]["surface"]))
    params = dict(block.params)
    if L is not None and not params.get("keep_pose"):
        rot = np.asarray(L.rotation)
        params["point"] = [float(v) for v in np.asarray(L.translation)]
        params["normal"] = [float(v) for v in rot[:, 2]]
    placed = Node(block.id, block.type_id, block.name, params)
    evaluator = ctx.evaluators.get(block.type_id)
    if evaluator is None:
        raise ExecutionError("exec.no_evaluator", f"no evaluator for {block.type_id}", block.id)
    rings = int(ctx.shared.get("rings", 30))
    waist = float(ctx.shared.get("gaussian_waist", 0.0) or 0.0)
    # --- part A: up to the block
    rays = _start_rays(ctx, compiled, path, rings, field, wl)
    pupil = np.stack([np.asarray(rays.pupil_x), np.asarray(rays.pupil_y)], -1).copy()
    apod = np.exp(-np.sum(pupil * pupil, axis=1) / waist**2) if waist > 0 else np.ones(len(pupil))
    _, rec_a = rt.trace.trace(
        compiled,
        rays,
        path=f"{path}#A",
        record_path=True,
        max_recorded_rays=max(rays.size, 1),
        polarization=True,
    )
    rows = np.flatnonzero(_passed(rec_a))
    if rows.size < 3:
        return []
    idx = np.asarray(rec_a.ray_indices)
    last = rec_a.slots - 1
    direction_a = np.asarray(rec_a.direction)
    prt_a = np.asarray(rec_a.prt)[rows, last]
    k_a = direction_a[rows, last]
    p_t = prt_a - np.einsum("ni,nj,njk->nik", k_a, k_a, prt_a)
    norm = 0.5 * np.sum(np.abs(p_t) ** 2, axis=(1, 2))
    weight = np.asarray(rec_a.weight)[rows, last]
    s_a = np.where(norm > 0, weight / np.where(norm > 0, norm, 1.0), 0.0)
    k0 = direction_a[rows, 0].mean(axis=0)
    k0 = k0 / np.linalg.norm(k0)
    seed = np.eye(3)[int(np.argmin(np.abs(k0)))]
    x = seed - k0 * float(k0 @ seed)
    x = x / np.linalg.norm(x)
    y = np.cross(k0, x)
    frame = _carried_frame(np.stack([x, y]), direction_a, rows, pupil, idx, last, rec_a, layouts)
    base = RayBundle(
        key=CoherenceKey("source", field, wl, 0),
        wavelength_um=float(compiled.wavelengths_um[wl]),
        position=np.asarray(rec_a.position)[rows, last],
        direction=k_a,
        opl=np.asarray(rec_a.opl)[rows, last]
        + _start_opl(rec_a, rows, _object_index(compiled, wl)),
        # the field is amplitude·prt·jones; the splitter's factors are in prt (power-normalised
        # P), the polarisation-independent ones (s) in the amplitude; the block that follows
        # splits the field again into scalar and prt
        amplitude=(np.sqrt(np.clip(s_a / 2.0, 0.0, None)) * apod[idx[rows]]).astype(complex),
        source_xy=pupil[idx[rows]],
        alive=np.ones(len(rows), dtype=bool),
        path=(path,),
        prt=prt_a.astype(complex),
        source_frame=frame,
        index=_segment_index(rec_a, rows, last - 1, last) if last >= 1 else 1.0,
        source_index=_segment_index(rec_a, rows, 0, 1),
    )
    sockets = node_sockets(placed, ctx.registry) if ctx.registry is not None else ()
    in_name = next((s.name for s in sockets if s.direction.value == "in"), "in")
    out_name = next((s.name for s in sockets if s.direction.value == "out"), "out")
    out: list[Any] = []
    for p, jones in enumerate((x, y)):
        incident = base.evolve(
            key=CoherenceKey("source", field, wl, p), jones=jones.astype(complex)
        )
        for b in evaluator(ctx, placed, out_name, {in_name: [incident]}):
            out += _trace_on(ctx, compiled, layouts, path, b, wl, field, fields_wanted)
    return out


def _trace_on(
    ctx: Context,
    compiled: Any,
    layouts: list[Any],
    path: str,
    b: RayBundle,
    wl: int,
    field: int,
    fields_wanted: bool,
) -> list[Any]:
    """Part B: trace a bundle that a block sent on, and combine OPL, field and polarisation."""
    return _trace_run(compiled, layouts, f"{path}#B", b, wl, field, fields_wanted)


def _trace_run(
    compiled: Any,
    layouts: list[Any],
    run: str,
    b: RayBundle,
    wl: int,
    field: int,
    fields_wanted: bool,
    first: bool = False,
    keep_all: bool = False,
) -> list[Any]:
    """Traces a bundle along the events of the path `run` of `compiled`, starting where the rays
    are, and combines OPL, field and polarisation (the field a·prt·jones goes on as a scalar
    along a real direction plus the rest in prt). `fields_wanted`: stop before the last event
    (the detector) and return fields. `first`: the bundle comes from the source – the medium
    there is taken from the first segment. `keep_all` (rays): the bundle keeps every ray of `b`,
    the ones lost on the way marked dead – the image triangulates over the source grid and must
    not bridge the gaps (as along a whole path)."""
    import raytatouille as rt

    alive = np.flatnonzero(np.asarray(b.alive) & (np.abs(np.asarray(b.amplitude)) > 0))
    if alive.size < 3:
        return []
    rays = rt.trace.RayBatch(int(alive.size))
    pos = np.asarray(b.position)[alive]
    dirs = np.asarray(b.direction)[alive]
    for name, values in (
        ("pos_x", pos[:, 0]),
        ("pos_y", pos[:, 1]),
        ("pos_z", pos[:, 2]),
        ("dir_x", dirs[:, 0]),
        ("dir_y", dirs[:, 1]),
        ("dir_z", dirs[:, 2]),
    ):
        np.asarray(getattr(rays, name))[:] = values
    np.asarray(rays.field)[:] = field
    np.asarray(rays.wl)[:] = wl
    _, rec = rt.trace.trace(
        compiled,
        rays,
        path=run,
        record_path=True,
        max_recorded_rays=max(rays.size, 1),
        polarization=True,
    )
    reached = _passed(rec)  # also a ray that starts on or behind the run's first surface: missed
    rows = np.flatnonzero(reached)
    if rows.size < 3:
        return []
    idx = np.asarray(rec.ray_indices)
    src = alive[idx[rows]]  # index into the incoming bundle
    slot = rec.slots - 2 if fields_wanted else rec.slots - 1
    slot = max(slot, 0)
    direction = np.asarray(rec.direction)
    p_b = np.asarray(rec.prt)[rows, slot]
    k_b = direction[rows, slot]
    p_t = p_b - np.einsum("ni,nj,njk->nik", k_b, k_b, p_b)
    norm = 0.5 * np.sum(np.abs(p_t) ** 2, axis=(1, 2))
    weight = np.asarray(rec.weight)[rows, slot]
    s_b = (
        np.where(norm > 0, weight / np.where(norm > 0, norm, 1.0), 0.0)
        if slot > 0
        else np.ones(len(rows))
    )
    prt_in = (
        np.asarray(b.prt)[src]
        if b.prt is not None
        else np.tile(np.eye(3, dtype=complex), (len(rows), 1, 1))
    )
    frame = (
        _carried_frame(
            np.asarray(b.source_frame),
            direction,
            rows,
            np.asarray(b.source_xy)[alive],
            idx,
            slot,
            rec,
            layouts,
        )
        if b.source_frame is not None
        else None
    )
    # the full field a·P_B·prt·jones·sqrt(s_B), split into a scalar along a real polarisation
    # direction (the ray image R uses only that) and the rest in prt (fields H use both)
    prt_total = np.einsum("nij,njk->nik", p_b, prt_in)
    jones = (
        np.asarray(b.jones, dtype=complex)
        if b.jones is not None
        else np.array([1.0, 0.0, 0.0], dtype=complex)
    )
    scale = np.asarray(b.amplitude)[src] * np.sqrt(np.clip(s_b, 0.0, None))
    f = scale[:, None] * np.einsum("nij,j->ni", prt_total, jones)
    c = _scalar_along(f)
    safe = np.where(np.abs(c) > 1e-300, c, 1.0)
    bundle = b.evolve(
        position=np.asarray(rec.position)[rows, slot],
        direction=k_b,
        opl=np.asarray(b.opl)[src] + np.asarray(rec.opl)[rows, slot],
        amplitude=c,
        source_xy=np.asarray(b.source_xy)[src],
        alive=np.abs(c) > 0,
        prt=(prt_total * (scale / safe)[:, None, None]).astype(complex),
        source_frame=frame,
        index=_segment_index(rec, rows, slot, slot + 1) if fields_wanted else b.index,
        source_index=_segment_index(rec, rows, 0, 1) if first else b.source_index,
        source_cos=(np.asarray(b.source_cos)[src] if np.ndim(b.source_cos) > 0 else b.source_cos),
    )
    if not fields_wanted:
        if keep_all:
            bundle = _full(b, bundle, src)
        return [bundle]
    try:
        return [rays_to_field(bundle)]
    except FieldError as e:
        raise ExecutionError("exec.field", f"{run}: {e}") from e


def _full(b: RayBundle, part: RayBundle, src: Any) -> RayBundle:
    """`part` (the rays `src` of `b` that came through) on all rays of `b`, the others dead."""
    n = b.size

    def spread(values: Any, fill: Any) -> Any:
        full = np.array(fill, copy=True)
        full[src] = values
        return full

    prt_fill = (
        np.asarray(b.prt) if b.prt is not None else np.tile(np.eye(3, dtype=complex), (n, 1, 1))
    )
    return part.evolve(
        position=spread(part.position, b.position),
        direction=spread(part.direction, b.direction),
        opl=spread(part.opl, b.opl),
        amplitude=spread(part.amplitude, np.zeros(n, dtype=complex)),
        source_xy=np.asarray(b.source_xy),
        alive=spread(part.alive, np.zeros(n, dtype=bool)),
        prt=spread(part.prt, prt_fill.astype(complex)),
        source_cos=(
            spread(part.source_cos, b.source_cos)
            if np.ndim(part.source_cos) > 0
            else part.source_cos
        ),
    )


def _scalar_along(f: Any) -> Any:
    """Complex amplitude of field vectors f (N, 3) along a real unit direction (exact for linear
    polarisation): the phase of the largest component is taken out, the real part normalised."""
    big = np.argmax(np.abs(f), axis=1)
    alpha = np.angle(f[np.arange(len(f)), big])
    e = np.real(f * np.exp(-1j * alpha)[:, None])
    norm = np.linalg.norm(e, axis=1)
    e = e / np.where(norm > 0, norm, 1.0)[:, None]
    return np.sum(f * e, axis=1)


def _input_representation(ctx: Context, node: Node) -> str:
    """The planned representation of the link into the detector (R without a plan)."""
    links = ctx.graph.links_to(node.id)
    if ctx.plan is None or not links:
        return "R"
    return ctx.plan.representations.get(link_key(links[0]), "R")


def _object_index(compiled: Any, wl: int) -> float:
    """Refractive index of object space (the system's ambient medium, the first of
    `compiled.media`) at wavelength `wl`; 1 without media."""
    media = list(getattr(compiled, "media", []) or [])
    if not media:
        return 1.0
    index = media[0].index
    return float(np.real(index[min(int(wl), len(index) - 1)]))


def _start_opl(rec: Any, rows: Any, n: float) -> Any:
    """The optical path the start rays carry already (issue #69): Raytatouille starts every ray with
    OPL 0 on its start plane; for a field point at infinity the start rays are parallel and that
    plane is not their wavefront unless the field point is on the axis. The plane wave's own phase
    is n·d̂·(p − p_ref) with the start positions p, the common direction d̂ and the index n of object
    space (not estimated from the first segment: an ideal lens as first surface adds its phase
    there); without it ψ = k₀·OPL contradicts the rays' directions (∇ψ = k n d̂⊥) and a field (H)
    focuses wrongly. The reference p_ref = (0, 0, z̄₀) – the axis at the mean start height of all
    recorded rays – depends on the start rays only, not on which of them survive, so that every
    path (e.g. both arms of an interferometer) has the same phase reference; on the axis (d̂ = ẑ,
    start rays on one plane) the term vanishes. Rays that do not start parallel (a point source)
    keep 0."""
    start = np.asarray(rec.position)[:, 0]
    pos = start[rows]
    d = np.asarray(rec.direction)[rows, 0]
    if len(rows) == 0 or float(np.max(np.linalg.norm(d - d[0], axis=1))) > 1e-12:
        return np.zeros(len(rows))
    ref = np.array([0.0, 0.0, float(np.mean(start[:, 2]))])
    return np.asarray(n * ((pos - ref) @ d[0]))


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
            f = f @ rotation_between(a, b).T
    return f


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
    idx = np.asarray(rec.ray_indices)
    reached = _passed(rec)
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
        opl=np.asarray(rec.opl)[rows, port] + _start_opl(rec, rows, _object_index(compiled, wl)),
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


# ---------------------------------------------------------------------- node-wise (R25)


def _nodewise(ctx: Context) -> bool:
    """Node-wise execution: asked for (shared["rtt_execution"] = "nodes"), or needed because a
    branch of the model is no path of the system; "paths" keeps whole paths (R15)."""
    mode = str(ctx.shared.get("rtt_execution", "auto"))
    if mode in ("nodes", "paths"):
        return mode == "nodes"
    system = ctx.shared.get("system")
    names = ctx.shared.get("system_path_names")
    stamp = (id(system), _generated(ctx))
    if names is None or names[0] != stamp:
        names = (stamp, {str(p.get("name")) for p in _base_json(ctx).get("paths", [])})
        ctx.shared["system_path_names"] = names
    return any(b not in names[1] for b in ctx.graph.meta.get("paths", {}))


def _base_id(node_id: str) -> str:
    return node_id.split("@")[0]


def _channel(events: list[Any]) -> str:
    """The channel of a run: the channels of its events (`importer.event_channel`)."""
    from rtt_fusion.graph.importer import event_channel

    return " ".join(c for c in (event_channel(e) for e in events) if c)


def _runs(
    ctx: Context,
) -> tuple[Any, list[Any], dict[Any, str], dict[tuple[str, int], str], dict[tuple[str, str], str]]:
    """The system compiled with one extra path per run (the events of one block in one system
    path, and every traversal a block holds), cached per offsets. Returns compiled, layouts, the
    runs by (block, in socket, out socket, channel), by (system path, position in it) and by
    (block, "<port in>><port out>|<channel>")."""
    import raytatouille as rt

    cache = ctx.shared.setdefault("run_systems", {})
    key = _variation(ctx)
    if key in cache:
        return cache[key]  # type: ignore[no-any-return]
    d = _system_json(ctx)
    by_key: dict[Any, str] = {}
    by_place: dict[tuple[str, int], str] = {}
    by_port: dict[tuple[str, str], str] = {}
    extra = []
    meta_paths = ctx.graph.meta.get("paths", {})
    # the blocks owning each surface (paths the model has no branch for are grouped by them)
    owner = {
        str(sid): n.id
        for n in ctx.graph.nodes.values()
        if n.type_id.startswith("rtt.") and "@" not in n.id
        for sid in n.params.get("surfaces", [])
    }
    for p in d.get("paths", []):
        name = str(p.get("name"))
        if name not in meta_paths:
            if not isinstance(p.get("events"), list):
                p["events"] = _auto_events(d)
            runs_of_path: list[tuple[str, list[Any]]] = []
            for e in p["events"]:
                block = owner.get(str(e["surface"]), "")
                if runs_of_path and runs_of_path[-1][0] == block:
                    runs_of_path[-1][1].append(e)
                else:
                    runs_of_path.append((block, [e]))
            for i, (block, events) in enumerate(runs_of_path):
                if not block:
                    continue
                run = f"{name}#{i}"
                extra.append({"name": run, "events": events})
                came = "source" if i == 0 else _base_id(runs_of_path[i - 1][0])
                goes = (
                    f"to_{_base_id(runs_of_path[i + 1][0])}" if i + 1 < len(runs_of_path) else None
                )
                by_key.setdefault((block, f"from_{came}", goes, _channel(events)), run)
            continue
        blocks = meta_paths[name]
        if not isinstance(p.get("events"), list):  # "auto": every surface in tree order,
            p["events"] = [  # mirrors reflect (a sequential path has no splitter)
                {"surface": str(sid)}
                | (
                    {"kind": "reflect"}
                    if ctx.graph.nodes[b].params.get("element_type") == "mirror"
                    else {}
                )
                for b in blocks
                for sid in ctx.graph.nodes[b].params.get("surfaces", [])
            ]
        groups = _events_of_blocks(ctx, name, d)
        for i, (block, events) in enumerate(zip(blocks, groups, strict=True)):
            if not events or not ctx.graph.nodes[block].type_id.startswith("rtt."):
                continue
            run = f"{name}#{i}"
            extra.append({"name": run, "events": events})
            came = "source" if i == 0 else _base_id(blocks[i - 1])
            goes = f"to_{_base_id(blocks[i + 1])}" if i + 1 < len(blocks) else None
            by_key.setdefault((_base_id(block), f"from_{came}", goes, _channel(events)), run)
            by_place[(name, i)] = run
    # the traversals the blocks hold themselves (S2: the workflow is the truth): a branch
    # through a block passes by its data, whatever path of the system there is
    for n in ctx.graph.nodes.values():
        if not n.type_id.startswith("rtt.") or "@" in n.id:
            continue
        for key, events in (n.params.get("runs") or {}).items():
            if not events:
                continue
            came, _, rest = str(key).partition(">")
            goes, _, channel = rest.partition("|")
            run = f"{n.id} {key}"
            extra.append({"name": run, "events": events})
            by_port[(n.id, str(key))] = run  # by the ports of the branch's step
            by_key.setdefault((n.id, f"from_{came}", f"to_{goes}" if goes else None, channel), run)
    # ghosts (R27): a step "ghost:s1,s2,…" reflects at the named surfaces of a section
    events_of: dict[str, Any] = {str(e["name"]): e["events"] for e in extra}
    for b in ctx.graph.meta.get("branches", []):
        steps = [s for s in b["steps"] if s[1] is not None]
        blocks = meta_paths.get(b["name"], [])
        for i, step in enumerate(steps):
            channel = str(step[3]) if len(step) > 3 else ""
            if not channel.startswith("ghost:") or i >= len(blocks):
                continue
            run = f"{b['name']}#{i}"
            extra.append(
                {"name": run, "events": _ghost_events(by_key, events_of, step, channel, blocks[i])}
            )
            by_place[(b["name"], i)] = run
    d["paths"] = [*d.get("paths", []), *extra]
    compiled = rt.compile(
        rt.System.from_json(json.dumps(d)), ctx.shared.get("materials", rt.MaterialLibrary())
    )
    cache[key] = (compiled, list(rt.layout.surfaces(compiled)), by_key, by_place, by_port)
    return cache[key]  # type: ignore[no-any-return]


def _auto_events(d: dict[str, Any]) -> list[Any]:
    """The events of an "auto" path: every surface in tree order, mirrors reflecting."""
    out: list[Any] = []

    def visit(node: dict[str, Any]) -> None:
        for s in node.get("surfaces", []):
            if isinstance(s, dict):
                event = {"surface": str(s.get("id"))}
                if node.get("type") == "mirror":
                    event["kind"] = "reflect"
                out.append(event)
        for c in node.get("children", []):
            if isinstance(c, dict):
                visit(c)

    visit(d["root"])
    return out


def _ghost_events(
    by_key: dict[Any, str], events_of: dict[str, Any], step: list[Any], channel: str, block: str
) -> list[Any]:
    """Events of a ghost through a section: along the order of its surfaces as a path of the
    system passes it from the port the light comes in, reflected at each named surface in turn,
    the rest transmitted; it leaves through the port it came in after an odd number of
    reflections, through the other after an even one."""
    node, came = _base_id(block), str(step[1])
    forward = next(
        (events_of[r] for k, r in by_key.items() if k[0] == node and k[1] == f"from_{came}"),
        None,
    )
    if forward is None:  # only passed the other way: reverse that run
        back = next(
            (events_of[r] for k, r in by_key.items() if k[0] == node and k[2] == f"to_{came}"),
            None,
        )
        if back is None:
            raise ExecutionError("exec.branch", f"no path of the system passes {node}", node)
        forward = list(reversed(back))
    order = [str(e["surface"]) for e in forward]
    # a surface passed on the way keeps the kind of the forward traversal (refraction at a lens
    # surface, a straight pass through a thin element); "transmit" would not refract
    passing = {str(e["surface"]): {k: v for k, v in e.items() if k != "order"} for e in forward}
    events: list[Any] = []
    pos, direction = -1, 1
    for target in [s for s in channel[6:].split(",") if s]:
        k = pos + direction
        while 0 <= k < len(order) and order[k] != target:
            events.append(dict(passing[order[k]]))
            k += direction
        if not 0 <= k < len(order):
            raise ExecutionError(
                "exec.branch", f"ghost {channel!r}: {target} not ahead of the light in {node}", node
            )
        events.append({"surface": target, "kind": "reflect"})
        pos, direction = k, -direction
    k = pos + direction
    while 0 <= k < len(order):
        events.append(dict(passing[order[k]]))
        k += direction
    return events


def _branch_steps(ctx: Context, name: str) -> list[Any]:
    """The steps of a branch after its source ([node, port in, port out, channel?])."""
    for b in ctx.graph.meta.get("branches", []):
        if b.get("name") == name:
            return [s for s in b["steps"] if s[1] is not None]
    return []


def _run_of(ctx: Context, node: Node, b: RayBundle, out_socket: str | None) -> tuple[str, int]:
    """The run a bundle takes at this block: the next visit of the block in its branch."""
    assert b.branch is not None
    name, last = b.branch
    blocks = ctx.graph.meta["paths"][name]
    i = next((k for k in range(last + 1, len(blocks)) if blocks[k] == node.id), -1)
    if i < 0:
        raise ExecutionError("exec.branch", f"branch {name!r} does not pass {node.id}", node.id)
    goes = f"to_{_base_id(blocks[i + 1])}" if i + 1 < len(blocks) else None
    if goes != out_socket:
        return "", i  # light of this branch leaves through another socket
    return _run_at(ctx, name, i), i


def _run_at(ctx: Context, name: str, i: int) -> str:
    """The run of the i-th block of a branch (`graph.meta["paths"][name][i]`)."""
    blocks = ctx.graph.meta["paths"][name]
    _, _, by_key, by_place, by_port = _runs(ctx)
    if (name, i) in by_place:
        return by_place[(name, i)]
    steps_here = _branch_steps(ctx, name)
    if i < len(steps_here):  # the traversal the block holds for the step's ports (S2, S4)
        s = steps_here[i]
        key = f"{s[1]}>{s[2] or ''}|{s[3] if len(s) > 3 else ''}"
        if (_base_id(blocks[i]), key) in by_port:
            return by_port[(_base_id(blocks[i]), key)]
    came = "source" if i == 0 else _base_id(blocks[i - 1])
    goes = f"to_{_base_id(blocks[i + 1])}" if i + 1 < len(blocks) else None
    steps = _branch_steps(ctx, name)
    channel = str(steps[i][3]) if i < len(steps) and len(steps[i]) > 3 else ""
    run = by_key.get((_base_id(blocks[i]), f"from_{came}", goes, channel))
    if run is None:
        raise ExecutionError(
            "exec.branch",
            f"branch {name!r}: no path of the system passes {blocks[i]} from {came} to "
            f"{goes or 'the detector'}",
            blocks[i],
        )
    return run


def _medium_after(ctx: Context, b: RayBundle, name: str, i: int) -> float:
    """Refractive index behind the last surface of the i-th block of a branch: of the first
    segment of the next block's run (fields H on the link need it, as along a whole path)."""
    import raytatouille as rt

    blocks = ctx.graph.meta["paths"][name]
    if i + 1 >= len(blocks) or not ctx.graph.nodes[blocks[i + 1]].type_id.startswith("rtt."):
        return b.index
    compiled, _, _, _, _ = _runs(ctx)
    alive = np.flatnonzero(np.asarray(b.alive))
    if alive.size < 3:
        return b.index
    rays = rt.trace.RayBatch(int(alive.size))
    for k, axis in enumerate("xyz"):
        np.asarray(getattr(rays, f"pos_{axis}"))[:] = np.asarray(b.position)[alive, k]
        np.asarray(getattr(rays, f"dir_{axis}"))[:] = np.asarray(b.direction)[alive, k]
    np.asarray(rays.field)[:] = b.key.field_index
    np.asarray(rays.wl)[:] = b.key.wavelength_index
    _, rec = rt.trace.trace(
        compiled,
        rays,
        path=_run_at(ctx, name, i + 1),
        record_path=True,
        max_recorded_rays=max(rays.size, 1),
    )
    rows = np.flatnonzero(np.asarray(rec.count) >= 2)
    return _segment_index(rec, rows, 0, 1) if rows.size else b.index


def _place_foreign(ctx: Context) -> None:
    """A block that replaced a Raytatouille block (a grating as the sample) sits where the
    replaced surface was – as in R19b – unless its parameter `keep_pose` is set."""
    import raytatouille as rt

    d = _system_json(ctx)
    compiled = _compiled(ctx)
    by_id = {L.id: L for L in rt.layout.surfaces(compiled)}
    for name, blocks in ctx.graph.meta.get("paths", {}).items():
        foreign = [
            k
            for k, b in enumerate(blocks)
            if b in ctx.graph.nodes and not ctx.graph.nodes[b].type_id.startswith("rtt.")
        ]
        if not foreign or not any(p.get("name") == name for p in d.get("paths", [])):
            continue
        groups = _events_of_blocks(ctx, name, d)
        for k in foreign:
            n = ctx.graph.nodes[blocks[k]]
            L = by_id.get(str(groups[k][0]["surface"])) if groups[k] else None
            if L is None or n.params.get("keep_pose"):
                continue
            rot = np.asarray(L.rotation)
            n.params["point"] = [float(v) for v in np.asarray(L.translation)]
            n.params["normal"] = [float(v) for v in rot[:, 2]]


def _aim_path(ctx: Context, first_block: str) -> str | None:
    """A path of the system that starts in `first_block` (rays are aimed along it)."""
    surfaces = {str(s) for s in ctx.graph.nodes[first_block].params.get("surfaces", [])}
    paths = _system_json(ctx).get("paths", [])
    for p in paths:
        events = p.get("events")
        if isinstance(events, list) and events and str(events[0]["surface"]) in surfaces:
            return str(p.get("name"))
    return str(paths[0].get("name")) if len(paths) == 1 else None


def _source_nodes(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
    """Start bundles of every branch leaving through `socket`: two input polarisations
    (unpolarised source) per field point and wavelength, aimed along a system path. Places the
    replaced blocks first (they are evaluated after the source)."""
    import raytatouille as rt

    if not _nodewise(ctx):
        return []
    _place_foreign(ctx)
    compiled = _compiled(ctx)
    system_paths = [str(p.get("name")) for p in _system_json(ctx).get("paths", [])]
    muted = {n.id for n in ctx.graph.nodes.values() if n.muted}
    rings = int(ctx.shared.get("rings", 30))
    waist = float(ctx.shared.get("gaussian_waist", 0.0) or 0.0)
    out: list[RayBundle] = []
    for name, blocks in ctx.graph.meta.get("paths", {}).items():
        if not blocks or socket != f"to_{_base_id(blocks[0])}":
            continue
        if muted & set(blocks):
            ctx.log.append(f"{node.id}: branch {name} left out (muted block on it)")
            continue
        aim = name if name in system_paths else _aim_path(ctx, blocks[0])
        if aim is None:
            raise ExecutionError(
                "exec.aiming", f"branch {name!r}: no system path starts like it", node.id
            )
        for field in range(compiled.field_count):
            for wl in range(len(compiled.wavelengths_um)):
                rays = _start_rays(ctx, compiled, aim, rings, field, wl)
                pupil = np.stack([np.asarray(rays.pupil_x), np.asarray(rays.pupil_y)], -1).copy()
                apod = (
                    np.exp(-np.sum(pupil * pupil, axis=1) / waist**2)
                    if waist > 0
                    else np.ones(len(pupil))
                )
                position = np.stack(
                    [np.asarray(rays.pos_x), np.asarray(rays.pos_y), np.asarray(rays.pos_z)], -1
                ).copy()
                direction = np.stack(
                    [np.asarray(rays.dir_x), np.asarray(rays.dir_y), np.asarray(rays.dir_z)], -1
                ).copy()
                # the start rays' own optical path (#69), with the index of the first segment; the
                # trace moves the rays along, so it comes after their start is copied
                _, rec0 = rt.trace.trace(
                    compiled, rays, path=aim, record_path=True, max_recorded_rays=max(rays.size, 1)
                )
                start_opl = np.zeros(len(pupil))
                n_obj = _object_index(compiled, wl)
                start_opl[np.asarray(rec0.ray_indices)] = _start_opl(
                    rec0, np.arange(len(np.asarray(rec0.ray_indices))), n_obj
                )
                k0 = direction.mean(axis=0)
                k0 = k0 / np.linalg.norm(k0)
                seed = np.eye(3)[int(np.argmin(np.abs(k0)))]
                x = seed - k0 * float(k0 @ seed)
                x = x / np.linalg.norm(x)
                y = np.cross(k0, x)
                n = len(pupil)
                base = RayBundle(
                    key=CoherenceKey("source", field, wl, 0),
                    wavelength_um=float(compiled.wavelengths_um[wl]),
                    position=position,
                    direction=direction,
                    opl=start_opl,
                    index=n_obj,
                    source_index=n_obj,
                    amplitude=(np.sqrt(0.5) * apod).astype(complex),
                    source_xy=pupil,
                    alive=np.ones(n, dtype=bool),
                    path=(name,),
                    prt=np.tile(np.eye(3, dtype=complex), (n, 1, 1)),
                    source_frame=np.stack([x, y]),
                    branch=(name, -1),
                )
                for p, jones in enumerate((x, y)):
                    out.append(
                        base.evolve(
                            key=CoherenceKey("source", field, wl, p), jones=jones.astype(complex)
                        )
                    )
    return out


def _block_nodes(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
    """A Raytatouille block passes each arriving bundle along its run in the bundle's branch."""
    if not _nodewise(ctx):
        return []
    compiled, layouts, _, _, _ = _runs(ctx)
    out: list[Any] = []
    for bundles in inputs.values():
        for b in bundles or []:
            if not isinstance(b, RayBundle) or b.branch is None:
                raise ExecutionError(
                    "exec.branch", "node-wise execution needs rays that know their branch", node.id
                )
            run, i = _run_of(ctx, node, b, socket)
            if not run:
                continue
            fields_next = ctx.plan is not None and any(
                ctx.plan.representations.get(link_key(lk)) == "H"
                for lk in ctx.graph.links_from(node.id, socket)
            )
            for traced in _trace_run(
                compiled,
                layouts,
                run,
                b,
                b.key.wavelength_index,
                b.key.field_index,
                False,
                first=b.branch[1] < 0,
                keep_all=True,
            ):
                traced = traced.evolve(branch=(b.branch[0], i), path=(*b.path, node.id))
                if fields_next:
                    traced = traced.evolve(index=_medium_after(ctx, traced, b.branch[0], i))
                out.append(traced)
    return out


def _detector_nodes(ctx: Context, node: Node, inputs: dict[str, Any]) -> Any:
    """Node-wise detector: traces the arriving bundles onto its surface (R) or turns them into
    fields after the last surface and propagates them (H)."""
    compiled, layouts, _, _, _ = _runs(ctx)
    fields_wanted = _input_representation(ctx, node) == "H"
    bundles: list[RayBundle] = []
    fields: list[FieldBundle] = []
    for arriving in inputs.values():
        for b in arriving or []:
            if isinstance(b, FieldBundle):  # turned into a field on the link (plan: H)
                fields.append(b)
                continue
            if not isinstance(b, RayBundle) or b.branch is None:
                continue
            run, _ = _run_of(ctx, node, b, None)
            traced = _trace_run(
                compiled,
                layouts,
                run,
                b,
                b.key.wavelength_index,
                b.key.field_index,
                fields_wanted,
                first=b.branch[1] < 0,
                keep_all=True,
            )
            if fields_wanted:
                fields += traced
            else:
                bundles += traced
    return _image(ctx, node, compiled, bundles, fields, fields_wanted, 0.3)


EVALUATORS: dict[str, Evaluator] = {
    "rtt.source": _source_nodes,
    "rtt.section": _block_nodes,
    "rtt.splitter": _block_nodes,
    "rtt.mirror": _block_nodes,
    "rtt.detector": detector,
    # S4: new Raytatouille elements from the library
    "rtt.part.mirror": _block_nodes,
    "rtt.part.lens": _block_nodes,
    "rtt.part.stop": _block_nodes,
    "rtt.part.splitter": _block_nodes,
    "rtt.part.detector": detector,
}
# node-wise, light on an H link into the detector arrives as fields (R25)
detector.representations = frozenset("RH")  # type: ignore[attr-defined]
