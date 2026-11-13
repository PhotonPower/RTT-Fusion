"""Layer stack block `fg.stack` (ADR grating as field operator and layer stacks, stage F4).

A laterally invariant stack of homogeneous layers – coating, filter, mirror layers – between two
media, computed exactly by FourierGrating (truncation N = 0: only the order 0 exists). Sockets and
ports as `optics.interface`: `front` is the side the normal points to (the configuration's
superstrate), `back` the other (its substrate); light in through one side leaves reflected through
the same side and transmitted through the other. Light from the back sees the stack turned round
(layers in reverse order, superstrate and substrate swapped).

Representation R (this stage): per ray the power-normalised Jones matrix of order 0 for its own
direction and wavelength, straight from FourierGrating (no table, no interpolation; directions
that occur several times are computed once), applied with its PRT matrix
(`couplings.jones_table.prt_matrix`) as in `fg.grating`. Thin element: rays enter and leave in the
plane `point`/`normal`; amplitude and phase are exact (FourierGrating refers the reflected wave to
the entry plane and the transmitted one to the exit plane), the lateral shift inside the stack
(at most its thickness · tan θ, µm) is left out. Fields (H) reach the block as rays through the
executor's conversion if a link carries R.

Representation H (stage F4b): the field is propagated into the stack's plane, turned into the
stack's frame (x = the plane's u, so that anisotropic layers see their axes), sampled on a grid
with its carrier removed and decomposed into plane waves (`graph.spectral`); every significant,
propagating plane wave goes through FourierGrating (`apply_many`, exact per direction), and the
outgoing spectrum is put back on the field's own samples with its phase ψ (the stack keeps k∥).
Plane waves below `field_threshold` (relative amplitude) and evanescent ones are left out; their
share of the power is logged.

The stack: a FourierGrating configuration (`fouriergrating`, a JSON file) with homogeneous layers
only, listed as FourierGrating lists them – from the substrate (back) upwards to the superstrate
(front), `docs/config_format.md` of FourierGrating –, materials as FourierGrating takes them (n,
n + ik, tensors). Its incidence, truncation and period are set by the block.
"""

from __future__ import annotations

import copy
import json
from typing import Any

import numpy as np

from rtt_fusion.couplings import jones_table as jt
from rtt_fusion.graph import spectral
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.field import (
    FieldBundle,
    FieldError,
    field_to_rays,
    propagate_to,
    rays_to_field,
)
from rtt_fusion.graph.grating import split_prt
from rtt_fusion.graph.light import Plane, RayBundle
from rtt_fusion.graph.model import (
    Category,
    Node,
    NodeType,
    ParamSpec,
    PortFrame,
    light_in,
    light_out,
)
from rtt_fusion.graph.registry import NodeTypeRegistry

TYPE_ID = "fg.stack"


def stack_config(config: dict[str, Any], from_back: bool = False) -> dict[str, Any]:
    """The FourierGrating configuration of the stack as the block computes it: incidence, output
    and sweep removed, N = 0, 1D; turned round for light from the back. ValueError for a layer
    that is not homogeneous."""
    cfg = copy.deepcopy(config)
    for key in ("incidence", "output", "sweep", "nearfield"):
        cfg.pop(key, None)
    for k, layer in enumerate(cfg.get("layers", [])):
        if layer.get("type") != "homogeneous":
            raise ValueError(
                f"layer {k} is {layer.get('type')!r}: a stack has homogeneous layers only "
                "(a grating belongs into fg.grating)"
            )
    cfg["dimension"] = "1D"
    cfg["truncation"] = {"Nx": 0}
    cfg["period_um"] = {"x": 1.0}
    if from_back:
        cfg["superstrate"], cfg["substrate"] = cfg["substrate"], cfg["superstrate"]
        cfg["layers"] = list(reversed(cfg.get("layers", [])))
    return cfg


def _handle(ctx: Context, node: Node, wavelength_um: float, from_back: bool) -> Any:
    path = str(node.params.get("fouriergrating", ""))
    if not path:
        raise ExecutionError(
            "exec.param", "the stack needs a FourierGrating configuration", node.id
        )
    cache = ctx.shared.setdefault("fg_stack_handles", {})
    key = (path, round(wavelength_um, 12), from_back)
    if key not in cache:
        try:
            from rtt_fusion.adapters import fouriergrating as fga

            with open(path, encoding="utf-8") as f:
                cfg = stack_config(json.load(f), from_back)
            cache[key] = fga.prepare(fga.Structure(cfg), (wavelength_um,))
        except (OSError, ValueError, KeyError) as e:
            raise ExecutionError("exec.param", f"stack {path}: {e}", node.id) from e
        except ImportError as e:
            raise ExecutionError(
                "exec.solver", f"FourierGrating is not installed: {e}", node.id
            ) from e
    return cache[key]


def _isotropic(ctx: Context, node: Node) -> bool:
    """True if every material of the stack is given by an index (no tensor)."""
    path = str(node.params.get("fouriergrating", ""))
    try:
        with open(path, encoding="utf-8") as f:
            materials = json.load(f).get("materials", {})
    except (OSError, ValueError):
        return False
    return all(set(m) == {"n"} for m in materials.values())


def _plane(node: Node) -> Plane:
    point = np.asarray(node.params.get("point", [0.0, 0.0, 0.0]), dtype=float)
    normal = np.asarray(node.params.get("normal", [0.0, 0.0, 1.0]), dtype=float)
    return Plane.make(point, normal / np.linalg.norm(normal))


def _apply(
    ctx: Context, node: Node, b: RayBundle, plane: Plane, from_back: bool, reflect: bool
) -> RayBundle | None:
    handle = _handle(ctx, node, b.wavelength_um, from_back)
    moved, _ = plane.intersect(b)
    alive = np.asarray(moved.alive)
    normal = np.asarray(plane.normal)
    z = -normal if not from_back else normal  # local z along the incident light
    x = np.asarray(plane.u) - z * float(np.asarray(plane.u) @ z)
    x = x / np.linalg.norm(x)
    rot = np.stack([x, np.cross(z, x), z])  # global -> local: v @ rot.T
    d_loc = np.asarray(moved.direction) @ rot.T
    n_in = handle.front_index.real
    n_out = n_in if reflect else handle.back_index.real
    if abs(float(b.index) - n_in) > 1e-9:
        ctx.log.append(
            f"{node.id}: rays arrive in index {b.index}, the stack expects {n_in} on that side"
        )
    sigma = n_in * d_loc[:, :2]
    side = "reflected" if reflect else "transmitted"
    keys = sorted(
        {
            (round(float(sigma[k, 0]), 12), round(float(sigma[k, 1]), 12))
            for k in np.flatnonzero(alive)
        }
    )
    try:
        jones = jt.jones_many(handle, b.wavelength_um, keys, side)
    except (ValueError, RuntimeError) as e:
        raise ExecutionError("exec.solver", f"FourierGrating: {e}", node.id) from e
    p_of: dict[tuple[float, float], Any] = {}
    for key, d in zip(keys, jones, strict=True):
        if 0 not in d or not np.all(np.isfinite(d[0])):
            p_of[key] = None  # total internal reflection: no transmitted order
            continue
        p_loc = jt.prt_matrix(d[0], key, key, n_in, n_out, reflect)
        p_of[key] = rot.T @ p_loc @ rot
    ok = alive.copy()
    p_glob = np.zeros((moved.size, 3, 3), dtype=complex)
    for k in np.flatnonzero(alive):
        p = p_of[(round(float(sigma[k, 0]), 12), round(float(sigma[k, 1]), 12))]
        if p is None:
            ok[k] = False
        else:
            p_glob[k] = p
    if not np.any(ok):
        return None
    sign = -1.0 if reflect else 1.0
    rho = np.hypot(sigma[:, 0], sigma[:, 1])
    cz = np.sqrt(np.clip(n_out**2 - rho**2, 0.0, None))
    d_out = np.stack([sigma[:, 0], sigma[:, 1], sign * cz], axis=-1) / n_out
    e_in = moved.field_direction()
    prt_in = (
        np.asarray(moved.prt)
        if moved.prt is not None
        else np.tile(np.eye(3, dtype=complex), (moved.size, 1, 1))
    )
    c, prt_out = split_prt(p_glob, e_in, prt_in)
    return moved.evolve(
        direction=d_out @ rot,
        amplitude=np.asarray(moved.amplitude) * np.where(ok, c, 0.0),
        alive=ok & (np.abs(c) > 0),
        prt=prt_out,
        index=n_out,
        path=(*b.path, f"{node.id}[{'r' if reflect else 't'}]"),
        source_frame=moved.source_frame,
    )


def _in_frame(field: Any, q: Any) -> Any:
    """A field with its in-plane coordinates turned by the 2x2 rotation q (new = q·old); a gridded
    field only by a multiple of 90° (its grid stays axis-parallel), else FieldError."""
    from fieldtracing.field import GridlessField

    if np.allclose(q, np.eye(2), atol=1e-12):
        return field
    if not isinstance(field, GridlessField):
        raise FieldError(
            "a gridded field whose axes are not those of the stack (anisotropic layers need "
            "the stack's frame): turn the stack's plane axes onto the field's"
        )
    u = np.asarray(field.amplitude)
    u = np.concatenate([u[:, :2] @ q.T, u[:, 2:]], axis=1)
    return GridlessField(
        np.asarray(field.points) @ q.T,
        np.asarray(field.phase),
        np.asarray(field.gradient) @ q.T,
        np.einsum("ij,njk,lk->nil", q, np.asarray(field.hessian), q),
        np.asarray(u, dtype=np.complex128),
        float(field.wavelength_um),
        float(field.refractive_index),
    )


def _to_plane(
    ctx: Context, node: Node, fb: FieldBundle, plane: Plane, from_back: bool
) -> FieldBundle:
    """The field in the stack's plane: geometrically (field -> rays -> field on the plane, as the
    ray path, keeping ψ) unless the rays have a caustic there – a focus on the stack –, then by
    Fieldtracing's wave propagation (`field.propagate_to`, possibly a gridded field)."""
    from fieldtracing.field import GridlessField

    light = np.asarray(plane.normal) * (1.0 if from_back else -1.0)
    ahead = float((np.asarray(plane.origin) - np.asarray(fb.plane.origin)) @ light)
    if float(np.asarray(fb.plane.normal) @ light) <= 0 or ahead < -1e-12:
        raise ExecutionError(
            "exec.field", "the field does not run into the stack (wrong side or direction)", node.id
        )
    if isinstance(fb.field, GridlessField):
        on = Plane.make(plane.origin, light, plane.v)
        try:
            return rays_to_field(field_to_rays(fb), on)
        except FieldError as e:
            ctx.log.append(f"{node.id}: {e} – wave propagation onto the stack instead")
    target, notes = propagate_to(fb, Plane.make(plane.origin, light, plane.v))
    ctx.log.extend(f"{node.id}: {n}" for n in notes)
    return target


def _apply_field(
    ctx: Context, node: Node, fb: FieldBundle, plane: Plane, from_back: bool, reflect: bool
) -> FieldBundle:
    from fieldtracing.field import GridlessField

    from rtt_fusion.adapters import fouriergrating as fga
    from rtt_fusion.ports import PlaneWaves

    target = _to_plane(ctx, node, fb, plane, from_back)
    z = np.asarray(target.plane.normal, dtype=float)
    light = np.asarray(plane.normal) * (1.0 if from_back else -1.0)
    if float(z @ light) <= 0:
        raise ExecutionError("exec.field", "the field does not run into the stack", node.id)
    x = np.asarray(plane.u) - z * float(np.asarray(plane.u) @ z)
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    q = np.array(
        [[x @ target.plane.u, x @ target.plane.v], [y @ target.plane.u, y @ target.plane.v]]
    )
    if not isinstance(target.field, GridlessField) and _isotropic(ctx, node):
        q = np.eye(2)  # isotropic layers: the frame's orientation in the plane does not matter
        x, y = np.asarray(target.plane.u), np.asarray(target.plane.v)
    field = _in_frame(target.field, q)
    lam = float(field.wavelength_um)
    handle = _handle(ctx, node, lam, from_back)
    n_in = handle.front_index.real
    n_out = n_in if reflect else handle.back_index.real
    if abs(float(field.refractive_index) - n_in) > 1e-9:
        ctx.log.append(
            f"{node.id}: field in index {field.refractive_index}, the stack expects {n_in}"
        )
    try:
        sampled = spectral.sample(field, int(node.params.get("field_grid_max", 1024)))
    except spectral.SpectralError as e:
        raise ExecutionError("exec.field", str(e), node.id) from e
    grid = sampled.grid
    spec = spectral.spectrum(sampled.residual)
    kxx, kyy = grid.frequencies()
    kx, ky = kxx + grid.carrier[0], kyy + grid.carrier[1]
    port = fga.front_port(handle, lam)
    propagating = np.hypot(kx, ky) < n_in * port.k0 * (1.0 - 1e-12)
    keep = propagating & spectral.significant(
        spec[..., :2], float(node.params.get("field_threshold", 1e-10))
    )
    total = float(np.sum(np.abs(spec[..., :2]) ** 2)) or 1.0
    dropped = float(np.sum(np.abs(spec[..., :2][~keep]) ** 2)) / total
    if dropped > 1e-12:
        ctx.log.append(
            f"{node.id}: {dropped:.2e} of the field's spectrum left out (evanescent or below "
            "the threshold)"
        )
    rows, cols = np.nonzero(keep)
    # reflection and transmission come from one solve per field (the two outputs are asked for
    # separately): the results are kept for the run
    store = ctx.shared.setdefault("fg_stack_field_results", {})
    run_key = (node.id, id(fb), from_back)
    if run_key not in store:
        c = grid.carrier
        units = [
            PlaneWaves(port, np.array([c[0]]), np.array([c[1]]), np.array([u], dtype=complex))
            for u in ((1.0, 0.0), (0.0, 1.0))
        ]
        incident = PlaneWaves(port, kx[rows, cols], ky[rows, cols], spec[rows, cols, :2])
        try:
            store[run_key] = (fb, fga.apply_many(handle, [incident, *units]))
        except (ValueError, RuntimeError) as e:
            raise ExecutionError("exec.solver", f"FourierGrating: {e}", node.id) from e
    results = store[run_key][1]
    side = "front" if reflect else "back"
    waves = results[0][side]
    # the carrier wave's transfer (3 x 2): columns for unit x and y input
    t0 = np.zeros((3, 2), dtype=complex)
    for j, r in enumerate(results[1:]):
        w = r[side]
        if len(w.kx):
            t0[:2, j] = np.asarray(w.field)[0]
            t0[2, j] = w.ez()[0]
    out_spec = np.zeros_like(spec)
    dk_x = 2.0 * np.pi / (grid.nx * grid.step)
    dk_y = 2.0 * np.pi / (grid.ny * grid.step)
    ix = np.rint((np.asarray(waves.kx) - grid.carrier[0]) / dk_x).astype(int) % grid.nx
    iy = np.rint((np.asarray(waves.ky) - grid.carrier[1]) / dk_y).astype(int) % grid.ny
    out_spec[iy, ix, :2] = np.asarray(waves.field)
    out_spec[iy, ix, 2] = waves.ez()
    new = sampled.back(spectral.inverse(out_spec), n_out, reflect, t0)
    origin = np.asarray(target.plane.origin)
    out_plane = Plane(origin, -z, x, -y) if reflect else Plane(origin, z, x, y)
    return FieldBundle(fb.key, out_plane, new, (*fb.path, f"{node.id}[{'r' if reflect else 't'}]"))


def stack(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
    """Light in through one side leaves reflected through it and transmitted through the other;
    rays as rays (R), fields as fields (H)."""
    plane = _plane(node)
    out_side = socket.split("_")[0]
    out: list[Any] = []
    for in_socket, bundles in inputs.items():
        in_side = in_socket.split("_")[0]
        for b in bundles or []:
            reflect = in_side == out_side
            if isinstance(b, FieldBundle):
                try:
                    out.append(_apply_field(ctx, node, b, plane, in_side == "back", reflect))
                except FieldError as e:
                    raise ExecutionError("exec.field", str(e), node.id) from e
                continue
            if not isinstance(b, RayBundle):
                continue
            result = _apply(ctx, node, b, plane, in_side == "back", reflect)
            if result is not None:
                out.append(result)
    return out


stack.representations = frozenset("RH")  # type: ignore[attr-defined]


def _frames(node: Node) -> dict[str, PortFrame]:
    plane = _plane(node)
    o = tuple(float(v) for v in plane.origin)
    n = tuple(float(v) for v in plane.normal)
    return {
        "front": PortFrame(o, n),  # type: ignore[arg-type]
        "back": PortFrame(o, tuple(-v for v in n)),  # type: ignore[arg-type]
    }


TYPE = NodeType(
    TYPE_ID,
    Category.DOMAIN,
    (
        light_in("front_in", "RH", "front"),
        light_out("front_out", "RH", "front"),
        light_in("back_in", "RH", "back"),
        light_out("back_out", "RH", "back"),
    ),
    (
        ParamSpec(
            "fouriergrating",
            "file",
            "",
            description="FourierGrating configuration of the stack: homogeneous layers, listed "
            "from the substrate (back) upwards to the superstrate (front)",
        ),
        ParamSpec("point", "point", [0.0, 0.0, 0.0], unit="mm"),
        ParamSpec(
            "normal", "point", [0.0, 0.0, 1.0], description="towards the front (superstrate)"
        ),
        ParamSpec(
            "field_threshold",
            "number",
            1e-10,
            group="Field path",
            description="plane waves below this share of the largest amplitude are left out (H)",
        ),
        ParamSpec(
            "field_grid_max",
            "integer",
            1024,
            group="Field path",
            description="at most this many grid points per axis for a field (H)",
        ),
    ),
    (
        ("front_out", ("front_in", "back_in")),
        ("back_out", ("front_in", "back_in")),
    ),
    scatter=tuple((i, o, "full", ()) for i in ("front", "back") for o in ("front", "back")),
    port_frames=_frames,
    description="Layer stack (FourierGrating, exact per ray): coating, filter, mirror layers",
)

EVALUATORS: dict[str, Evaluator] = {TYPE_ID: stack}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    registry.add(TYPE)
    return dict(EVALUATORS)
