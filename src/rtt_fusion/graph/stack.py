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

Representation S (stage R30c): a plane-wave spectrum is moved onto the stack's plane in the
spectrum (`spectrum.propagate`, exact; resampled through H only from a tilted plane), and every
wave goes through FourierGrating on its own (`couplings.grating_field.spectrum_through`): the
stack is diagonal in k∥, so reflected and transmitted light leave as spectra on the same k∥ with
the sampling kept.

The stack: a FourierGrating configuration (`fouriergrating`, a JSON file) with homogeneous layers
only, listed as FourierGrating lists them – from the substrate (back) upwards to the superstrate
(front), `docs/config_format.md` of FourierGrating –, materials as FourierGrating takes them (n,
n + ik, tensors). Its incidence, truncation and period are set by the block.
"""

from __future__ import annotations

import copy
import json
from dataclasses import replace
from typing import Any

import numpy as np

from rtt_fusion.couplings import grating_field as gf
from rtt_fusion.couplings import jones_table as jt
from rtt_fusion.graph import embedded, spectral
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.field import FieldBundle, FieldError, field_tolerance
from rtt_fusion.graph.grating import CONFIG, spectrum_onto, split_prt
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
from rtt_fusion.graph.spectrum import SpectrumBundle

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


def _config(node: Node) -> dict[str, Any]:
    """The stack's FourierGrating configuration: held (`config`) or read from its file."""
    cfg = embedded.content(node, *CONFIG)
    if cfg is None:
        raise ExecutionError(
            "exec.param", "the stack needs a FourierGrating configuration", node.id
        )
    return cfg


def _handle(ctx: Context, node: Node, wavelength_um: float, from_back: bool) -> Any:
    raw = _config(node)
    cache = ctx.shared.setdefault("fg_stack_handles", {})
    key = (json.dumps(raw, sort_keys=True), round(wavelength_um, 12), from_back)
    if key not in cache:
        try:
            from rtt_fusion.adapters import fouriergrating as fga

            cfg = stack_config(raw, from_back)
            cache[key] = fga.prepare(fga.Structure(cfg), (wavelength_um,))
        except (ValueError, KeyError) as e:
            raise ExecutionError("exec.param", f"stack configuration: {e}", node.id) from e
        except ImportError as e:
            raise ExecutionError(
                "exec.solver", f"FourierGrating is not installed: {e}", node.id
            ) from e
    return cache[key]


def _isotropic(ctx: Context, node: Node) -> bool:
    """True if every material of the stack is given by an index (no tensor)."""
    try:
        materials = _config(node).get("materials", {})
    except ExecutionError:
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


def _to_plane(
    ctx: Context, node: Node, fb: FieldBundle, plane: Plane, from_back: bool
) -> FieldBundle:
    """The field in the stack's plane (`couplings.grating_field.onto_plane`)."""
    log: list[str] = []
    try:
        target = gf.onto_plane(fb, plane, from_back, log)
    except FieldError as e:
        raise ExecutionError("exec.field", str(e), node.id) from e
    finally:
        ctx.log.extend(f"{node.id}: {line}" for line in log)
    return target


def _apply_field(
    ctx: Context, node: Node, fb: FieldBundle, plane: Plane, from_back: bool, reflect: bool
) -> FieldBundle:
    from fieldtracing.field import GridlessField

    target = _to_plane(ctx, node, fb, plane, from_back)
    keep_axes = not isinstance(target.field, GridlessField) and _isotropic(ctx, node)
    fr = gf.frame(target, plane.u, keep_axes)
    try:
        field = gf.in_frame(target.field, fr.q)
    except FieldError as e:
        raise ExecutionError("exec.field", str(e), node.id) from e
    lam = float(field.wavelength_um)
    handle = _handle(ctx, node, lam, from_back)
    n_in = handle.front_index.real
    if abs(float(field.refractive_index) - n_in) > 1e-9:
        ctx.log.append(
            f"{node.id}: field in index {field.refractive_index}, the stack expects {n_in}"
        )
    # reflection and transmission come from one solve per field (the two outputs are asked for
    # separately): the results are kept for the run
    store = ctx.shared.setdefault("fg_stack_field_results", {})
    run_key = (node.id, id(fb), from_back)
    if run_key not in store:
        log_lines: list[str] = []
        try:
            applied = gf.apply(
                handle,
                field,
                threshold=float(node.params.get("field_threshold", 1e-10)),
                grid_max=int(node.params.get("field_grid_max", 1024)),
                transfer=str(node.params.get("field_transfer", "exact")),
                log=log_lines,
                cache=ctx.shared.setdefault("fg_stack_transfers", {}).setdefault(
                    handle.content_hash, {}
                ),
            )
        except spectral.SpectralError as e:
            raise ExecutionError("exec.field", str(e), node.id) from e
        except (ValueError, RuntimeError) as e:
            raise ExecutionError("exec.solver", f"FourierGrating: {e}", node.id) from e
        ctx.log.extend(f"{node.id}: {line}" for line in log_lines)
        if applied.dropped > 1e-12:
            ctx.log.append(
                f"{node.id}: {applied.dropped:.2e} of the field's spectrum left out (evanescent "
                "or below the threshold)"
            )
        store[run_key] = (fb, applied)
    applied = store[run_key][1]
    new = applied.fields.get(("front" if reflect else "back", (0, 0)))
    if new is None:
        raise ExecutionError("exec.field", "no light leaves on this side", node.id)
    origin = np.asarray(target.plane.origin)
    out_plane = Plane(origin, -fr.z, fr.x, -fr.y) if reflect else Plane(origin, fr.z, fr.x, fr.y)
    label = f"{node.id}[{'r' if reflect else 't'}]"
    return FieldBundle(fb.key, out_plane, new, (*fb.path, label), fb.branch)


def _apply_spectrum(
    ctx: Context, node: Node, sb: SpectrumBundle, plane: Plane, from_back: bool, reflect: bool
) -> SpectrumBundle:
    """A plane-wave spectrum through the stack (stage R30c): moved onto the stack's plane in the
    spectrum (S4), every wave through FourierGrating on its own (diagonal in k∥), in the frame x =
    the plane's u (anisotropic layers see their axes); the spectrum leaves with the incident
    axes, turned to (u, −v, −n) when reflected (`couplings.grating_field.spectrum_through`)."""
    light = np.asarray(plane.normal) * (1.0 if from_back else -1.0)
    if float(np.asarray(sb.plane.normal) @ light) <= 0.0:
        raise ExecutionError(
            "exec.field", "the spectrum does not run into the block (wrong side)", node.id
        )
    log: list[str] = []
    try:
        moved = spectrum_onto(sb, plane.origin, light, log, field_tolerance(ctx.shared))
        handle = _handle(ctx, node, sb.wavelength_um, from_back)
        side = "front" if reflect else "back"
        by_order = gf.spectrum_through(
            handle,
            moved,
            plane.u,
            side,
            cache=ctx.shared.setdefault("fg_stack_transfers", {}).setdefault(
                handle.content_hash, {}
            ),
            log=log,
            transfer=str(node.params.get("field_transfer", "exact")),
        )
    except FieldError as e:
        raise ExecutionError("exec.field", str(e), node.id) from e
    except (ValueError, RuntimeError) as e:
        raise ExecutionError("exec.solver", f"FourierGrating: {e}", node.id) from e
    finally:
        ctx.log.extend(f"{node.id}: {line}" for line in log)
    out: SpectrumBundle | None = by_order.get((0, 0))
    if out is None:
        raise ExecutionError("exec.field", "no light leaves on this side", node.id)
    return replace(out, path=(*sb.path, f"{node.id}[{'r' if reflect else 't'}]"))


def stack(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
    """Light in through one side leaves reflected through it and transmitted through the other;
    rays as rays (R), fields as fields (H), spectra as spectra (S)."""
    plane = _plane(node)
    out_side = socket.split("_")[0]
    out: list[Any] = []
    for in_socket, bundles in inputs.items():
        in_side = in_socket.split("_")[0]
        for b in bundles or []:
            reflect = in_side == out_side
            if isinstance(b, SpectrumBundle):
                out.append(_apply_spectrum(ctx, node, b, plane, in_side == "back", reflect))
                continue
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


stack.representations = frozenset("RHS")  # type: ignore[attr-defined]


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
        light_in("front_in", "RHS", "front"),
        light_out("front_out", "RHS", "front"),
        light_in("back_in", "RHS", "back"),
        light_out("back_out", "RHS", "back"),
    ),
    (
        ParamSpec(
            "fouriergrating",
            "file",
            "",
            description="FourierGrating configuration of the stack: homogeneous layers, listed "
            "from the substrate (back) upwards to the superstrate (front); Check workflow stores "
            "it in the block (config) and clears the path",
        ),
        ParamSpec(
            "config",
            "data",
            {},
            description="the FourierGrating configuration held by the block (before the file)",
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
            "field_transfer",
            "choice",
            "exact",
            options=("exact", "interpolated"),
            group="Field path",
            description="exact: FourierGrating per plane wave; interpolated: from a verified "
            "sub-grid of the spectrum (1e-8), much faster for wide beams",
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
