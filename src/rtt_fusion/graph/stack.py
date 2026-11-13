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
executor's conversion; the spectral field path follows in stage F4b.

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
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
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


def stack(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[RayBundle]:
    """Light in through one side leaves reflected through it and transmitted through the other."""
    plane = _plane(node)
    out_side = socket.split("_")[0]
    out: list[RayBundle] = []
    for in_socket, bundles in inputs.items():
        in_side = in_socket.split("_")[0]
        for b in bundles or []:
            if not isinstance(b, RayBundle):
                continue
            result = _apply(ctx, node, b, plane, in_side == "back", in_side == out_side)
            if result is not None:
                out.append(result)
    return out


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
