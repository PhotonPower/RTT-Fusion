"""Rigorous grating block in the working graph (stage R19, concept 4: sample -> FourierGrating).

`fusion.grating` takes the place of a mirror or a plate (same sockets "in" and "out", so that
`Graph.replace_node` keeps the links) and sends every requested diffraction order on: a
splitter whose channels are the orders (concept 3.1, fusion plan 7.1). The interaction is the
tabulated one of stage R6 (`couplings.jones_table`): power-normalised Jones matrices J_m over
(λ, σ) from FourierGrating – either read from a stored table (`raytatouille-result` file, runs
without the solver) or, with `fouriergrating` (a FourierGrating configuration file, as in
`docs/konventionen.md`), built by the block itself around the incident rays: λ ± 1 nm and the σ
range of the rays plus a margin, 4 nodes per axis at least (cubic interpolation), cached per
configuration and grid. Building needs the solver in the environment.

Per ray and order m (local frame: x = grating vector, y = grooves, z = surface normal towards the
light's side of exit for "transmitted", towards the incident side for "reflected"):

- direction from the grating equation σ_m = σ + m·λ/Λ·x̂ (σ = n·k̂⊥, direction cosines),
- field E_m = P_m·E with the order's PRT matrix (`jones_table.prt_matrix`, P_m k_in = k_m),
- phase: the table refers to x = 0 (ADR 0025 point 1); a hit at x adds the phase layer's
  m·2π·x/Λ, i.e. OPL + m·λ·x/Λ,
- all orders keep the coherence key of the incident bundle (they interfere where they overlap).

Rays whose order is evaporating or whose (λ, σ) lies outside the table or near a Rayleigh anomaly
(`jones_table.near_rayleigh_anomaly`) are dropped with a log line.
"""

from __future__ import annotations

import json
from typing import Any

import numpy as np

from rtt_fusion._types import FloatArray
from rtt_fusion.couplings import jones_table as jt
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.light import Plane, RayBundle, rotation_between
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
from rtt_fusion.results import load_json


def load_table(path: str) -> jt.JonesTable:
    """A Jones table from a `raytatouille-result` file (parameters: orders, side, period_um,
    n_in, n_out; inputs: wavelengths_um, sigma_x, sigma_y, jones), as written in stage R7."""
    with open(path, encoding="utf-8") as f:
        _, d = load_json(f.read())
    p, i = d["parameters"], d["inputs"]
    return jt.JonesTable(
        tuple(int(m) for m in p["orders"]),
        str(p["side"]),
        np.asarray(i["wavelengths_um"], dtype=float),
        np.asarray(i["sigma_x"], dtype=float),
        np.asarray(i["sigma_y"], dtype=float),
        float(p["period_um"]),
        float(p["n_in"]),
        float(p["n_out"]),
        np.asarray(i["jones"], dtype=complex),
    )


def _grid(values: FloatArray, margin: float, step: float) -> FloatArray:
    """At least 4 equidistant nodes around the values (cubic stencils inside)."""
    lo, hi = float(np.min(values)) - margin, float(np.max(values)) + margin
    n = max(4, int(np.ceil((hi - lo) / step)) + 1)
    return np.asarray(np.linspace(lo, hi, n), dtype=np.float64)


def _built_table(
    ctx: Context, node: Node, wavelength_um: float, sigma: FloatArray
) -> jt.JonesTable:
    """A Jones table from FourierGrating around the incident directions (cached)."""
    path = str(node.params.get("fouriergrating", ""))
    try:
        from rtt_fusion.adapters import fouriergrating as fga
    except ImportError as e:  # pragma: no cover - the adapter itself imports lazily
        raise ExecutionError("exec.solver", f"FourierGrating adapter: {e}", node.id) from e
    try:
        with open(path, encoding="utf-8") as f:
            config = json.load(f)
    except (OSError, ValueError) as e:
        raise ExecutionError("exec.param", f"cannot read {path}: {e}", node.id) from e
    for key in ("incidence", "wavelength_um", "output"):
        config.pop(key, None)  # the block sets the incidence itself
    side = str(node.params.get("side", "reflected"))
    orders = tuple(int(m) for m in node.params.get("orders", []) or (-1, 0, 1))
    margin = float(node.params.get("sigma_margin", 0.01))
    sx = _grid(sigma[:, 0], margin, 0.01)
    sy = _grid(sigma[:, 1], margin, 0.01)
    lam = np.linspace(wavelength_um - 1.5e-3, wavelength_um + 1.5e-3, 4)
    key = json.dumps([config, side, orders, lam.tolist(), sx.tolist(), sy.tolist()], sort_keys=True)
    cache = ctx.shared.setdefault("built_jones_tables", {})
    if key not in cache:
        try:
            handle = fga.prepare(fga.Structure(config), tuple(float(v) for v in lam))
            cache[key] = jt.build_table(handle, lam, sx, sy, orders, side)
        except ImportError as e:
            raise ExecutionError(
                "exec.solver", f"FourierGrating is not installed: {e}", node.id
            ) from e
        except (ValueError, RuntimeError, KeyError) as e:
            raise ExecutionError("exec.solver", f"FourierGrating: {e}", node.id) from e
        ctx.log.append(
            f"{node.id}: Jones table from FourierGrating, {len(lam)}×{len(sx)}×{len(sy)} nodes, "
            f"orders {list(orders)}, {side}"
        )
    table: jt.JonesTable = cache[key]
    return table


def _table(ctx: Context, node: Node) -> jt.JonesTable:
    path = str(node.params.get("table", ""))
    if not path:
        raise ExecutionError(
            "exec.param", "the grating needs a table or a FourierGrating configuration", node.id
        )
    cache = ctx.shared.setdefault("jones_tables", {})
    if path not in cache:
        try:
            cache[path] = load_table(path)
        except (OSError, KeyError, ValueError) as e:
            raise ExecutionError("exec.param", f"cannot read the table {path}: {e}", node.id) from e
    table: jt.JonesTable = cache[path]
    return table


def _inside(axis: Any, x: float) -> bool:
    a = np.asarray(axis)
    return bool(a.min() - 1e-12 <= x <= a.max() + 1e-12)


def grating(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[RayBundle]:
    built = bool(node.params.get("fouriergrating")) and not node.params.get("table")
    table = None if built else _table(ctx, node)
    point = np.asarray(node.params.get("point", [0.0, 0.0, 0.0]), dtype=float)
    normal = np.asarray(node.params.get("normal", [0.0, 0.0, 1.0]), dtype=float)
    normal = normal / np.linalg.norm(normal)
    gvec = np.asarray(node.params.get("grating_vector", [1.0, 0.0, 0.0]), dtype=float)
    out: list[RayBundle] = []
    for b in inputs.get("in", []):
        if not isinstance(b, RayBundle):
            continue
        d_mean = np.asarray(b.direction)[np.asarray(b.alive)].mean(axis=0)
        z = normal if float(d_mean @ normal) > 0 else -normal  # local z along the incident light
        x = gvec - z * float(gvec @ z)
        x = x / np.linalg.norm(x)
        y = np.cross(z, x)
        rot = np.stack([x, y, z])  # global -> local: v @ rot.T
        plane = Plane(point, z, x, y)
        moved, uv = plane.intersect(b)
        alive = np.asarray(moved.alive)
        d_loc = np.asarray(moved.direction) @ rot.T
        lam = b.wavelength_um
        if built:
            n_in0 = float(node.params.get("n_in", 1.0))
            table = _built_table(ctx, node, lam, n_in0 * d_loc[alive][:, :2])
        assert table is not None
        reflected = table.side == "reflected"
        orders = [int(m) for m in node.params.get("orders", []) or table.orders]
        unknown = [m for m in orders if m not in table.orders]
        if unknown:
            raise ExecutionError("exec.param", f"orders {unknown} not in the table", node.id)
        sigma = table.n_in * d_loc[:, :2]
        if not _inside(table.wavelengths_um, lam):
            raise ExecutionError(
                "exec.grating_range", f"wavelength {lam} µm outside the table", node.id
            )
        e_in = moved.field_direction()
        prt_in = (
            np.asarray(moved.prt)
            if moved.prt is not None
            else np.tile(np.eye(3, dtype=complex), (moved.size, 1, 1))
        )
        for m in orders:
            sig_m = sigma + np.array([m * lam / table.period_um, 0.0])
            rho = np.hypot(sig_m[:, 0], sig_m[:, 1])
            ok = alive & (rho < table.n_out)
            cz = np.sqrt(np.clip(table.n_out**2 - rho**2, 0.0, None))
            sign = -1.0 if reflected else 1.0
            d_out = np.stack([sig_m[:, 0], sig_m[:, 1], sign * cz], axis=-1) / table.n_out
            p_glob = np.zeros((moved.size, 3, 3), dtype=complex)
            cache: dict[tuple[float, float], Any] = {}
            dropped = 0
            for k in np.flatnonzero(ok):
                key = (round(float(sigma[k, 0]), 12), round(float(sigma[k, 1]), 12))
                if key not in cache:
                    if not (
                        _inside(table.sigma_x, key[0]) and _inside(table.sigma_y, key[1])
                    ) or jt.near_rayleigh_anomaly(table, lam, *key):
                        cache[key] = None
                    else:
                        jm = jt.interpolate(table, m, lam, *key)
                        if not np.all(np.isfinite(jm)):
                            cache[key] = None
                        else:
                            p_loc = jt.prt_matrix(
                                jm,
                                key,
                                (key[0] + m * lam / table.period_um, key[1]),
                                table.n_in,
                                table.n_out,
                                reflected,
                            )
                            cache[key] = rot.T @ p_loc @ rot
                p = cache[key]
                if p is None:
                    ok[k] = False
                    dropped += 1
                    continue
                p_glob[k] = p
            if dropped:
                ctx.log.append(
                    f"{node.id}: order {m}: {dropped} rays outside the table or near a Rayleigh "
                    "anomaly"
                )
            if not np.any(ok):
                continue
            f_out = np.einsum("nij,nj->ni", p_glob, e_in)
            # scalar amplitude along a real polarisation direction (for the ray image R), the rest
            # in prt so that a·prt·jones is the full field E_m (for fields H)
            big = np.argmax(np.abs(f_out), axis=1)
            alpha = np.angle(f_out[np.arange(moved.size), big])
            e_ref = np.real(f_out * np.exp(-1j * alpha)[:, None])
            norm = np.linalg.norm(e_ref, axis=1)
            e_ref = e_ref / np.where(norm > 0, norm, 1.0)[:, None]
            c = np.sum(f_out * e_ref, axis=1)
            safe = np.where(np.abs(c) > 1e-300, c, 1.0)
            prt_out = np.einsum("nij,njk->nik", p_glob, prt_in) / safe[:, None, None]
            opl = np.asarray(moved.opl) + m * (lam * 1e-3) * uv[:, 0] / (table.period_um * 1e-3)
            out.append(
                moved.evolve(
                    direction=d_out @ rot,
                    opl=opl,
                    amplitude=np.asarray(moved.amplitude) * np.where(ok, c, 0.0),
                    alive=ok & (np.abs(c) > 0),
                    prt=prt_out,
                    path=(*b.path, f"{node.id}[{m}]"),
                    source_frame=_carry(moved.source_frame, d_mean, (d_out @ rot)[ok].mean(axis=0)),
                )
            )
    return out


def _carry(frame: Any, a: Any, b: Any) -> Any:
    """Source axes carried from direction a to b: reflection if the order runs back, else the
    smallest rotation (as for refraction)."""
    if frame is None:
        return None
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    f = np.asarray(frame, dtype=float)
    if float(a @ b) < 0:
        n = (b - a) / np.linalg.norm(b - a)
        return f - 2.0 * np.outer(f @ n, n)
    return f @ rotation_between(a, b).T


TYPE = NodeType(
    "fusion.grating",
    Category.SPLITTER,
    (light_in("in", "RH"), light_out("out", "RH", label="orders")),
    (
        ParamSpec("table", "file", "", description="Jones table (raytatouille-result file)"),
        ParamSpec(
            "orders", "numbers", [], description="orders to pass on; empty: all of the table"
        ),
        ParamSpec("point", "point", [0.0, 0.0, 0.0], unit="mm"),
        ParamSpec("normal", "point", [0.0, 0.0, 1.0]),
        ParamSpec("grating_vector", "point", [1.0, 0.0, 0.0], description="local x, ⟂ grooves"),
    ),
    (("out", ("in",)),),
    ports=(PortSpec("front", "in", "out", "front"),),
    scatter=(("front", "front", "full", ()),),
    description="Periodic sample (FourierGrating table): diffraction orders as channels",
)

EVALUATORS: dict[str, Evaluator] = {"fusion.grating": grating}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    registry.add(TYPE)
    return dict(EVALUATORS)
