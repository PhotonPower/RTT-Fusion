"""Condensing a group (ADR ports, poses, branches, section 7; stage R29).

A group with `condense` set is not resolved into its blocks: it becomes one block
(`fusion.condensed`) that applies, per inner branch, a table over incident plane waves computed
once from the group's own blocks and kept in the operator store (key: the content of the group,
the branch and the wavelength – not its pose, the table lives in the group's frame).

Table: for each direction σ = (σx, σy) of the incident wave on a square grid (`sigma_max`,
`count`, odd) the coefficients for s and p (s ⟂ plane of incidence; p_in = s × k_in, p_out =
s × k_out), the direction out, the optical path from the entry plane to the exit plane and the
lateral shift between them. That describes a laterally invariant group (layer stacks, plates,
plane splitters) exactly at the grid points; between them bilinear interpolation (the tolerance of
the table). The probes check it: a group whose result depends on where the light enters (a lens)
or that mixes s and p (anisotropic) cannot be condensed this way (`exec.condense`); its operator
needs the adapters' operators (L4). The probes run in air (index 1 at the entry port).

Plane-wave spectra (S, stage R30c, `apply_spectrum`): the table per wave at its direction, the
phase from the probe's optical path and lateral shift; the block takes rays and spectra natively
(fields reach it as rays).

An outer branch through a condensed group still names its inner branch; the block applies exactly
the inner branches the branches of the model use (decision P6), one output bundle per inner
branch.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from rtt_fusion import ports
from rtt_fusion._types import FloatArray
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError, Executor, _hash
from rtt_fusion.graph.light import Plane, RayBundle, physical_split, rotation_between
from rtt_fusion.graph.model import Node
from rtt_fusion.graph.pose import placement
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.graph.spectrum import SpectrumBundle, propagate

TYPE_ID = "fusion.condensed"
INVARIANCE_TOL = 1e-9


@dataclass(frozen=True)
class Table:
    """The plane-wave table of one inner branch at one wavelength (group frame)."""

    sigma: FloatArray  # (n,) grid of σx and σy
    c_s: Any  # (n, n) complex
    c_p: Any
    k_out: FloatArray  # (n, n, 3)
    d_opl: FloatArray  # (n, n)
    shift: FloatArray  # (n, n, 3)
    entry: tuple[
        FloatArray, FloatArray, FloatArray, FloatArray
    ]  # origin, normal (along light), u, v


def _frame(origin: Any, along: Any) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    p = Plane.make(origin, along)
    return (
        np.asarray(p.origin, dtype=float),
        np.asarray(p.normal, dtype=float),
        np.asarray(p.u, dtype=float),
        np.asarray(p.v, dtype=float),
    )


def _variant(ctx: Context, node: Node) -> tuple[Any, dict[str, tuple[str, str]], dict[str, str]]:
    from rtt_fusion.graph import io

    data = node.params["_variant"]
    model = io.model_from_dict(
        {"format": io.FORMAT, "version": io.MODEL_VERSION, **data["model"]}, ctx.registry
    )
    ports = {str(p): (str(e[0]), str(e[1])) for p, e in data["ports"].items()}
    params = {str(k): str(v) for k, v in data.get("params", {}).items()}
    for exported, address in params.items():  # the group's values of its exported parameters
        if exported not in node.params:
            continue
        block, _, param = address.partition(".")
        name, _, index = param.partition("[")
        target = model.nodes[block].params
        if index:
            values = list(target.get(name) or [])
            k = int(index.rstrip("]"))
            while len(values) <= k:
                values.append(0.0)
            values[k] = node.params[exported]
            target[name] = values
        else:
            target[name] = node.params[exported]
    return model, ports, params


def _probe(
    ctx: Context,
    model: Any,
    branch: str,
    entry: tuple[str, str],
    exit_: tuple[str, str],
    frames: dict[tuple[str, str], Any],
    k_in: FloatArray,
    pol: FloatArray,
    offset: FloatArray,
    wavelength_um: float,
) -> tuple[Any, FloatArray, float, FloatArray]:
    """One plane wave through the inner branch: field out for unit field in, direction out,
    optical path and shift from the entry plane to where the light leaves."""
    from rtt_fusion.graph.topology import Model, unroll

    fe, fx = frames[entry], frames[exit_]
    o_e, n_e = np.asarray(fe.origin), np.asarray(fe.normal)
    o_x, n_x = np.asarray(fx.origin), np.asarray(fx.normal)
    probe = Model(
        nodes=dict(model.nodes), connections=list(model.connections), links=list(model.links)
    )
    centre = o_e + offset - 2.0 * k_in
    probe.add(
        Node(
            "_probe",
            "optics.collimated_source",
            params={
                "center": [float(x) for x in centre],
                "direction": [float(x) for x in k_in],
                "radius": 1e-3,
                "rings": 1,
                "wavelengths_um": [wavelength_um],
                "polarisation": [float(x) for x in pol],
            },
        )
    )
    probe.add(
        Node(
            "_sink",
            "optics.detector",
            params={
                "point": [float(x) for x in o_x + 5.0 * n_x],
                "normal": [float(x) for x in -n_x],
            },
        )
    )
    probe.connect(("_probe", "out"), entry)
    probe.connect(exit_, ("_sink", "in"))
    inner = next(b for b in model.branches if b.name == branch)
    probe.branch("probe", ("_probe", None, "out"), *inner.steps, ("_sink", "in", None))
    g = unroll(probe, ctx.registry)
    source = g.links_to("_sink")[0].source
    out = Executor(g, ctx.registry, ctx.evaluators).evaluate(*source)
    if not out:
        raise ExecutionError("exec.condense", f"no light through the inner branch {branch!r}")
    b: RayBundle = out[0]
    i = int(np.argmin(np.linalg.norm(np.asarray(b.source_xy), axis=1)))
    if not bool(np.asarray(b.alive)[i]):
        raise ExecutionError("exec.condense", f"the probe is lost in the inner branch {branch!r}")
    t = float(((o_e - centre) @ n_e) / (k_in @ n_e))
    r_entry = centre + t * k_in
    field = np.asarray(b.amplitude)[i] * b.field_direction()[i]
    return (
        field,
        np.asarray(b.direction)[i],
        float(np.asarray(b.opl)[i] - t),
        np.asarray(b.position)[i] - r_entry,
    )


def build_table(ctx: Context, node: Node, branch: str, wavelength_um: float) -> Table:
    """The table of one inner branch from probe waves through the group's blocks."""
    from rtt_fusion.graph.topology import port_frames

    model, ports, _ = _variant(ctx, node)
    entry_port, exit_port = node.params["_branches"][branch]
    entry, exit_ = ports[entry_port], ports[exit_port]
    frames = port_frames(model, ctx.registry)
    if entry not in frames or exit_ not in frames:
        raise ExecutionError(
            "exec.condense", "the blocks at the group's ports give no port frames", node.id
        )
    settings = node.params.get("condense") or {}
    sigma_max = float(settings.get("sigma_max", 0.2))
    count = int(settings.get("count", 5))
    if count < 2 or count % 2 == 0:
        raise ExecutionError("exec.condense", "count: an odd number ≥ 3 (0 on the grid)", node.id)
    sigma = np.linspace(-sigma_max, sigma_max, count)
    o_e, along, u, v = _frame(frames[entry].origin, -np.asarray(frames[entry].normal))
    c_s = np.zeros((count, count), dtype=complex)
    c_p = np.zeros((count, count), dtype=complex)
    k_out = np.zeros((count, count, 3))
    d_opl = np.zeros((count, count))
    shift = np.zeros((count, count, 3))
    for a, sx in enumerate(sigma):
        for c, sy in enumerate(sigma):
            sz2 = 1.0 - sx * sx - sy * sy
            if sz2 <= 0:
                raise ExecutionError(
                    "exec.condense", "sigma_max reaches grazing incidence", node.id
                )
            k = sx * u + sy * v + np.sqrt(sz2) * along
            s = np.cross(k, along)
            s = s / np.linalg.norm(s) if np.linalg.norm(s) > 1e-12 else u
            p_in = np.cross(s, k)
            f_s, ko, dl, sh = _probe(
                ctx, model, branch, entry, exit_, frames, k, s, np.zeros(3), wavelength_um
            )
            f_p, _, _, _ = _probe(
                ctx, model, branch, entry, exit_, frames, k, p_in, np.zeros(3), wavelength_um
            )
            p_out = np.cross(s, ko)
            if max(abs(p_out @ f_s), abs(s @ f_p)) > INVARIANCE_TOL * max(
                abs(s @ f_s), abs(p_out @ f_p), 1e-300
            ):
                raise ExecutionError(
                    "exec.condense", "the group mixes s and p: no table of this kind", node.id
                )
            c_s[a, c], c_p[a, c] = s @ f_s, p_out @ f_p
            k_out[a, c], d_opl[a, c], shift[a, c] = ko, dl, sh
            if (a, c) in ((count // 2, count // 2), (0, count // 2)):  # lateral invariance
                g_s, ko2, dl2, sh2 = _probe(
                    ctx, model, branch, entry, exit_, frames, k, s, 0.2 * u + 0.1 * v, wavelength_um
                )
                if (
                    abs(s @ g_s - c_s[a, c]) > INVARIANCE_TOL * max(abs(c_s[a, c]), 1e-300)
                    or np.abs(ko2 - ko).max() > INVARIANCE_TOL
                    or abs(dl2 - dl) > INVARIANCE_TOL
                    or np.abs(sh2 - sh).max() > INVARIANCE_TOL
                ):
                    raise ExecutionError(
                        "exec.condense",
                        "the group is not laterally invariant (e.g. a lens): no plane-wave table",
                        node.id,
                    )
    return Table(sigma, c_s, c_p, k_out, d_opl, shift, (o_e, along, u, v))


def _bilinear(table: Table, sx: FloatArray, sy: FloatArray, values: Any) -> Any:
    n = len(table.sigma)
    step = table.sigma[1] - table.sigma[0]
    fx = np.clip((sx - table.sigma[0]) / step, 0.0, n - 1.0)
    fy = np.clip((sy - table.sigma[0]) / step, 0.0, n - 1.0)
    i = np.minimum(np.floor(fx).astype(int), n - 2)
    j = np.minimum(np.floor(fy).astype(int), n - 2)
    wx, wy = fx - i, fy - j
    shape = (-1,) + (1,) * (values.ndim - 2)
    wx, wy = wx.reshape(shape), wy.reshape(shape)
    return (
        values[i, j] * (1 - wx) * (1 - wy)
        + values[i + 1, j] * wx * (1 - wy)
        + values[i, j + 1] * (1 - wx) * wy
        + values[i + 1, j + 1] * wx * wy
    )


def _table(ctx: Context, node: Node, branch: str, wavelength_um: float) -> Table:
    store = ctx.operators if ctx.operators is not None else ctx.shared.setdefault("operators", {})
    local = {
        k: v
        for k, v in node.params.items()
        if k not in ("pose", "pose_offset", "_channels") and not k.startswith("pose_by_")
    }
    key = f"condensed|{_hash(local)}|{branch}|{wavelength_um!r}"
    if key not in store:
        store[key] = build_table(ctx, node, branch, wavelength_um)
        ctx.log.append(f"{node.id}: table of {branch!r} at {wavelength_um} µm built")
    table: Table = store[key]
    return table


def apply(table: Table, b: RayBundle, iso: Any, name: str) -> RayBundle:
    """The table applied to a bundle (world coordinates; `iso` the group's placement)."""
    o_e = np.asarray(iso.point(table.entry[0]))
    along, u, v = (np.asarray(iso.vector(x)) for x in table.entry[1:])
    moved, _ = Plane(o_e, along, u, v).intersect(b)
    k = np.asarray(moved.direction)
    sx, sy = k @ u, k @ v
    inside = (np.abs(sx) <= table.sigma[-1] + 1e-12) & (np.abs(sy) <= table.sigma[-1] + 1e-12)
    c_s = _bilinear(table, sx, sy, table.c_s)
    c_p = _bilinear(table, sx, sy, table.c_p)
    k_out = _bilinear(table, sx, sy, table.k_out) @ iso.rotation.T
    k_out = k_out / np.linalg.norm(k_out, axis=1)[:, None]
    d_opl = _bilinear(table, sx, sy, table.d_opl)
    shift = _bilinear(table, sx, sy, table.shift) @ iso.rotation.T
    s = np.cross(k, along)
    norm = np.linalg.norm(s, axis=1)
    s = np.where(norm[:, None] > 1e-12, s / np.where(norm > 0, norm, 1.0)[:, None], u)
    p_in, p_out = np.cross(s, k), np.cross(s, k_out)
    geometric = np.einsum("ni,nj->nij", s, s) + np.einsum("ni,nj->nij", p_out, p_in)
    physical = c_s[:, None, None] * np.einsum("ni,nj->nij", s, s) + c_p[:, None, None] * np.einsum(
        "ni,nj->nij", p_out, p_in
    )
    # the s/p coefficients of the table in prt (R31: mixed s/p light exact), the scalar factor
    # along the geometrically carried polarisation as for an interface
    factor, prt_out = physical_split(physical, geometric, moved)
    alive = np.asarray(moved.alive) & inside
    return moved.evolve(
        position=np.asarray(moved.position) + shift,
        direction=k_out,
        opl=np.asarray(moved.opl) + d_opl,
        amplitude=np.asarray(moved.amplitude) * np.where(alive, factor, 0.0),
        alive=alive,
        prt=prt_out,
        index=1.0,
        path=(*b.path, name),
    )


def _exit_frame(along: FloatArray, u: FloatArray, v: FloatArray, k_out: FloatArray) -> Any:
    """Axes (u', v') of the plane the spectrum leaves on, normal `k_out`: the entry axes carried
    from `along` to `k_out` by the smallest rotation, or mirrored (v flipped, to keep the frame
    right-handed) if the light runs back."""
    if float(along @ k_out) > 0:
        rot = rotation_between(along, k_out)
        return rot @ u, rot @ v
    n = (k_out - along) / np.linalg.norm(k_out - along)

    def mirror(a: FloatArray) -> FloatArray:
        return np.asarray(a - 2.0 * float(a @ n) * n)

    return mirror(u), -mirror(v)


def apply_spectrum(
    table: Table, sb: SpectrumBundle, iso: Any, name: str, log: list[str]
) -> SpectrumBundle | None:
    """The table applied to a plane-wave spectrum, wave by wave (stage R30c): moved onto the
    entry plane in the spectrum (`spectrum.propagate`, S4), per wave the table at its direction σ
    (bilinear between the grid points, the tolerance of the table): the field c_s·(s·E)·s +
    c_p·(p_in·E)·p_out, the direction k_out and the phase k₀·OPL of the probe ray that entered at
    the entry origin o_e and left at o_e + shift. On the exit plane (origin o_e + shift(σ = 0),
    normal k_out(σ = 0)) the wave has the phase k₀ OPL + k_out·(o_x − o_e − shift): for a plate
    this is k_z′d, the plane-wave transfer. Δ²k is kept (a laterally invariant group maps k∥
    rigidly), dΩ follows in the medium behind (index 1, as the probes). The spectrum keeps no
    sampling (the beam's lateral shift in the group is in the phases, not in a sampling);
    waves outside the table, evanescent or not leaving forward are left out
    (`exec.spectrum_cut`)."""
    o_e = np.asarray(iso.point(table.entry[0]))
    along, u, v = (np.asarray(iso.vector(x)) for x in table.entry[1:])
    moved = propagate(sb, Plane(o_e, along, u, v), log)
    w = moved.waves
    n = moved.refractive_index
    k0 = moved.k0
    if abs(n - 1.0) > 1e-9:
        log.append(f"{name}: spectrum in index {n}, the table was probed in air (index 1)")
    k_loc = np.asarray(w.port.kz(w.kx, w.ky))
    rot = np.stack([u, v, along])
    e_glob = np.column_stack([np.asarray(w.field), w.ez()]) @ rot
    k_dir = np.column_stack([w.kx, w.ky, np.real(k_loc)]) / (n * k0)
    sx, sy = k_dir[:, 0], k_dir[:, 1]
    limit = table.sigma[-1] + 1e-12
    ok = (np.abs(np.imag(k_loc)) == 0.0) & (np.abs(sx) <= limit) & (np.abs(sy) <= limit)
    k = k_dir @ rot
    zero = np.zeros(1)
    k_c = (_bilinear(table, zero, zero, table.k_out) @ iso.rotation.T)[0]
    k_c = k_c / np.linalg.norm(k_c)
    shift_c = (_bilinear(table, zero, zero, table.shift) @ iso.rotation.T)[0]
    o_x = o_e + shift_c
    u_x, v_x = _exit_frame(along, u, v, k_c)
    c_s = _bilinear(table, sx, sy, table.c_s)
    c_p = _bilinear(table, sx, sy, table.c_p)
    k_out = _bilinear(table, sx, sy, table.k_out) @ iso.rotation.T
    k_out = k_out / np.linalg.norm(k_out, axis=1)[:, None]
    d_opl = _bilinear(table, sx, sy, table.d_opl)
    shift = _bilinear(table, sx, sy, table.shift) @ iso.rotation.T
    ok &= k_out @ k_c > 1e-9
    weight = np.sum(np.abs(e_glob) ** 2, axis=1)
    if w.solid_angle is not None:
        weight = weight / np.asarray(w.solid_angle)
    share = float(weight[~ok].sum() / (weight.sum() or 1.0))
    if share > 0.0:
        log.append(
            f"exec.spectrum_cut: {share:.2e} of the power left out of the spectrum (outside the "
            "table, evanescent or not leaving forward)"
        )
    s = np.cross(k, along)
    norm = np.linalg.norm(s, axis=1)
    s = np.where(norm[:, None] > 1e-12, s / np.where(norm > 0, norm, 1.0)[:, None], u)
    p_in, p_out = np.cross(s, k), np.cross(s, k_out)
    e_out = (
        c_s[:, None] * np.sum(s * e_glob, axis=1)[:, None] * s
        + c_p[:, None] * np.sum(p_in * e_glob, axis=1)[:, None] * p_out
    )
    phase = k0 * d_opl + k0 * np.sum(k_out * (o_x - o_e - shift), axis=1)
    e_out = e_out * np.exp(1j * phase)[:, None]
    keep = np.flatnonzero(ok)
    if len(keep) == 0:
        return None
    kvec = k0 * np.column_stack([k_out @ u_x, k_out @ v_x])[keep]
    solid = None
    if w.solid_angle is not None:
        dk2 = np.asarray(w.solid_angle)[keep] * n * k0 * np.real(k_loc[keep])
        solid = dk2 / (k0 * k0 * (k_out[keep] @ np.cross(u_x, v_x)))  # Δ²k/(k k_z), k = k₀
    waves = ports.PlaneWaves(
        ports.Port(0.0, 1, 1.0 + 0j, moved.wavelength_um),
        np.ascontiguousarray(kvec[:, 0]),
        np.ascontiguousarray(kvec[:, 1]),
        np.ascontiguousarray(np.column_stack([e_out @ u_x, e_out @ v_x])[keep]),
        orders=None if w.orders is None else np.asarray(w.orders)[keep],
        solid_angle=solid,
    )
    if moved.support is not None:
        log.append(f"{name}: the spectrum behind the condensed group keeps no sampling for S -> H")
    exit_plane = Plane(o_x, np.cross(u_x, v_x), u_x, v_x)
    return SpectrumBundle(sb.key, exit_plane, waves, (*sb.path, name), sb.branch)


def condensed(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
    """The condensed group: per arriving bundle, the tables of the inner branches the model's
    branches take from its port in to this port out; rays as rays, spectra as spectra (R30c)."""
    exit_port = socket[: -len("_out")]
    channels = node.params.get("_channels") or {}
    iso = placement(node)
    out: list[Any] = []
    for in_socket, bundles in inputs.items():
        entry_port = in_socket[: -len("_in")]
        for name in channels.get(f"{entry_port}>{exit_port}", []):
            for b in bundles or []:
                if isinstance(b, RayBundle):
                    out.append(
                        apply(
                            _table(ctx, node, name, b.wavelength_um), b, iso, f"{node.id}[{name}]"
                        )
                    )
                elif isinstance(b, SpectrumBundle):
                    log: list[str] = []
                    try:
                        got = apply_spectrum(
                            _table(ctx, node, name, b.wavelength_um),
                            b,
                            iso,
                            f"{node.id}[{name}]",
                            log,
                        )
                    finally:
                        ctx.log.extend(f"{node.id}: {line}" for line in log)
                    if got is not None:
                        out.append(got)
    return out


#: rays and spectra natively; fields reach the block as rays (it has no field path)
condensed.representations = frozenset("RS")  # type: ignore[attr-defined]


def condensed_node(gid: str, node: Node, variant: Any, settings: dict[str, Any]) -> Node:
    """The block standing for a condensed group: one port per exported port, the inner branches
    and the variant as data."""
    from rtt_fusion.graph import io
    from rtt_fusion.graph.model import PortSpec, light_in, light_out

    body = io.model_to_dict(variant.model)
    del body["format"], body["version"]
    branches: dict[str, list[str]] = {}
    for b in variant.model.branches:
        first, last = b.steps[0], b.steps[-1]
        entry = next(p for p, e in variant.ports.items() if e == (first.node, first.port_in))
        exit_ = next(p for p, e in variant.ports.items() if e == (last.node, last.port_out))
        branches[b.name] = [entry, exit_]
    through: dict[str, list[str]] = {}
    for entry, exit_ in branches.values():
        through.setdefault(f"{exit_}_out", [])
        if f"{entry}_in" not in through[f"{exit_}_out"]:
            through[f"{exit_}_out"].append(f"{entry}_in")
    names = list(variant.ports)
    return replace(
        node,
        type_id=TYPE_ID,
        params={
            **node.params,
            "_variant": {
                "model": json.loads(json.dumps(body)),
                "ports": {p: list(e) for p, e in variant.ports.items()},
                "params": dict(variant.params),
            },
            "_branches": branches,
            "condense": dict(settings),
        },
        sockets=tuple(
            s for p in names for s in (light_in(f"{p}_in", "RHS"), light_out(f"{p}_out", "RHS"))
        ),
        ports=tuple(PortSpec(p, f"{p}_in", f"{p}_out") for p in names),
        through=tuple((o, tuple(i)) for o, i in through.items()),
    )


EVALUATORS: dict[str, Evaluator] = {TYPE_ID: condensed}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    return dict(EVALUATORS)


__all__ = ["TYPE_ID", "Table", "apply", "build_table", "condensed", "condensed_node", "register"]
