"""A small analytic ray backend with plane elements (stage R15).

Collimated source, plane mirror, plane beam splitter and plane detector – enough for
interferometers such as the Michelson or a Linnik without objectives, and independent of
Raytatouille, so the execution, the coherence rules and the coherent image (pilot W9) can be
checked against the analytic two-beam interference. Elements intersect the incoming rays with
their own plane (free propagation in air, n = 1, in between).

Beam splitter channels follow the Raytatouille splitter block (`rtt.splitter`): light from `in`
leaves reflected (`r`) and transmitted (`t`); light returning on the transmitted arm (`t_in`) is
reflected towards `merged`, light returning on the reflected arm (`r_in`) is transmitted towards
`merged` – the Michelson/Linnik arrangement.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from rtt_fusion._types import FloatArray
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.field import FieldBundle, FieldError, field_image
from rtt_fusion.graph.imaging import coherent_image
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle
from rtt_fusion.graph.model import (
    Category,
    Node,
    NodeType,
    ParamSpec,
    PortFrame,
    PortSpec,
    data_in,
    data_out,
    light_in,
    light_out,
)
from rtt_fusion.graph.registry import NodeTypeRegistry

SQRT_HALF = float(np.sqrt(0.5))


def _vec(node: Node, name: str, default: tuple[float, float, float]) -> FloatArray:
    return np.asarray(node.params.get(name, default), dtype=float)


def hexapolar(rings: int) -> tuple[FloatArray, FloatArray]:
    """Unit-disc points as rt.trace.HexapolarPupil: centre, then 6k points on ring k, from +y."""
    xs, ys = [0.0], [0.0]
    for k in range(1, rings + 1):
        for j in range(6 * k):
            a = 2 * np.pi * j / (6 * k)
            xs.append(k / rings * np.sin(a))
            ys.append(k / rings * np.cos(a))
    return np.asarray(xs), np.asarray(ys)


# ------------------------------------------------------------------------------ evaluators


def source(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[RayBundle]:
    """Collimated beam: disc of `radius` around `center` perpendicular to `direction`; one bundle
    per field angle (deg, about the frame's u and v axes) and wavelength. With `gaussian_waist`
    w > 0 the amplitude is exp(−ρ²/w²) (a Gaussian beam at its waist, cut at `radius`);
    `irradiance` is then the irradiance on the axis."""
    centre = _vec(node, "center", (0.0, 0.0, 0.0))
    frame = Plane.make(
        centre, _vec(node, "direction", (0.0, 0.0, 1.0)), _vec(node, "up", (0.0, 1.0, 0.0))
    )
    radius = float(node.params.get("radius", 1.0))
    px, py = hexapolar(int(node.params.get("rings", 20)))
    xy = np.stack([px * radius, py * radius], axis=-1)
    pos = centre + xy[:, :1] * frame.u + xy[:, 1:] * frame.v
    out = []
    wavelengths = [float(w) for w in node.params.get("wavelengths_um", [0.6328])]
    fields = node.params.get("field_angles_deg", [[0.0, 0.0]])
    for fi, (ax, ay) in enumerate(fields):
        d = np.tan(np.radians(ax)) * frame.u + np.tan(np.radians(ay)) * frame.v + frame.normal
        d = d / np.linalg.norm(d)
        opl0 = (pos - centre) @ d  # plane wave: equal phase on planes ⟂ d
        waist = float(node.params.get("gaussian_waist", 0.0))
        profile = (
            np.exp(-np.sum(xy * xy, axis=1) / waist**2) if waist > 0 else np.ones(len(xy))
        ).astype(complex)
        pol = np.asarray(node.params.get("polarisation", [1.0, 0.0, 0.0]), dtype=complex)
        pol = pol - d * (d @ pol)
        if np.linalg.norm(pol) < 1e-12:
            raise ExecutionError("exec.param", "polarisation along the beam", node.id)
        pol = pol / np.linalg.norm(pol)
        for wi, wl in enumerate(wavelengths):
            out.append(
                RayBundle(
                    key=CoherenceKey(node.id, fi, wi, 0),
                    wavelength_um=wl,
                    position=pos.copy(),
                    direction=np.tile(d, (len(pos), 1)),
                    opl=opl0.copy(),
                    amplitude=profile.copy(),
                    source_xy=xy.copy(),
                    alive=np.ones(len(pos), dtype=bool),
                    irradiance=float(node.params.get("irradiance", 1.0)),
                    path=(node.id,),
                    jones=pol,
                    source_cos=float(abs(d @ frame.normal)),
                    source_frame=np.stack([frame.u, frame.v]),
                )
            )
    return out


def _plane(node: Node) -> Plane:
    """The element's plane; `drive` (mm, set by a phase shifter) shifts it along its normal."""
    normal = _vec(node, "normal", (0.0, 0.0, -1.0))
    normal = normal / np.linalg.norm(normal)
    point = _vec(node, "point", (0.0, 0.0, 0.0)) + float(node.params.get("drive", 0.0)) * normal
    return Plane.make(point, normal, _vec(node, "up", (0.0, 1.0, 0.0)))


def _reflect(b: RayBundle, plane: Plane, factor: complex, name: str) -> RayBundle:
    moved, _ = plane.intersect(b)
    d = np.asarray(moved.direction)
    r = d - 2.0 * (d @ plane.normal)[:, None] * plane.normal
    householder = np.eye(3) - 2.0 * np.outer(plane.normal, plane.normal)
    prt_new = (
        np.tile(householder, (moved.size, 1, 1)).astype(complex)
        if moved.prt is None
        else np.einsum("ij,njk->nik", householder, np.asarray(moved.prt))
    )
    return moved.evolve(
        direction=r,
        amplitude=np.asarray(moved.amplitude) * factor,
        path=(*b.path, name),
        prt=prt_new,
        source_frame=(
            None if moved.source_frame is None else np.asarray(moved.source_frame) @ householder.T
        ),
    )


def _transmit(b: RayBundle, plane: Plane, factor: complex, name: str) -> RayBundle:
    moved, _ = plane.intersect(b)
    return moved.evolve(amplitude=np.asarray(moved.amplitude) * factor, path=(*b.path, name))


def mirror(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[RayBundle]:
    plane = _plane(node)
    r = complex(node.params.get("reflectance_amplitude", -1.0))
    return [_reflect(b, plane, r, node.id) for b in inputs.get("in", [])]


def beam_splitter(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[RayBundle]:
    plane = _plane(node)
    r = complex(node.params.get("r", SQRT_HALF))
    t = complex(node.params.get("t", SQRT_HALF))
    if socket == "r":
        return [_reflect(b, plane, r, node.id) for b in inputs.get("in", [])]
    if socket == "t":
        return [_transmit(b, plane, t, node.id) for b in inputs.get("in", [])]
    if socket == "merged":
        return [_reflect(b, plane, r, node.id) for b in inputs.get("t_in", [])] + [
            _transmit(b, plane, t, node.id) for b in inputs.get("r_in", [])
        ]
    raise ExecutionError("exec.socket_unknown", f"beam splitter has no output {socket!r}", node.id)


def aplanatic_lens(
    ctx: Context, node: Node, socket: str, inputs: dict[str, Any]
) -> list[RayBundle]:
    """Ideal aplanatic lens for collimated light along `axis`: a ray at height h (|h| < f, inside
    the aperture) leaves the Gaussian reference sphere (radius f about the focus F = point + f·axis)
    at sin θ = h/f towards F (sine condition); the OPL is set so that all rays of a plane wave
    along the axis meet F in phase; the polarisation turns in the meridional plane (the
    illumination of Richards and Wolf 1959)."""
    o = _vec(node, "point", (0.0, 0.0, 0.0))
    axis = _vec(node, "axis", (0.0, 0.0, 1.0))
    axis = axis / np.linalg.norm(axis)
    f = float(node.params.get("focal_length", 10.0))
    a = float(node.params.get("aperture_radius", 0.5 * f))
    plane = Plane.make(o, axis)
    focus = o + f * axis
    out = []
    for b in inputs.get("in", []):
        moved, uv = plane.intersect(b)
        h_vec = np.asarray(moved.position) - o
        h_vec = h_vec - (h_vec @ axis)[:, None] * axis
        h = np.linalg.norm(h_vec, axis=1)
        alive = np.asarray(moved.alive) & (h <= a) & (h < f)
        rho_hat = np.where(h[:, None] > 0, h_vec / np.where(h > 0, h, 1.0)[:, None], plane.u)
        sin_t = np.clip(h / f, 0.0, 1.0 - 1e-12)
        cos_t = np.sqrt(1.0 - sin_t**2)
        k = cos_t[:, None] * axis - sin_t[:, None] * rho_hat
        start = focus - f * k
        s_hat = np.cross(axis, rho_hat)
        p_out = cos_t[:, None] * rho_hat + sin_t[:, None] * axis
        rot = (
            np.einsum("ni,nj->nij", s_hat, s_hat)
            + np.einsum("ni,nj->nij", p_out, rho_hat)
            + np.einsum("ni,j->nij", k, axis)
        )
        prt = rot if moved.prt is None else np.einsum("nij,njk->nik", rot, np.asarray(moved.prt))
        out.append(
            moved.evolve(
                position=start,
                direction=k,
                opl=np.asarray(moved.opl),  # every ray then reaches F after f: in phase
                alive=alive,
                prt=prt.astype(complex),
                path=(*b.path, node.id),
            )
        )
    return out


def thin_lens(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[RayBundle]:
    """Paraxial thin lens (R28): a ray at height h with direction tangents τ (transverse
    direction over the axial component) leaves with τ − h/f, for any incidence; the OPL gains
    −h²/(2f) (the paraxial phase of a thin lens). Exact in the sense of matrix optics."""
    o = _vec(node, "point", (0.0, 0.0, 0.0))
    axis = _vec(node, "axis", (0.0, 0.0, 1.0))
    axis = axis / np.linalg.norm(axis)
    f = float(node.params.get("focal_length", 10.0))
    a = float(node.params.get("aperture_radius", 1e9))
    plane = Plane.make(o, axis)
    out = []
    for b in inputs.get("in", []):
        moved, _ = plane.intersect(b)
        h_vec = np.asarray(moved.position) - o
        h_vec = h_vec - (h_vec @ axis)[:, None] * axis
        h2 = np.sum(h_vec * h_vec, axis=1)
        d = np.asarray(moved.direction)
        dz = d @ axis
        tau = (d - dz[:, None] * axis) / np.where(np.abs(dz) > 1e-15, dz, 1.0)[:, None]
        tau_out = tau - h_vec / f
        k = axis[None, :] + tau_out
        k = k / np.linalg.norm(k, axis=1)[:, None]
        out.append(
            moved.evolve(
                direction=k,
                opl=np.asarray(moved.opl) - h2 / (2.0 * f),
                alive=np.asarray(moved.alive) & (h2 <= a * a) & (dz > 0),
                path=(*b.path, node.id),
            )
        )
    return out


def detector(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
    """Coherent image of rays (R, `imaging.coherent_image`) and fields (H, `field.field_image`)."""
    plane = _plane(node)
    hu, hv = (float(x) for x in node.params.get("half_size", [1.0, 1.0]))
    nx, ny = (int(x) for x in node.params.get("pixels", [128, 128]))
    items = list(inputs.get("in", []))
    rays = [b for b in items if isinstance(b, RayBundle)]
    fields = [b for b in items if isinstance(b, FieldBundle)]
    if not fields:
        return coherent_image(rays, plane, hu, hv, nx, ny)
    if rays:
        raise ExecutionError(
            "exec.mixed_representations", "rays and fields on one detector input", node.id
        )
    try:
        return field_image(fields, plane, hu, hv, nx, ny, ctx.log)
    except FieldError as e:
        raise ExecutionError("exec.field", str(e), node.id) from e


detector.representations = frozenset("RH")  # type: ignore[attr-defined]


# ------------------------------------------------------------------------------ registration


def _fresnel(
    b: RayBundle, plane: Plane, n1: float, n2: float, reflect: bool, name: str
) -> RayBundle:
    """Light of `b` (in medium n1) at an interface to n2: reflected or transmitted, Fresnel
    coefficients for s and p (R27). The scalar amplitude is power-normalised (transmission:
    t·sqrt(n2 cos θt / (n1 cos θi))) and taken along the geometrically carried polarisation – exact
    for s or p light and at normal incidence; prt carries the geometry. Total internal reflection:
    the transmitted rays are lost."""
    moved, _ = plane.intersect(b)
    k = np.asarray(moved.direction)
    n = plane.normal
    cos_i = k @ n
    side = np.sign(np.where(cos_i == 0, 1.0, cos_i))[:, None]
    nn = n[None, :] * side  # the normal along the light
    ci = np.abs(cos_i)
    sin_t2 = (n1 / n2) ** 2 * np.clip(1.0 - ci**2, 0.0, None)
    tir = sin_t2 > 1.0
    ct = np.sqrt(np.clip(1.0 - sin_t2, 0.0, None))
    if reflect:
        k_out = k - 2.0 * (k * nn).sum(axis=1)[:, None] * nn
    else:
        k_out = (n1 / n2) * k + ((ct - (n1 / n2) * ci)[:, None]) * nn
        k_out = k_out / np.linalg.norm(k_out, axis=1)[:, None]
    # s ⟂ plane of incidence; p_in = s × k; p_out = k_out × s (reflected) or s × k_out
    # (transmitted), so that at normal incidence the geometric map is the identity transversally
    s_vec = np.cross(k, nn)
    norm = np.linalg.norm(s_vec, axis=1)
    fallback = np.cross(np.tile(plane.u, (len(k), 1)), k)  # normal incidence: any transverse axis
    s_vec = np.where(
        norm[:, None] > 1e-12, s_vec / np.where(norm > 0, norm, 1.0)[:, None], fallback
    )
    s_vec = s_vec / np.linalg.norm(s_vec, axis=1)[:, None]
    p_in = np.cross(s_vec, k)
    p_out = np.cross(k_out, s_vec) if reflect else np.cross(s_vec, k_out)
    with np.errstate(invalid="ignore", divide="ignore"):
        if reflect:
            ct_c = np.sqrt((1.0 - sin_t2).astype(complex))  # imaginary under TIR
            c_s = (n1 * ci - n2 * ct_c) / (n1 * ci + n2 * ct_c)
            c_p = (n1 * ct_c - n2 * ci) / (n1 * ct_c + n2 * ci)
        else:
            gain = np.sqrt(np.where(tir, 0.0, n2 * ct / np.where(ci > 0, n1 * ci, 1.0)))
            c_s = (2 * n1 * ci / (n1 * ci + n2 * ct)) * gain
            c_p = (2 * n1 * ci / (n2 * ci + n1 * ct)) * gain
    geometric = np.einsum("ni,nj->nij", s_vec, s_vec) + np.einsum("ni,nj->nij", p_out, p_in)
    physical = c_s[:, None, None] * np.einsum("ni,nj->nij", s_vec, s_vec) + c_p[
        :, None, None
    ] * np.einsum("ni,nj->nij", p_out, p_in)
    field = moved.field_direction()
    carried = np.einsum("nij,nj->ni", geometric, field)
    out_field = np.einsum("nij,nj->ni", physical, field)
    norm2 = np.sum(np.abs(carried) ** 2, axis=1)
    factor = np.sum(out_field * np.conj(carried), axis=1) / np.where(norm2 > 0, norm2, 1.0)
    prt = np.asarray(moved.prt) if moved.prt is not None else np.tile(np.eye(3), (len(k), 1, 1))
    alive = np.asarray(moved.alive) & (~tir if not reflect else True)
    return moved.evolve(
        direction=k_out,
        amplitude=np.asarray(moved.amplitude) * np.where(alive, factor, 0.0),
        alive=alive,
        prt=np.einsum("nij,njk->nik", geometric, prt).astype(complex),
        index=n2 if not reflect else n1,
        path=(*b.path, name),
        source_frame=moved.source_frame,
    )


def interface(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[RayBundle]:
    """Plane interface between `n_front` (the side `normal` points to) and `n_back` (R27): light
    in through one side leaves reflected through the same side, transmitted through the other."""
    plane = _plane(node)
    n_side = {
        "front": float(node.params.get("n_front", 1.0)),
        "back": float(node.params.get("n_back", 1.5)),
    }
    out_side = socket.split("_")[0]
    out: list[RayBundle] = []
    for in_socket, bundles in inputs.items():
        in_side = in_socket.split("_")[0]
        other = "back" if in_side == "front" else "front"
        for b in bundles or []:
            out.append(
                _fresnel(b, plane, n_side[in_side], n_side[other], in_side == out_side, node.id)
            )
    return out


def splitter(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[RayBundle]:
    """Lossless thin splitter with four ports (R27): front_a, front_b on the side `normal`
    points to, back_a, back_b behind; front_a and back_a lie on one line, as front_b and back_b.
    Light in at a front port leaves reflected (r) at the other front port and transmitted (t) at
    the back port on its line; from the back r' = −r*·t/t* (Stokes, lossless and reciprocal:
    S = Sᵀ unitary)."""
    plane = _plane(node)
    r = complex(node.params.get("r", SQRT_HALF))
    t = complex(node.params.get("t", SQRT_HALF))
    r_back = -np.conj(r) * t / np.conj(t)
    out_side, out_line = socket.split("_")[0], socket.split("_")[1]
    out: list[RayBundle] = []
    for in_socket, bundles in inputs.items():
        in_side, in_line = in_socket.split("_")[0], in_socket.split("_")[1]
        for b in bundles or []:
            if in_side == out_side and in_line != out_line:
                out.append(_reflect(b, plane, r if in_side == "front" else r_back, node.id))
            elif in_side != out_side and in_line == out_line:
                out.append(_transmit(b, plane, t, node.id))
    return out


def _interface_frames(node: Node) -> dict[str, PortFrame]:
    plane = _plane(node)
    o, n = _point(plane.origin), np.asarray(plane.normal)
    return {
        "front": PortFrame(o, _unit(n), None, float(node.params.get("n_front", 1.0))),
        "back": PortFrame(o, _unit(-n), None, float(node.params.get("n_back", 1.5))),
    }


def _unit(v: Any) -> tuple[float, float, float]:
    a = np.asarray(v, dtype=float)
    a = a / np.linalg.norm(a)
    return (float(a[0]), float(a[1]), float(a[2]))


def _point(v: Any) -> tuple[float, float, float]:
    a = np.asarray(v, dtype=float)
    return (float(a[0]), float(a[1]), float(a[2]))


def _source_frames(node: Node) -> dict[str, PortFrame]:
    return {
        "out": PortFrame(
            _point(_vec(node, "center", (0.0, 0.0, 0.0))),
            _unit(_vec(node, "direction", (0.0, 0.0, 1.0))),
            float(node.params.get("radius", 1.0)),
        )
    }


def _mirror_frames(node: Node) -> dict[str, PortFrame]:
    plane = _plane(node)
    return {"front": PortFrame(_point(plane.origin), _unit(plane.normal))}


def _detector_frames(node: Node) -> dict[str, PortFrame]:
    plane = _plane(node)
    hu, hv = (float(x) for x in node.params.get("half_size", [1.0, 1.0]))
    return {"in": PortFrame(_point(plane.origin), _unit(plane.normal), float(np.hypot(hu, hv)))}


def _lens_frames(node: Node) -> dict[str, PortFrame]:
    point = _point(_vec(node, "point", (0.0, 0.0, 0.0)))
    axis = np.asarray(_unit(_vec(node, "axis", (0.0, 0.0, 1.0))))
    f = float(node.params.get("focal_length", 10.0))
    radius = float(node.params.get("aperture_radius", 0.5 * f))
    return {
        "in": PortFrame(point, _unit(-axis), radius),
        "out": PortFrame(point, _unit(axis), radius),
    }


_PLANE_PARAMS = (
    ParamSpec("point", "point", [0.0, 0.0, 0.0], unit="mm"),
    ParamSpec("normal", "point", [0.0, 0.0, -1.0]),
    ParamSpec("up", "point", [0.0, 1.0, 0.0], description="direction of the local v axis"),
    ParamSpec("drive", "number", 0.0, unit="mm", description="shift along the normal"),
)

TYPES = (
    NodeType(
        "optics.collimated_source",
        Category.SOURCE,
        (light_out("out", "RH"),),
        (
            ParamSpec("center", "point", [0.0, 0.0, 0.0], unit="mm"),
            ParamSpec("direction", "point", [0.0, 0.0, 1.0]),
            ParamSpec("radius", "number", 1.0, unit="mm", minimum=0.0),
            ParamSpec("rings", "integer", 20, minimum=1),
            ParamSpec("wavelengths_um", "numbers", [0.6328], unit="µm"),
            ParamSpec("field_angles_deg", "numbers", [[0.0, 0.0]], unit="deg"),
            ParamSpec("irradiance", "number", 1.0, unit="W/mm²"),
            ParamSpec("polarisation", "point", [1.0, 0.0, 0.0], description="field direction"),
            ParamSpec("up", "point", [0.0, 1.0, 0.0], description="direction of the local v axis"),
            ParamSpec(
                "gaussian_waist",
                "number",
                0.0,
                unit="mm",
                minimum=0.0,
                description="0: uniform; > 0: Gaussian amplitude exp(−ρ²/w²)",
            ),
        ),
        port_frames=_source_frames,
        description="Collimated beam (plane wave per field angle and wavelength)",
    ),
    NodeType(
        "optics.plane_mirror",
        Category.DOMAIN,
        (light_in("in", "RH"), light_out("out", "RH"), data_in("drive", "drive", "drive")),
        (*_PLANE_PARAMS, ParamSpec("reflectance_amplitude", "number", -1.0)),
        (("out", ("in",)),),
        ports=(PortSpec("front", "in", "out", "front"),),
        port_frames=_mirror_frames,
        description="Plane mirror",
    ),
    NodeType(
        "optics.beam_splitter",
        Category.SPLITTER,
        (
            light_in("in", "RH"),
            light_out("r", "RH"),
            light_out("t", "RH"),
            light_in("r_in", "RH"),
            light_in("t_in", "RH"),
            light_out("merged", "RH"),
        ),
        (*_PLANE_PARAMS, ParamSpec("r", "number", SQRT_HALF), ParamSpec("t", "number", SQRT_HALF)),
        (("r", ("in",)), ("t", ("in",)), ("merged", ("r_in", "t_in"))),
        description="Plane beam splitter (Michelson / Linnik arrangement)",
    ),
    NodeType(
        "optics.thin_lens",
        Category.DOMAIN,
        (light_in("in", "RH"), light_out("out", "RH")),
        (
            ParamSpec("point", "point", [0.0, 0.0, 0.0], unit="mm"),
            ParamSpec("axis", "point", [0.0, 0.0, 1.0]),
            ParamSpec("focal_length", "number", 10.0, unit="mm"),
            ParamSpec("aperture_radius", "number", 1e9, unit="mm"),
        ),
        (("out", ("in",)),),
        port_frames=_lens_frames,
        description="Paraxial thin lens (matrix optics, any incidence)",
    ),
    NodeType(
        "optics.interface",
        Category.DOMAIN,
        (
            light_in("front_in", "RH", "front"),
            light_out("front_out", "RH", "front"),
            light_in("back_in", "RH", "back"),
            light_out("back_out", "RH", "back"),
        ),
        (
            *_PLANE_PARAMS,
            ParamSpec(
                "n_front", "number", 1.0, description="index on the side the normal points to"
            ),
            ParamSpec("n_back", "number", 1.5, description="index behind"),
        ),
        (
            ("front_out", ("front_in", "back_in")),
            ("back_out", ("front_in", "back_in")),
        ),
        scatter=tuple((i, o, "full", ()) for i in ("front", "back") for o in ("front", "back")),
        port_frames=_interface_frames,
        description="Plane interface of two media (Fresnel, s and p)",
    ),
    NodeType(
        "optics.splitter",
        Category.SPLITTER,
        tuple(
            spec
            for side in ("front", "back")
            for line in ("a", "b")
            for spec in (
                light_in(f"{side}_{line}_in", "RH", side),
                light_out(f"{side}_{line}_out", "RH", side),
            )
        ),
        (*_PLANE_PARAMS, ParamSpec("r", "number", SQRT_HALF), ParamSpec("t", "number", SQRT_HALF)),
        (
            ("front_a_out", ("front_b_in", "back_a_in")),
            ("front_b_out", ("front_a_in", "back_b_in")),
            ("back_a_out", ("back_b_in", "front_a_in")),
            ("back_b_out", ("back_a_in", "front_b_in")),
        ),
        description="Lossless splitter with four ports (Stokes relations)",
    ),
    NodeType(
        "optics.aplanatic_lens",
        Category.DOMAIN,
        (light_in("in", "RH"), light_out("out", "RH")),
        (
            ParamSpec("point", "point", [0.0, 0.0, 0.0], unit="mm"),
            ParamSpec("axis", "point", [0.0, 0.0, 1.0]),
            ParamSpec("focal_length", "number", 10.0, unit="mm"),
            ParamSpec("aperture_radius", "number", 5.0, unit="mm"),
        ),
        (("out", ("in",)),),
        port_frames=_lens_frames,
        description="Ideal aplanatic lens (sine condition, Richards-Wolf illumination)",
    ),
    NodeType(
        "optics.detector",
        Category.DETECTOR,
        (light_in("in", "RH"), data_out("image", "result")),
        (
            *_PLANE_PARAMS,
            ParamSpec("half_size", "numbers", [1.0, 1.0], unit="mm"),
            ParamSpec("pixels", "numbers", [128, 128]),
        ),
        (("image", ("in",)),),
        port_frames=_detector_frames,
        description="Plane detector: coherent image (amplitudes per coherence group)",
    ),
)

EVALUATORS: dict[str, Evaluator] = {
    "optics.collimated_source": source,
    "optics.plane_mirror": mirror,
    "optics.interface": interface,
    "optics.thin_lens": thin_lens,
    "optics.splitter": splitter,
    "optics.beam_splitter": beam_splitter,
    "optics.aplanatic_lens": aplanatic_lens,
    "optics.detector": detector,
}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    """Adds the block types to a registry and returns their evaluators."""
    for t in TYPES:
        registry.add(t)
    return dict(EVALUATORS)
