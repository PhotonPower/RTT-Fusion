"""Emitters at a sensor surface as a light source (pilot K5, stage K5b): `emitter.dipole_layer`.

Fluorophores of a receptor layer (Q from `biooptomd.layer`, or one dipole given by parameters)
emit through a planar stack – substrate, optional layers below the surface (e.g. a gold film),
the adlayer n(z) from N above it, the solution – into the half-space the block collects from.
Each emitter is solved with the stack (`couplings.dipole_layers`, near field and supercritical
light included) and leaves as one continuous spectrum S: emitters are mutually incoherent, so
each gets its own coherence group (`CoherenceKey(node, emitter, 0, 0)`).

Frames (D13, DP8): the stack is computed in BioOptoMD's sensor frame (surface z = 0, solution
z > 0, lengths nm for Q). The block places it with its pose: the collected half-space lies along
`direction`; for "substrate" the sensor frame is mirrored (x, y, z) → (u, v, −n) about the surface,
for "solution" it is (u, v, n) – the rule of `DipoleSeries.in_rtt_frame` (light_from "solution"
resp. "substrate"), with `center` the point on the surface under the anchor (the mean lateral
position of the samples, as BindOptiX's `from_tracker`) and `up` fixing u, v (`Plane.make`).

Each spectrum carries the sampling of a point source at the emitter's foot point on the surface
(ψ̃ = 0, −ρ = 0, Hessian 0, as the particle's far field, R30d) and is moved onto the block's port
plane `distance` along `direction` (`spectrum.propagate`), so S → H → R works behind the block.
Waves: Gauss–Legendre in θ on [0, θ_c] and [θ_c, θ_max] (θ_c the critical angle, where the
spectrum has a kink), a uniform φ grid. Power: a free dipole in the emitter's medium would emit
`power_w`; each spectrum carries `power_w` times its share `down` (resp. `up`) of `rates`.
"""

from __future__ import annotations

from itertools import pairwise
from typing import Any

import numpy as np

from rtt_fusion.adapters import biooptomd as bmd
from rtt_fusion.couplings import dipole_layers as dl
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.light import CoherenceKey, Plane
from rtt_fusion.graph.model import (
    Category,
    Node,
    NodeType,
    ParamSpec,
    PortFrame,
    data_in,
    light_out,
)
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.ports import PlaneWaves, Port

NM = 1e-6  # mm per nm


def _vec(node: Node, name: str, default: tuple[float, float, float]) -> Any:
    return np.asarray(node.params.get(name, default), dtype=float)


def stack_of(node: Node, profile: bmd.Profile | None) -> dl.Stack:
    """The stack in the sensor frame: substrate, `layers` below the surface (rows thickness nm,
    Re n, Im n, bottom-up, ending at z = 0), the adlayer of N above the surface, the solution."""
    rows = [list(map(float, r)) for r in node.params.get("layers", []) or []]
    below: list[tuple[complex, float]] = []
    for r in rows:
        if len(r) not in (2, 3):
            raise ExecutionError("exec.param", "layers: rows (thickness nm, Re n[, Im n])", node.id)
        below.append((complex(r[1], r[2] if len(r) == 3 else 0.0), r[0] * NM))
    n_solution = complex(float(node.params.get("n_solution", 1.333)))
    above: list[tuple[complex, float]] = []
    if profile is not None:
        thick, n = bmd.profile_layers(profile)
        above = [(complex(v), float(t) * NM) for t, v in zip(thick, n, strict=True)]
        n_solution = complex(profile.n_solvent)
    z_bottom = -sum(t for _, t in below)
    return dl.Stack(
        complex(float(node.params.get("n_substrate", 1.52))),
        tuple(below + above),
        n_solution,
        z_bottom,
    )


def emitters(node: Node, q: bmd.DipoleSeries | None) -> tuple[Any, Any]:
    """(positions (N, 3) in nm in the sensor frame relative to the anchor, unit dipoles (N, 3))."""
    if q is None:
        height = float(node.params.get("height_nm", 5.0))
        mu = _vec(node, "dipole", (0.0, 0.0, 1.0))
        if np.linalg.norm(mu) == 0:
            raise ExecutionError("exec.param", "dipole: zero orientation", node.id)
        return np.array([[0.0, 0.0, height]]), (mu / np.linalg.norm(mu))[None, :]
    state = str(node.params.get("state", "") or "") or q.states[0]
    if state not in q.series:
        raise ExecutionError("exec.data", f"state {state!r} not in Q {q.states}", node.id)
    every: Any = np.asarray(q.series[state], dtype=float).reshape(-1, 7)
    a = every
    count = int(node.params.get("samples", 20))
    if count < len(every):
        rng = np.random.default_rng(int(node.params.get("seed", 0)))
        a = every[np.sort(rng.choice(len(every), count, replace=False))]
    xy = a[:, 0:2] - every[:, 0:2].mean(axis=0)
    pos = np.column_stack([xy, a[:, 6]])
    mu = a[:, 3:6] / np.linalg.norm(a[:, 3:6], axis=1, keepdims=True)
    return pos, mu


def _angles(n_out: float, n_other: float, na: float, n_theta: int) -> tuple[Any, Any]:
    """θ nodes and weights sin θ dθ: Gauss–Legendre on [0, θ_c] and [θ_c, θ_max]."""
    theta_max = float(np.arcsin(min(na / n_out, 0.999)))
    edges = [0.0, theta_max]
    if n_other < n_out:
        crit = float(np.arcsin(n_other / n_out))
        if crit < theta_max:
            edges = [0.0, crit, theta_max]
    x, w = np.polynomial.legendre.leggauss(n_theta)
    nodes, weights = [], []
    for a, b in pairwise(edges):
        t = a + (b - a) * (x + 1) / 2
        nodes.append(t)
        weights.append(w * (b - a) / 2 * np.sin(t))
    return np.concatenate(nodes), np.concatenate(weights)


def dipole_layer(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
    """One continuous spectrum per emitter on the port plane (module docstring)."""
    from fieldtracing.field import GridlessField

    from rtt_fusion.graph.spectrum import SpectrumBundle, propagate

    q = inputs.get("dipoles") or None  # an input without a link arrives empty
    if isinstance(q, list):
        q = q[0] if len(q) == 1 else None
    if q is not None and not isinstance(q, bmd.DipoleSeries):
        raise ExecutionError("exec.data", "dipoles: not a BioOptoMD dipole series (Q)", node.id)
    profile = inputs.get("profile") or None
    if isinstance(profile, list):
        profile = profile[0] if len(profile) == 1 else None
    if profile is not None and not isinstance(profile, bmd.Profile):
        raise ExecutionError("exec.data", "profile: not a BioOptoMD index profile (N)", node.id)
    stack = stack_of(node, profile)
    pos_nm, mus = emitters(node, q)
    lam = float(node.params.get("wavelength_um", 0.67))
    side = "bottom" if str(node.params.get("collect", "substrate")) == "substrate" else "top"
    n_out = float(np.real(stack.n_bottom if side == "bottom" else stack.n_top))
    n_other = float(np.real(stack.n_top if side == "bottom" else stack.n_bottom))
    na = float(node.params.get("na", n_out))
    theta, wt = _angles(n_out, n_other, na, int(node.params.get("n_theta", 32)))
    n_phi = int(node.params.get("n_phi", 48))
    phi = np.arange(n_phi) * 2 * np.pi / n_phi
    weight = (wt[:, None] * np.full(n_phi, 2 * np.pi / n_phi)[None, :]).ravel()
    # local frame: n along `direction`; sensor (x, y, ±z) → (u, v, n)
    plane = Plane.make(
        _vec(node, "center", (0.0, 0.0, 0.0)),
        _vec(node, "direction", (0.0, 0.0, 1.0)),
        _vec(node, "up", (0.0, 1.0, 0.0)),
    )
    zsign = -1.0 if side == "bottom" else 1.0
    k0 = 2 * np.pi / (lam * 1e-3)
    k = n_out * k0
    tt, pp = np.meshgrid(theta, phi, indexing="ij")
    kx = (k * np.sin(tt) * np.cos(pp)).ravel()
    ky = (k * np.sin(tt) * np.sin(pp)).ravel()
    kz = (k * np.cos(tt)).ravel()
    port = Port(0.0, 1, complex(n_out), lam)
    power_w = float(node.params.get("power_w", 1.0))
    distance = float(node.params.get("distance", 0.01))
    target = Plane(plane.origin + distance * plane.normal, plane.normal, plane.u, plane.v)
    surface_z = stack.interfaces[0] if side == "bottom" else stack.interfaces[-1]
    out = []
    notes: list[str] = []
    cache: dict[float, Any] = {}
    for i, (r, mu) in enumerate(zip(pos_nm, mus, strict=True)):
        z0 = float(r[2]) * NM
        if z0 not in cache:
            try:
                cache[z0] = dl.angular_spectrum_grid(stack, lam, z0, theta, phi, side)
            except ValueError as e:
                raise ExecutionError("exec.param", str(e), node.id) from e
        a_sensor = np.einsum("j,jtpi->tpi", mu, cache[z0]).reshape(-1, 3)
        scale = np.sqrt(power_w / (2 * dl.free_power(stack, lam, z0)))
        a = a_sensor * scale * np.array([1.0, 1.0, zsign])  # sensor → local (u, v, n)
        # referred to r_ref of the module: the outer interface of the stack under the emitter
        foot = plane.origin + r[0] * NM * plane.u + r[1] * NM * plane.v
        at = Plane(foot + zsign * surface_z * plane.normal, plane.normal, plane.u, plane.v)
        field = a * weight[:, None]
        support = GridlessField(
            np.column_stack([kx, ky]),
            np.zeros(len(kx)),
            np.zeros((len(kx), 2)),
            np.zeros((len(kx), 2, 2)),
            np.asarray(a / (k * kz)[:, None], dtype=np.complex128),
            lam,
            n_out,
            "frequency",
        )
        sb = SpectrumBundle(
            CoherenceKey(node.id, i, 0, 0),
            at,
            PlaneWaves(port, kx, ky, np.ascontiguousarray(field[:, :2]), solid_angle=weight),
            (node.id,),
            None,
            support,
        )
        out.append(propagate(sb, target, notes))
    ctx.log.extend(f"{node.id}: {line}" for line in dict.fromkeys(notes))
    return out


dipole_layer.representations = frozenset("S")  # type: ignore[attr-defined]


def _frames(node: Node) -> dict[str, PortFrame]:
    centre = _vec(node, "center", (0.0, 0.0, 0.0))
    d = _vec(node, "direction", (0.0, 0.0, 1.0))
    d = d / np.linalg.norm(d)
    dist = float(node.params.get("distance", 0.01))
    o = centre + dist * d
    return {"out": PortFrame((float(o[0]), float(o[1]), float(o[2])), (d[0], d[1], d[2]))}


TYPES = (
    NodeType(
        "emitter.dipole_layer",
        Category.SOURCE,
        (data_in("dipoles", "Q"), data_in("profile", "N"), light_out("out", "S")),
        (
            ParamSpec(
                "center",
                "point",
                [0.0, 0.0, 0.0],
                unit="mm",
                description="point on the sensor surface under the emitters",
            ),
            ParamSpec(
                "direction",
                "point",
                [0.0, 0.0, 1.0],
                description="direction of the collected light (away from the surface)",
            ),
            ParamSpec("up", "point", [0.0, 1.0, 0.0], description="fixes the in-plane axes"),
            ParamSpec(
                "collect",
                "choice",
                "substrate",
                options=("substrate", "solution"),
                description="half-space the light is collected from",
            ),
            ParamSpec(
                "wavelength_um",
                "number",
                0.67,
                unit="µm",
                minimum=0.0,
                description="emission wavelength",
            ),
            ParamSpec("n_substrate", "number", 1.52, minimum=1.0),
            ParamSpec(
                "n_solution",
                "number",
                1.333,
                minimum=1.0,
                description="index above the surface (N's solvent if N is connected)",
            ),
            ParamSpec(
                "layers",
                "numbers",
                [],
                description="rows (thickness nm, Re n, Im n) "
                "between substrate and surface, bottom-up (e.g. 50, 0.197, 3.44 for gold)",
            ),
            ParamSpec("state", "text", "", description="receptor state of Q (empty: the first)"),
            ParamSpec(
                "samples",
                "integer",
                20,
                minimum=1,
                description="emitters drawn from Q (each a coherence group)",
            ),
            ParamSpec("seed", "integer", 0),
            ParamSpec(
                "height_nm",
                "number",
                5.0,
                unit="nm",
                minimum=0.0,
                description="without Q: height of the one dipole above the surface",
            ),
            ParamSpec(
                "dipole",
                "point",
                [0.0, 0.0, 1.0],
                description="without Q: orientation of the one dipole (sensor frame)",
            ),
            ParamSpec(
                "power_w",
                "number",
                1.0,
                unit="W",
                minimum=0.0,
                description="power a free dipole in the emitter's medium would emit",
            ),
            ParamSpec(
                "na",
                "number",
                1.52,
                minimum=0.0,
                description="largest n sin θ of the waves (default: up to grazing)",
            ),
            ParamSpec("n_theta", "integer", 32, minimum=4),
            ParamSpec("n_phi", "integer", 48, minimum=8),
            ParamSpec(
                "distance",
                "number",
                0.01,
                unit="mm",
                minimum=0.0,
                description="port plane in front of the surface along direction",
            ),
        ),
        port_frames=_frames,
        description="Fluorophores at a sensor surface (Q, N) emitting through a planar stack, "
        "near field and supercritical light included (S, one coherence group per emitter)",
    ),
)

EVALUATORS: dict[str, Evaluator] = {"emitter.dipole_layer": dipole_layer}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    for t in TYPES:
        registry.add(t)
    return dict(EVALUATORS)
