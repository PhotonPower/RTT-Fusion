"""Ideal aplanatic imaging of spectra (pilot K5, stage K5c): `optics.ideal_imaging`.

An infinity-corrected microscope as one block: objective (object-side index n, the medium of the
incoming spectrum, numerical aperture `na`) and tube lens (image-side index n'), lateral
magnification M (the image is inverted). It works on the angular spectrum S of light from near
its object point (an emitter, a scatterer), exactly in the sense of Richards and Wolf:

- the spectrum is first moved onto the object plane (`point`, normal `axis`; `spectrum.propagate`);
- every wave k̂ = (sin θ cos φ, sin θ sin φ, cos θ) with n sin θ ≤ NA becomes the image-side wave
  k̂' = (−sin θ' cos φ, −sin θ' sin φ, cos θ'), n' sin θ' = n sin θ / M (k'∥ = −k∥/M: a point at
  Δ in the object plane goes to −M·Δ in the image plane);
- its field turns in the meridional plane as through an aplanatic objective and tube lens
  (θ̂ → ρ̂ in the pupil → p' = cos θ' ρ̂ + sin θ' n̂; φ̂ stays):
  A' = c(θ) [(A·φ̂) φ̂ + (A·θ̂) p'], c = (n'M/n)·√(n'/n)·√(cos θ'/cos θ), so that every wave keeps
  its power ((2π)²/k ∫|A|² dΩ = (2π)²/k' ∫|A'|² dΩ');
- the phase is kept: the object point on the axis focuses at `image_point` in phase.

The output spectrum lies on the image plane (`image_point`, normal `axis`) with the sampling of a
point source there (as the particle's far field, R30d), so a detector near the image plane sums
it directly (the Debye integral without its approximation, `spectrum.spectrum_image`) and S → H
works. Waves beyond the NA are left out (logged with their share of the power). Valid near the
object point (aplanatic: the sine condition, no field aberrations modelled).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.light import Plane
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
from rtt_fusion.ports import PlaneWaves, Port


def _vec(node: Node, name: str, default: tuple[float, float, float]) -> Any:
    return np.asarray(node.params.get(name, default), dtype=float)


def image_spectrum(sb: Any, node: Node, notes: list[str]) -> Any:
    """One spectrum through the ideal microscope (module docstring)."""
    from fieldtracing.field import GridlessField

    from rtt_fusion.graph.spectrum import SpectrumBundle, propagate

    if not sb.continuous:
        raise ExecutionError("exec.param", "ideal imaging needs a continuous spectrum", node.id)
    axis = _vec(node, "axis", (0.0, 0.0, 1.0))
    axis = axis / np.linalg.norm(axis)
    up = _vec(node, "up", (0.0, 1.0, 0.0))
    obj = Plane.make(_vec(node, "point", (0.0, 0.0, 0.0)), axis, up)
    at = propagate(sb, obj, notes)
    w = at.waves
    n = float(np.real(w.port.index))
    n_img = float(node.params.get("n_image", 1.0))
    mag = float(node.params.get("magnification", 40.0))
    na = float(node.params.get("na", 1.4))
    if mag <= 0 or n_img <= 0:
        raise ExecutionError("exec.param", "magnification and n_image must be positive", node.id)
    k0 = w.port.k0
    k, k_img = n * k0, n_img * k0
    kx, ky = np.asarray(w.kx, dtype=float), np.asarray(w.ky, dtype=float)
    kz = np.real(w.port.kz(kx, ky))
    e = np.column_stack([np.asarray(w.field), w.ez()]) / np.asarray(w.solid_angle)[:, None]  # A
    sin_t = np.hypot(kx, ky) / k
    keep = (sin_t * n <= na) & (kz > 0)
    total = float(np.sum(np.sum(np.abs(e) ** 2, 1) * np.asarray(w.solid_angle)))
    lost = float(np.sum((np.sum(np.abs(e) ** 2, 1) * np.asarray(w.solid_angle))[~keep]))
    if total > 0 and lost > 1e-9 * total:
        notes.append(f"{lost / total:.2e} of the spectrum's power lies beyond NA {na:g}")
    kx, ky, kz, e, sin_t = kx[keep], ky[keep], kz[keep], e[keep], sin_t[keep]
    omega = np.asarray(w.solid_angle)[keep]
    cos_t = kz / k
    phi = np.arctan2(ky, kx)
    c, s = np.cos(phi), np.sin(phi)
    zeros = np.zeros_like(c)
    th_hat = np.stack([cos_t * c, cos_t * s, -sin_t], -1)
    ph_hat = np.stack([-s, c, zeros], -1)
    sin_i = n * sin_t / (n_img * mag)
    if np.any(sin_i >= 1):
        raise ExecutionError("exec.param", "the image-side NA exceeds n_image", node.id)
    cos_i = np.sqrt(1 - sin_i**2)
    p_img = np.stack([cos_i * c, cos_i * s, sin_i], -1)
    scale = (n_img * mag / n) * np.sqrt(n_img / n) * np.sqrt(cos_i / cos_t)
    a_img = scale[:, None] * (
        np.einsum("ni,ni->n", e, ph_hat)[:, None] * ph_hat
        + np.einsum("ni,ni->n", e, th_hat)[:, None] * p_img
    )
    kx_i, ky_i, kz_i = -k_img * sin_i * c, -k_img * sin_i * s, k_img * cos_i
    # dΩ' = (n/(n'M))² cos θ / cos θ' dΩ
    omega_i = omega * (n / (n_img * mag)) ** 2 * cos_t / cos_i
    image = Plane(_vec(node, "image_point", (0.0, 0.0, 200.0)), obj.normal, obj.u, obj.v)
    port = Port(0.0, 1, complex(n_img), w.port.wavelength_um)
    support = GridlessField(
        np.column_stack([kx_i, ky_i]),
        np.zeros(len(kx_i)),
        np.zeros((len(kx_i), 2)),
        np.zeros((len(kx_i), 2, 2)),
        np.asarray(a_img / (k_img * kz_i)[:, None], dtype=np.complex128),
        float(w.port.wavelength_um),
        n_img,
        "frequency",
    )
    field = a_img * omega_i[:, None]
    return SpectrumBundle(
        sb.key,
        image,
        PlaneWaves(port, kx_i, ky_i, np.ascontiguousarray(field[:, :2]), solid_angle=omega_i),
        (*sb.path, node.id),
        sb.branch,
        support,
    )


def ideal_imaging(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
    from rtt_fusion.graph.spectrum import SpectrumBundle

    notes: list[str] = []
    out = []
    for sb in inputs.get("in", []):
        if not isinstance(sb, SpectrumBundle):
            raise ExecutionError("exec.param", "ideal imaging takes spectra (S)", node.id)
        out.append(image_spectrum(sb, node, notes))
    ctx.log.extend(f"{node.id}: {line}" for line in dict.fromkeys(notes))
    return out


ideal_imaging.representations = frozenset("S")  # type: ignore[attr-defined]


def _frames(node: Node) -> dict[str, PortFrame]:
    p = _vec(node, "point", (0.0, 0.0, 0.0))
    q = _vec(node, "image_point", (0.0, 0.0, 200.0))
    a = _vec(node, "axis", (0.0, 0.0, 1.0))
    a = a / np.linalg.norm(a)
    return {
        "in": PortFrame((p[0], p[1], p[2]), (-a[0], -a[1], -a[2])),
        "out": PortFrame((q[0], q[1], q[2]), (a[0], a[1], a[2])),
    }


TYPES = (
    NodeType(
        "optics.ideal_imaging",
        Category.DOMAIN,
        (light_in("in", "S"), light_out("out", "S")),
        (
            ParamSpec(
                "point",
                "point",
                [0.0, 0.0, 0.0],
                unit="mm",
                description="object point (on the axis, in the object plane)",
            ),
            ParamSpec("axis", "point", [0.0, 0.0, 1.0]),
            ParamSpec("up", "point", [0.0, 1.0, 0.0], description="fixes the in-plane axes"),
            ParamSpec(
                "image_point",
                "point",
                [0.0, 0.0, 200.0],
                unit="mm",
                description="where the object point is imaged",
            ),
            ParamSpec(
                "magnification",
                "number",
                40.0,
                minimum=0.0,
                description="lateral magnification |M| (the image is inverted)",
            ),
            ParamSpec("na", "number", 1.4, minimum=0.0, description="object-side NA"),
            ParamSpec("n_image", "number", 1.0, minimum=1.0, description="image-side index"),
        ),
        port_frames=_frames,
        description="Ideal aplanatic microscope (objective + tube lens, magnification M) for "
        "angular spectra (S), Richards–Wolf",
    ),
)

EVALUATORS: dict[str, Evaluator] = {"optics.ideal_imaging": ideal_imaging}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    for t in TYPES:
        registry.add(t)
    return dict(EVALUATORS)
