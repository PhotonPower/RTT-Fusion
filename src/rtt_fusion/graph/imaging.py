"""Coherent image of ray bundles on a detector plane (stage R15, concept section 5).

Each bundle (one coherence group) becomes a scalar field on the pixel grid,
E(u, v) = A(u, v) · exp(i (k₀ · OPL(u, v) + φ)), with
- OPL interpolated linearly between the ray hits (it is smooth where the bundle has no caustic;
  interpolating the phase itself would wrap),
- A = |a| · sqrt(I₀ · dA_source / dA_detector) from energy conservation along the ray tube: the
  ray density is taken from a triangulation of the start coordinates and its image on the detector
  (the geometric part of Fieldtracing ADR 0004; Jacobian per triangle, averaged at the vertices),
- φ the phase of the collected amplitude factor a (beam splitter, mirror).
Fields of the same coherence key add; intensities of different keys add. Outside the convex hull of
a bundle's hits its field is 0.

Validity (R15): bundles without caustics between source and detector, and an OPL that is close to
linear between neighbouring rays (exact for plane and gently curved waves; the error of the linear
interpolation grows with the wavefront curvature times the square of the ray spacing).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay, QhullError

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle


@dataclass
class Image:
    """Detector image: values[j, i] at (u_i, v_j), W/mm²; extent (u0, u1, v0, v1) in mm."""

    values: FloatArray
    extent: tuple[float, float, float, float]
    ray_power: float = 0.0  # Σ of the power carried by the rays that reached the plane (W)
    groups: list[CoherenceKey] = field(default_factory=list)

    @property
    def pixel_area(self) -> float:
        ny, nx = self.values.shape
        u0, u1, v0, v1 = self.extent
        return float((u1 - u0) * (v1 - v0) / (nx * ny))

    @property
    def power(self) -> float:
        """∫ I dA over the detector (W)."""
        return float(self.values.sum() * self.pixel_area)


def _vertex_jacobian(src: FloatArray, det: FloatArray) -> tuple[FloatArray, FloatArray, float]:
    """dA_source/dA_detector at each ray, the triangles of the source triangulation and the
    largest relative change of a triangle edge (for diagnostics)."""
    tri = Delaunay(src).simplices

    def area(p: FloatArray) -> FloatArray:
        a, b, c = p[tri[:, 0]], p[tri[:, 1]], p[tri[:, 2]]
        return np.asarray(
            0.5
            * np.abs(
                (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
            )
        )

    a_src, a_det = area(src), area(det)
    ok = a_det > 1e-30
    jac_t = np.where(ok, a_src / np.where(ok, a_det, 1.0), 0.0)
    acc = np.zeros(len(src))
    cnt = np.zeros(len(src))
    for k in range(3):
        np.add.at(acc, tri[:, k], jac_t)
        np.add.at(cnt, tri[:, k], 1.0)
    return np.where(cnt > 0, acc / np.maximum(cnt, 1.0), 0.0), np.asarray(tri), float(a_src.sum())


def coherent_image(
    bundles: list[RayBundle],
    plane: Plane,
    half_u: float,
    half_v: float,
    nx: int,
    ny: int,
) -> Image:
    """Intensity on the pixel grid of `plane` (u ∈ [−half_u, half_u], v likewise)."""
    u = (np.arange(nx) + 0.5) / nx * 2 * half_u - half_u
    v = (np.arange(ny) + 0.5) / ny * 2 * half_v - half_v
    uu, vv = np.meshgrid(u, v)
    fields: dict[CoherenceKey, ComplexArray] = {}
    ray_power = 0.0
    for b in bundles:
        moved, uv = plane.intersect(b)
        alive = np.asarray(moved.alive)
        if alive.sum() < 3:
            continue
        src = np.asarray(moved.source_xy)[alive]
        det = uv[alive]
        try:
            jac, _, area_src = _vertex_jacobian(src, det)
        except QhullError:
            continue
        amp = np.asarray(moved.amplitude)[alive]
        a = np.abs(amp) * np.sqrt(b.irradiance * jac)
        opl = np.asarray(moved.opl)[alive]
        phi0 = float(np.angle(np.mean(amp))) if np.any(np.abs(amp) > 0) else 0.0
        ray_power += float(b.irradiance * area_src * np.mean(np.abs(amp) ** 2))
        interp_a = LinearNDInterpolator(det, a, fill_value=0.0)
        interp_opl = LinearNDInterpolator(det, opl, fill_value=np.nan)
        a_px = interp_a(uu, vv)
        opl_px = interp_opl(uu, vv)
        e = np.where(
            np.isfinite(opl_px), a_px * np.exp(1j * (b.k0 * np.nan_to_num(opl_px) + phi0)), 0.0
        )
        e = np.asarray(e, dtype=complex)
        fields[b.key] = fields[b.key] + e if b.key in fields else e
    values = np.zeros((ny, nx))
    for e in fields.values():
        values = values + np.abs(e) ** 2
    return Image(values, (-half_u, half_u, -half_v, half_v), ray_power, list(fields))
