"""Coherent image of ray bundles on a detector plane (stage R15, concept section 5; vectorial
since R32).

Each bundle (one coherence group) becomes a vector field on the pixel grid,
E(u, v) = A(u, v) · exp(i k₀ · OPL(u, v)), with
- OPL interpolated linearly between the ray hits (it is smooth where the bundle has no caustic;
  interpolating the phase itself would wrap),
- A the complex field vector a·prt·jones of the rays (`RayBundle.field_vectors`, R31: amplitude
  factors of beam splitters and mirrors, Fresnel and grating matrices, elliptical polarisation)
  times sqrt(I₀ · dA_source / dA_detector) from energy conservation along the ray tube: the ray
  density is taken from a triangulation of the start coordinates and its image on the detector
  (the geometric part of Fieldtracing ADR 0004; Jacobian per triangle, averaged at the vertices);
  A is interpolated linearly as well (it holds no propagation phase, so it is smooth).
Fields of the same coherence key add as vectors, the image is Σ_keys |E|²: light of orthogonal
polarisations does not interfere, as in the field (H) and spectrum (S) images. Before R32 the
image was scalar (|a| with the mean phase of a), so crossed polarisations made fringes. Outside the
convex hull of a bundle's hits its field is 0.

Validity (R15): bundles without caustics between source and detector, and an OPL that is close to
linear between neighbouring rays (exact for plane and gently curved waves; the error of the linear
interpolation grows with the wavefront curvature times the square of the ray spacing). The field
vector is interpolated linearly as well: where the polarisation turns between neighbouring rays by
Δ (curved surfaces, high NA), |E| between them is low by about Δ²/8 (1.7e-5 of the peak for the
ghost of a singlet, W52, with 0.5° between rays); denser rays reduce it.
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
    #: a warning for the user, empty if none: the light covers too few pixel centres – the field
    #: is sampled there, so a spot smaller than a few pixels is missed (partly or wholly)
    note: str = ""

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
    covered: list[int] = []  # pixel centres inside the footprint of each bundle
    spot = [np.inf, -np.inf, np.inf, -np.inf]  # extent of all hits (u0, u1, v0, v1)
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
        try:
            mesh = Delaunay(det)
        except QhullError:
            continue
        vec = np.asarray(moved.field_vectors())[alive]  # a·prt·jones (R31)
        a = vec * np.sqrt(b.irradiance * jac)[:, None]
        opl = np.asarray(moved.opl)[alive]
        ray_power += float(b.irradiance * area_src * np.mean(np.sum(np.abs(vec) ** 2, axis=1)))
        a_px = LinearNDInterpolator(mesh, a, fill_value=0.0)(uu, vv)  # (ny, nx, 3)
        opl_px = LinearNDInterpolator(mesh, opl, fill_value=np.nan)(uu, vv)
        covered.append(int(np.isfinite(opl_px).sum()))
        spot = [
            min(spot[0], float(det[:, 0].min())),
            max(spot[1], float(det[:, 0].max())),
            min(spot[2], float(det[:, 1].min())),
            max(spot[3], float(det[:, 1].max())),
        ]
        inside = np.isfinite(opl_px)
        phase = np.exp(1j * b.k0 * np.nan_to_num(opl_px))
        e = np.asarray(np.where(inside[..., None], a_px * phase[..., None], 0.0), dtype=complex)
        fields[b.key] = fields[b.key] + e if b.key in fields else e
    values = np.zeros((ny, nx))
    for e in fields.values():
        values = values + np.sum(np.abs(e) ** 2, axis=-1)
    note = _undersampled(covered, spot, half_u, half_v, nx, ny)
    return Image(values, (-half_u, half_u, -half_v, half_v), ray_power, list(fields), note)


#: fewer pixel centres than this inside a bundle's footprint: the image misses it (partly)
MIN_COVERED_PIXELS = 4


def _undersampled(
    covered: list[int], spot: list[float], half_u: float, half_v: float, nx: int, ny: int
) -> str:
    """The warning for light that covers too few pixel centres, with a window that fits it.
    Geometric on purpose: the power in the image may differ from the arriving power by
    interference (a dark port), which is no reason to warn."""
    if not covered or min(covered) >= MIN_COVERED_PIXELS:
        return ""
    width, height = spot[1] - spot[0], spot[3] - spot[2]
    reach = max(abs(spot[0]), abs(spot[1]), abs(spot[2]), abs(spot[3]))
    pixel = max(2 * half_u / nx, 2 * half_v / ny)
    window = 1.2 * reach  # centred on the detector, as the image is
    message = (
        f"the light covers only {min(covered)} pixel centre(s) of the {nx}×{ny} image "
        f"({pixel:.3g} mm per pixel; the spot is {width:.3g} × {height:.3g} mm), so the image "
        "misses (part of) it"
    )
    if reach > max(half_u, half_v):
        return (
            message
            + " – it lies (partly) outside the window: set the detector's half size to about "
            + f"{window:.2g} mm"
        )
    if window < max(half_u, half_v):
        return message + f" – set the detector's half size to about {window:.2g} mm"
    return message + " – use more pixels"
