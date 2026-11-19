"""Light on the links of the executed graph (stage R15): representation R as ray bundles.

A `RayBundle` is one coherence group (concept section 5): all its rays come from one source point
or mode, one wavelength and one polarisation state of the source, identified by `key`. Rays carry
position and direction (global, mm), the optical path length OPL (mm, vacuum-equivalent, as
Raytatouille's `opl`), a complex amplitude factor collected at interactions (beam splitter r/t,
mirror reflectance; scalar in R15, the vector form comes with the polarised tracer) and their
start coordinates in the source plane, from which the ray density – and with it the amplitude –
follows at a detector (energy conservation along the ray tube, Fieldtracing ADR 0004 for the
geometric part).

Amplitudes are normalised to the source: `irradiance` is the irradiance at the source plane in
W/mm² (uniform over the bundle), so |E|² of the detector field is W/mm².

Polarisation (R18): `jones` is the field direction of the source (global, transverse to the start
direction) and `prt` the accumulated polarisation matrix of each ray (global, complex 3×3;
None = identity): mirrors multiply it with I − 2nnᵀ, an aplanatic lens with its meridional
rotation, interfaces, gratings and Raytatouille paths with their Jones/PRT matrices. The field of a
ray is a·prt·jones, a the scalar `amplitude` (beam splitter r/t, mirror reflectance).

Complex polarisation (R31, docs/plan-r31-polarisation.md): ê = prt·jones is a complex vector, so
elliptical polarisation survives every conversion. Where a block splits the field f of a ray into
a scalar and the rest, it takes f = a·ê with a complex unit vector ê and |a| = |f|
(`split_field`); for linear polarisation ê is real and a is what it was before R31. The ray
image (`imaging.coherent_image`, vectorial since R32) and the field representation H
(`graph.field`) use a·prt·jones.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np
from numpy.typing import NDArray

from rtt_fusion._types import ComplexArray, FloatArray

BoolArray = NDArray[np.bool_]


@dataclass(frozen=True)
class CoherenceKey:
    """Identity of a coherence group: amplitudes add within one key, intensities across keys."""

    source: str
    field_index: int = 0
    wavelength_index: int = 0
    polarisation: int = 0


@dataclass(frozen=True)
class RayBundle:
    """Rays of one coherence group.

    Attributes:
        key: coherence group.
        wavelength_um: vacuum wavelength.
        position: (N, 3) mm, global.
        direction: (N, 3) unit vectors.
        opl: (N,) optical path length from the source plane, mm.
        amplitude: (N,) complex amplitude factor (1 at the source).
        source_xy: (N, 2) start coordinates in the source plane, mm (for the ray density).
        alive: (N,) bool, False for rays lost on the way (outside an aperture, missed).
        irradiance: source irradiance in W/mm².
        path: names of the blocks passed (for diagnostics).
        jones: (3,) complex field direction of the source (global; None: x̂ projected).
        prt: (N, 3, 3) complex polarisation matrix (global; None: identity); prt·jones is the
            complex polarisation of each ray (R31). Real arrays (port data written before R31)
            are taken as complex with imaginary part 0.
        source_cos: cos γ of the rays on the source plane (1 for a plane wave ⟂ the plane).
        source_frame: (2, 3) the axes of the `source_xy` coordinates carried along the rays
            geometrically (mirror: reflected, refraction: smallest rotation k_in -> k_out), global,
            at the present position. It removes the rotation and mirroring of the ray map that
            come from the geometry of the path, so that the Maslov index counts caustics only.
        index, source_index: refractive indices of the present medium and at the source plane
            (absolute; the OPL is vacuum-equivalent).
    """

    key: CoherenceKey
    wavelength_um: float
    position: FloatArray
    direction: FloatArray
    opl: FloatArray
    amplitude: ComplexArray
    source_xy: FloatArray
    alive: BoolArray
    irradiance: float = 1.0
    path: tuple[str, ...] = field(default=())
    jones: ComplexArray | None = None
    prt: ComplexArray | None = None
    source_cos: FloatArray | float = 1.0  # cos γ between the rays and the source plane
    source_frame: FloatArray | None = None  # (2, 3) source axes carried along (global, now)
    index: float = 1.0  # refractive index of the medium the rays are in now
    source_index: float = 1.0  # refractive index at the source plane
    #: the branch the light follows and the position of its last Raytatouille block in
    #: `graph.meta["paths"][branch]` (-1: just left the source) – node-wise execution (R25)
    branch: tuple[str, int] | None = None

    def __post_init__(self) -> None:
        # real polarisation data (older port data, callers with float arrays) as complex128
        for name in ("jones", "prt"):
            value = getattr(self, name)
            if value is not None and not np.iscomplexobj(value):
                object.__setattr__(self, name, np.asarray(value, dtype=np.complex128))

    @property
    def size(self) -> int:
        return int(np.asarray(self.opl).shape[0])

    @property
    def k0(self) -> float:
        """Vacuum wavenumber, 1/mm."""
        return float(2.0 * np.pi / (self.wavelength_um * 1e-3))

    def field_direction(self) -> ComplexArray:
        """(N, 3) prt·jones: the polarisation of each ray (unit for unit jones)."""
        j = (
            np.asarray(self.jones, dtype=np.complex128)
            if self.jones is not None
            else np.array([1.0, 0.0, 0.0], dtype=np.complex128)
        )
        if self.prt is None:
            return np.asarray(np.tile(j, (self.size, 1)))
        return np.asarray(np.einsum("nij,j->ni", np.asarray(self.prt), j))

    def field_vectors(self) -> ComplexArray:
        """(N, 3) a·prt·jones: the complex field of each ray relative to the source (R31)."""
        return np.asarray(np.asarray(self.amplitude)[:, None] * self.field_direction())

    def evolve(self, **changes: object) -> RayBundle:
        return replace(self, **changes)  # type: ignore[arg-type]


@dataclass(frozen=True)
class Plane:
    """A plane with a local frame: origin, unit normal, in-plane axes u and v (mm)."""

    origin: FloatArray
    normal: FloatArray
    u: FloatArray
    v: FloatArray

    @staticmethod
    def make(origin: object, normal: object, up: object = (0.0, 1.0, 0.0)) -> Plane:
        o = np.asarray(origin, dtype=float)
        n = np.asarray(normal, dtype=float)
        n = n / np.linalg.norm(n)
        upv = np.asarray(up, dtype=float)
        if abs(float(np.dot(upv, n))) > 0.99:
            upv = np.array([1.0, 0.0, 0.0])
        u = np.cross(upv, n)
        u = u / np.linalg.norm(u)
        v = np.cross(n, u)
        return Plane(o, n, u, v)

    def intersect(
        self, bundle: RayBundle, index: float | None = None
    ) -> tuple[RayBundle, FloatArray]:
        """Moves the rays to the plane (OPL grows by n·s, n the medium the rays are in:
        `bundle.index` unless given); rays that do not reach it forward are lost. Returns the
        moved bundle and the local (u, v) coordinates of the hits."""
        if index is None:
            index = bundle.index
        d = np.asarray(bundle.direction)
        p = np.asarray(bundle.position)
        denom = d @ self.normal
        with np.errstate(divide="ignore", invalid="ignore"):
            s = ((self.origin - p) @ self.normal) / denom
        ok = np.asarray(bundle.alive) & np.isfinite(s) & (s > -1e-12) & (np.abs(denom) > 1e-14)
        s = np.where(ok, s, 0.0)
        hit = p + s[:, None] * d
        moved = bundle.evolve(position=hit, opl=np.asarray(bundle.opl) + index * s, alive=ok)
        rel = hit - self.origin
        return moved, np.stack([rel @ self.u, rel @ self.v], axis=-1)


#: |f·f| below this share of |f|² counts as circular polarisation in `split_field`
CIRCULAR = 1e-9


def split_field(f: Any) -> tuple[ComplexArray, ComplexArray]:
    """Field vectors f (N, 3) as f = a·ê (R31): a (N,) complex with |a| = |f|, ê (N, 3) a complex
    unit vector (ê†ê = 1, ê = 0 where f = 0). The phase of a is that of the field along the major
    axis of its polarisation ellipse, ½ arg(f·f), with the sign that makes the largest component
    of Re ê positive (taken at the largest component of f where Re ê has one); for (nearly)
    circular polarisation the phase of the largest component of f. For linear polarisation ê is
    real and a the complex amplitude along it, as before R31."""
    f = np.asarray(f, dtype=np.complex128)
    mag2 = f.real**2 + f.imag**2
    size = np.sqrt(mag2.sum(axis=1))
    rows = np.arange(len(f))
    big = np.argmax(mag2, axis=1)
    square = np.sum(f * f, axis=1)
    at_big = f[rows, big]
    # the phase factor e^(iφ): sqrt of the unit f·f (½ arg), or of the largest component
    circular = np.abs(square) <= CIRCULAR * size**2
    unit = np.where(
        circular,
        at_big / np.where(circular & (size > 0), np.abs(at_big), 1.0),
        np.sqrt(square / np.where(circular, 1.0, np.abs(square))),
    )
    unit = np.where(size > 0, unit, 1.0)
    major = np.real(f * np.conj(unit)[:, None])
    # the sign from the largest component of f (as before R31, so that linear polarisation keeps
    # its sign also where two components are equally large), else from the largest of Re ê
    lead = major[rows, big]
    weak = np.abs(lead) <= 1e-6 * np.sqrt(np.sum(major * major, axis=1))
    if np.any(weak):
        lead = np.where(weak, major[rows, np.argmax(np.abs(major), axis=1)], lead)
    unit = np.where(lead < 0, -unit, unit)
    a = size * unit
    e = f / np.where(size > 0, a, 1.0)[:, None]
    return np.asarray(a, dtype=np.complex128), np.asarray(e, dtype=np.complex128)


def split_along(f: Any, carried: Any) -> ComplexArray:
    """The scalar factor c of new field vectors f (N, 3) against the geometrically carried
    polarisation `carried` (N, 3) (R31): |c| = |f|/|carried|, the phase that of carried†·f; where
    f is (nearly) orthogonal to `carried` the phase of `split_field`. For f = λ·carried it is λ,
    the factor of the interfaces before R31."""
    f = np.asarray(f, dtype=np.complex128)
    g = np.asarray(carried, dtype=np.complex128)
    size_f = np.linalg.norm(f, axis=1)
    size_g = np.linalg.norm(g, axis=1)
    overlap = np.sum(np.conj(g) * f, axis=1)
    a, _ = split_field(f)
    use = np.abs(overlap) > 1e-12 * size_f * size_g
    unit = np.where(
        use, overlap / np.where(use, np.abs(overlap), 1.0), a / np.where(size_f > 0, size_f, 1.0)
    )
    c = np.where(size_g > 0, size_f / np.where(size_g > 0, size_g, 1.0), 0.0) * unit
    return np.asarray(c, dtype=np.complex128)


def physical_split(
    physical: ComplexArray, geometric: FloatArray, b: RayBundle
) -> tuple[ComplexArray, ComplexArray]:
    """The new field P·prt·jones of the rays of `b` behind an element with the physical matrices
    P (N, 3, 3) and their geometric part G (R31): (c, prt_out) with the scalar factor c
    (`split_along`: |P·ê|/|G·ê| with the phase of the projection onto the geometrically
    carried polarisation G·ê) and prt_out = P·prt/c, so that a·c·prt_out·jones is the full field
    – elliptical polarisation and mixed s/p light included. For s or p light, and at normal
    incidence, c and prt_out·jones are the factor and G·prt·jones of R27."""
    field = b.field_direction()
    carried = np.einsum("nij,nj->ni", geometric, field)
    c = split_along(np.einsum("nij,nj->ni", physical, field), carried)
    prt = (
        np.asarray(b.prt)
        if b.prt is not None
        else np.tile(np.eye(3, dtype=complex), (len(field), 1, 1))
    )
    safe = np.where(np.abs(c) > 1e-300, c, 1.0)
    prt_out = np.einsum("nij,njk->nik", physical, prt) / safe[:, None, None]
    # rays without a physical matrix (grazing incidence, 0/0): no field, as before R31
    bad = ~np.isfinite(c) | ~np.all(np.isfinite(prt_out), axis=(1, 2))
    if np.any(bad):
        c = np.where(bad, 0.0, c)
        prt_out[bad] = 0.0
    return c, np.asarray(prt_out, dtype=np.complex128)


def rotation_between(a: FloatArray, b: FloatArray) -> FloatArray:
    """Smallest rotation taking the unit vector a to b (Rodrigues); −I for opposite vectors."""
    v = np.cross(a, b)
    c = float(a @ b)
    s = float(np.linalg.norm(v))
    if s < 1e-15:
        return np.eye(3) if c > 0 else -np.eye(3)
    k = v / s
    kx = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    return np.asarray(np.eye(3) + s * kx + (1 - c) * (kx @ kx))


def concat(bundles: list[RayBundle]) -> list[RayBundle]:
    """Bundles on a link: kept separate per coherence group (a link may carry several)."""
    return list(bundles)
