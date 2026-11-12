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
direction) and `prt` the accumulated geometric polarisation matrix of each ray (global, 3×3;
None = identity): mirrors multiply it with I − 2nnᵀ, an aplanatic lens with its meridional
rotation. The field of a ray is a·prt·jones, a the scalar `amplitude` (beam splitter r/t, mirror
reflectance). The scalar image of R15 (`imaging.coherent_image`) does not use them; the field
representation H does (`graph.field`).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

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
        prt: (N, 3, 3) geometric polarisation matrix (global; None: identity).
        source_cos: cos γ of the rays on the source plane (1 for a plane wave ⟂ the plane).
        source_frame: (2, 3) global axes of the `source_xy` coordinates; with `prt` it gives the
            orientation of the ray map that comes from reflections (not from caustics).
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
    source_frame: FloatArray | None = None  # (2, 3) axes u, v of `source_xy` (global)

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

    def intersect(self, bundle: RayBundle, index: float = 1.0) -> tuple[RayBundle, FloatArray]:
        """Moves the rays to the plane (OPL grows by n·s); rays that do not reach it forward are
        lost. Returns the moved bundle and the local (u, v) coordinates of the hits."""
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


def concat(bundles: list[RayBundle]) -> list[RayBundle]:
    """Bundles on a link: kept separate per coherence group (a link may carry several)."""
    return list(bundles)
