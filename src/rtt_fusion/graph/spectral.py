"""Fields (H) on a grid and as plane waves inside a block (grating ADR, stages F4b and F1).

A block that acts on the plane-wave spectrum of a field – a layer stack (diagonal in k∥), a grating
(k∥ → k∥ + n·K) – samples the incoming hybrid field in its own plane on a regular grid, applies its
operator to the spectrum and puts the result back on the field's own sample points. The
representations G and S stay inside the block (ADR section 2).

**Carrier.** The grid holds the residual R(x) = E(x)·e^(−i k_c·x) with the power-weighted mean
tangential wave vector k_c of the field: a tilted beam needs the grid only for its envelope and its
curvature, not for the tilt. Spectral point κ of the grid is the plane wave k∥ = k_c + κ.

**Sampling.** Step ≤ π / (2·max|∇ψ − k_c|) over the samples (the phase curvature left after the
carrier, with a factor 2 against aliasing) and ≤ half the median distance of neighbouring samples
(the envelope); the window is the samples' bounding box enlarged by `PAD` and rounded up to a power
of two points per axis, at most `max_points`.

**Gridded fields** (Fieldtracing's `GriddedField`, e.g. after a propagation into a focus): their own
square grid is used as it is, with the carrier k_c = 0 and the full field values·e^(i/2·pᵀCp); the
outgoing field is a `GriddedField` on the same grid with the same carrier C and new values, no
interpolation. Only plain lattices (L = I) in the space domain.

**Back to the samples.** Blocks that keep k∥ (stacks) and blocks that add n·K (gratings, per order)
keep the tangential phase: ψ, ∇ψ and its Hessian of the incident field carry over to the outgoing
one, only the residual U changes. The outgoing field is therefore a `GridlessField` on the same
points with the same ψ and the new U = E_out·e^(−iψ), E_out interpolated (cubic) from the outgoing
residual grid. Reflected light leaves into the half space it came from: its frame is (u, −v, −n).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree

from rtt_fusion._types import ComplexArray, FloatArray

if TYPE_CHECKING:
    from fieldtracing.field import GridlessField

#: the window is the samples' bounding box times this (zero padding against the FFT's wrap-around)
PAD = 1.5


class SpectralError(ValueError):
    """The field cannot be put on an admissible grid (too many points, too few samples)."""


@dataclass(frozen=True)
class Grid:
    """A centred grid in the plane of the field: centre (mm), step (mm), points per axis, carrier
    k_c (1/mm)."""

    centre: FloatArray
    step: float
    nx: int
    ny: int
    carrier: FloatArray

    def axes(self) -> tuple[FloatArray, FloatArray]:
        """x and y of the grid nodes (mm, absolute in the plane)."""
        x = self.centre[0] + (np.arange(self.nx) - self.nx // 2) * self.step
        y = self.centre[1] + (np.arange(self.ny) - self.ny // 2) * self.step
        return np.asarray(x), np.asarray(y)

    def frequencies(self) -> tuple[FloatArray, FloatArray]:
        """κ of the FFT bins (1/mm), numpy order, shape (ny, nx) each."""
        kx = 2.0 * np.pi * np.fft.fftfreq(self.nx, d=self.step)
        ky = 2.0 * np.pi * np.fft.fftfreq(self.ny, d=self.step)
        kxx, kyy = np.meshgrid(kx, ky)
        return np.asarray(kxx), np.asarray(kyy)


def _pow2(n: float) -> int:
    return int(2 ** max(3, int(np.ceil(np.log2(max(n, 8.0))))))


def choose_grid(field: GridlessField, max_points: int = 1024) -> Grid:
    """The grid for a hybrid field (module docstring); SpectralError if it would need more than
    `max_points` per axis."""
    pts = np.asarray(field.points)
    if len(pts) < 6:
        raise SpectralError("too few samples for a grid")
    u = np.asarray(field.amplitude)
    weight = np.sum(np.abs(u) ** 2, axis=1)
    weight = weight / (weight.sum() or 1.0)
    grad = np.asarray(field.gradient)
    carrier = np.asarray(weight @ grad)
    spread = float(np.max(np.linalg.norm(grad - carrier, axis=1)))
    near, _ = cKDTree(pts).query(pts, k=2)
    spacing = float(np.median(near[:, 1]))
    step = spacing / 2.0
    if spread > 0:
        step = min(step, np.pi / (2.0 * spread))
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    centre = 0.5 * (lo + hi)
    extent = np.maximum(hi - lo, step) * PAD
    nx, ny = _pow2(extent[0] / step), _pow2(extent[1] / step)
    if max(nx, ny) > max_points:
        raise SpectralError(
            f"the field needs {max(nx, ny)} grid points per axis (at most {max_points}): its "
            f"phase curvature or size is too large for the spectral path – use rays (R)"
        )
    return Grid(np.asarray(centre, dtype=float), float(step), nx, ny, carrier)


def residual_grid(field: GridlessField, grid: Grid) -> ComplexArray:
    """R = E·e^(−i k_c·x) on the grid nodes, (ny, nx, 3) local components; zero outside the
    samples."""
    x, y = grid.axes()
    xx, yy = np.meshgrid(x, y)
    q = np.stack([xx.ravel(), yy.ravel()], axis=-1)
    e = field.residual_at(q) * np.exp(1j * (field.phase_at(q) - q @ grid.carrier))[:, None]
    return np.asarray(e.reshape(grid.ny, grid.nx, 3), dtype=np.complex128)


def spectrum(residual: ComplexArray) -> ComplexArray:
    """FFT of a residual grid with the phase referred to the grid centre, (ny, nx, c)."""
    ny, nx = residual.shape[:2]
    s = np.fft.fft2(residual, axes=(0, 1)) / (nx * ny)
    return np.asarray(s * _centre_shift(nx, ny)[..., None])


def inverse(spec: ComplexArray) -> ComplexArray:
    """The inverse of `spectrum`."""
    ny, nx = spec.shape[:2]
    s = spec / _centre_shift(nx, ny)[..., None]
    return np.asarray(np.fft.ifft2(s, axes=(0, 1)) * (nx * ny))


def _centre_shift(nx: int, ny: int) -> ComplexArray:
    """e^(i κ·x₀) with x₀ the offset of node (0, 0) from the centre node, in index units."""
    fx = np.fft.fftfreq(nx) * (nx // 2)
    fy = np.fft.fftfreq(ny) * (ny // 2)
    return np.asarray(np.exp(2j * np.pi * (fx[None, :] + fy[:, None])))


def at_samples(residual_out: ComplexArray, grid: Grid, field: GridlessField) -> ComplexArray:
    """E_out at the field's samples from an outgoing residual grid (cubic), (N, 3)."""
    x, y = grid.axes()
    pts = np.asarray(field.points)
    out = np.zeros((len(pts), 3), dtype=np.complex128)
    for c in range(3):
        for part, unit in ((np.real, 1.0), (np.imag, 1j)):
            interp = RegularGridInterpolator(
                (y, x),
                part(residual_out[:, :, c]),
                method="cubic",
                bounds_error=False,
                fill_value=0.0,
            )
            out[:, c] += unit * interp(np.stack([pts[:, 1], pts[:, 0]], axis=-1))
    return np.asarray(out * np.exp(1j * (pts @ grid.carrier))[:, None])


def outgoing(
    field: GridlessField, e_out: ComplexArray, refractive_index: float, reflected: bool
) -> GridlessField:
    """The outgoing field on the incident field's samples with its ψ (module docstring); for
    reflected light in the mirrored frame (u, −v, −n)."""
    from fieldtracing.field import GridlessField

    pts = np.asarray(field.points, dtype=float)
    grad = np.asarray(field.gradient, dtype=float)
    hess = np.asarray(field.hessian, dtype=float)
    u = np.asarray(e_out) * np.exp(-1j * np.asarray(field.phase))[:, None]
    if reflected:
        flip = np.array([1.0, -1.0])
        pts, grad = pts * flip, grad * flip
        hess = hess * np.array([[1.0, -1.0], [-1.0, 1.0]])
        u = u * np.array([1.0, -1.0, -1.0])
    return GridlessField(
        pts,
        np.asarray(field.phase, dtype=float),
        grad,
        hess,
        np.asarray(u, dtype=np.complex128),
        float(field.wavelength_um),
        float(refractive_index),
    )


@dataclass(frozen=True)
class Sampled:
    """A field on its grid: the grid, the residual (ny, nx, 3) and how to make the outgoing field
    from an outgoing residual grid (`back(residual_out, refractive_index, reflected)`)."""

    grid: Grid
    residual: ComplexArray
    back: Any


def sample(field: Any, max_points: int = 1024) -> Sampled:
    """A hybrid or gridded field on its grid (module docstring)."""
    from fieldtracing.field import GridlessField

    if isinstance(field, GridlessField):
        grid = choose_grid(field, max_points)

        res_in = residual_grid(field, grid)

        def back_gridless(
            res: ComplexArray, n: float, reflected: bool, t0: ComplexArray | None = None
        ) -> Any:
            if t0 is None:
                return outgoing(field, at_samples(res, grid, field), n, reflected)
            # the carrier wave's transfer t0 (3 x 2) pointwise on the samples – exact for a plane
            # wave and free of the reconstruction's ripple at hard edges –, only the rest
            # (T(κ) − t0)·Ẽ(κ) through the grid
            e_in = (
                np.asarray(field.amplitude)[:, :2] * np.exp(1j * np.asarray(field.phase))[:, None]
            )
            rest = res - np.einsum("ij,yxj->yxi", t0, res_in[..., :2])
            e_out = e_in @ t0.T + at_samples(rest, grid, field)
            return outgoing(field, e_out, n, reflected)

        return Sampled(grid, res_in, back_gridless)
    return _sample_gridded(field, max_points)


def _sample_gridded(field: Any, max_points: int) -> Sampled:
    from fieldtracing.field import GriddedField

    if field.domain != "space" or not np.allclose(np.asarray(field.lattice), np.eye(2)):
        raise SpectralError("only gridded fields in space on a plain lattice can be sampled")
    n = int(np.asarray(field.values).shape[0])
    if n > max_points:
        raise SpectralError(f"the gridded field has {n} points per axis (at most {max_points})")
    try:
        full = np.asarray(field.physical(), dtype=np.complex128)
    except ValueError as e:
        raise SpectralError(f"the gridded field's carrier is not resolved: {e}") from e
    grid = Grid(np.zeros(2), float(field.step), n, n, np.zeros(2))
    x, y = grid.axes()
    xx, yy = np.meshgrid(x, y)
    c = np.asarray(field.curvature, dtype=float)
    carrier = np.exp(0.5j * (c[0, 0] * xx**2 + (c[0, 1] + c[1, 0]) * xx * yy + c[1, 1] * yy**2))

    def back_gridded(
        res: ComplexArray, refractive_index: float, reflected: bool, t0: Any = None
    ) -> Any:  # on its own nodes the inverse FFT is exact: t0 is not needed
        values = res / carrier[..., None]
        curvature = c
        if reflected:  # frame (u, −v, −n): y -> −y on the centred grid
            values = np.roll(values[::-1], 1, axis=0) * np.array([1.0, -1.0, -1.0])
            if n % 2:
                values = np.roll(values, -1, axis=0)
            curvature = c * np.array([[1.0, -1.0], [-1.0, 1.0]])
        return GriddedField(
            float(field.step),
            np.asarray(values, dtype=np.complex128),
            float(field.wavelength_um),
            float(refractive_index),
            np.asarray(curvature, dtype=float),
        )

    return Sampled(grid, full, back_gridded)


def significant(spec: ComplexArray, threshold: float) -> Any:
    """Mask of the spectral points whose amplitude exceeds `threshold` times the largest."""
    mag = np.max(np.abs(spec), axis=-1)
    return mag > threshold * float(mag.max() or 1.0)
