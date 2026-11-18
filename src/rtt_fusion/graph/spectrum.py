"""Representation S on links: plane-wave spectra (stage R30a, docs/plan-s-spektren.md).

A `SpectrumBundle` is one coherence group as plane waves (`ports.PlaneWaves`) on a port plane:
wave vectors k∥ = (k_u, k_v) and transverse amplitudes (E_u, E_v) in the plane's frame (u, v, n),
n the propagation axis, every wave with phase 1 at the plane's origin (docs/konventionen.md);
E_n follows from transversality. Two kinds, told apart by the fields of the waves (decision S1):

- **continuous** (`solid_angle` set): a sampled angular spectrum of a bounded beam,
  E(r) = ∫ A(k̂) e^(ik k̂·r) dΩ with field_j = A(k̂_j)·dΩ_j. With the Fourier spectrum of the
  field in the plane, E(ρ) = ∫ Ẽ(k∥) e^(ik∥·ρ) d²k∥, this is field_j = Ẽ(k_j)·Δ²k_j and
  dΩ_j = Δ²k_j/(k·k_z,j) (d²k∥ = k k_z dΩ).
- **discrete** (`solid_angle` None, e.g. diffraction orders): plane waves without extent.

Operators (fusion plan, table 6.3; plan S3):

- `field_to_spectrum` (H → S): the pointwise Fourier transform of Fieldtracing
  (`GridlessField.pft`, stationary phase: k_j = ∇ψ(ρ_j)) where the local Fresnel numbers are
  large, otherwise an FFT of the field on a grid (the grid of `graph.spectral`: residual with the
  carrier k_c, interpolated cubically, Ẽ(k_c + κ) from the FFT); a gridded field of Fieldtracing
  (after a propagation into a focus) is representation G and goes through `ports.to_plane_waves`
  (G → S). The PFT's quadrature weight Δ²k_j is a third of the adjacent Delaunay triangles of
  the wave vectors (the quadrature `fusion.particle` uses on the direction cosines): the spectrum
  is integrated where it lives, and a focused beam whose rays are uniform in the pupil, hence in
  k∥, gets the trapezoid rule there (the field's own quadrature in the port plane,
  `field.field_power`, is coarse at the rim of a converging bundle: about 1 % low for NA 0.8).
  Evanescent waves (|k∥| ≥ n k₀) are left out; a dropped share above `CUT_WARNING` of the power
  is reported as `exec.spectrum_cut`.
- `spectrum_to_field` (S → H): the inverse PFT onto the port plane, or the inverse FFT back on
  the field's own samples (with its phase ψ, as `graph.spectral` does for stacks and gratings),
  or back onto the gridded field's grid. The bundle keeps the sampling it came from (`support`).
- `rays_to_spectrum` (R → S): always through H (decision 2; `field.rays_to_field` then
  `field_to_spectrum`). The Debye way per ray stays inside `fusion.particle` for now.
- `spectrum_to_rays` (S → R): S → H → R (`field.field_to_rays`); `exec.spectrum_rays` if the
  bundle is not large against λ (rays from a beam of a few wavelengths say little).

Free propagation and image (stage R30b, plan S4):

- `propagate` (S4): the spectrum on another port plane. Parallel planes: every wave gets the phase
  e^(i k·Δr) of the shift (evanescent waves decay; behind the plane they would grow and are left
  out), an in-plane turn of the frame turns k∥ and E. Tilted planes: every wave is turned into the
  new frame (u, v, n of the target, n along the light); evanescent waves and waves that do not
  cross the new plane forward are left out with `exec.spectrum_cut`. The weights dΩ belong to the
  directions and stay. The spectral phase of a PFT spectrum moves with it (ψ̃ + k·Δr, its gradient
  −ρ by the ray shift −(k∥/k_n)Δn, its Hessian by −Δn ∂²k_n/∂k∥²), so the inverse PFT still
  gives the field on the new plane. A grid spectrum keeps its bins on a grid that moves with the
  light (R30c), with samples carried there geometrically from the old plane, as long as the beam
  stays inside the grid's window; tilted planes lose the sampling (`support` None). Like fields,
  spectra are not moved by the executor: a block that takes S natively moves them to its port
  plane (the detector evaluates them in space directly).

Blocks with S (stage R30c): gratings, layer stacks and condensed groups take spectra natively
(`graph.grating`, `graph.stack`, `graph.condense`); `order_support` carries the sampling through
a map k∥ -> F·k∥ + n·K, `overlapping_orders` is the `exec.spectrum_rays` check of orders that
are not separated (S -> R). `GridSupport` is plain data, so spectra travel as port data to
blocks of other processes.
- `spectrum_image`: the detector's coherent image from spectra (and fields). The field of the
  waves at the pixel centres in space, E(r) = Σ_j E_j e^(i k_j·(r − o)), with no plane in between:
  exact for discrete waves (diffraction orders, interference of plane waves) and for a tilted
  detector, and for a continuous spectrum the quadrature of the angular-spectrum integral – exact
  free propagation without Fresnel or Debye approximation, as long as the waves are dense enough
  for the window. That is checked per wave j from where its ray lands, ρ_j (−∇ψ̃ moved along k_j,
  or the field's extent on the grid of an FFT spectrum): the phase of the integrand changes
  between neighbouring waves by |r − ρ_j|·δk_j, which must stay below `GRID_STEP_LIMIT` for a
  grid spectrum (beyond it the replicas of the periodic sum enter the window) and
  `SCATTERED_STEP_LIMIT` for the scattered waves of a PFT. Near a focus this always holds (all
  rays land at the focus: the sum is the Debye integral without its approximation); far from a
  focus – a large beam on a window far away – it does not, and the spectrum goes back to a field
  on its plane (S → H, the inverse transform) and through Fieldtracing's propagation planner as
  H does (`field.field_sums`), with a line in the log. Cost: one complex exponential per pixel
  and wave, in blocks of `SUM_BLOCK` (a PFT spectrum of 2800 waves on 256² pixels: 1.8·10⁸
  terms, a few seconds); a grid spectrum on a detector parallel to its plane with the same axes
  (up to sign and order) factorises into two matrix products over the grid's rows and columns
  (`grid_sum`), and a spectrum whose sum would exceed `SUM_BUDGET` terms goes through H as well.

Power: `spectrum_power` in the units of `field.field_power` (W): the flux of `ports` divided by
k₀ (`ports.spectrum_power` counts Re k_z |E|², `field_power` n cos γ |U|²).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.spatial import Delaunay, QhullError

from rtt_fusion import ports
from rtt_fusion._types import ComplexArray, FloatArray, IntArray
from rtt_fusion.graph import spectral
from rtt_fusion.graph.field import (
    FieldBundle,
    FieldError,
    field_power,
    field_sums,
    field_to_rays,
    image_of,
    pixel_centres,
    rays_to_field,
)
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle

if TYPE_CHECKING:
    from fieldtracing.field import GridlessField

#: waves of a grid spectrum whose amplitude is below this share of the largest are left out
SPECTRUM_THRESHOLD = 1e-8
#: share of the power dropped from a spectrum above which `exec.spectrum_cut` is reported
CUT_WARNING = 1e-6
#: waves closer to grazing than this (|k∥|/k ≥ 1 − GRAZING) count as evanescent
GRAZING = 1e-9
#: Fieldtracing's default ε of its selection rules: the PFT is used where 1/F ≤ ε
DEFAULT_TOLERANCE = 0.02
#: below this local Fresnel number the PFT is not used even when no grid fits
MIN_FRESNEL = 10.0
#: S -> R: a bundle narrower than this many wavelengths (in the medium) gets `exec.spectrum_rays`
BUNDLE_MIN_WAVELENGTHS = 10.0
#: image by the sum of the waves: largest phase step (rad) of the integrand between neighbouring
#: waves over the window – a grid spectrum is periodic with the window L of its grid (Poisson:
#: replicas L apart, 2π per bin), the scattered waves of a PFT get half of that
GRID_STEP_LIMIT = 2.0 * np.pi
SCATTERED_STEP_LIMIT = np.pi
#: waves whose Fourier amplitude |Ẽ| is below this share of the largest do not count for the
#: sampling check of the sum (their replicas stay below it)
SIGNIFICANT = 1e-4
#: complex numbers per block of the direct sum (pixels × waves), about 64 MB
SUM_BLOCK = 4_000_000
#: waves × pixels above which a spectrum with a sampling goes through H instead of the direct
#: sum (a few seconds); a grid spectrum on a parallel detector factorises and has no limit
SUM_BUDGET = 500_000_000


@dataclass(frozen=True)
class GridSupport:
    """How a grid spectrum was sampled, as plain data (port data carries it across a process
    boundary, R30c): the grid in the bundle's plane (centre, step, size, carrier k_c), the grid
    bin (row, column) of every wave, so that S → H puts the waves back on the same grid, and
    the samples the field comes back to (`field`): a `GridlessField` whose points, phase ψ,
    gradient and Hessian give the geometry of the field (its amplitude is not used; the field is
    E·e^(−iψ) interpolated from the grid, as `spectral.outgoing`), or a `GriddedField` whose grid
    and curvature C do (values E·e^(−i/2 xᵀCx) on its own nodes)."""

    grid: spectral.Grid
    bins: IntArray
    field: Any

    @property
    def sampled(self) -> spectral.Sampled:
        """A view as `spectral.Sampled` (grid and the way back; the residual is not kept)."""
        empty = np.zeros((0, 0, 3), dtype=np.complex128)
        return spectral.Sampled(self.grid, empty, lambda r, n, _: self.to_field(r, n))

    def to_field(self, residual: ComplexArray, refractive_index: float) -> Any:
        """The field object from a residual grid (ny, nx, 3) of this support. SpectralError if a
        gridded field's grid does not resolve the carrier."""
        from fieldtracing.field import GridlessField

        g, f = self.grid, self.field
        if isinstance(f, GridlessField):
            e = spectral.at_samples(residual, g, f)
            return spectral.outgoing(f, e, refractive_index, False)
        from fieldtracing.field import GriddedField

        if float(np.max(np.abs(g.carrier))) * g.step > np.pi:
            raise spectral.SpectralError(
                "the gridded field's grid does not resolve the spectrum's carrier "
                f"(|k_c|·step = {float(np.max(np.abs(g.carrier))) * g.step:.2f} > π)"
            )
        x, y = g.axes()
        xx, yy = np.meshgrid(x, y)
        c = np.asarray(f.curvature, dtype=float)
        quad = 0.5 * (c[0, 0] * xx**2 + (c[0, 1] + c[1, 0]) * xx * yy + c[1, 1] * yy**2)
        values = residual * np.exp(1j * (g.carrier[0] * xx + g.carrier[1] * yy - quad))[..., None]
        return GriddedField(
            float(f.step),
            np.asarray(values, dtype=np.complex128),
            float(f.wavelength_um),
            float(refractive_index),
            c,
        )


@dataclass(frozen=True)
class SpectrumBundle:
    """One coherence group as plane waves on a port plane (module docstring).

    Attributes:
        key: coherence group.
        plane: the port plane; phase reference at its origin, frame (u, v, n).
        waves: the plane waves (port z = 0, direction +1, medium and wavelength of the light).
        path: names of the blocks passed (diagnostics).
        branch: branch and position of the last Raytatouille block, as for rays and fields.
        support: the sampling the spectrum came from, for S → H: the frequency-domain
            `GridlessField` of the PFT (spectral phase ψ̃, its gradient −ρ and Hessian; one
            sample per wave) or a `GridSupport`; None for spectra without one (discrete waves).
    """

    key: CoherenceKey
    plane: Plane
    waves: ports.PlaneWaves
    path: tuple[str, ...] = ()
    branch: tuple[str, int] | None = None
    support: Any = None

    @property
    def wavelength_um(self) -> float:
        return float(self.waves.port.wavelength_um)

    @property
    def k0(self) -> float:
        return float(self.waves.port.k0)

    @property
    def refractive_index(self) -> float:
        return float(np.real(self.waves.port.index))

    @property
    def continuous(self) -> bool:
        """True for a sampled angular spectrum (weights dΩ), False for discrete waves."""
        return self.waves.solid_angle is not None


def spectrum_power(sb: SpectrumBundle) -> float:
    """Power of a continuous spectrum in the units of `field.field_power` (W)."""
    return ports.spectrum_power(sb.waves) / sb.k0


def _port(wavelength_um: float, n: float) -> ports.Port:
    return ports.Port(0.0, 1, complex(n), float(wavelength_um))


def _vertex_areas(points: FloatArray) -> FloatArray:
    """A third of the areas of the adjacent Delaunay triangles per point."""
    try:
        tri = Delaunay(points).simplices
    except QhullError as e:
        raise FieldError(f"the field's samples span no area: {e}") from e
    a, b, c = points[tri[:, 0]], points[tri[:, 1]], points[tri[:, 2]]
    area = 0.5 * np.abs(
        (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    )
    out = np.zeros(len(points))
    for j in range(3):
        np.add.at(out, tri[:, j], area / 3.0)
    return out


def beam_size(points: FloatArray, weight: FloatArray) -> float:
    """√2 times the power-weighted RMS distance of the samples from their centroid (mm): the
    radius of a uniform disc, the waist w of a Gaussian beam – the length L_U over which the
    amplitude varies."""
    w = np.asarray(weight, dtype=float)
    total = float(w.sum())
    if total <= 0:
        return 0.0
    p = np.asarray(points, dtype=float)
    centre = (w @ p) / total
    return float(np.sqrt(2.0 * (w @ np.sum((p - centre) ** 2, axis=1)) / total))


def fresnel_number(f: GridlessField) -> float:
    """Smallest local Fresnel number |h|·L_U² of a hybrid field (Fieldtracing, paper 7.1): h the
    smallest |eigenvalue| of the phase Hessians, L_U the beam size (`beam_size`); 0 where a
    Hessian is singular or not finite."""
    hessian = np.asarray(f.hessian, dtype=float)
    if not np.all(np.isfinite(hessian)):
        return 0.0
    eigen = np.abs(np.linalg.eigvalsh(0.5 * (hessian + np.swapaxes(hessian, 1, 2))))
    weight = np.sum(np.abs(np.asarray(f.amplitude)) ** 2, axis=1)
    size = beam_size(np.asarray(f.points), weight)
    return float(np.min(eigen)) * size * size


def _report_cut(notes: list[str] | None, share: float, what: str) -> None:
    if notes is not None and share > CUT_WARNING:
        notes.append(
            f"exec.spectrum_cut: {share:.2e} of the power left out of the spectrum ({what})"
        )


def field_to_spectrum(
    fb: FieldBundle,
    method: str = "auto",
    tolerance: float | None = None,
    notes: list[str] | None = None,
    max_points: int = 1024,
) -> SpectrumBundle:
    """H → S on the field's plane.

    Args:
        fb: the field (a space-domain `GridlessField`, or a `GriddedField` = representation G).
        method: "pft", "fft" or "auto": the PFT where the smallest local Fresnel number F is at
            least 1/ε (its error is O(1/F)), else the FFT on a grid; if no grid fits
            (`spectral.SpectralError`) the PFT down to F = MIN_FRESNEL, with a note of its error.
        tolerance: ε (None: Fieldtracing's default 0.02).
        notes: receives the diagnostics (`exec.spectrum_cut`) and the method chosen.
        max_points: grid points per axis at most (FFT).
    """
    from fieldtracing.field import GridlessField

    f = fb.field
    if not isinstance(f, GridlessField):
        return _gridded_to_spectrum(fb, notes, max_points)
    if f.domain != "space":
        raise FieldError("H -> S needs a field in the space domain")
    eps = float(tolerance) if tolerance else DEFAULT_TOLERANCE
    fresnel = fresnel_number(f)
    if method == "auto":
        method = "pft" if fresnel * eps >= 1.0 else "fft"
        if method == "fft":
            try:
                return _fft_spectrum(fb, f, notes, max_points)
            except spectral.SpectralError as e:
                if fresnel < MIN_FRESNEL:
                    raise FieldError(f"no spectrum for the field: {e}") from e
                if notes is not None:
                    notes.append(
                        f"{fb.key.source}: H -> S by the PFT at a local Fresnel number "
                        f"{fresnel:.3g} (estimated error {1.0 / fresnel:.2g}; no grid fits: {e})"
                    )
                method = "pft"
    if method == "pft":
        return _pft_spectrum(fb, f, notes)
    if method == "fft":
        try:
            return _fft_spectrum(fb, f, notes, max_points)
        except spectral.SpectralError as e:
            raise FieldError(f"no grid spectrum for the field: {e}") from e
    raise ValueError(f"unknown method {method!r} (pft, fft or auto)")


def _pft_spectrum(fb: FieldBundle, f: GridlessField, notes: list[str] | None) -> SpectrumBundle:
    """H → S by Fieldtracing's PFT (stationary phase), weights Δ²k from the wave vectors."""
    det = np.linalg.det(np.asarray(f.hessian, dtype=float))
    if not np.all(np.isfinite(det)) or np.any(det == 0.0):
        raise FieldError("the PFT needs a regular phase Hessian at every sample (use the FFT)")
    spec = f.pft()
    k = float(f.k)
    kvec = np.asarray(spec.points, dtype=float)
    q = k * k - np.sum(kvec * kvec, axis=1)
    keep = q > 2.0 * GRAZING * k * k
    # share of the power in evanescent waves, by the field's own quadrature (area in the plane)
    density = np.sum(np.abs(np.asarray(f.amplitude)) ** 2, axis=1) * _vertex_areas(
        np.asarray(f.points, dtype=float)
    )
    _report_cut(notes, float(density[~keep].sum() / (density.sum() or 1.0)), "evanescent")
    if keep.sum() < 3:
        raise FieldError("the spectrum is evanescent")
    kvec, q = kvec[keep], q[keep]
    area_k = _vertex_areas(kvec)
    field = (
        np.asarray(spec.amplitude)[keep]
        * np.exp(1j * np.asarray(spec.phase)[keep])[:, None]
        * area_k[:, None]
    )
    kz = np.sqrt(q)
    waves = ports.PlaneWaves(
        _port(f.wavelength_um, f.refractive_index),
        kvec[:, 0].copy(),
        kvec[:, 1].copy(),
        np.asarray(field[:, :2], dtype=np.complex128),
        solid_angle=area_k / (k * kz),
    )
    support = type(spec)(
        kvec,
        np.asarray(spec.phase)[keep],
        np.asarray(spec.gradient)[keep],
        np.asarray(spec.hessian)[keep],
        np.asarray(spec.amplitude)[keep],
        spec.wavelength_um,
        spec.refractive_index,
        "frequency",
    )
    if notes is not None:
        notes.append(f"{fb.key.source}: H -> S by the PFT, {len(kvec)} waves")
    return SpectrumBundle(fb.key, fb.plane, waves, fb.path, fb.branch, support)


def _grid_waves(
    spec: Any,
    grid: spectral.Grid,
    wavelength_um: float,
    n: float,
    notes: list[str] | None,
) -> tuple[ports.PlaneWaves, IntArray]:
    """Plane waves from a grid spectrum (ny, nx, 3) with the phase referred to the grid centre:
    k = k_c + κ, field = Ẽ(κ)·e^(−iκ·c) (phase 1 at the plane origin)."""
    kap_x, kap_y = grid.frequencies()
    kx = grid.carrier[0] + kap_x
    ky = grid.carrier[1] + kap_y
    k = n * 2.0 * np.pi / (wavelength_um * 1e-3)
    q = k * k - (kx * kx + ky * ky)
    magnitude = np.sum(np.abs(spec) ** 2, axis=-1)
    propagating = q > 2.0 * GRAZING * k * k
    total = float(magnitude.sum()) or 1.0
    _report_cut(notes, float(magnitude[~propagating].sum()) / total, "evanescent")
    amp = np.sqrt(magnitude)
    keep = propagating & (amp > SPECTRUM_THRESHOLD * float(amp.max() or 1.0))
    rows, cols = np.nonzero(keep)
    shift = np.exp(-1j * (kap_x * grid.centre[0] + kap_y * grid.centre[1]))
    field = spec[rows, cols] * shift[rows, cols][:, None]
    dk2 = (2.0 * np.pi) ** 2 / (grid.nx * grid.ny * grid.step * grid.step)
    kz = np.sqrt(q[rows, cols])
    waves = ports.PlaneWaves(
        _port(wavelength_um, n),
        np.asarray(kx[rows, cols], dtype=float),
        np.asarray(ky[rows, cols], dtype=float),
        np.asarray(field[:, :2], dtype=np.complex128),
        solid_angle=dk2 / (k * kz),
    )
    return waves, np.stack([rows, cols], axis=-1)


def _sample(f: GridlessField, max_points: int) -> tuple[spectral.Grid, ComplexArray]:
    """The field on the grid of `spectral.choose_grid`, as `spectral.sample` but with the
    residual interpolated cubically (Clough–Tocher) onto the grid nodes: the linear interpolation
    there loses about 2e−3 of a Gaussian beam's power at 10 samples per waist, the cubic one
    5e−5. Back on the samples bicubically with the field's phase ψ (`GridSupport.to_field`)."""
    grid = spectral.choose_grid(f, max_points)
    x, y = grid.axes()
    xx, yy = np.meshgrid(x, y)
    q = np.stack([xx.ravel(), yy.ravel()], axis=-1)
    e = f.residual_at(q, method="cubic") * np.exp(1j * (f.phase_at(q) - q @ grid.carrier))[:, None]
    return grid, np.asarray(e.reshape(grid.ny, grid.nx, 3), dtype=np.complex128)


def _fft_spectrum(
    fb: FieldBundle, f: GridlessField, notes: list[str] | None, max_points: int
) -> SpectrumBundle:
    """H → S by an FFT of the field's residual on a grid (`graph.spectral`)."""
    g, residual = _sample(f, max_points)
    spec = spectral.spectrum(residual)
    waves, bins = _grid_waves(spec, g, float(f.wavelength_um), float(f.refractive_index), notes)
    if notes is not None:
        notes.append(
            f"{fb.key.source}: H -> S by an FFT on {g.nx}×{g.ny} points, {len(waves.kx)} waves"
        )
    return SpectrumBundle(
        fb.key, fb.plane, waves, fb.path, fb.branch, GridSupport(g, np.asarray(bins), f)
    )


def _gridded_to_spectrum(
    fb: FieldBundle, notes: list[str] | None, max_points: int
) -> SpectrumBundle:
    """G → S: a gridded field of Fieldtracing (space, plain lattice) by `ports.to_plane_waves`."""
    f = fb.field
    try:
        sampled = spectral.sample(f, max_points)
    except spectral.SpectralError as e:
        raise FieldError(f"no spectrum for the gridded field: {e}") from e
    grid = sampled.grid
    n = float(f.refractive_index)
    port = _port(float(f.wavelength_um), n)
    values = np.asarray(sampled.residual)  # the physical field on the centred grid
    g = ports.GriddedField(port, grid.step, grid.step, values[..., :2])
    plain = ports.to_plane_waves(g)  # all bins, phase 1 at the centre node = plane origin
    spec = np.zeros((grid.ny, grid.nx, 3), dtype=np.complex128)
    rows = np.rint(plain.ky * grid.ny * grid.step / (2 * np.pi)).astype(int) % grid.ny
    cols = np.rint(plain.kx * grid.nx * grid.step / (2 * np.pi)).astype(int) % grid.nx
    spec[rows, cols, :2] = plain.field
    # E_n of the grid enters only the dropped share (the waves take it from transversality)
    spec[..., 2] = spectral.spectrum(values[..., 2:])[..., 0]
    waves, bins = _grid_waves(spec, grid, float(f.wavelength_um), n, notes)
    if notes is not None:
        notes.append(f"{fb.key.source}: G -> S on {grid.nx}×{grid.ny} points")
    return SpectrumBundle(
        fb.key, fb.plane, waves, fb.path, fb.branch, GridSupport(grid, np.asarray(bins), f)
    )


def _field3(w: ports.PlaneWaves) -> Any:
    """(N, 3) local field of the waves, E_n from transversality."""
    e = np.asarray(w.field)
    return np.column_stack([e[:, 0], e[:, 1], w.ez()])


def _grid_residual(sb: SpectrumBundle) -> ComplexArray:
    """The residual (ny, nx, 3) of a grid spectrum on its grid (the inverse FFT of its waves)."""
    support = sb.support
    assert isinstance(support, GridSupport)
    grid = support.grid
    rows, cols = support.bins[:, 0], support.bins[:, 1]
    kap_x, kap_y = grid.frequencies()
    shift = np.exp(1j * (kap_x * grid.centre[0] + kap_y * grid.centre[1]))
    spec = np.zeros((grid.ny, grid.nx, 3), dtype=np.complex128)
    spec[rows, cols] = _field3(sb.waves) * shift[rows, cols][:, None]
    return spectral.inverse(spec)


def spectrum_to_field(sb: SpectrumBundle, notes: list[str] | None = None) -> FieldBundle:
    """S → H onto the bundle's plane: the inverse PFT for a PFT spectrum, the inverse FFT back on
    the samples (or the grid) for a grid spectrum."""
    from fieldtracing.field import GridlessField

    w = sb.waves
    support = sb.support
    if w.solid_angle is None or support is None:
        raise FieldError(
            "discrete plane waves (no solid angles or no sampling) do not become a field on a "
            "plane; a detector takes them as they are (the sum of the waves)"
        )
    n = sb.refractive_index
    k = n * sb.k0
    if isinstance(support, GridlessField):
        kvec = np.column_stack([w.kx, w.ky])
        if kvec.shape != np.asarray(support.points).shape or not np.allclose(
            kvec, support.points, rtol=0.0, atol=1e-9 * k
        ):
            raise FieldError("the plane waves no longer match the spectrum's sampling")
        kz = np.real(w.port.kz(w.kx, w.ky))
        area_k = np.asarray(w.solid_angle) * k * kz
        residual = _field3(w) / area_k[:, None] * np.exp(-1j * np.asarray(support.phase))[:, None]
        freq = GridlessField(
            np.asarray(support.points),
            np.asarray(support.phase),
            np.asarray(support.gradient),
            np.asarray(support.hessian),
            np.asarray(residual, dtype=np.complex128),
            sb.wavelength_um,
            n,
            "frequency",
        )
        space = freq.inverse_pft()
        if notes is not None:
            notes.append(f"{sb.key.source}: S -> H by the inverse PFT, {len(w.kx)} samples")
        return FieldBundle(sb.key, sb.plane, space, sb.path, sb.branch)
    if isinstance(support, GridSupport):
        grid = support.grid
        if len(support.bins) != len(w.kx):
            raise FieldError("the plane waves no longer match the spectrum's grid")
        residual = _grid_residual(sb)
        try:
            out = support.to_field(residual, n)
        except spectral.SpectralError as e:
            raise FieldError(str(e)) from e
        if notes is not None:
            notes.append(f"{sb.key.source}: S -> H by an inverse FFT on {grid.nx}×{grid.ny}")
        return FieldBundle(sb.key, sb.plane, out, sb.path, sb.branch)
    raise FieldError(f"unknown sampling of the spectrum: {type(support).__name__}")


def rays_to_spectrum(
    bundle: RayBundle,
    plane: Plane | None = None,
    tolerance: float | None = None,
    notes: list[str] | None = None,
) -> SpectrumBundle:
    """R → S, always through H (decision 2): `field.rays_to_field` on `plane` (default: the
    port plane of the rays), then `field_to_spectrum`."""
    return field_to_spectrum(rays_to_field(bundle, plane), tolerance=tolerance, notes=notes)


def spectrum_to_rays(sb: SpectrumBundle, notes: list[str] | None = None) -> RayBundle:
    """S → R through H: the field on the port plane (`spectrum_to_field`), then rays along its
    local wave vectors (`field.field_to_rays`); `exec.spectrum_rays` if the bundle is not large
    against the wavelength in the medium."""
    from fieldtracing.field import GridlessField

    fb = spectrum_to_field(sb, notes)
    f = fb.field
    if isinstance(f, GridlessField) and notes is not None:
        weight = np.sum(np.abs(np.asarray(f.amplitude)) ** 2, axis=1)
        size = beam_size(np.asarray(f.points), weight)
        wavelength = sb.wavelength_um * 1e-3 / sb.refractive_index
        if size < BUNDLE_MIN_WAVELENGTHS * wavelength:
            notes.append(
                f"exec.spectrum_rays: the bundle ({size * 1e3:.3g} µm) is not large against the "
                f"wavelength ({wavelength * 1e3:.3g} µm): rays describe it poorly"
            )
    return field_to_rays(fb)


# ------------------------------------------------------------------ R30b: propagation and image


def _frame(plane: Plane) -> FloatArray:
    """Rows u, v, n: global -> local is x @ R.T, local -> global is x @ R."""
    return np.asarray(np.stack([plane.u, plane.v, plane.normal]))


def _local_waves(sb: SpectrumBundle) -> tuple[ComplexArray, ComplexArray]:
    """Wave vectors (k_u, k_v, k_n) (complex k_n for evanescent waves) and fields (E_u, E_v, E_n)
    in the frame of the bundle's plane."""
    w = sb.waves
    kz = np.asarray(w.port.kz(w.kx, w.ky), dtype=np.complex128)
    k = np.column_stack([np.asarray(w.kx, dtype=np.complex128), w.ky, kz])
    return np.asarray(k), np.asarray(_field3(w), dtype=np.complex128)


def _weights(sb: SpectrumBundle, e: ComplexArray) -> FloatArray:
    """The share of the power per wave: |E|²/dΩ (continuous) or |E|² (discrete waves)."""
    magnitude = np.sum(np.abs(e) ** 2, axis=1)
    if sb.waves.solid_angle is None:
        return np.asarray(magnitude)
    return np.asarray(magnitude / np.asarray(sb.waves.solid_angle))


def propagate(sb: SpectrumBundle, plane: Plane, notes: list[str] | None = None) -> SpectrumBundle:
    """S4: the spectrum on `plane` (module docstring). The new frame is the plane's (u, v, n),
    turned to (u, −v, −n) where its normal points against the light."""
    rot_s = _frame(sb.plane)
    cos = float(np.asarray(sb.plane.normal) @ np.asarray(plane.normal))
    if cos >= 0.0:
        target = Plane(np.asarray(plane.origin), plane.normal, plane.u, plane.v)
    else:
        target = Plane(np.asarray(plane.origin), -plane.normal, plane.u, -plane.v)
    rot_t = _frame(target)
    k_loc, e_loc = _local_waves(sb)
    w = sb.waves
    n = sb.refractive_index
    k = n * sb.k0
    delta = rot_s @ (np.asarray(target.origin) - np.asarray(sb.plane.origin))  # source frame
    evanescent = np.abs(np.imag(k_loc[:, 2])) > 0.0
    parallel = abs(cos) > 1.0 - 1e-12
    turn = rot_t @ rot_s.T  # source local -> target local (block diagonal for parallel planes)
    if parallel and np.allclose(turn, np.eye(3), rtol=0.0, atol=1e-12) and not delta.any():
        return SpectrumBundle(sb.key, target, w, sb.path, sb.branch, sb.support)  # same plane
    weight = _weights(sb, e_loc)
    total = float(weight.sum()) or 1.0
    if parallel:
        keep = ~evanescent | (delta[2] >= 0.0)
        _report_cut(notes, float(weight[~keep].sum()) / total, "evanescent waves behind the plane")
        k_new = k_loc[keep] @ turn.T
        k_new[:, 2] = k_loc[keep, 2]
    else:
        k_real = np.real(k_loc) @ rot_s @ rot_t.T  # global, then the target's frame
        forward = k_real[:, 2] > np.sqrt(2.0 * GRAZING) * k
        keep = ~evanescent & forward
        _report_cut(
            notes,
            float(weight[~keep].sum()) / total,
            "evanescent or not forward through the tilted plane",
        )
        k_new = k_real[keep].astype(np.complex128)
    phase = np.exp(1j * (k_loc[keep] @ delta))
    e_new = (e_loc[keep] * phase[:, None]) @ rot_s @ rot_t.T
    port = _port(sb.wavelength_um, n)
    waves = ports.PlaneWaves(
        port,
        np.ascontiguousarray(np.real(k_new[:, 0])),
        np.ascontiguousarray(np.real(k_new[:, 1])),
        np.ascontiguousarray(e_new[:, :2]),
        orders=None if w.orders is None else np.asarray(w.orders)[keep],
        solid_angle=None if w.solid_angle is None else np.asarray(w.solid_angle)[keep],
    )
    support: Any = None
    why = "tilted plane" if not parallel else "waves left out"
    if parallel and _is_pft(sb.support) and bool(keep.all()):
        support = _moved_support(sb.support, k_loc, delta, turn[:2, :2], k_new)
    elif parallel and isinstance(sb.support, GridSupport) and bool(keep.all()):
        support, why = _moved_grid(sb, target, turn, delta)
    if support is None and sb.support is not None and notes is not None:
        notes.append(
            f"{sb.key.source}: the spectrum moved in S keeps no sampling for S -> H ({why})"
        )
    return SpectrumBundle(sb.key, target, waves, sb.path, sb.branch, support)


def _moved_grid(
    sb: SpectrumBundle, target: Plane, turn: FloatArray, delta: FloatArray
) -> tuple[GridSupport | None, str]:
    """The sampling of a grid spectrum on a parallel plane shifted by `delta` (source frame), or
    None and why. The waves keep their bins (same carrier, step and size); the grid moves with
    the light and the samples come from the field on the old plane carried there geometrically
    (H -> R -> H, as `couplings.grating_field.onto_plane`): S -> H on the new plane is the inverse
    FFT of the moved waves – exact as long as the beam stays inside the grid's window, which is
    checked against the rule the grid was chosen by (samples' extent · `spectral.PAD`)."""
    from fieldtracing.field import GridlessField

    support = sb.support
    assert isinstance(support, GridSupport)
    if not np.allclose(turn, np.eye(3), rtol=0.0, atol=1e-12):
        return None, "frame turned in the plane"
    grid, f = support.grid, support.field
    shift = np.asarray(delta[:2], dtype=float)
    if abs(float(delta[2])) <= 1e-12:  # in the plane: new coordinates, the same field
        if not isinstance(f, GridlessField):
            if shift.any():
                return None, "a gridded field shifted in its plane"
            return support, ""
        moved = GridlessField(
            np.asarray(f.points) - shift,
            np.asarray(f.phase),
            np.asarray(f.gradient),
            np.asarray(f.hessian),
            np.asarray(f.amplitude),
            float(f.wavelength_um),
            float(f.refractive_index),
        )
        centre = np.asarray(grid.centre) - shift
        return GridSupport(_regrid(grid, centre), support.bins, moved), ""
    if not isinstance(f, GridlessField):
        return None, "a gridded field moved along the normal"
    try:
        fb = spectrum_to_field(sb)
        geometric: Any = rays_to_field(field_to_rays(fb), target).field
    except FieldError as e:
        return None, f"no geometric field on the new plane: {e}"
    pts = np.asarray(geometric.points, dtype=float)
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    window = np.array([grid.nx, grid.ny]) * grid.step
    if np.any((hi - lo) * spectral.PAD > window):
        return None, "the beam outgrows the grid of its spectrum"
    return GridSupport(_regrid(grid, 0.5 * (lo + hi)), support.bins, geometric), ""


def _regrid(grid: spectral.Grid, centre: Any, carrier: Any = None) -> spectral.Grid:
    return spectral.Grid(
        np.asarray(centre, dtype=float),
        grid.step,
        grid.nx,
        grid.ny,
        np.asarray(grid.carrier if carrier is None else carrier, dtype=float),
    )


def order_support(
    support: Any, take: IntArray, flip: bool, shift: FloatArray, kvec: FloatArray, n_out: float
) -> Any:
    """The sampling of a spectrum after a block that maps the wave k∥ to F·k∥ + `shift` on the
    same plane (R30c: an order of a grating, a stack with shift 0), F = diag(1, −1) for light
    that leaves into the half space it came from (`flip`, frame (u, −v, −n)), else the identity;
    `take` the incident waves that have an outgoing one, `kvec` (M, 2) their new k∥, `n_out` the
    medium the light leaves into. A PFT spectrum keeps ψ̃, its gradient (−ρ: the light leaves
    where it arrived) and Hessian on the shifted points; a grid spectrum keeps its bins on a grid
    with the carrier F·k_c + shift and samples with ψ + shift·x (the residual on the same bins,
    as `couplings.grating_field` does for fields: no window commensurate with the period is
    needed). None stays None."""
    if support is None:
        return None
    from fieldtracing.field import GridlessField

    sign = np.array([1.0, -1.0]) if flip else np.ones(2)
    sign3 = np.array([1.0, -1.0, -1.0]) if flip else np.ones(3)
    hess_sign = np.outer(sign, sign)
    shift = np.asarray(shift, dtype=float)
    if isinstance(support, GridlessField):  # PFT: frequency domain, one sample per wave
        return GridlessField(
            np.ascontiguousarray(kvec, dtype=float),
            np.asarray(support.phase)[take],
            np.asarray(support.gradient)[take] * sign,
            np.asarray(support.hessian)[take] * hess_sign,
            np.asarray(support.amplitude)[take] * sign3,
            float(support.wavelength_um),
            float(n_out),
            "frequency",
        )
    if not isinstance(support, GridSupport):
        return None
    grid = support.grid
    bins = np.asarray(support.bins)[take].copy()
    if flip:  # κ_v -> −κ_v
        bins[:, 0] = (-bins[:, 0]) % grid.ny
    moved = _regrid(grid, np.asarray(grid.centre) * sign, np.asarray(grid.carrier) * sign + shift)
    f = support.field
    if isinstance(f, GridlessField):
        pts = np.asarray(f.points) * sign
        sample: Any = GridlessField(
            pts,
            np.asarray(f.phase) + pts @ shift,
            np.asarray(f.gradient) * sign + shift,
            np.asarray(f.hessian) * hess_sign,
            np.asarray(f.amplitude) * sign3,
            float(f.wavelength_um),
            float(n_out),
        )
    else:
        from fieldtracing.field import GriddedField

        sample = GriddedField(
            float(f.step),
            np.asarray(f.values),
            float(f.wavelength_um),
            float(n_out),
            np.asarray(f.curvature, dtype=float) * hess_sign,
        )
    return GridSupport(moved, bins, sample)


def overlapping_orders(spectra: list[SpectrumBundle]) -> list[str]:
    """`exec.spectrum_rays` for orders that are not separated (plan S3): two spectra of one
    coherence group with different diffraction orders (`waves.orders`) on parallel planes whose
    significant waves overlap in k∥. Rays of either order then describe directions the other one
    has as well – the orders' interference within a direction is lost in R."""
    out: list[str] = []
    marked = [
        sb
        for sb in spectra
        if sb.continuous and sb.waves.orders is not None and len(sb.waves.kx) > 0
    ]
    boxes = []
    for sb in marked:
        w = sb.waves
        assert w.solid_angle is not None
        kz = np.real(w.port.kz(w.kx, w.ky))
        amp = np.linalg.norm(np.asarray(w.field), axis=1) / (np.asarray(w.solid_angle) * kz)
        big = amp >= SIGNIFICANT * float(amp.max() or 1.0)
        k = np.column_stack([w.kx, w.ky])[big]
        glob = k[:, 0:1] * sb.plane.u + k[:, 1:2] * sb.plane.v
        boxes.append(glob)
    for i in range(len(marked)):
        for j in range(i + 1, len(marked)):
            a, b = marked[i], marked[j]
            if a.key != b.key or float(np.asarray(a.plane.normal) @ b.plane.normal) < 1 - 1e-9:
                continue
            oa = {tuple(o) for o in np.asarray(a.waves.orders).tolist()}
            ob = {tuple(o) for o in np.asarray(b.waves.orders).tolist()}
            if oa == ob:
                continue
            frame = np.stack([a.plane.u, a.plane.v])
            pa, pb = boxes[i] @ frame.T, boxes[j] @ frame.T
            if np.all(pa.max(axis=0) >= pb.min(axis=0)) and np.all(
                pb.max(axis=0) >= pa.min(axis=0)
            ):
                out.append(
                    f"exec.spectrum_rays: the orders {a.path[-1] if a.path else '?'} and "
                    f"{b.path[-1] if b.path else '?'} overlap in their directions (the beam's "
                    "angular spread reaches the next order): as rays their interference within a "
                    "direction is lost"
                )
    return out


def _is_pft(support: Any) -> bool:
    try:
        from fieldtracing.field import GridlessField
    except ImportError:  # pragma: no cover - fieldtracing is needed for supports
        return False
    return isinstance(support, GridlessField)


def _moved_support(
    support: Any,
    k_loc: ComplexArray,
    delta: FloatArray,
    turn: FloatArray,
    k_new: ComplexArray,
) -> Any:
    """The frequency-domain field of a PFT spectrum on a parallel plane shifted by `delta` (source
    frame) and turned in the plane by `turn` (2 × 2)."""
    kpar = np.real(k_loc[:, :2])
    kn = np.real(k_loc[:, 2])
    dn = float(delta[2])
    phase = np.asarray(support.phase) + kpar @ delta[:2] + kn * dn
    gradient = np.asarray(support.gradient) + delta[:2] - kpar / kn[:, None] * dn
    eye = np.eye(2)
    d2kn = -(eye / kn[:, None, None] + np.einsum("ni,nj->nij", kpar, kpar) / kn[:, None, None] ** 3)
    hessian = np.asarray(support.hessian) + dn * d2kn
    turn3 = np.eye(3)
    turn3[:2, :2] = turn
    return type(support)(
        np.ascontiguousarray(np.real(k_new[:, :2])),
        phase,
        gradient @ turn.T,
        np.einsum("ij,njk,lk->nil", turn, hessian, turn),
        np.asarray(support.amplitude) @ turn3.T,
        support.wavelength_um,
        support.refractive_index,
        "frequency",
    )


def _landing(sb: SpectrumBundle) -> tuple[FloatArray, FloatArray] | None:
    """Where the significant waves start in the bundle's plane: per wave (N, 2) for a PFT
    spectrum (ρ = −∇ψ̃), or the box (2, 2) [lo, hi] of the field on the grid of an FFT
    spectrum; and the mask of the significant waves. None without a sampling."""
    w = sb.waves
    assert w.solid_angle is not None
    k = sb.refractive_index * sb.k0
    kz = np.real(w.port.kz(w.kx, w.ky))
    dk2 = np.asarray(w.solid_angle) * k * kz
    amp = np.linalg.norm(np.asarray(w.field), axis=1) / dk2
    significant = amp >= SIGNIFICANT * float(amp.max() or 1.0)
    support = sb.support
    if _is_pft(support):
        return -np.asarray(support.gradient, dtype=float), significant
    if isinstance(support, GridSupport):
        grid = support.grid
        magnitude = np.sum(np.abs(_grid_residual(sb)) ** 2, axis=-1)
        rows, cols = np.nonzero(magnitude >= SIGNIFICANT**2 * float(magnitude.max() or 1.0))
        x, y = grid.axes()
        box = np.array([[x[cols].min(), y[rows].min()], [x[cols].max(), y[rows].max()]])
        return box, significant
    return None


def sum_phase_step(sb: SpectrumBundle, corners: FloatArray) -> tuple[float, float]:
    """The largest phase step of the sum of the waves between neighbouring waves over a window
    with these (global) corners, and its limit (module docstring); (0, inf) for discrete waves,
    (nan, inf) for a continuous spectrum without sampling."""
    w = sb.waves
    if w.solid_angle is None:
        return 0.0, np.inf
    found = _landing(sb)
    if found is None:
        return float("nan"), np.inf
    start, significant = found
    q = (np.asarray(corners, dtype=float) - sb.plane.origin) @ _frame(sb.plane).T  # (C, 3)
    kz = np.real(w.port.kz(w.kx, w.ky))
    tan = np.column_stack([w.kx, w.ky]) / kz[:, None]
    tan = tan[significant]
    shift = q[:, None, :2] - tan[None, :, :] * q[:, None, 2:]  # (C, N, 2): r − Δn·tan θ
    if isinstance(sb.support, GridSupport):
        grid = sb.support.grid
        dk = 2.0 * np.pi / (np.array([grid.nx, grid.ny]) * grid.step)
        reach = np.maximum(np.abs(shift - start[0]), np.abs(shift - start[1]))
        return float(np.max(reach * dk)), GRID_STEP_LIMIT
    k = sb.refractive_index * sb.k0
    dk = np.sqrt(np.asarray(w.solid_angle)[significant] * k * kz[significant])
    reach = np.max(np.abs(shift - start[significant][None]), axis=2)  # ∞-norm
    return float(np.max(reach * dk[None, :])), SCATTERED_STEP_LIMIT


def wave_sum(sb: SpectrumBundle, points: FloatArray) -> ComplexArray:
    """The field Σ_j E_j e^(i k_j·(r − o)) of the waves at global points (M, 3), as global
    components (M, 3). Evanescent waves count only in front of the plane."""
    rot = _frame(sb.plane)
    k_loc, e_loc = _local_waves(sb)
    evanescent = np.abs(np.imag(k_loc[:, 2])) > 0.0
    q = (np.asarray(points, dtype=float) - sb.plane.origin) @ rot.T
    out = np.zeros((len(q), 3), dtype=np.complex128)
    block = max(1, SUM_BLOCK // max(len(k_loc), 1))
    for a in range(0, len(q), block):
        part = q[a : a + block]
        expo = np.exp(1j * (part @ k_loc.T))
        if evanescent.any():
            expo[:, evanescent] *= (part[:, 2:3] >= 0.0).astype(float)
        out[a : a + block] = expo @ e_loc
    return np.asarray(out @ rot)


def grid_sum(
    sb: SpectrumBundle, plane: Plane, half_u: float, half_v: float, nx: int, ny: int
) -> ComplexArray | None:
    """`wave_sum` on the pixel centres of a detector (`field.pixel_centres`) for a grid spectrum on
    a plane parallel to the detector, whose axes are those of the detector up to sign and order:
    the sum then factorises into two matrix products over the grid's columns and rows (a partial
    inverse DFT onto the pixels, exact as the sum). None where it does not apply."""
    support = sb.support
    if not isinstance(support, GridSupport) or len(sb.waves.kx) == 0:
        return None
    rot = _frame(sb.plane)
    if abs(float(np.asarray(plane.normal) @ rot[2])) < 1.0 - 1e-12:
        return None
    axes = np.array([[plane.u @ rot[0], plane.u @ rot[1]], [plane.v @ rot[0], plane.v @ rot[1]]])
    swapped = abs(axes[0, 1]) > 0.5
    signs = (axes[0, 1], axes[1, 0]) if swapped else (axes[0, 0], axes[1, 1])
    if not all(abs(abs(x) - 1.0) < 1e-12 for x in signs):
        return None
    grid = support.grid
    rows, cols = support.bins[:, 0], support.bins[:, 1]
    kap_x = 2.0 * np.pi * np.fft.fftfreq(grid.nx, d=grid.step)
    kap_y = 2.0 * np.pi * np.fft.fftfreq(grid.ny, d=grid.step)
    kx_col = grid.carrier[0] + kap_x  # k_u of every column of the grid
    ky_row = grid.carrier[1] + kap_y  # k_v of every row
    k_loc, e_loc = _local_waves(sb)
    q0 = rot @ (np.asarray(plane.origin) - np.asarray(sb.plane.origin))
    coeff = np.zeros((grid.ny, grid.nx, 3), dtype=np.complex128)
    coeff[rows, cols] = e_loc * np.exp(1j * k_loc[:, 2] * q0[2])[:, None]
    s = (np.arange(nx) + 0.5) / nx * 2 * half_u - half_u  # pixel coordinates along u, v
    t = (np.arange(ny) + 0.5) / ny * 2 * half_v - half_v
    if not swapped:  # q_u = q0_u + σ_u s, q_v = q0_v + σ_v t
        x = np.exp(1j * np.outer(q0[0] + signs[0] * s, kx_col))  # (nx, columns)
        y = np.exp(1j * np.outer(q0[1] + signs[1] * t, ky_row))  # (ny, rows)
        out = np.stack([y @ coeff[..., c] @ x.T for c in range(3)], axis=-1)
    else:  # q_u = q0_u + σ t (v of the detector along u), q_v = q0_v + σ' s
        x = np.exp(1j * np.outer(q0[1] + signs[0] * s, ky_row))  # (nx, rows)
        y = np.exp(1j * np.outer(q0[0] + signs[1] * t, kx_col))  # (ny, columns)
        out = np.stack([y @ coeff[..., c].T @ x.T for c in range(3)], axis=-1)
    return np.asarray(out.reshape(ny * nx, 3) @ rot)


def _on_pixels(
    sb: SpectrumBundle,
    plane: Plane,
    half_u: float,
    half_v: float,
    nx: int,
    ny: int,
    pixels: FloatArray,
    corners: FloatArray,
    notes: list[str],
) -> ComplexArray | FieldBundle:
    """The global field of one spectrum on the pixel centres, or – where the sum would alias or
    cost too much – the spectrum as a field on its plane (S -> H) for the field propagation."""
    count = len(sb.waves.kx)
    step, limit = sum_phase_step(sb, corners)
    if step > limit:
        notes.append(
            f"{sb.key.source}: the sum of {count} waves would alias on the window (phase step "
            f"{step:.3g} rad between neighbouring waves, at most {limit:.3g}) – S -> H on its "
            "plane and propagated as a field"
        )
        return spectrum_to_field(sb, notes)
    e = grid_sum(sb, plane, half_u, half_v, nx, ny)
    cost = count * len(pixels)
    if e is None and cost > SUM_BUDGET and sb.support is not None:
        notes.append(
            f"{sb.key.source}: the sum of {count} waves on {len(pixels)} pixels would cost "
            f"{cost:.2g} terms (at most {SUM_BUDGET:.2g}) – S -> H on its plane and propagated "
            "as a field"
        )
        return spectrum_to_field(sb, notes)
    if np.isnan(step):
        notes.append(
            f"{sb.key.source}: a continuous spectrum without its sampling – the sum of its "
            f"{count} waves is not checked for aliasing"
        )
    detail = "" if np.isnan(step) or step == 0.0 else f" (phase step {step:.2g} rad)"
    notes.append(f"{sb.key.source}: image from the sum of {count} plane waves{detail}")
    return e if e is not None else wave_sum(sb, pixels)


def spectrum_image(
    spectra: list[SpectrumBundle],
    fields: list[FieldBundle],
    plane: Plane,
    half_u: float,
    half_v: float,
    nx: int,
    ny: int,
    log: list[str] | None = None,
    progress: Callable[[float, str], None] | None = None,
    tolerance: float | None = None,
) -> Image:
    """Intensity on the pixel grid of `plane` from spectra and fields, amplitudes per coherence
    group (module docstring). `ray_power`: the power of the continuous spectra
    (`spectrum_power`) and of the fields (`field.field_power`); discrete waves carry none."""
    pixels = pixel_centres(plane, half_u, half_v, nx, ny)
    corners = np.asarray(
        [
            plane.origin + a * half_u * plane.u + b * half_v * plane.v
            for a in (-1, 1)
            for b in (-1, 1)
        ]
    )
    total = max(len(spectra) + len(fields), 1)
    sums: dict[CoherenceKey, ComplexArray] = {}
    power = sum(field_power(fb) for fb in fields)
    rest: list[FieldBundle] = []
    notes: list[str] = []
    try:
        for i, sb in enumerate(spectra):
            if progress is not None:
                progress(i / total, f"plane waves {i + 1}/{len(spectra)}")
            if sb.continuous:
                power += spectrum_power(sb)
            e = _on_pixels(sb, plane, half_u, half_v, nx, ny, pixels, corners, notes)
            if isinstance(e, FieldBundle):
                rest.append(e)
            else:
                sums[sb.key] = sums[sb.key] + e if sb.key in sums else e
    finally:
        if log is not None:
            log.extend(notes)
    rest += fields
    if rest:
        done = len(spectra) + len(fields) - len(rest)

        def part(p: float, text: str) -> None:
            if progress is not None:
                progress((done + p * len(rest)) / total, text)

        sums, _ = field_sums(rest, plane, half_u, half_v, nx, ny, log, part, tolerance, sums)
    return image_of(sums, half_u, half_v, nx, ny, power)


__all__ = [
    "BUNDLE_MIN_WAVELENGTHS",
    "CUT_WARNING",
    "SPECTRUM_THRESHOLD",
    "GridSupport",
    "SpectrumBundle",
    "beam_size",
    "field_to_spectrum",
    "fresnel_number",
    "grid_sum",
    "propagate",
    "rays_to_spectrum",
    "spectrum_image",
    "spectrum_power",
    "spectrum_to_field",
    "spectrum_to_rays",
    "sum_phase_step",
    "wave_sum",
]
