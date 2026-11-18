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

Power: `spectrum_power` in the units of `field.field_power` (W): the flux of `ports` divided by
k₀ (`ports.spectrum_power` counts Re k_z |E|², `field_power` n cos γ |U|²).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.spatial import Delaunay, QhullError

from rtt_fusion import ports
from rtt_fusion._types import FloatArray, IntArray
from rtt_fusion.graph import spectral
from rtt_fusion.graph.field import FieldBundle, FieldError, field_to_rays, rays_to_field
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


@dataclass(frozen=True)
class GridSupport:
    """How a grid spectrum was sampled: the field on its grid (`spectral.Sampled`) and the grid
    bin (row, column) of every wave, so that S → H puts the waves back on the same grid."""

    sampled: spectral.Sampled
    bins: IntArray


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


def _sample(f: GridlessField, max_points: int) -> spectral.Sampled:
    """The field on the grid of `spectral.choose_grid`, as `spectral.sample` but with the
    residual interpolated cubically (Clough–Tocher) onto the grid nodes: the linear interpolation
    there loses about 2e−3 of a Gaussian beam's power at 10 samples per waist, the cubic one
    5e−5. Back on the samples bicubically with the field's phase ψ (`spectral.outgoing`)."""
    grid = spectral.choose_grid(f, max_points)
    x, y = grid.axes()
    xx, yy = np.meshgrid(x, y)
    q = np.stack([xx.ravel(), yy.ravel()], axis=-1)
    e = f.residual_at(q, method="cubic") * np.exp(1j * (f.phase_at(q) - q @ grid.carrier))[:, None]
    residual = np.asarray(e.reshape(grid.ny, grid.nx, 3), dtype=np.complex128)

    def back(res: Any, n: float, reflected: bool) -> Any:
        return spectral.outgoing(f, spectral.at_samples(res, grid, f), n, reflected)

    return spectral.Sampled(grid, residual, back)


def _fft_spectrum(
    fb: FieldBundle, f: GridlessField, notes: list[str] | None, max_points: int
) -> SpectrumBundle:
    """H → S by an FFT of the field's residual on a grid (`graph.spectral`)."""
    sampled = _sample(f, max_points)
    spec = spectral.spectrum(sampled.residual)
    waves, bins = _grid_waves(
        spec, sampled.grid, float(f.wavelength_um), float(f.refractive_index), notes
    )
    if notes is not None:
        g = sampled.grid
        notes.append(
            f"{fb.key.source}: H -> S by an FFT on {g.nx}×{g.ny} points, {len(waves.kx)} waves"
        )
    return SpectrumBundle(
        fb.key, fb.plane, waves, fb.path, fb.branch, GridSupport(sampled, np.asarray(bins))
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
        fb.key, fb.plane, waves, fb.path, fb.branch, GridSupport(sampled, np.asarray(bins))
    )


def _field3(w: ports.PlaneWaves) -> Any:
    """(N, 3) local field of the waves, E_n from transversality."""
    e = np.asarray(w.field)
    return np.column_stack([e[:, 0], e[:, 1], w.ez()])


def spectrum_to_field(sb: SpectrumBundle, notes: list[str] | None = None) -> FieldBundle:
    """S → H onto the bundle's plane: the inverse PFT for a PFT spectrum, the inverse FFT back on
    the samples (or the grid) for a grid spectrum."""
    from fieldtracing.field import GridlessField

    w = sb.waves
    support = sb.support
    if w.solid_angle is None or support is None:
        raise FieldError(
            "discrete plane waves (no solid angles or no sampling) cannot become a field yet "
            "(detector and propagation in the spectrum come with R30b)"
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
        grid = support.sampled.grid
        rows, cols = support.bins[:, 0], support.bins[:, 1]
        if len(rows) != len(w.kx):
            raise FieldError("the plane waves no longer match the spectrum's grid")
        kap_x, kap_y = grid.frequencies()
        shift = np.exp(1j * (kap_x * grid.centre[0] + kap_y * grid.centre[1]))
        spec = np.zeros((grid.ny, grid.nx, 3), dtype=np.complex128)
        spec[rows, cols] = _field3(w) * shift[rows, cols][:, None]
        residual = spectral.inverse(spec)
        try:
            out = support.sampled.back(residual, n, False)
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


__all__ = [
    "BUNDLE_MIN_WAVELENGTHS",
    "CUT_WARNING",
    "SPECTRUM_THRESHOLD",
    "GridSupport",
    "SpectrumBundle",
    "beam_size",
    "field_to_spectrum",
    "fresnel_number",
    "rays_to_spectrum",
    "spectrum_power",
    "spectrum_to_field",
    "spectrum_to_rays",
]
