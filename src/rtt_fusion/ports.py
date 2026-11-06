"""Port representations G (gridded field) and S (plane waves), conversion and power.

Conventions: docs/konventionen.md. A port is a plane z = const with a direction (+1: wave runs in
+z), the medium at the port and the vacuum wavelength. Fields are the transverse components
(E_x, E_y); E_z follows from transversality, E_z = −direction·(k_x E_x + k_y E_y)/k_z.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rtt_fusion._types import ComplexArray, FloatArray

GRAZING_TOLERANCE = 1e-12  # relative to |ñ²k₀²|, same rule as BWPM


@dataclass(frozen=True)
class Port:
    """Port plane: z (mm), direction ±1, medium index ñ, vacuum wavelength (µm)."""

    z: float
    direction: int
    index: complex
    wavelength_um: float

    def __post_init__(self) -> None:
        if self.direction not in (1, -1):
            raise ValueError("direction must be +1 or -1")
        if self.wavelength_um <= 0:
            raise ValueError("wavelength must be positive")

    @property
    def k0(self) -> float:
        """Vacuum wavenumber in 1/mm."""
        return float(2.0 * np.pi / (self.wavelength_um * 1e-3))

    def kz(self, kx: FloatArray, ky: FloatArray) -> ComplexArray:
        """k_z = sqrt(ñ²k₀² − k⊥²), Im ≥ 0 (Re ≥ 0 if real), grazing rule (1/mm)."""
        n2k2 = complex(self.index) ** 2 * self.k0 * self.k0
        q = np.asarray(n2k2 - (np.asarray(kx) ** 2 + np.asarray(ky) ** 2), dtype=np.complex128)
        q = np.where(np.abs(q) <= GRAZING_TOLERANCE * abs(n2k2), 0.0, q)
        kz = np.sqrt(q)
        flip = (kz.imag < 0) | ((kz.imag == 0) & (kz.real < 0))
        return np.asarray(np.where(flip, -kz, kz))


@dataclass(frozen=True)
class GriddedField:
    """Representation G: transverse field (ny, nx, 2) on the centred grid of the port."""

    port: Port
    dx: float
    dy: float
    field: ComplexArray

    def __post_init__(self) -> None:
        f = np.asarray(self.field)
        if f.ndim != 3 or f.shape[2] != 2:
            raise ValueError("field must have shape (ny, nx, 2)")

    @property
    def shape(self) -> tuple[int, int]:
        f = np.asarray(self.field)
        return int(f.shape[0]), int(f.shape[1])

    def frequencies(self) -> tuple[FloatArray, FloatArray]:
        """(k_x, k_y) of the grid, shape (ny, nx) each, numpy FFT order (1/mm)."""
        ny, nx = self.shape
        kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=self.dx)
        ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=self.dy)
        kxx, kyy = np.meshgrid(kx, ky)
        return np.asarray(kxx), np.asarray(kyy)


@dataclass(frozen=True)
class PlaneWaves:
    """Representation S: plane waves (k_x, k_y, E_x, E_y), phase 1 at (0, 0) in the port plane.

    Attributes:
        port: the port.
        kx, ky: transverse wave vectors (1/mm), shape (N,).
        field: transverse amplitudes (N, 2).
        orders: optional diffraction orders (N, 2) as integers (n, m).
        solid_angle: optional quadrature weights dΩ (N,) when the waves discretise a continuous
            angular spectrum E = ∫ A(k̂) e^(ik k̂·r) dΩ, field_j = A(k̂_j)·dΩ_j (then `power` applies).
    """

    port: Port
    kx: FloatArray
    ky: FloatArray
    field: ComplexArray
    orders: np.ndarray | None = None  # type: ignore[type-arg]
    solid_angle: FloatArray | None = None

    def __post_init__(self) -> None:
        n = np.asarray(self.kx).shape
        if np.asarray(self.ky).shape != n or np.asarray(self.field).shape != (*n, 2):
            raise ValueError("kx, ky must have shape (N,), field (N, 2)")

    def ez(self) -> ComplexArray:
        """E_z of every wave from transversality."""
        kz = self.port.kz(self.kx, self.ky)
        e = np.asarray(self.field)
        with np.errstate(divide="ignore", invalid="ignore"):
            out = -self.port.direction * (self.kx * e[:, 0] + self.ky * e[:, 1]) / kz
        return np.asarray(np.where(kz != 0, out, 0.0))

    def power_density(self) -> FloatArray:
        """Poynting flux density in z per wave, Re(k_z)·|E|² (arbitrary common factor)."""
        kz = self.port.kz(self.kx, self.ky)
        e = np.asarray(self.field)
        magnitude = np.sum(np.abs(e) ** 2, axis=1) + np.abs(self.ez()) ** 2
        return np.asarray(np.real(kz) * magnitude)


def spectrum_power(s: PlaneWaves, cone_na: float | None = None) -> float:
    """Power of a continuous angular spectrum (waves with `solid_angle`), optionally only of the
    waves within a cone sin θ ≤ NA/n around the port normal.

    E = ∫ A e^(ik k̂·r) dΩ has the flux (2π)²/k ∫ |A|² dΩ in the units of `power_density`
    (Re k_z |E|² per unit area; |A|² with the longitudinal part), k = n k₀.
    """
    if s.solid_angle is None:
        raise ValueError("plane waves without solid-angle weights")
    k = float(np.real(s.port.index)) * s.port.k0
    w = np.asarray(s.solid_angle)
    e = np.asarray(s.field)
    magnitude = np.sum(np.abs(e) ** 2, axis=1) + np.abs(s.ez()) ** 2
    terms = magnitude / w
    if cone_na is not None:
        sin_theta = np.hypot(s.kx, s.ky) / k
        terms = np.where(sin_theta <= cone_na / float(np.real(s.port.index)), terms, 0.0)
    return float((2 * np.pi) ** 2 / k * np.sum(terms))


def _centre_phase(g: GriddedField) -> ComplexArray:
    """exp(i k·r₀) with r₀ the position of index (0, 0) relative to the grid centre."""
    ny, nx = g.shape
    kx, ky = g.frequencies()
    return np.asarray(np.exp(-1j * (kx * (nx // 2) * g.dx + ky * (ny // 2) * g.dy)))


def to_plane_waves(g: GriddedField, threshold: float = 0.0) -> PlaneWaves:
    """G → S: Ẽ(k) = (1/N) Σ E(r) e^(−i k·r) on the centred grid; waves with |Ẽ| ≤ threshold
    (relative to the largest) are left out."""
    ny, nx = g.shape
    spectrum: ComplexArray = np.fft.fft2(np.asarray(g.field), axes=(0, 1)) / (nx * ny)
    spectrum = np.asarray(spectrum / _centre_phase(g)[..., None])
    kx, ky = g.frequencies()
    magnitude = np.max(np.abs(spectrum), axis=-1)
    keep = magnitude > threshold * float(np.max(magnitude) or 1.0)
    return PlaneWaves(g.port, kx[keep], ky[keep], spectrum[keep])


def to_gridded(s: PlaneWaves, nx: int, ny: int, dx: float, dy: float) -> GriddedField:
    """S → G: E(r) = Σ Ẽ e^(i k·r) on the centred grid; every wave must lie on a grid frequency."""
    mx = np.asarray(s.kx) * nx * dx / (2 * np.pi)
    my = np.asarray(s.ky) * ny * dy / (2 * np.pi)
    if (
        np.max(np.abs(mx - np.rint(mx)), initial=0) > 1e-9
        or np.max(np.abs(my - np.rint(my)), initial=0) > 1e-9
    ):
        raise ValueError("plane waves must lie on the grid frequencies (periodic window)")
    x = (np.arange(nx) - nx // 2) * dx
    y = (np.arange(ny) - ny // 2) * dy
    field = np.zeros((ny, nx, 2), dtype=np.complex128)
    amplitudes = np.asarray(s.field)
    for j in range(amplitudes.shape[0]):
        phase = np.exp(1j * (s.kx[j] * x[None, :] + s.ky[j] * y[:, None]))
        field += phase[..., None] * amplitudes[j][None, None, :]
    return GriddedField(s.port, dx, dy, field)


def power(g: GriddedField) -> float:
    """Poynting flux through the window of a gridded port field, Σ Re(k_z)|E|²·dx·dy/N (Parseval)
    with Ẽ as in `to_plane_waves`; times the window area this is the power."""
    s = to_plane_waves(g)
    ny, nx = g.shape
    return float(np.sum(s.power_density()) * nx * ny * g.dx * g.dy)
