"""A dipole emitter in a planar multilayer (pilot K5, stage K5a): its angular spectrum in the two
half-spaces – supercritical light included – and its decay rates.

A fluorophore a few nm above a sensor surface couples its near field into the substrate: light
beyond the critical angle (supercritical angle fluorescence, SPCE through a gold film) and a
changed total decay rate (Chance–Prock–Silbey). A far field passed through a stack misses both,
because the evanescent part of the dipole's spectrum never leaves the solution as a plane wave.
Here the dipole is solved together with the stack, plane wave by plane wave (Weyl).

Geometry (along z, lengths in mm): a bottom half-space (index `n_bottom`, e.g. the glass
substrate), layers from the bottom up starting at `z_bottom` (index, thickness; e.g. a gold film,
the adlayer of BioOptoMD's n(z)), a top half-space (`n_top`, e.g. water). The dipole p sits at
(0, 0, z0) in any of them, in a lossless medium. Indices may be complex (Im > 0 absorbs, e^(−iωt)).

Method: for each transverse wave number κ the field is TE (amplitude ψ = E_y') or TM (ψ = Z₀H_y')
in the frame x' ∥ κ̂, y' = ẑ × x'. In layer l, ψ = u e^(ik_z(z − z_b)) + d e^(−ik_z(z − z_t))
(up-going referred to the layer's bottom, down-going to its top: only decaying exponentials), the
second tangential component is Y(u-part − d-part) with Y = k_z/k₀ (TE) resp. k_z/(k₀n²) (TM).
Tangential E and H are continuous at every interface; the half-spaces carry outgoing waves only.
The dipole enters as its primary field in its own layer, from the Weyl expansion

    E_prim(r) = ∫∫ i k₀²/(8π² k_z) (I − k k/k²) p e^(i(k∥·r∥ + k_z|z − z0|)) dk_x dk_y

(units ε₀ = μ₀ = c = 1, so ω = k₀ and Z₀ = 1; Novotny & Hecht, Principles of Nano-Optics, 2.12
and 10.6). The rates are normalised to the free dipole in the medium of the emitter,
P₀ = k₀⁴ n |p|²/(12π):

- `total`: P/P₀ = 1 + 6π/(n k₀³|p|²) Im(p*·E_s(r0)), E_s the field of the stack at the dipole
  (Chance–Prock–Silbey; Novotny & Hecht 8.140);
- `up`, `down`: the power leaving through the half-spaces (Poynting flux of the propagating
  waves, (2π)² ∫∫ n/2 |Ê|² cos θ dk∥), `down` also split by the numerical aperture of a cone in the
  bottom half-space (`collected`);
- `total − up − down` is what the layers absorb (0 for lossless stacks: a check of the method).

The angular spectrum in a half-space is given per solid angle, E(r) = ∫ A(k̂) e^(ik k̂·(r − r_ref)) dΩ
(the form of `ports.PlaneWaves` with `solid_angle`), A = k² |cos θ| Ê, with r_ref on the
interface below resp. above the dipole.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from typing import Any

import numpy as np

from rtt_fusion._types import ComplexArray, FloatArray

#: Gauss–Legendre points per κ interval (between branch points); 384 converges the total rate of
#: a dipole 10 nm above a gold film to ~1e-6 (96: 3e-3)
N_KAPPA = 384


@dataclass(frozen=True)
class Stack:
    """Bottom half-space, layers (index, thickness mm) from the bottom up starting at `z_bottom`,
    top half-space."""

    n_bottom: complex
    layers: tuple[tuple[complex, float], ...]
    n_top: complex
    z_bottom: float = 0.0

    @property
    def indices(self) -> ComplexArray:
        """Index per medium: bottom, layers, top."""
        return np.array(
            [self.n_bottom, *[n for n, _ in self.layers], self.n_top], dtype=np.complex128
        )

    @property
    def interfaces(self) -> FloatArray:
        """z of the L + 1 interfaces, bottom up."""
        z = [self.z_bottom]
        for _, t in self.layers:
            z.append(z[-1] + float(t))
        return np.asarray(z)

    def medium_at(self, z: float) -> int:
        """Index of the medium containing z (0: bottom half-space, L + 1: top)."""
        zi = self.interfaces
        return int(np.searchsorted(zi, z, side="right"))


def _kz(n: ComplexArray, k0: float, kappa: FloatArray) -> ComplexArray:
    """k_z with Im ≥ 0 (decaying or outgoing), shape (len(kappa), len(n))."""
    kz = np.sqrt(
        (np.asarray(n)[None, :] * k0) ** 2 - np.asarray(kappa, dtype=complex)[:, None] ** 2
    )
    return np.asarray(np.where(np.imag(kz) < 0, -kz, kz), dtype=np.complex128)


@dataclass(frozen=True)
class _Response:
    """Per κ and polarisation (TE, TM), the response to unit primary waves (up, down) at z0:
    amplitude in the bottom half-space (down-going, at its top), in the top half-space
    (up-going, at its bottom), and the stack's own field at z0 split into up- and down-going."""

    bottom: ComplexArray  # (K, 2 pol, 2 source)
    top: ComplexArray
    at_up: ComplexArray
    at_down: ComplexArray


def _solve(stack: Stack, k0: float, z0: float, kappa: FloatArray) -> _Response:
    n = stack.indices
    m = len(n)  # media: 0 bottom, 1..L layers, m-1 top
    zi = stack.interfaces
    s = stack.medium_at(z0)
    kz = _kz(n, k0, kappa)  # (K, m)
    big = len(kappa)
    thick = np.diff(zi)  # (L,)
    ph = np.ones((big, m), dtype=complex)  # e^(ik_z t) per finite layer
    if m > 2:
        ph[:, 1:-1] = np.exp(1j * kz[:, 1:-1] * thick[None, :])
    # unknowns: d_0, (u_l, d_l) for l = 1..L, u_top
    nu = 2 * (m - 1)

    def col_u(layer: int) -> int:
        return 2 * layer - 1 if layer < m - 1 else nu - 1

    def col_d(layer: int) -> int:
        return 0 if layer == 0 else 2 * layer

    out: dict[str, ComplexArray] = {k: np.zeros((big, 2, 2), complex) for k in ("b", "t", "u", "d")}
    for pol in (0, 1):  # TE, TM
        y = kz / k0 if pol == 0 else kz / (k0 * n[None, :] ** 2)
        a = np.zeros((big, nu, nu), dtype=complex)
        for i in range(m - 1):  # interface between medium i (below) and i + 1 (above)
            lo, hi = i, i + 1
            r0, r1 = 2 * i, 2 * i + 1  # rows: ψ, flux
            # below: its top; up-going e^(ik t) (half-space bottom: none), down-going 1
            if lo > 0:
                a[:, r0, col_u(lo)] += ph[:, lo]
                a[:, r1, col_u(lo)] += y[:, lo] * ph[:, lo]
            a[:, r0, col_d(lo)] += 1.0
            a[:, r1, col_d(lo)] += -y[:, lo]
            # above: its bottom; up-going 1, down-going e^(ik t) (half-space top: none)
            a[:, r0, col_u(hi)] -= 1.0
            a[:, r1, col_u(hi)] -= y[:, hi]
            if hi < m - 1:
                a[:, r0, col_d(hi)] -= ph[:, hi]
                a[:, r1, col_d(hi)] -= -y[:, hi] * ph[:, hi]
        for src in (0, 1):  # unit primary wave going up (0) or down (1) from z0
            rhs = np.zeros((big, nu), dtype=complex)
            if src == 0 and s < m - 1:  # reaches the top of the source medium (interface s)
                e = np.exp(1j * kz[:, s] * (zi[s] - z0))
                rhs[:, 2 * s] -= e  # the primary is on the "below" side of interface s
                rhs[:, 2 * s + 1] -= y[:, s] * e
            if src == 1 and s > 0:  # reaches the bottom of the source medium (interface s − 1)
                e = np.exp(1j * kz[:, s] * (z0 - zi[s - 1]))
                rhs[:, 2 * (s - 1)] += e  # on the "above" side of interface s − 1
                rhs[:, 2 * (s - 1) + 1] += -y[:, s] * e
            x = np.linalg.solve(a, rhs[..., None])[..., 0]
            out["b"][:, pol, src] = x[:, col_d(0)]
            out["t"][:, pol, src] = x[:, col_u(m - 1)]
            # the stack's field at z0 in the source medium
            if s == 0:
                up = np.zeros(big, complex)
                down = x[:, col_d(0)] * np.exp(-1j * kz[:, 0] * (z0 - zi[0]))
            elif s == m - 1:
                up = x[:, col_u(s)] * np.exp(1j * kz[:, s] * (z0 - zi[-1]))
                down = np.zeros(big, complex)
            else:
                up = x[:, col_u(s)] * np.exp(1j * kz[:, s] * (z0 - zi[s - 1]))
                down = x[:, col_d(s)] * np.exp(-1j * kz[:, s] * (z0 - zi[s]))
            out["u"][:, pol, src] = up
            out["d"][:, pol, src] = down
    return _Response(out["b"], out["t"], out["u"], out["d"])


def _frame(phi: FloatArray) -> tuple[FloatArray, FloatArray]:
    """x' = κ̂ and y' = ẑ × x' for the azimuths phi, (P, 3) each."""
    c, s = np.cos(phi), np.sin(phi)
    z = np.zeros_like(phi)
    return np.stack([c, s, z], -1), np.stack([-s, c, z], -1)


def _primary(
    p: ComplexArray, n: complex, k0: float, kappa: FloatArray, phi: FloatArray
) -> tuple[ComplexArray, ComplexArray]:
    """TE and TM amplitudes ψ of the primary waves going up and down, (K, P, 2 pol, 2 dir)."""
    kz = _kz(np.array([n]), k0, kappa)[:, 0]  # (K,)
    xp, yp = _frame(phi)
    out = np.zeros((len(kappa), len(phi), 2, 2), dtype=complex)
    k2 = (n * k0) ** 2
    for j, sign in enumerate((1.0, -1.0)):
        kv = (
            kappa[:, None, None] * xp[None, :, :]
            + sign * kz[:, None, None] * np.array([0.0, 0.0, 1.0])[None, None, :]
        )  # (K, P, 3)
        kp = np.einsum("kpi,i->kp", kv, p)
        e = (1j * k0**2 / (8 * np.pi**2 * kz))[:, None, None] * (
            p[None, None, :] - kv * (kp / k2)[..., None]
        )
        h = np.cross(kv, e) / k0  # Z₀H
        out[:, :, 0, j] = np.einsum("kpi,pi->kp", e, yp)
        out[:, :, 1, j] = np.einsum("kpi,pi->kp", h, yp)
    return out, kz


def _vector(
    psi: ComplexArray, n: complex, k0: float, kappa: FloatArray, phi: FloatArray, direction: float
) -> ComplexArray:
    """E of TE/TM waves psi (K, P, 2) going up (+1) or down (−1) in medium n, (K, P, 3)."""
    kz = _kz(np.array([n]), k0, kappa)[:, 0]
    xp, yp = _frame(phi)
    te, tm = psi[..., 0], psi[..., 1]
    ez = np.array([0.0, 0.0, 1.0])
    e_tm = (
        direction * kz[:, None, None] * xp[None, :, :] - kappa[:, None, None] * ez[None, None, :]
    ) / (k0 * n**2)
    return np.asarray(te[..., None] * yp[None, :, :] + tm[..., None] * e_tm)


def _kappa_nodes(
    stack: Stack, k0: float, z0: float, n_per: int, extra: tuple[float, ...] = ()
) -> tuple[FloatArray, FloatArray]:
    """κ nodes and weights for ∫ dκ: Gauss–Legendre between the branch points n k₀ of all media
    with the map κ = a + (b − a)(1 − cos t)/2 (dense at both ends, where k_z has square-root
    points; `extra`: further edges, e.g. the NA of a collection cone), then the evanescent tail
    κ = κ_max + v²/(2h), v ∈ [0, √46] (h: distance of the dipole to the nearest interface;
    e^(−2κh) has fallen by e^(−46); the square root removes the 1/k_z point at κ_max)."""
    n = np.real(stack.indices)
    zi = stack.interfaces
    h = max(float(np.min(np.abs(zi - z0))), 1e-9)
    kmax = float(n.max()) * k0
    edges = sorted({0.0, *[float(v) * k0 for v in n if v > 0], *[e for e in extra if 0 < e < kmax]})
    x, w = np.polynomial.legendre.leggauss(n_per)
    t = (x + 1) * np.pi / 2
    nodes, weights = [], []
    for a, b in pairwise(edges):
        nodes.append(a + (b - a) * (1 - np.cos(t)) / 2)
        weights.append(w * (np.pi / 2) * (b - a) * np.sin(t) / 2)
    vmax = np.sqrt(46.0)
    v = (x + 1) * vmax / 2
    nodes.append(kmax + v**2 / (2 * h))
    weights.append(w * (vmax / 2) * v / h)
    return np.concatenate(nodes), np.concatenate(weights)


@dataclass(frozen=True)
class Rates:
    """Powers normalised to the free dipole in the emitter's medium."""

    total: float
    up: float
    down: float
    collected: float  # down within the numerical aperture `na`
    na: float

    @property
    def absorbed(self) -> float:
        return self.total - self.up - self.down


def rates(
    stack: Stack,
    wavelength_um: float,
    z0: float,
    p: Any,
    na: float = 0.0,
    n_per: int = N_KAPPA,
    n_phi: int = 8,
) -> Rates:
    """Decay rates of the dipole p (3,) at height z0 (mm) in the stack."""
    k0 = 2 * np.pi / (wavelength_um * 1e-3)
    p = np.asarray(p, dtype=complex)
    s = stack.medium_at(z0)
    n_s = stack.indices[s]
    if abs(np.imag(n_s)) > 0:
        raise ValueError("the emitter must sit in a lossless medium")
    kappa, wk = _kappa_nodes(stack, k0, z0, n_per, (na * k0,) if na > 0 else ())
    phi = np.arange(n_phi) * 2 * np.pi / n_phi
    wphi = 2 * np.pi / n_phi
    resp = _solve(stack, k0, z0, kappa)
    prim, _ = _primary(p, n_s, k0, kappa, phi)  # (K, P, pol, dir)
    # stack field at z0: Σ over sources of the response × primary amplitude
    up = np.einsum("kad,kpad->kpa", resp.at_up, prim)
    down = np.einsum("kad,kpad->kpa", resp.at_down, prim)
    e_s = _vector(up, n_s, k0, kappa, phi, 1.0) + _vector(down, n_s, k0, kappa, phi, -1.0)
    integrand = np.einsum("kpi,i->kp", e_s, np.conj(p))
    e_dot = np.sum(integrand * (kappa * wk)[:, None]) * wphi
    p0_norm = 6 * np.pi / (float(np.real(n_s)) * k0**3 * float(np.vdot(p, p).real))
    total = 1.0 + p0_norm * float(np.imag(e_dot))
    # radiated powers through the half-spaces (propagating waves only)
    p0 = k0**4 * float(np.real(n_s)) * float(np.vdot(p, p).real) / (12 * np.pi)

    def flux(amp: ComplexArray, medium: int, direction: float, cone: float | None) -> float:
        n_m = stack.indices[medium]
        if abs(np.imag(n_m)) > 0:
            return 0.0
        nm = float(np.real(n_m))
        prop = kappa < nm * k0
        if cone is not None:
            prop &= kappa <= cone * k0
        psi = np.einsum("kad,kpad->kpa", amp, prim)
        e = _vector(psi, n_m, k0, kappa, phi, direction)
        kz = np.real(_kz(np.array([n_m]), k0, kappa)[:, 0])
        density = (2 * np.pi) ** 2 * nm / 2 * np.sum(np.abs(e) ** 2, -1) * (kz / (nm * k0))[:, None]
        return float(np.sum((density * (kappa * wk)[:, None])[prop]) * wphi) / p0

    m = len(stack.indices)
    zi = stack.interfaces
    kz_s = _kz(np.array([n_s]), k0, kappa)[:, 0]
    bottom, top = resp.bottom.copy(), resp.top.copy()
    if s == 0:  # the primary wave going down leaves directly, referred to the interface
        bottom[:, :, 1] += np.exp(1j * kz_s * (z0 - zi[0]))[:, None]
    if s == m - 1:
        top[:, :, 0] += np.exp(1j * kz_s * (zi[-1] - z0))[:, None]
    p_down = flux(bottom, 0, -1.0, None)
    p_up = flux(top, m - 1, 1.0, None)
    collected = flux(bottom, 0, -1.0, na) if na > 0 else 0.0
    return Rates(total, p_up, p_down, collected, na)


def angular_spectrum(
    stack: Stack,
    wavelength_um: float,
    z0: float,
    p: Any,
    theta: FloatArray,
    phi: FloatArray,
    side: str = "bottom",
) -> ComplexArray:
    """A(k̂) per solid angle (E = ∫ A e^(ik k̂·(r − r_ref)) dΩ) in the bottom (side "bottom",
    k̂_z < 0, r_ref = (0, 0, z_bottom)) or top half-space (k̂_z > 0, r_ref on the top interface),
    for polar angles theta (from −ẑ resp. +ẑ) and azimuths phi of the same shape; (N, 3)."""
    k0 = 2 * np.pi / (wavelength_um * 1e-3)
    p = np.asarray(p, dtype=complex)
    s = stack.medium_at(z0)
    m = len(stack.indices)
    medium = 0 if side == "bottom" else m - 1
    n_m = complex(stack.indices[medium])
    if abs(np.imag(n_m)) > 0:
        raise ValueError("the half-space must be lossless")
    theta = np.asarray(theta, dtype=float).ravel()
    phi = np.asarray(phi, dtype=float).ravel()
    kappa = np.real(n_m) * k0 * np.sin(theta)
    resp = _solve(stack, k0, z0, kappa)
    out = np.zeros((len(theta), 3), dtype=complex)
    for i in range(len(theta)):  # one azimuth per direction
        prim, _ = _primary(p, stack.indices[s], k0, kappa[i : i + 1], phi[i : i + 1])
        amp = (resp.bottom if side == "bottom" else resp.top)[i : i + 1].copy()
        if s == medium:  # the primary wave leaves directly (referred to the interface)
            zi = stack.interfaces
            kz_s = _kz(np.array([stack.indices[s]]), k0, kappa[i : i + 1])[0, 0]
            if side == "bottom":
                amp[:, :, 1] += np.exp(1j * kz_s * (z0 - zi[0]))
            else:
                amp[:, :, 0] += np.exp(1j * kz_s * (zi[-1] - z0))
        psi = np.einsum("kad,kpad->kpa", amp, prim)
        e = _vector(
            psi, n_m, k0, kappa[i : i + 1], phi[i : i + 1], -1.0 if side == "bottom" else 1.0
        )
        k = np.real(n_m) * k0
        out[i] = k**2 * abs(np.cos(theta[i])) * e[0, 0]
    return out


def angular_spectrum_grid(
    stack: Stack,
    wavelength_um: float,
    z0: float,
    theta: FloatArray,
    phi: FloatArray,
    side: str = "bottom",
) -> ComplexArray:
    """A(k̂) as `angular_spectrum` on the grid theta × phi for the three unit dipoles x̂, ŷ, ẑ:
    shape (3 dipoles, Nθ, Nφ, 3); the stack is solved once per θ. A(p) = Σ_j p_j A[j]."""
    k0 = 2 * np.pi / (wavelength_um * 1e-3)
    s = stack.medium_at(z0)
    m = len(stack.indices)
    medium = 0 if side == "bottom" else m - 1
    n_m = complex(stack.indices[medium])
    if abs(np.imag(n_m)) > 0:
        raise ValueError("the half-space must be lossless")
    theta = np.asarray(theta, dtype=float).ravel()
    phi = np.asarray(phi, dtype=float).ravel()
    kappa = np.real(n_m) * k0 * np.sin(theta)
    resp = _solve(stack, k0, z0, kappa)
    amp = (resp.bottom if side == "bottom" else resp.top).copy()
    if s == medium:
        zi = stack.interfaces
        kz_s = _kz(np.array([stack.indices[s]]), k0, kappa)[:, 0]
        if side == "bottom":
            amp[:, :, 1] += np.exp(1j * kz_s * (z0 - zi[0]))[:, None]
        else:
            amp[:, :, 0] += np.exp(1j * kz_s * (zi[-1] - z0))[:, None]
    direction = -1.0 if side == "bottom" else 1.0
    k = np.real(n_m) * k0
    out = np.zeros((3, len(theta), len(phi), 3), dtype=complex)
    for j in range(3):
        p = np.zeros(3, dtype=complex)
        p[j] = 1.0
        prim, _ = _primary(p, stack.indices[s], k0, kappa, phi)
        psi = np.einsum("kad,kpad->kpa", amp, prim)
        e = _vector(psi, n_m, k0, kappa, phi, direction)
        out[j] = k**2 * np.abs(np.cos(theta))[:, None, None] * e
    return out


def free_power(stack: Stack, wavelength_um: float, z0: float) -> float:
    """P₀ of a unit dipole in the emitter's medium in the units of this module (k₀⁴ n/(12π))."""
    k0 = 2 * np.pi / (wavelength_um * 1e-3)
    n = float(np.real(stack.indices[stack.medium_at(z0)]))
    return float(k0**4 * n / (12 * np.pi))


def band_power(
    stack: Stack,
    wavelength_um: float,
    z0: float,
    p: Any,
    side: str,
    na_lo: float,
    na_hi: float,
    n_per: int = N_KAPPA,
    n_phi: int = 8,
) -> float:
    """Power leaving through one half-space ("bottom" or "top") in the band of directions
    na_lo ≤ n sin θ < na_hi of that half-space, normalised to the free dipole in the emitter's
    medium (as `rates`); the primary wave is included where the dipole sits in that half-space."""
    k0 = 2 * np.pi / (wavelength_um * 1e-3)
    p = np.asarray(p, dtype=complex)
    s = stack.medium_at(z0)
    m = len(stack.indices)
    medium = 0 if side == "bottom" else m - 1
    n_m = complex(stack.indices[medium])
    if abs(np.imag(n_m)) > 0:
        raise ValueError("the half-space must be lossless")
    nm = float(np.real(n_m))
    lo, hi = max(na_lo, 0.0), min(na_hi, nm)
    if hi <= lo:
        return 0.0
    kappa, wk = _kappa_nodes(stack, k0, z0, n_per, (lo * k0, hi * k0))
    phi = np.arange(n_phi) * 2 * np.pi / n_phi
    resp = _solve(stack, k0, z0, kappa)
    prim, _ = _primary(p, stack.indices[s], k0, kappa, phi)
    amp = (resp.bottom if side == "bottom" else resp.top).copy()
    if s == medium:
        zi = stack.interfaces
        kz_s = _kz(np.array([stack.indices[s]]), k0, kappa)[:, 0]
        if side == "bottom":
            amp[:, :, 1] += np.exp(1j * kz_s * (z0 - zi[0]))[:, None]
        else:
            amp[:, :, 0] += np.exp(1j * kz_s * (zi[-1] - z0))[:, None]
    psi = np.einsum("kad,kpad->kpa", amp, prim)
    e = _vector(psi, n_m, k0, kappa, phi, -1.0 if side == "bottom" else 1.0)
    kz = np.real(_kz(np.array([n_m]), k0, kappa)[:, 0])
    density = (2 * np.pi) ** 2 * nm / 2 * np.sum(np.abs(e) ** 2, -1) * (kz / (nm * k0))[:, None]
    inside = (kappa >= lo * k0) & (kappa < hi * k0)
    p0 = free_power(stack, wavelength_um, z0) * float(np.vdot(p, p).real)
    return float(np.sum((density * (kappa * wk)[:, None])[inside]) * 2 * np.pi / n_phi) / p0
