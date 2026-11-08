"""Phase-space packets: field (port G) → photon packets (port P) and coupling into a mode
(Scattuccino S9b, decisions D8 and D10).

**Field → packets (Husimi).** For a transverse field E(x, y) the Husimi distribution
Q(r, k) = Σ_{c=x,y} |∫ E_c(r′) φ*(r′ − r) e^(−ik·r′) d²r′|², with the Gaussian window
φ(r) = exp(−|r|²/(4σ²)) (|φ|² has variance σ² per axis, its spectrum 1/(4σ²)), is the Wigner
function of E smoothed with that of φ: always ≥ 0, ∫Q ∝ power. Packets are drawn from Q:
window centre on the grid (stride), then k from the window's spectrum, both jittered uniformly
within their cell; the direction is (k⊥, k_z)/(n k₀) with k_z = sqrt(n²k₀² − k⊥²) along the port
direction; evanescent parts are dropped and booked as `evanescent`. Second moments are those of
the Wigner function plus σ² (position) and 1/(4σ²) (wavevector) per axis – the smoothing.

**Coupling into a Gaussian mode (D10).** For a mode of waist w₀ (|E|² ∝ exp(−2r²/w₀²)) the
Wigner density per axis is Gaussian with variances w₀²/4 (position) and 1/w₀² (wavevector,
1/mm²). The power coupled from light with phase-space density p is
η = (2π)² ∫ p_m p dr dk; with packets: η ≈ (2π)² Σᵢ (wᵢ/P) p_m(rᵢ, kᵢ). Exact for light whose
phase-space density is known exactly (incoherent sources, multiply scattered light); for packets
drawn from the Husimi distribution of a coherent field the smoothing biases η: for the mode itself
η = Π_axes [ (v_x,m + v_x,Q)(v_k,m + v_k,Q) ]^(−1/2) instead of 1 (4/9 for a matched window
σ = w₀/2) – coherent (ballistic) light is coupled through the field overlap instead (S9b-3).

**Resolution.** A packet's position at a distance z behind the hand-off spreads by
σ² + z²/(4σ²n²k₀²) on top of the wave's own second moment; the smoothing is smallest for
σ_opt = sqrt(z/(2nk₀)). A focus far behind the hand-off plane is therefore blurred by the
Husimi window – the hand-off should lie at or behind the focus, or the focus inside the
wave-optical region (consequence for D9).

Units: mm, 1/mm.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rtt_fusion._types import FloatArray
from rtt_fusion.ports import GriddedField, power


@dataclass
class Packets:
    x: FloatArray
    y: FloatArray
    ux: FloatArray
    uy: FloatArray
    uz: FloatArray
    weight: FloatArray  # power per packet
    stokes: FloatArray  # (N, 4) normalised, frame (x̂, ŷ)
    evanescent: float  # power share dropped as evanescent

    def k_perp(self, n_k0: float) -> tuple[FloatArray, FloatArray]:
        return np.asarray(self.ux * n_k0), np.asarray(self.uy * n_k0)


def _window(g: GriddedField, sigma: float, pad: int) -> tuple[int, int, FloatArray, FloatArray]:
    half = int(np.ceil(5 * sigma / g.dx))
    m = 2 * half + 1
    size = int(2 ** np.ceil(np.log2(m * pad)))
    o = np.arange(-half, half + 1)
    win = np.exp(-((o[None, :] * g.dx) ** 2 + (o[:, None] * g.dy) ** 2) / (4 * sigma**2))
    k = 2 * np.pi * np.fft.fftfreq(size, d=g.dx)
    return half, size, np.asarray(win), np.asarray(k)


def field_to_packets(
    g: GriddedField,
    n: int,
    sigma: float,
    rng: np.random.Generator,
    stride: int = 1,
    pad: int = 4,
    periodic: bool = False,
) -> Packets:
    """n packets from the Husimi distribution of g (window σ in mm; window centres every
    `stride` grid points; spectrum zero-padded by `pad`). periodic: the field is periodic on its
    grid (e.g. BWPM without absorbing boundary) – windows wrap around and every grid point can be
    a centre; otherwise only centres whose window lies inside the grid."""
    e = np.asarray(g.field)
    ny, nx = g.shape
    half, size, win, k = _window(g, sigma, pad)
    if 2 * half + 1 > min(nx, ny):
        raise ValueError("the window (±5σ) must be smaller than the grid")
    # power per window centre: Σ |E|² |φ|² (Parseval), on the centres
    inten = np.sum(np.abs(e) ** 2, axis=2)
    pad_ker = np.zeros((ny, nx))
    pad_ker[: 2 * half + 1, : 2 * half + 1] = win**2
    ker = np.roll(pad_ker, (-half, -half), axis=(0, 1))
    conv = np.real(np.fft.ifft2(np.fft.fft2(inten) * np.conj(np.fft.fft2(ker))))
    if periodic:
        ic, jc = np.arange(0, nx, stride), np.arange(0, ny, stride)
    else:
        ic = np.arange(half, nx - half, stride)
        jc = np.arange(half, ny - half, stride)
    wsum = np.maximum(conv[np.ix_(jc, ic)], 0.0)
    counts = rng.multinomial(n, (wsum / wsum.sum()).ravel()).reshape(wsum.shape)
    xs = (np.arange(nx) - nx // 2) * g.dx
    ys = (np.arange(ny) - ny // 2) * g.dy
    dk = abs(float(k[1] - k[0]))
    n_k0 = float(np.real(g.port.index)) * g.port.k0
    cols: dict[str, list[FloatArray]] = {
        key: [] for key in ("x", "y", "kx", "ky", "s1", "s2", "s3")
    }
    evanescent = 0
    for a, b in zip(*np.nonzero(counts), strict=True):
        c = int(counts[a, b])
        j, i = jc[a], ic[b]
        if periodic:
            rows = np.arange(j - half, j + half + 1)
            cols_ = np.arange(i - half, i + half + 1)
            patch = (
                e.take(rows, axis=0, mode="wrap").take(cols_, axis=1, mode="wrap") * win[:, :, None]
            )
        else:
            patch = e[j - half : j + half + 1, i - half : i + half + 1, :] * win[:, :, None]
        fx = np.fft.fft2(patch[:, :, 0], s=(size, size))
        fy = np.fft.fft2(patch[:, :, 1], s=(size, size))
        q = np.abs(fx) ** 2 + np.abs(fy) ** 2
        pick = rng.choice(q.size, size=c, p=(q / q.sum()).ravel())
        kyi, kxi = np.unravel_index(pick, q.shape)
        kxs = k[kxi] + (rng.random(c) - 0.5) * dk
        kys = k[kyi] + (rng.random(c) - 0.5) * dk
        ok = kxs**2 + kys**2 < n_k0**2
        evanescent += int((~ok).sum())
        fxs, fys = fx[kyi, kxi][ok], fy[kyi, kxi][ok]
        it = np.abs(fxs) ** 2 + np.abs(fys) ** 2
        cc = fxs * np.conj(fys)
        cols["x"].append(xs[i] + (rng.random(int(ok.sum())) - 0.5) * stride * g.dx)
        cols["y"].append(ys[j] + (rng.random(int(ok.sum())) - 0.5) * stride * g.dy)
        cols["kx"].append(kxs[ok])
        cols["ky"].append(kys[ok])
        cols["s1"].append((np.abs(fxs) ** 2 - np.abs(fys) ** 2) / it)
        cols["s2"].append(2 * cc.real / it)
        cols["s3"].append(-2 * cc.imag / it)
    arr = {key: np.concatenate(v) if v else np.zeros(0) for key, v in cols.items()}
    m = arr["x"].size
    kz = np.sqrt(np.maximum(n_k0**2 - arr["kx"] ** 2 - arr["ky"] ** 2, 0.0)) * g.port.direction
    ev = evanescent / n
    p_total = power(g) * (1.0 - ev)
    stokes = np.column_stack([np.ones(m), arr["s1"], arr["s2"], arr["s3"]])
    return Packets(
        arr["x"],
        arr["y"],
        arr["kx"] / n_k0,
        arr["ky"] / n_k0,
        kz / n_k0,
        np.full(m, p_total / max(m, 1)),
        stokes,
        ev,
    )


def gaussian_mode_density(
    x: FloatArray, y: FloatArray, kx: FloatArray, ky: FloatArray, w0: float
) -> FloatArray:
    """Wigner density (normalised to 1) of a Gaussian mode of waist w0 at its waist, centred."""
    vx, vk = w0**2 / 4, 1 / w0**2
    return np.asarray(
        np.exp(-(x**2 + y**2) / (2 * vx) - (kx**2 + ky**2) / (2 * vk))
        / ((2 * np.pi) ** 2 * vx * vk)
    )


def coupling_into_gaussian_mode(
    x: FloatArray, y: FloatArray, kx: FloatArray, ky: FloatArray, weight: FloatArray, w0: float
) -> float:
    """η = (2π)² Σ (wᵢ/P) p_m(rᵢ, kᵢ) – fraction of the packets' power coupled into the mode."""
    p = weight / weight.sum()
    return float((2 * np.pi) ** 2 * np.sum(p * gaussian_mode_density(x, y, kx, ky, w0)))
