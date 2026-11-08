"""Phase-space hand-off (Scattuccino S9b-1, D8/D10): field → packets → mode coupling."""

from __future__ import annotations

import numpy as np

from rtt_fusion.phasespace import _window, coupling_into_gaussian_mode, field_to_packets
from rtt_fusion.ports import GriddedField, Port, power

LAM = 0.8
K0 = 2 * np.pi / (LAM * 1e-3)


def _beam(nx: int, dx: float, w0: float, f: float | None) -> GriddedField:
    x = (np.arange(nx) - nx // 2) * dx
    xx, yy = np.meshgrid(x, x)
    phase = 0.0 if f is None else -K0 * (xx**2 + yy**2) / (2 * f)
    e = np.zeros((nx, nx, 2), dtype=complex)
    e[:, :, 0] = np.exp(-(xx**2 + yy**2) / w0**2 + 1j * phase)
    return GriddedField(Port(0.0, 1, 1.0 + 0j, LAM), dx, dx, e)


def test_collimated_beam_moments_power_and_coupling() -> None:
    dx, w0, sig = 0.002, 0.03, 0.01
    g = _beam(160, dx, w0, None)
    pk = field_to_packets(g, 200_000, sig, np.random.default_rng(1), pad=2)
    assert abs(pk.weight.sum() - power(g)) < 1e-12 * power(g) and pk.evanescent == 0.0
    _, size, _, _ = _window(g, sig, 2)
    dk = 2 * np.pi / (size * dx)
    vx = w0**2 / 4 + sig**2 + dx**2 / 12  # Wigner + window + jitter
    vk = 1 / w0**2 + 1 / (4 * sig**2) + dk**2 / 12
    n = pk.x.size
    for v, ref in (
        (np.var(pk.x), vx),
        (np.var(pk.y), vx),
        (np.var(pk.ux * K0), vk),
        (np.var(pk.uy * K0), vk),
    ):
        assert abs(v / ref - 1) < 4 * np.sqrt(2 / n) + 0.005
    # coupling of the Husimi packets into the same mode: product of Gaussian overlaps
    eta = coupling_into_gaussian_mode(pk.x, pk.y, pk.ux * K0, pk.uy * K0, pk.weight, w0)
    ref_eta = 1 / ((w0**2 / 4 + vx) * (1 / w0**2 + vk))
    assert abs(eta / ref_eta - 1) < 0.02
    # linear x polarisation is carried as Stokes (1, 1, 0, 0)
    assert np.allclose(pk.stokes, [1.0, 1.0, 0.0, 0.0])


def test_focused_beam_packets_follow_the_wave() -> None:
    nx, dx, w0, f, sig = 256, 0.002, 0.06, 1.0, 0.01
    g = _beam(nx, dx, w0, f)
    pk = field_to_packets(g, 200_000, sig, np.random.default_rng(2), pad=2)
    e = np.asarray(g.field)[:, :, 0]
    spec = np.fft.fft2(e)
    k = 2 * np.pi * np.fft.fftfreq(nx, dx)
    kx, ky = np.meshgrid(k, k)
    kz = np.sqrt(K0**2 - kx**2 - ky**2 + 0j)
    x = (np.arange(nx) - nx // 2) * dx
    xx, _ = np.meshgrid(x, x)
    for z in (0.5, 1.0, 1.5):
        inten = np.abs(np.fft.ifft2(spec * np.exp(1j * kz * z))) ** 2
        wave = np.sum(inten * xx**2) / inten.sum()
        pred = wave + sig**2 + z**2 / (4 * sig**2) / K0**2 + dx**2 / 12
        packets = np.var(pk.x + z * pk.ux / pk.uz)
        assert abs(packets / pred - 1) < 0.02  # measured ≤ 0.8 %


def test_periodic_windows_keep_the_power_of_a_wide_field() -> None:
    # a beam filling much of a periodic grid: with wrapping windows no power is lost at the edges
    g = _beam(96, 0.002, 0.08, None)
    lim = field_to_packets(g, 50_000, 0.008, np.random.default_rng(3), pad=2)
    per = field_to_packets(g, 50_000, 0.008, np.random.default_rng(3), pad=2, periodic=True)
    assert abs(per.weight.sum() - power(g)) < 1e-12 * power(g)
    # the restricted centres misplace the edge power onto the inner centres (narrower spread)
    assert np.var(per.x) > np.var(lim.x)
