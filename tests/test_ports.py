"""P1 (R0): port representations G and S against analytic plane waves."""

from __future__ import annotations

import numpy as np
import pytest

from rtt_fusion.ports import GriddedField, PlaneWaves, Port, power, to_gridded, to_plane_waves

PORT = Port(z=0.0, direction=1, index=1.5, wavelength_um=0.5)
NX, NY, DX, DY = 24, 16, 0.5e-3 / 8, 0.5e-3 / 6


def _plane_wave(mx: int, my: int, e: tuple[complex, complex]) -> GriddedField:
    x = (np.arange(NX) - NX // 2) * DX
    y = (np.arange(NY) - NY // 2) * DY
    kx, ky = 2 * np.pi * mx / (NX * DX), 2 * np.pi * my / (NY * DY)
    phase = np.exp(1j * (kx * x[None, :] + ky * y[:, None]))
    return GriddedField(PORT, DX, DY, np.stack([e[0] * phase, e[1] * phase], axis=-1))


@pytest.mark.parametrize(("mx", "my"), [(0, 0), (3, -2), (-5, 1)])
def test_plane_wave_to_spectrum_and_back(mx: int, my: int) -> None:
    g = _plane_wave(mx, my, (0.3 - 0.2j, 0.7 + 0.1j))
    s = to_plane_waves(g, threshold=1e-9)
    assert len(s.kx) == 1
    assert abs(s.kx[0] - 2 * np.pi * mx / (NX * DX)) < 1e-9
    assert abs(s.ky[0] - 2 * np.pi * my / (NY * DY)) < 1e-9
    # phase referred to (0, 0): the amplitude is exactly the one put in
    assert np.max(np.abs(s.field[0] - np.array([0.3 - 0.2j, 0.7 + 0.1j]))) < 1e-13
    back = to_gridded(s, NX, NY, DX, DY)
    assert np.max(np.abs(back.field - g.field)) < 1e-13


def test_off_grid_wave_is_rejected() -> None:
    s = PlaneWaves(PORT, np.array([1.234]), np.array([0.0]), np.array([[1.0, 0.0]]))
    with pytest.raises(ValueError):
        to_gridded(s, NX, NY, DX, DY)


def test_longitudinal_component_and_power_of_tm_wave() -> None:
    # TM in the x-z plane: E ∥ (cos θ, 0, −sin θ) for k ∥ (sin θ, 0, cos θ), E_x = cos θ
    k = 1.5 * PORT.k0
    theta = 0.4
    s = PlaneWaves(
        PORT, np.array([k * np.sin(theta)]), np.array([0.0]), np.array([[np.cos(theta), 0.0]])
    )
    assert abs(s.ez()[0] - (-np.sin(theta))) < 1e-13
    assert abs(s.power_density()[0] - k * np.cos(theta) * 1.0) < 1e-9 * k
    # backward wave: E_z changes sign
    back = PlaneWaves(Port(0.0, -1, 1.5, 0.5), s.kx, s.ky, s.field)
    assert abs(back.ez()[0] - np.sin(theta)) < 1e-13


def test_power_of_gridded_plane_wave_and_evanescent() -> None:
    g = _plane_wave(3, -2, (0.3, 0.4j))
    s = to_plane_waves(g, threshold=1e-9)
    area = NX * DX * NY * DY
    assert abs(power(g) - float(s.power_density()[0]) * area) < 1e-12 * power(g)
    # an evanescent wave carries no power
    kx_ev = 2 * np.pi * 11 / (NX * DX)  # > 1.5 k0
    assert kx_ev > 1.5 * PORT.k0
    ev = PlaneWaves(PORT, np.array([kx_ev]), np.array([0.0]), np.array([[0.0, 1.0]]))
    assert ev.power_density()[0] == 0.0


def test_power_equals_bwpm_port_power() -> None:
    bwpm_adapter = pytest.importorskip("bwpm.adapter")
    rng = np.random.default_rng(1)
    field = rng.normal(size=(NY, NX, 2)) + 1j * rng.normal(size=(NY, NX, 2))
    g = GriddedField(PORT, DX, DY, field)
    p = bwpm_adapter.Port(0.0, 1, 1.5, 0.5, NX, NY, DX, DY)
    assert abs(
        power(g) - bwpm_adapter.port_power(bwpm_adapter.PortField(p, field))
    ) < 1e-12 * power(g)
