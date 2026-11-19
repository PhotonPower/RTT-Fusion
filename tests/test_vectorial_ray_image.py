"""The vectorial ray image (stage R32): fields of one coherence key add as vectors.

Reference, independent of the code under test (W91): two plane waves of one coherence group,
E_m = e_m exp(i k (±x sin θ + z cos θ)) with transverse unit polarisations e_m, on the plane z = 0
give I(x) = cos θ · (|e₁|² + |e₂|² + 2 Re(e₁·ē₂ exp(2ik x sin θ))) (two-beam interference, Born/Wolf
7.2; cos θ the projection of the irradiance onto the plane). Parallel polarisations make fringes of
visibility 1, crossed ones none, opposite circular ones the visibility (1 − cos 2θ)/2. The scalar
image before R32 gave visibility 1 in every case.
"""

from __future__ import annotations

import numpy as np
import pytest

from rtt_fusion.graph.imaging import coherent_image
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle
from rtt_fusion.graph.optics import hexapolar

LAM = 0.6  # µm
K = 2 * np.pi / (LAM * 1e-3)  # 1/mm
THETA = 0.01  # rad: fringe period λ/(2 sin θ) = 0.03 mm
HALF, PIXELS = 0.15, 120
KEY = CoherenceKey("src", 0, 0, 0)


def plane_wave(sign: float, jones: list[complex], key: CoherenceKey = KEY) -> RayBundle:
    """A plane wave along (±sin θ, 0, cos θ), started on a disc 1 mm before z = 0, phase 0 at the
    origin of z = 0."""
    d = np.array([sign * np.sin(THETA), 0.0, np.cos(THETA)])
    frame = Plane.make(-1.0 * d, d, np.array([0.0, 1.0, 0.0]))
    px, py = hexapolar(40)
    xy = np.stack([px, py], axis=-1) * 0.5  # radius 0.5 mm, spacing 12.5 µm
    pos = -1.0 * d + xy[:, :1] * frame.u + xy[:, 1:] * frame.v
    return RayBundle(
        key=key,
        wavelength_um=LAM,
        position=pos,
        direction=np.tile(d, (len(pos), 1)),
        opl=pos @ d,  # phase k·r of the plane wave
        amplitude=np.ones(len(pos), dtype=complex),
        source_xy=xy,
        alive=np.ones(len(pos), dtype=bool),
        irradiance=1.0,
        path=("src",),
        jones=np.asarray(jones, dtype=complex),
        source_cos=1.0,
        source_frame=np.stack([frame.u, frame.v]),
    )


def detector() -> Plane:
    return Plane.make(np.zeros(3), np.array([0.0, 0.0, -1.0]), np.array([0.0, 1.0, 0.0]))


def reference(e1: np.ndarray, e2: np.ndarray, plane: Plane) -> np.ndarray:
    u = (np.arange(PIXELS) + 0.5) / PIXELS * 2 * HALF - HALF
    uu, vv = np.meshgrid(u, u)
    x = uu * plane.u[0] + vv * plane.v[0]
    cross = np.sum(e1 * np.conj(e2))
    return np.asarray(
        np.cos(THETA)
        * (
            np.sum(np.abs(e1) ** 2)
            + np.sum(np.abs(e2) ** 2)
            + 2 * np.real(cross * np.exp(2j * K * x * np.sin(THETA)))
        )
    )


Y = np.array([0.0, 1.0, 0.0], dtype=complex)


def p_axis(sign: float) -> np.ndarray:
    """The polarisation in the plane of incidence, transverse to (±sin θ, 0, cos θ)."""
    return np.array([np.cos(THETA), 0.0, -sign * np.sin(THETA)], dtype=complex)


CASES = {
    "parallel": (Y, Y),
    "crossed": (Y, p_axis(-1.0)),
    "circular, opposite": (
        (Y + 1j * p_axis(1.0)) / np.sqrt(2),
        (Y - 1j * p_axis(-1.0)) / np.sqrt(2),
    ),
    "circular, same": ((Y + 1j * p_axis(1.0)) / np.sqrt(2), (Y + 1j * p_axis(-1.0)) / np.sqrt(2)),
    "elliptical against linear": ((0.8 * Y + 0.6j * p_axis(1.0)), p_axis(-1.0)),
}


@pytest.mark.parametrize("case", list(CASES))
def test_w91_two_beam_interference_follows_the_polarisations(case: str) -> None:
    e1, e2 = CASES[case]
    plane = detector()
    image = coherent_image(
        [plane_wave(1.0, list(e1)), plane_wave(-1.0, list(e2))], plane, HALF, HALF, PIXELS, PIXELS
    )
    want = reference(e1, e2, plane)
    assert np.abs(image.values - want).max() < 1e-9 * want.max()
    # the visibility of the reference (the pixel centres miss the exact extrema)
    visibility = 2 * abs(np.sum(e1 * np.conj(e2))) / np.sum(np.abs(e1) ** 2 + np.abs(e2) ** 2)
    # fringes exactly where the reference has them; with crossed light none at all
    assert np.ptp(image.values) == pytest.approx(np.ptp(want), abs=1e-9)
    expected = {
        "parallel": 1.0,
        "crossed": 0.0,
        "circular, opposite": (1 - np.cos(2 * THETA)) / 2,
        "circular, same": (1 + np.cos(2 * THETA)) / 2,
    }
    if case in expected:
        assert visibility == pytest.approx(expected[case], abs=1e-12)


def test_w91_different_coherence_keys_add_intensities() -> None:
    """Two keys never interfere, whatever their polarisation."""
    plane = detector()
    other = CoherenceKey("src", 0, 0, 1)
    image = coherent_image(
        [plane_wave(1.0, list(Y)), plane_wave(-1.0, list(Y), other)],
        plane,
        HALF,
        HALF,
        PIXELS,
        PIXELS,
    )
    assert np.abs(image.values - 2 * np.cos(THETA)).max() < 1e-9


def test_w91_the_field_of_the_rays_counts_not_only_the_amplitude() -> None:
    """prt·jones enters the image: a ray field rotated by prt from y to x interferes with an x wave
    but not with a y wave, and |prt·jones| scales the power (ray_power likewise)."""
    plane = detector()
    rot = np.array([[0.0, 0.5, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.5]], dtype=complex)
    a = plane_wave(1.0, list(Y))
    a = a.evolve(prt=np.tile(rot, (a.size, 1, 1)))
    x_wave = plane_wave(-1.0, list(p_axis(-1.0)))
    y_wave = plane_wave(-1.0, list(Y))
    with_y = coherent_image([a, y_wave], plane, HALF, HALF, PIXELS, PIXELS)
    assert np.abs(with_y.values - 1.25 * np.cos(THETA)).max() < 1e-9
    with_x = coherent_image([a, x_wave], plane, HALF, HALF, PIXELS, PIXELS)
    e1 = rot @ Y
    assert np.abs(with_x.values - reference(e1, p_axis(-1.0), plane)).max() < 1e-9
    alone = coherent_image([a], plane, HALF, HALF, PIXELS, PIXELS)
    assert alone.ray_power == pytest.approx(0.25 * np.pi * 0.25, rel=0.02)  # |f|² · disc area
