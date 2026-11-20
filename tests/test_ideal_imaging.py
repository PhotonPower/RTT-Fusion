"""Ideal aplanatic imaging of spectra (pilot K5, stage K5c): emitter → `optics.ideal_imaging` →
detector.

References (W96), independent of the block:
- (a) energy: the image spectrum carries the power of the object spectrum within the NA
  (`ports.spectrum_power` with a cone), for an emitter on glass with supercritical light;
- (b) low NA: the image of a dipole parallel to the surface in a homogeneous medium is the Airy
  pattern (2 J₁(v)/v)², v = k₀ (NA/M) ρ (Born & Wolf 8.5.2; vectorial corrections O(NA²));
- (c) magnification: an emitter shifted by Δ in the object plane gives the same image shifted by
  −M·Δ (inverted), whole pixels here;
- (d) symmetry: a dipole along the axis images to a ring, a minimum on the axis.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from scipy.special import j1

from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import Executor

pytest.importorskip("fieldtracing")
from rtt_fusion.graph.spectrum import spectrum_power  # noqa: E402
from rtt_fusion.ports import spectrum_power as cone_power  # noqa: E402

LAM, M, IMG = 0.6, 10.0, [0.0, 0.0, 50.0]


def _run(
    dipole: list[float],
    na: float,
    n_solution: float = 1.52,
    centre: tuple[float, float] = (0.0, 0.0),
    half: float = 0.06,
    pixels: int = 64,
    image: bool = True,
) -> tuple[Any, Any, Any]:
    reg, ev = workflow.registry()
    g = Graph()
    g.add(
        Node(
            "e",
            "emitter.dipole_layer",
            params={
                "dipole": dipole,
                "height_nm": 1.0,
                "wavelength_um": LAM,
                "n_substrate": 1.52,
                "n_solution": n_solution,
                "na": na,  # the pupil's edge on the emitter's θ quadrature
                "n_theta": 48,
                "distance": 0.005,
                "center": [centre[0], centre[1], 0.0],
            },
        )
    )
    g.add(
        Node(
            "m",
            "optics.ideal_imaging",
            params={"point": [0.0, 0.0, 1e-6], "image_point": IMG, "magnification": M, "na": na},
        )
    )
    g.add(
        Node(
            "d",
            "optics.detector",
            params={
                "point": IMG,
                "normal": [0.0, 0.0, -1.0],
                "half_size": [half, half],
                "pixels": [pixels, pixels],
            },
        )
    )
    g.connect(("e", "out"), ("m", "in"))
    g.connect(("m", "out"), ("d", "in"))
    for link in g.links:
        g.set_transport(link, "S")
    ex = Executor(g, reg, ev, available=workflow.available())
    img = ex.evaluate("d", "image") if image else None
    return ex.evaluate("e", "out"), ex.evaluate("m", "out"), img


def test_w96_power_within_the_na_is_kept() -> None:
    (s_in,), (s_out,), _ = _run([0.6, 0.0, 0.8], na=1.3, n_solution=1.333, image=False)
    want = cone_power(s_in.waves, cone_na=1.3) / s_in.k0
    assert want == pytest.approx(spectrum_power(s_in), rel=1e-12)  # the emitter stops at the NA
    assert spectrum_power(s_out) == pytest.approx(want, rel=1e-9)


def test_w96_low_na_image_is_airy() -> None:
    na = 0.1
    _, _, img = _run([1.0, 0.0, 0.0], na=na, half=0.09, pixels=72)
    v = np.asarray(img.values)
    u0, u1, _, _ = img.extent
    x = u0 + (np.arange(v.shape[1]) + 0.5) * (u1 - u0) / v.shape[1]
    xx, yy = np.meshgrid(x, x)
    arg = 2 * np.pi / (LAM * 1e-3) * (na / M) * np.hypot(xx, yy)
    airy = (2 * j1(arg) / arg) ** 2
    airy = airy / airy.max()  # both on the same pixels (none on the axis)
    assert np.abs(v / v.max() - airy).max() < 1e-3  # measured 1.8e-4 (vectorial: O(NA²))


def test_w96_shift_is_magnified_and_inverted() -> None:
    half, pixels = 0.06, 60
    step = 2 * half / pixels  # 2 µm per pixel
    delta = 3 * step / M  # three pixels in the image
    _, _, a = _run([1.0, 0.0, 0.0], na=0.3, half=half, pixels=pixels)
    _, _, b = _run([1.0, 0.0, 0.0], na=0.3, centre=(delta, 0.0), half=half, pixels=pixels)
    va, vb = np.asarray(a.values), np.asarray(b.values)
    # the image moved by −MΔ along world x; the detector faces −z, so its u axis is −x: +3 pixels
    assert np.abs(vb[:, 10:-10] - va[:, 7:-13]).max() < 1e-9 * va.max()


def test_w96_axial_dipole_is_a_ring() -> None:
    """Rotationally symmetric, a minimum on the axis (only the small longitudinal field E_z of the
    image side, sin θ' ≤ NA/(nM), is left there)."""
    _, _, img = _run([0.0, 0.0, 1.0], na=0.9, half=0.03, pixels=61)
    v = np.asarray(img.values)
    assert np.abs(v - v.T).max() < 1e-9 * v.max()
    assert np.abs(v - v[::-1, :]).max() < 1e-9 * v.max()
    assert v[30, 30] < 0.03 * v.max()
    assert v[30, 30] == v[28:33, 28:33].min()
