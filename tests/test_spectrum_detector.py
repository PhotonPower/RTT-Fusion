"""Representation S at the detector and in free propagation (R30b, docs/plan-s-spektren.md).

References, independent of the code under test:
- W79 two plane waves: I(x) = a₁² + a₂² + 2a₁a₂ cos(2k sin θ·x − φ) for waves at ±θ in the x–z
  plane, polarised along y (fringe period λ/(2 sin θ), contrast 2a₁a₂/(a₁² + a₂²), phase φ);
  waves of different coherence groups add as intensities. The focus of a Gaussian beam through
  an aplanatic lens (NA 0.5): the Richards–Wolf integrals (Bessel form, Gauss–Legendre in θ) with
  the Gaussian pupil apodisation written out here, and the image of the same link with H
  (Fieldtracing's Debye focus) – "focus via S = focus via H" to 1e−3. The top-hat focus of W13
  against `fieldtracing.focus.richards_wolf`.
- W80 free propagation in the spectrum: a Gaussian beam w₀ moved by z in the spectrum against the
  Fresnel integral of a Gaussian, the paraxial Gaussian beam
  E = (w₀/w)·exp(−ρ²/w²)·exp(i(kz + kρ²/(2R) − ζ)) (non-paraxial terms 1/(k w₀)² ≈ 2e−5 here);
  a PFT spectrum moved by Δz and taken back by the inverse PFT against Fieldtracing's geometric
  propagation of the field (`GridlessField.propagate`); plane waves are the same in any frame:
  a spectrum moved onto a tilted or turned plane gives the field of the original at the same
  points in space; an evanescent wave decays as exp(−κΔz), κ = sqrt(k∥² − k²).
- The Raytatouille detector (paraboloid of W17): S on the link into the focus equals H (whole
  paths and node-wise).
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.special import jv

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.graph import Category, Graph, Node, NodeType, builtin_registry, optics
from rtt_fusion.graph.execute import WITH_SPECTRA, Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.light import CoherenceKey, Plane
from rtt_fusion.graph.model import light_out
from rtt_fusion.ports import PlaneWaves, Port

pytest.importorskip("fieldtracing")
from fieldtracing.field import GriddedField, GridlessField  # noqa: E402

from rtt_fusion.graph import spectrum as sp  # noqa: E402
from rtt_fusion.graph.field import FieldBundle, pixel_centres  # noqa: E402

DATA = Path(__file__).parent / "data"
LAM = 0.6  # µm
K = 2 * np.pi / (LAM * 1e-3)  # 1/mm
KEY = CoherenceKey("waves")


def pixel_points(plane: Plane, half: float, n: int) -> FloatArray:
    """Global pixel centres of a square detector, row by row (rows along v)."""
    u = (np.arange(n) + 0.5) / n * 2 * half - half
    uu, vv = np.meshgrid(u, u)
    return np.asarray(plane.origin + uu.reshape(-1, 1) * plane.u + vv.reshape(-1, 1) * plane.v)


# ------------------------------------------------------------------------------ W79 detector


def wave_source(bundles: list[sp.SpectrumBundle]) -> tuple[Executor, Graph]:
    """A source block that emits the given spectra (S only) into an optics detector."""
    reg = builtin_registry()
    ev: dict[str, Any] = optics.register(reg)
    reg.add(NodeType("t.waves", Category.SOURCE, (light_out("out", "S"),)))

    def emit(ctx: Any, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
        return list(bundles)

    emit.representations = frozenset("S")  # type: ignore[attr-defined]
    ev["t.waves"] = emit
    g = Graph()
    g.add(Node("src", "t.waves"))
    return Executor(g, reg, ev, available=WITH_SPECTRA), g


def two_waves(
    theta: float, a1: float, a2: float, phi: float, keys: tuple[str, str] = ("w", "w")
) -> list[sp.SpectrumBundle]:
    """Waves at ±θ in the x–z plane, polarised along y, phase φ of the second at the origin."""
    plane = Plane.make([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])  # u = x, v = y
    port = Port(0.0, 1, 1.0, LAM)
    kx = np.array([K * np.sin(theta), -K * np.sin(theta)])
    field = np.array([[0.0, a1], [0.0, a2 * np.exp(1j * phi)]], dtype=complex)
    if keys[0] == keys[1]:
        waves = PlaneWaves(port, kx, np.zeros(2), field, orders=np.array([[1, 0], [-1, 0]]))
        return [sp.SpectrumBundle(CoherenceKey(keys[0]), plane, waves)]
    return [
        sp.SpectrumBundle(
            CoherenceKey(key), plane, PlaneWaves(port, kx[i : i + 1], np.zeros(1), field[i : i + 1])
        )
        for i, key in enumerate(keys)
    ]


@pytest.mark.parametrize("tilt", [0.0, 0.3])
def test_w79_two_plane_waves_give_analytic_fringes(tilt: float) -> None:
    theta, a1, a2, phi = 0.01, 1.0, 0.5, 0.7
    ex, g = wave_source(two_waves(theta, a1, a2, phi))
    half, n = 0.1, 64
    normal = [np.sin(tilt), 0.0, -np.cos(tilt)]  # tilted about y: the fringes stretch by 1/cos
    params = {"point": [0, 0, 5.0], "normal": normal, "half_size": [half, half], "pixels": [n, n]}
    g.add(Node("cam", "optics.detector", params=params))
    g.connect(("src", "out"), ("cam", "in"))  # auto: both blocks take S natively
    im: Image = ex.evaluate("cam", "image")
    assert ex.plan is not None and ex.plan.representations["src.out -> cam.in"] == "S"
    assert any("image from the sum of 2 plane waves" in line for line in ex.context.log)
    x = pixel_points(Plane.make([0, 0, 5.0], normal), half, n)[:, 0].reshape(n, n)
    period = LAM * 1e-3 / (2 * np.sin(theta))
    mean, contrast = a1**2 + a2**2, 2 * a1 * a2 / (a1**2 + a2**2)
    ref = mean * (1 + contrast * np.cos(2 * np.pi * x / period - phi))
    assert np.abs(im.values - ref).max() < 1e-10 * ref.max()  # phases k·r ≈ 5e4 rad
    assert im.ray_power == 0.0  # plane waves without extent carry no finite power
    # the fringes read off the image: the period from the strongest FFT bin of a row (along u,
    # λ/(2 sin θ cos(tilt)) on the tilted detector), contrast and phase by least squares in x
    row = im.values[n // 2]
    spectrum = np.abs(np.fft.rfft(row - row.mean()))
    measured = 1.0 / np.fft.rfftfreq(n, 2 * half / n)[np.argmax(spectrum)]
    assert measured == pytest.approx(period / np.cos(tilt), rel=0.1)  # one FFT bin
    q = 2 * np.pi / period * x.ravel()
    basis = np.column_stack([np.ones_like(q), np.cos(q), np.sin(q)])
    a, b, c = np.linalg.lstsq(basis, im.values.ravel(), rcond=None)[0]
    assert np.hypot(b, c) / a == pytest.approx(contrast, rel=1e-9)
    assert np.arctan2(c, b) == pytest.approx(phi, abs=1e-9)


def test_w79_waves_of_different_groups_add_as_intensities() -> None:
    ex, g = wave_source(two_waves(0.01, 1.0, 0.5, 0.7, keys=("a", "b")))
    params = {"point": [0, 0, 5.0], "half_size": [0.1, 0.1], "pixels": [16, 16]}
    g.add(Node("cam", "optics.detector", params=params))
    g.connect(("src", "out"), ("cam", "in"))
    im: Image = ex.evaluate("cam", "image")
    assert len(im.groups) == 2
    assert np.abs(im.values - 1.25).max() < 1e-12


def focus_graph(
    transport: str, waist: float = 0.0, rings: int = 40, n: int = 41, half: float = 2.5e-3
) -> Executor:
    """W13: collimated beam (radius 1 mm, Gaussian waist if given) through an aplanatic lens
    f = 2 mm (NA 0.5), detector in the focus."""
    f, a = 2.0, 1.0
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    src: dict[str, Any] = {"radius": a, "rings": rings}
    if waist:
        src["gaussian_waist"] = waist
    g.add(Node("src", "optics.collimated_source", params=src))
    g.add(
        Node(
            "lens",
            "optics.aplanatic_lens",
            params={"focal_length": f, "aperture_radius": a * (1 + 1e-6), "point": [0, 0, 1]},
        )
    )
    cam = {"point": [0, 0, 1 + f], "normal": [0, 0, -1], "half_size": [half, half]}
    g.add(Node("cam", "optics.detector", params={**cam, "pixels": [n, n]}))
    g.connect(("src", "out"), ("lens", "in"))
    g.connect(("lens", "out"), ("cam", "in"), transport)
    return Executor(g, reg, ev, available=WITH_SPECTRA)


def richards_wolf_gaussian(
    points: FloatArray, f: float, k: float, sin_max: float, waist: float
) -> ComplexArray:
    """Focal field of an aplanatic lens for x-polarised light with the pupil amplitude
    exp(−h²/w²), h = f sin θ (Richards and Wolf 1959: I₀, I₁, I₂ with the apodisation)."""
    rho = np.hypot(points[:, 0], points[:, 1])
    phi = np.arctan2(points[:, 1], points[:, 0])
    xg, wg = np.polynomial.legendre.leggauss(400)
    t_max = np.arcsin(sin_max)
    t, w = 0.5 * t_max * (xg + 1), 0.5 * t_max * wg
    c, s = np.cos(t), np.sin(t)
    weight = w * np.sqrt(c) * np.exp(-((f * s) ** 2) / waist**2) * s
    x = k * np.outer(rho, s)
    i0 = (weight * (1 + c) * jv(0, x)).sum(axis=1)
    i1 = (weight * s * jv(1, x)).sum(axis=1)
    i2 = (weight * (1 - c) * jv(2, x)).sum(axis=1)
    pre = -0.5j * k * f
    return np.stack(
        [
            pre * (i0 + i2 * np.cos(2 * phi)),
            pre * i2 * np.sin(2 * phi),
            -2j * pre * i1 * np.cos(phi),
        ],
        axis=1,
    )


def test_w79_gaussian_focus_via_s_equals_via_h() -> None:
    waist, n, half = 0.5, 41, 2.5e-3
    ex = focus_graph("S", waist)
    via_s: Image = ex.evaluate("cam", "image")
    key = "lens.out -> cam.in"
    assert any(key in line and "H -> S by the PFT" in line for line in ex.context.log)
    assert any("image from the sum of" in line for line in ex.context.log)
    assert not any("S -> H" in line for line in ex.context.log)  # the detector takes S
    via_h: Image = focus_graph("H", waist).evaluate("cam", "image")
    assert np.abs(via_s.values - via_h.values).max() < 1e-3 * via_h.values.max()  # 7.5e-4
    k = 2 * np.pi / 0.6328e-3  # the source's default wavelength
    plane = Plane.make([0, 0, 3.0], [0, 0, -1])
    pts = pixel_points(plane, half, n)[:, :2]
    ref = np.sum(np.abs(richards_wolf_gaussian(pts, 2.0, k, 0.5, waist)) ** 2, axis=1)
    ref = ref.reshape(n, n)
    # the sum of the waves is the Debye integral without approximation: closer to Richards–Wolf
    # than the H image (measured 6.8e-5 against 7.8e-4)
    assert np.abs(via_s.values - ref).max() < 2e-4 * ref.max()
    assert via_s.ray_power == pytest.approx(via_h.ray_power, rel=2e-3)


def test_w79_top_hat_focus_against_richards_wolf() -> None:
    # W13 with S on the link (R30a gave the same image as H by S -> H; now the sum of the waves)
    fieldtracing_focus = pytest.importorskip("fieldtracing.focus")
    n, half = 41, 2.5e-3
    im: Image = focus_graph("S", rings=30).evaluate("cam", "image")
    pts = pixel_points(Plane.make([0, 0, 3.0], [0, 0, -1]), half, n)[:, :2]
    k = 2 * np.pi / 0.6328e-3
    ref = np.sum(np.abs(fieldtracing_focus.richards_wolf(pts, 2.0, k, 0.5)) ** 2, axis=1)
    ref = ref.reshape(n, n)
    assert im.values.max() == pytest.approx(ref.max(), rel=1e-3)
    assert np.abs(im.values - ref).max() / ref.max() < 1e-3  # measured 4e-4 (H: 1.2e-3)


def test_w79_mixed_fields_and_spectra_interfere() -> None:
    # the same field once as H and once as S in one coherence group: four times the intensity
    ex = focus_graph("S", rings=12, n=9)
    lens = ex.evaluate("lens", "out")[0]
    from rtt_fusion.graph.field import rays_to_field

    fb = rays_to_field(lens)
    sb = sp.field_to_spectrum(fb)
    plane = Plane.make([0, 0, 3.0], [0, 0, -1])
    one = sp.spectrum_image([sb], [], plane, 2e-3, 2e-3, 9, 9)
    both = sp.spectrum_image([sb], [fb], plane, 2e-3, 2e-3, 9, 9)
    assert len(both.groups) == 1
    assert np.abs(both.values - 4 * one.values).max() < 2e-2 * 4 * one.values.max()


def test_w79_far_from_the_focus_the_detector_goes_through_h() -> None:
    # 0.5 mm behind the focus the cone is 0.6 mm wide: the waves (Δk ≈ k·NA/12) would alias, the
    # spectrum goes back to a field and is propagated as H – the same image as with H
    ex = focus_graph("S", rings=12, n=21, half=0.4)
    ex.graph.nodes["cam"].params["point"] = [0, 0, 3.5]
    via_s: Image = ex.evaluate("cam", "image")
    assert any("would alias" in line for line in ex.context.log)
    hx = focus_graph("H", rings=12, n=21, half=0.4)
    hx.graph.nodes["cam"].params["point"] = [0, 0, 3.5]
    via_h: Image = hx.evaluate("cam", "image")
    assert np.abs(via_s.values - via_h.values).max() < 1e-9 * via_h.values.max()


def test_auto_plans_with_native_detectors() -> None:
    # a detector takes S natively, but "auto" takes S only where the block in front takes it too
    # and R is not cheaper: the existing plans stay R; a spectrum from fusion.convert reaches the
    # detector as it is (as fields did on R links since R30a) and is summed there
    from rtt_fusion.graph import conversion

    ex = focus_graph("auto", rings=12, n=9)
    ex.check()
    assert ex.plan is not None and set(ex.plan.representations.values()) == {"R"}
    ev = {**ex.evaluators, **conversion.register(ex.registry)}
    g = ex.graph
    g.add(Node("toS", "fusion.convert", params={"to": "S"}))
    g.disconnect(g.links[1])
    g.connect(("lens", "out"), ("toS", "in"))
    g.connect(("toS", "out"), ("cam", "in"))
    ex = Executor(g, ex.registry, ev, available=WITH_SPECTRA)
    im: Image = ex.evaluate("cam", "image")
    assert ex.plan is not None and set(ex.plan.representations.values()) == {"R"}
    assert any("image from the sum of" in line for line in ex.context.log)
    assert not any("S -> H" in line for line in ex.context.log)
    assert im.values.max() > 0


def test_a_sum_beyond_the_budget_goes_through_h(monkeypatch: Any) -> None:
    monkeypatch.setattr(sp, "SUM_BUDGET", 1000)
    ex = focus_graph("S", rings=12, n=9)
    via_s: Image = ex.evaluate("cam", "image")
    assert any("would cost" in line for line in ex.context.log)
    via_h: Image = focus_graph("H", rings=12, n=9).evaluate("cam", "image")
    assert np.abs(via_s.values - via_h.values).max() < 1e-9 * via_h.values.max()


# ------------------------------------------------------------------------------ W80 propagation


def gridded_gaussian(w0: float, n: int, step: float) -> FieldBundle:
    """A Gaussian exp(−ρ²/w₀²) (x-polarised) at its waist on a centred grid (representation G)."""
    x = (np.arange(n) - n // 2) * step
    xx, yy = np.meshgrid(x, x)
    values = np.zeros((n, n, 3), dtype=np.complex128)
    values[..., 0] = np.exp(-(xx**2 + yy**2) / w0**2)
    plane = Plane.make([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
    return FieldBundle(KEY, plane, GriddedField(step, values, LAM, 1.0))


def gaussian_beam(points: FloatArray, w0: float, z: float) -> ComplexArray:
    """The Fresnel integral of exp(−ρ²/w₀²) over z: the paraxial Gaussian beam with e^(ikz)."""
    zr = K * w0**2 / 2
    w = w0 * np.sqrt(1 + (z / zr) ** 2)
    r2 = points[:, 0] ** 2 + points[:, 1] ** 2
    inv_r = z / (z * z + zr * zr)
    gouy = np.arctan(z / zr)
    return np.asarray(
        (w0 / w) * np.exp(-r2 / w**2) * np.exp(1j * (K * z + 0.5 * K * r2 * inv_r - gouy))
    )


@pytest.mark.parametrize("z_over_zr", [0.5, 1.0, 2.0])
def test_w80_gaussian_beam_propagated_in_the_spectrum(z_over_zr: float) -> None:
    w0 = 0.02
    zr = K * w0**2 / 2
    z = z_over_zr * zr
    sb = sp.field_to_spectrum(gridded_gaussian(w0, 256, w0 / 8))  # window 32 w0
    notes: list[str] = []
    target = Plane.make([0.0, 0.0, z], [0.0, 0.0, 1.0])
    moved = sp.propagate(sb, target, notes)
    assert moved.plane.origin[2] == z and moved.support is None
    assert any("keeps no sampling" in line for line in notes)  # a grid spectrum
    assert sp.spectrum_power(moved) == pytest.approx(sp.spectrum_power(sb), rel=1e-12)
    w = w0 * np.sqrt(1 + z_over_zr**2)
    pts = pixel_points(target, 2 * w, 33)
    ours = sp.wave_sum(moved, pts)[:, 0]
    ref = gaussian_beam(pts, w0, z)
    assert np.abs(ours - ref).max() < 1e-4  # measured 1.6e-5 (non-paraxial terms)
    # moving is exact: the waves on the old plane give the same field at the same points
    assert np.abs(sp.wave_sum(sb, pts)[:, 0] - ours).max() < 1e-11


@pytest.mark.parametrize(("z_over_zr", "aliased"), [(1.0, False), (2.0, True)])
def test_w80_gaussian_beam_detector_behind_a_spectrum_link(z_over_zr: float, aliased: bool) -> None:
    # W14 with S: H -> S by an FFT (Gaussian at its waist), the detector sums the waves on its
    # pixels by two matrix products (grid spectrum, parallel planes); at 2 z_R the beam would
    # meet its periodic replica of the grid window, and the detector goes through H instead
    wl = 0.6328
    w0 = 0.05
    zr = np.pi * w0**2 / (wl * 1e-3)
    z = z_over_zr * zr
    w = w0 * np.sqrt(1 + z_over_zr**2)
    half, n = 2.0 * w, 41
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    src = {"radius": 2.5 * w0, "rings": 25, "gaussian_waist": w0}
    g.add(Node("src", "optics.collimated_source", params=src))
    cam = {"point": [0, 0, z], "normal": [0, 0, -1], "half_size": [half, half], "pixels": [n, n]}
    g.add(Node("cam", "optics.detector", params=cam))
    g.connect(("src", "out"), ("cam", "in"), "S")
    ex = Executor(g, reg, ev, available=WITH_SPECTRA)
    im: Image = ex.evaluate("cam", "image")
    assert any("H -> S by an FFT" in line for line in ex.context.log)
    assert any("would alias" in line for line in ex.context.log) == aliased
    u = (np.arange(n) + 0.5) / n * 2 * half - half
    uu, vv = np.meshgrid(u, u)
    ref = (w0 / w) ** 2 * np.exp(-2 * (uu**2 + vv**2) / w**2)
    # the source is cut at 2.5 w0 (amplitude 2e−3 there); measured 3.8e-3 at z_R (H: 7.8e-3)
    assert np.abs(im.values - ref).max() / ref.max() < 1e-2
    assert im.ray_power == pytest.approx(np.pi * w0**2 / 2, rel=1e-3)


def test_w80_grid_sum_equals_the_sum_of_the_waves() -> None:
    sb = sp.field_to_spectrum(gridded_gaussian(0.02, 64, 0.02 / 4))
    for normal, up in (([0, 0, -1], [0, 1, 0]), ([0, 0, 1], [1, 0, 0]), ([0, 0, -1], [-1, 0, 0])):
        plane = Plane.make([0.003, -0.002, 0.4], normal, up)
        fast = sp.grid_sum(sb, plane, 0.04, 0.03, 12, 10)
        assert fast is not None
        pts = pixel_centres(plane, 0.04, 0.03, 12, 10)
        direct = sp.wave_sum(sb, pts)
        assert np.abs(fast - direct).max() < 1e-12 * np.abs(direct).max()
    turned = Plane.make([0.0, 0.0, 0.4], [0, 0, 1], [1, 1, 0])  # axes at 45°: no product form
    assert sp.grid_sum(sb, turned, 0.04, 0.03, 12, 10) is None


def gaussian_far(z: float, w0: float = 0.005, samples: int = 61) -> FieldBundle:
    """The hybrid field of a Gaussian beam z behind its waist (as in test_spectrum_links)."""
    zr = K * w0**2 / 2
    q = z - 1j * zr
    w = w0 * np.sqrt(1 + (z / zr) ** 2)
    x = np.linspace(-3.5 * w, 3.5 * w, samples)
    xx, yy = np.meshgrid(x, x)
    p = np.column_stack([xx.ravel(), yy.ravel()])
    r2 = np.sum(p * p, axis=1)
    e = (-1j * zr / q) * np.exp(1j * K * r2 / (2 * q))
    curvature = K * z / (z * z + zr * zr)
    psi = 0.5 * curvature * r2
    grad = curvature * p
    kz = np.sqrt(K * K - np.sum(grad * grad, axis=1))
    u = e * np.exp(-1j * psi)
    amp = np.column_stack([u, np.zeros_like(u), -u * grad[:, 0] / kz]).astype(np.complex128)
    hess = np.tile(curvature * np.eye(2), (len(p), 1, 1))
    plane = Plane.make([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
    return FieldBundle(KEY, plane, GridlessField(p, psi, grad, hess, amp, LAM, 1.0))


def test_w80_pft_spectrum_moved_and_back_is_the_geometric_field() -> None:
    zr = K * 0.005**2 / 2
    fb = gaussian_far(100 * zr)
    sb = sp.field_to_spectrum(fb)
    assert sp._is_pft(sb.support)
    dz = 50 * zr
    moved = sp.propagate(sb, Plane.make([0.0, 0.0, dz], [0.0, 0.0, 1.0]))
    back = sp.spectrum_to_field(moved).field  # the inverse PFT on the new plane
    assert isinstance(back, GridlessField) and isinstance(fb.field, GridlessField)
    geo = fb.field.propagate(dz)  # Fieldtracing: rays along ∇ψ, amplitude by the Jacobian
    assert np.abs(np.asarray(back.points) - np.asarray(geo.points)).max() < 1e-12
    assert np.abs(np.asarray(back.phase) - np.asarray(geo.phase)).max() < 1e-8
    amp, ref = np.asarray(back.amplitude), np.asarray(geo.amplitude)
    assert np.abs(amp - ref).max() < 1e-10 * np.abs(ref).max()


def focus_spectrum() -> sp.SpectrumBundle:
    """The spectrum of the converging beam behind an aplanatic lens (NA 0.5)."""
    from rtt_fusion.graph.field import rays_to_field

    lens = focus_graph("S", rings=12, n=9).evaluate("lens", "out")[0]
    return sp.field_to_spectrum(rays_to_field(lens))


def test_w80_moved_spectrum_is_the_same_field_in_space() -> None:
    sb = focus_spectrum()
    rng = np.random.default_rng(1)
    pts = rng.normal(size=(30, 3)) * np.array([2e-3, 2e-3, 1e-2]) + np.array([0.0, 0.0, 3.0])
    field = sp.wave_sum(sb, pts)
    for plane in (
        Plane.make([0.01, -0.02, 2.9], [0, 0, 1], [1, 1, 0]),  # parallel, turned in the plane
        Plane.make([0.01, -0.02, 2.9], [0, 0, -1]),  # parallel, normal against the light
        Plane.make([0.01, 0.02, 2.9], [np.sin(0.2), 0, np.cos(0.2)]),  # tilted, all forward
    ):
        notes: list[str] = []
        moved = sp.propagate(sb, plane, notes)
        assert not any("exec.spectrum_cut" in line for line in notes)
        assert float(moved.plane.normal @ np.array([0, 0, 1.0])) > 0  # along the light
        assert np.abs(sp.wave_sum(moved, pts) - field).max() < 1e-10 * np.abs(field).max()
        assert sp.spectrum_power(moved) == pytest.approx(sp.spectrum_power(sb), rel=1e-12)


def test_w80_tilted_plane_drops_waves_that_do_not_cross_it() -> None:
    sb = focus_spectrum()  # directions up to 30° off the axis
    tilt = np.radians(70.0)
    normal = np.array([np.sin(tilt), 0.0, np.cos(tilt)])
    notes: list[str] = []
    moved = sp.propagate(sb, Plane.make([0, 0, 2.5], normal), notes)
    w = sb.waves
    kz = np.real(w.port.kz(w.kx, w.ky))
    k_global = np.column_stack([w.kx, w.ky, kz]) @ np.stack(
        [sb.plane.u, sb.plane.v, sb.plane.normal]
    )
    behind = k_global @ normal <= 0
    assert 0 < behind.sum() < len(behind) and len(moved.waves.kx) == int((~behind).sum())
    terms = np.sum(np.abs(sp._field3(w)) ** 2, axis=1) / np.asarray(w.solid_angle)
    cut = [line for line in notes if line.startswith("exec.spectrum_cut")]
    assert len(cut) == 1
    assert float(cut[0].split()[1]) == pytest.approx(terms[behind].sum() / terms.sum(), rel=1e-2)


def test_w80_evanescent_wave_decays_and_is_dropped_where_it_cannot_go() -> None:
    port = Port(0.0, 1, 1.0, LAM)
    kx = np.array([0.5 * K, 1.2 * K])
    field = np.array([[0.0, 1.0], [0.0, 1.0]], dtype=complex)
    plane = Plane.make([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
    sb = sp.SpectrumBundle(KEY, plane, PlaneWaves(port, kx, np.zeros(2), field))
    dz = 2e-4
    moved = sp.propagate(sb, Plane.make([0.0, 0.0, dz], [0.0, 0.0, 1.0]))
    kappa = K * np.sqrt(1.2**2 - 1)
    assert abs(moved.waves.field[1, 1]) == pytest.approx(np.exp(-kappa * dz), rel=1e-12)
    assert moved.waves.field[0, 1] == pytest.approx(np.exp(1j * K * np.sqrt(0.75) * dz))
    for target in (
        Plane.make([0.0, 0.0, -dz], [0.0, 0.0, 1.0]),  # behind: it would grow
        Plane.make([0.0, 0.0, dz], [0.1, 0.0, 1.0]),  # tilted: no real direction to turn
    ):
        notes: list[str] = []
        out = sp.propagate(sb, target, notes)
        assert len(out.waves.kx) == 1 and any("exec.spectrum_cut" in n for n in notes)


# ------------------------------------------------------------------------------ Raytatouille


def rtt_image(transport: str, mode: str) -> tuple[Image, Executor]:
    from rtt_fusion.graph import backend_rtt, graph_from_raytatouille

    rt = pytest.importorskip("raytatouille")
    path = DATA / "paraboloid_f50.rtt.json"
    g = graph_from_raytatouille(str(path))
    g.nodes["focus"].params["half_size"] = [0.1, 0.1]
    for link in g.links_to("focus"):
        g.set_transport(link, transport)
    shared = {"system": rt.load(str(path)), "rings": 20, "pixels": [31, 31], "rtt_execution": mode}
    ex = Executor(g, builtin_registry(), backend_rtt.EVALUATORS, shared, available=WITH_SPECTRA)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        im: Image = ex.evaluate("focus", "image")
    return im, ex


def test_w79_raytatouille_detector_takes_the_spectrum() -> None:
    # R30a: with S on the link the detector of whole paths computed with rays; now the fields
    # after the last surface become spectra and their waves are summed in the focus. Reference:
    # the Airy pattern of W17 (the H image of the same link differs by 5e-4, at 150 times the
    # cost: Fieldtracing's Debye focus on a large grid)
    via_s, ex = rtt_image("S", "paths")
    assert any("H -> S by the PFT" in line for line in ex.context.log)
    assert any("image from the sum of" in line for line in ex.context.log)
    assert len(via_s.groups) == 2  # unpolarised: two input polarisations
    n, half, wavelength_mm = 31, 0.1, 0.5876e-3
    u = (np.arange(n) + 0.5) / n * 2 * half - half
    uu, vv = np.meshgrid(u, u)
    na = np.sin(2 * np.arctan(2.0 / (4 * 100.0)))
    v = 2 * np.pi / wavelength_mm * na * np.hypot(uu, vv)
    airy = np.where(v > 0, (2 * jv(1, v) / np.where(v > 0, v, 1.0)) ** 2, 1.0)
    assert np.abs(via_s.values / via_s.values.max() - airy).max() < 1e-3  # measured 5.6e-4
    peak = np.pi * (np.pi * na**2 / wavelength_mm**2)
    assert via_s.values.max() == pytest.approx(peak, rel=1e-2)
    nodes, ex = rtt_image("S", "nodes")  # node-wise: the executor converts on the link
    assert any("H -> S by the PFT" in line for line in ex.context.log)
    assert np.abs(nodes.values - via_s.values).max() < 1e-9 * via_s.values.max()


def test_sockets_of_imported_systems_take_spectra() -> None:
    from rtt_fusion.graph import graph_from_raytatouille

    d = json.loads((DATA / "paraboloid_f50.rtt.json").read_text(encoding="utf-8"))
    g = graph_from_raytatouille(d)
    reg = builtin_registry()
    from rtt_fusion.graph.validate import node_sockets

    for node in g.nodes.values():
        for s in node_sockets(node, reg):
            if s.kind.name == "LIGHT":
                assert "S" in s.accepts, (node.id, s.name)
