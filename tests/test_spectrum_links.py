"""Representation S on links (R30a): plane-wave spectra between blocks.

References, independent of the code under test:
- W76 H -> S: the angular spectrum of a paraxial Gaussian beam at distance z from its waist,
  Ẽ(κ) = w₀²/(4π)·exp(−κ²w₀²/4)·exp(−iκ²z/(2k)) for E(ρ) = ∫ Ẽ e^(iκ·ρ) d²κ (Fourier transform
  of exp(−ρ²/w₀²), Fresnel propagator). The PFT is a stationary-phase evaluation: for the
  Gaussian its prefactor is 1/(1 + i z_R/z) off the exact one (the Gaussian integral with the
  complex width 1/w² − ik/(2R) against its imaginary part alone), so its error is z_R/z – the
  test checks that value, i.e. the order O(1/F) and its constant. At the waist (flat phase) the
  FFT; a gridded field (G -> S) with the discrete Fourier transform of a well-sampled Gaussian.
  Power π w₀²/2 (unit amplitude). The evanescent share of a Gaussian narrower than λ,
  exp(−k²w₀²/2) (radial integral of |Ẽ|²), against `exec.spectrum_cut`.
- W77 round trips: S -> H -> S and H -> S -> H give back what went in (the inverse PFT is the
  Legendre inverse of the PFT; the FFT path interpolates); R -> S -> R equals R -> H -> R and the
  spectrum carries the rays' power; a link with S into a detector gives the image of a link with
  H (since R30b the detector sums the waves itself: to 1e−3, see test_spectrum_detector.py).
- W78 focus of an aplanatic lens (NA 0.8) from the spectrum of its rays (R -> H -> S): the
  Richards–Wolf spectrum of `couplings.focus` (Gauss–Legendre, no rays) as in W29, and the
  detector image of a link with S against Fieldtracing's `focus.richards_wolf` (Bessel
  integrals) as in W13.
- Planner (decision 3 of docs/plan-s-spektren.md): "auto" takes S only between two blocks that
  both take S natively; existing plans stay as they were.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.graph import (
    Category,
    Graph,
    Node,
    NodeType,
    NodeTypeRegistry,
    Requirements,
    builtin_registry,
    conversion,
    graph_from_raytatouille,
    optics,
    plan,
)
from rtt_fusion.graph.execute import WITH_FIELDS, WITH_SPECTRA, ExecutionError, Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle
from rtt_fusion.graph.model import data_out, light_in, light_out
from rtt_fusion.graph.planner import link_key
from rtt_fusion.ports import PlaneWaves, Port

pytest.importorskip("fieldtracing")
from fieldtracing.field import GriddedField, GridlessField  # noqa: E402

from rtt_fusion.couplings.focus import richards_wolf  # noqa: E402
from rtt_fusion.graph import spectrum as sp  # noqa: E402
from rtt_fusion.graph.field import (  # noqa: E402
    FieldBundle,
    FieldError,
    field_to_rays,
    rays_to_field,
)
from rtt_fusion.graph.particle import ray_powers  # noqa: E402

DATA = Path(__file__).parent / "data"
LAM = 0.6  # µm
K = 2 * np.pi / (LAM * 1e-3)  # 1/mm
W0 = 0.005  # waist, mm
ZR = K * W0**2 / 2  # Rayleigh length, mm
KEY = CoherenceKey("beam")
PLANE = Plane.make([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])


def gaussian(z: float, w0: float = W0, samples: int = 81, extent: float = 3.5) -> FieldBundle:
    """Paraxial Gaussian beam (x-polarised, unit amplitude at the waist) z mm behind its waist as a
    hybrid field on a square of ±extent·w(z): E = (−i z_R/q)·exp(ikρ²/(2q)), q = z − i z_R (the
    factor e^(ikz) left out), split into ψ = kρ²/(2R) and U = E·e^(−iψ); E_z from transversality."""
    zr = K * w0**2 / 2
    q = z - 1j * zr
    w = w0 * np.sqrt(1 + (z / zr) ** 2)
    x = np.linspace(-extent * w, extent * w, samples)
    xx, yy = np.meshgrid(x, x)
    p = np.column_stack([xx.ravel(), yy.ravel()])
    r2 = np.sum(p * p, axis=1)
    e = (-1j * zr / q) * np.exp(1j * K * r2 / (2 * q))
    curvature = 0.0 if z == 0 else K * z / (z * z + zr * zr)
    psi = 0.5 * curvature * r2
    grad = curvature * p
    kz = np.sqrt(K * K - np.sum(grad * grad, axis=1))
    u = e * np.exp(-1j * psi)
    amp = np.column_stack([u, np.zeros_like(u), -u * grad[:, 0] / kz]).astype(np.complex128)
    hess = np.tile(curvature * np.eye(2), (len(p), 1, 1))
    return FieldBundle(KEY, PLANE, GridlessField(p, psi, grad, hess, amp, LAM, 1.0))


def angular_spectrum(kx: FloatArray, ky: FloatArray, z: float, w0: float = W0) -> ComplexArray:
    kap2 = kx * kx + ky * ky
    return np.asarray(
        w0**2 / (4 * np.pi) * np.exp(-kap2 * w0**2 / 4) * np.exp(-1j * kap2 * z / (2 * K))
    )


def fourier(sb: sp.SpectrumBundle) -> ComplexArray:
    """Ẽ_x(k_j) of the waves: field_j / Δ²k_j with Δ²k = dΩ·k·k_z."""
    w = sb.waves
    assert w.solid_angle is not None
    dk2 = np.asarray(w.solid_angle) * K * np.real(w.port.kz(w.kx, w.ky))
    return np.asarray(w.field[:, 0] / dk2)


def relative(a: ComplexArray, b: ComplexArray) -> float:
    return float(np.abs(a - b).max() / np.abs(b).max())


# ------------------------------------------------------------------------------ W76 H -> S


@pytest.mark.parametrize("z_over_zr", [100.0, 1000.0])
def test_w76_gaussian_by_the_pft(z_over_zr: float) -> None:
    fb = gaussian(z_over_zr * ZR)
    notes: list[str] = []
    sb = sp.field_to_spectrum(fb, notes=notes)  # "auto": F = 2z/z_R ≥ 1/ε
    assert any("by the PFT" in n for n in notes)
    assert sb.continuous and sb.plane is PLANE and sb.key == KEY
    err = relative(fourier(sb), angular_spectrum(sb.waves.kx, sb.waves.ky, z_over_zr * ZR))
    # stationary phase: the prefactor is 1/(1 + i z_R/z) off, |1 − that| = z_R/z (to O((z_R/z)²))
    assert err == pytest.approx(1.0 / z_over_zr, rel=2e-2)
    assert sp.spectrum_power(sb) == pytest.approx(np.pi * W0**2 / 2, rel=1e-3)


def test_w76_gaussian_at_the_waist_by_the_fft() -> None:
    fb = gaussian(0.0)
    notes: list[str] = []
    sb = sp.field_to_spectrum(fb, notes=notes)  # flat phase: no PFT, an FFT on a grid
    assert any("by an FFT" in n for n in notes)
    assert not any("exec.spectrum_cut" in n for n in notes)
    err = relative(fourier(sb), angular_spectrum(sb.waves.kx, sb.waves.ky, 0.0))
    assert err < 1e-4  # measured 2.2e-5 (cubic interpolation of the samples onto the grid)
    assert sp.spectrum_power(sb) == pytest.approx(np.pi * W0**2 / 2, rel=1e-3)


def test_w76_gridded_field_by_ports() -> None:
    # G -> S: a Gaussian on a centred grid (16 w0 wide, 8 points per w0) through
    # ports.to_plane_waves; the discrete transform of a well-sampled Gaussian is its continuous
    # one up to aliasing exp(−(π/step)²w0²/4) and the window, both far below rounding
    n, step = 128, W0 / 8
    x = (np.arange(n) - n // 2) * step
    xx, yy = np.meshgrid(x, x)
    values = np.zeros((n, n, 3), dtype=np.complex128)
    values[..., 0] = np.exp(-(xx**2 + yy**2) / W0**2)
    fb = FieldBundle(KEY, PLANE, GriddedField(step, values, LAM, 1.0))
    sb = sp.field_to_spectrum(fb)
    assert relative(fourier(sb), angular_spectrum(sb.waves.kx, sb.waves.ky, 0.0)) < 1e-12
    assert sp.spectrum_power(sb) == pytest.approx(np.pi * W0**2 / 2, rel=1e-6)
    back = sp.spectrum_to_field(sb).field
    assert isinstance(back, GriddedField)
    # back on the grid: up to the waves below SPECTRUM_THRESHOLD (1e−8 of the largest) left out
    assert np.abs(np.asarray(back.values)[..., :2] - values[..., :2]).max() < 1e-7


def test_w76_evanescent_share_is_reported() -> None:
    # a waist narrower than λ: a share of the spectrum is evanescent and left out. Reference: the
    # analytic spectrum on the frequency grid of the FFT (the window is 15 w0, Δk = k/5, so the
    # sum over the grid differs from the continuous share exp(−k²w0²/2) = 0.112 by the bins at
    # the circle |k∥| = k)
    w0 = 0.2e-3
    fb = gaussian(0.0, w0=w0, samples=61)
    notes: list[str] = []
    sb = sp.field_to_spectrum(fb, notes=notes)
    cut = [n for n in notes if n.startswith("exec.spectrum_cut")]
    assert len(cut) == 1
    share = float(cut[0].split()[1])
    assert isinstance(sb.support, sp.GridSupport)
    kx, ky = sb.support.sampled.grid.frequencies()
    k2 = kx * kx + ky * ky
    density = np.abs(angular_spectrum(kx, ky, 0.0, w0)) ** 2
    assert share == pytest.approx(density[k2 >= K * K].sum() / density.sum(), rel=1e-2)
    assert share == pytest.approx(np.exp(-(K**2) * w0**2 / 2), rel=0.3)
    assert np.all(np.hypot(sb.waves.kx, sb.waves.ky) < K)


# ------------------------------------------------------------------------------ W77 round trips


@pytest.mark.parametrize(("z_over_zr", "tol"), [(100.0, 1e-12), (0.0, 1e-4)])
def test_w77_spectrum_field_spectrum(z_over_zr: float, tol: float) -> None:
    # PFT: the inverse PFT is the Legendre inverse, exact; FFT: interpolation (measured 1.7e-5;
    # H -> S -> H below 2.8e-4)
    sb = sp.field_to_spectrum(gaussian(z_over_zr * ZR))
    fb = sp.spectrum_to_field(sb)
    assert isinstance(fb.field, GridlessField) and fb.field.domain == "space"
    again = sp.field_to_spectrum(fb)
    assert np.allclose(again.waves.kx, sb.waves.kx, rtol=0, atol=1e-9 * K)
    assert relative(again.waves.field, sb.waves.field) < tol
    assert sp.spectrum_power(again) == pytest.approx(sp.spectrum_power(sb), rel=1e-3)


@pytest.mark.parametrize(("z_over_zr", "tol"), [(100.0, 1e-12), (0.0, 1e-3)])
def test_w77_field_spectrum_field(z_over_zr: float, tol: float) -> None:
    fb = gaussian(z_over_zr * ZR)
    back = sp.spectrum_to_field(sp.field_to_spectrum(fb)).field
    f = fb.field
    assert isinstance(back, GridlessField) and isinstance(f, GridlessField)
    assert np.abs(np.asarray(back.points) - np.asarray(f.points)).max() < 1e-12
    assert np.abs(np.asarray(back.phase) - np.asarray(f.phase)).max() < 1e-9
    assert relative(np.asarray(back.amplitude)[:, 0], np.asarray(f.amplitude)[:, 0]) < tol


def lens_rays(na: float = 0.8, f: float = 2.0, rings: int = 30) -> RayBundle:
    """Rays behind an aplanatic lens (unit irradiance on a pupil of radius NA·f)."""
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    a = na * f
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": a, "rings": rings, "center": [0, 0, -1], "wavelengths_um": [LAM]},
        )
    )
    # the aperture a hair larger than the beam: no ray of the rim lost by rounding
    g.add(
        Node(
            "lens",
            "optics.aplanatic_lens",
            params={"focal_length": f, "aperture_radius": a * (1 + 1e-6)},
        )
    )
    g.connect(("src", "out"), ("lens", "in"))
    rays: RayBundle = Executor(g, reg, ev).evaluate("lens", "out")[0]
    return rays


def test_w77_rays_spectrum_rays() -> None:
    rays = lens_rays()
    notes: list[str] = []
    sb = sp.rays_to_spectrum(rays, notes=notes)
    assert any("by the PFT" in n for n in notes)  # a converging bundle: large Fresnel numbers
    power = float(ray_powers(rays, np.asarray(rays.alive)).sum())
    assert power == pytest.approx(np.pi * 1.6**2, rel=1e-3)  # the source power
    assert sp.spectrum_power(sb) == pytest.approx(power, rel=1e-3)
    out = sp.spectrum_to_rays(sb, notes)
    ref = field_to_rays(rays_to_field(rays))
    for name in ("position", "direction", "opl", "amplitude", "source_xy"):
        a, b = np.asarray(getattr(out, name)), np.asarray(getattr(ref, name))
        assert np.abs(a - b).max() <= 1e-12 * max(1.0, float(np.abs(b).max())), name
    assert not any("exec.spectrum_rays" in n for n in notes)
    # the directions are those of the rays that went in
    alive = np.asarray(rays.alive)
    assert np.abs(np.asarray(out.direction) - np.asarray(rays.direction)[alive]).max() < 1e-12


def test_w77_narrow_bundle_as_rays_is_reported() -> None:
    fb = gaussian(0.0, w0=0.002)  # 2 µm: not large against λ = 0.6 µm
    notes: list[str] = []
    rays = sp.spectrum_to_rays(sp.field_to_spectrum(fb), notes)
    assert isinstance(fb.field, GridlessField)
    assert rays.size == len(np.asarray(fb.field.points))
    assert any(n.startswith("exec.spectrum_rays") for n in notes)


def test_w77_discrete_waves_do_not_become_a_field() -> None:
    port = Port(0.0, 1, 1.0, LAM)
    waves = PlaneWaves(port, np.zeros(1), np.zeros(1), np.array([[1.0, 0.0]], dtype=complex))
    with pytest.raises(FieldError, match="discrete"):
        sp.spectrum_to_field(sp.SpectrumBundle(KEY, PLANE, waves))


def focus_graph(transport: str, n: int = 41, half: float = 2.5e-3) -> Executor:
    """W13 with a link of the given transport into the detector in the focus (NA 0.5)."""
    f, a = 2.0, 1.0
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(Node("src", "optics.collimated_source", params={"radius": a, "rings": 30}))
    g.add(
        Node(
            "lens",
            "optics.aplanatic_lens",
            params={"focal_length": f, "aperture_radius": a * (1 + 1e-6), "point": [0, 0, 1]},
        )
    )
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, 1 + f],
                "normal": [0, 0, -1],
                "half_size": [half, half],
                "pixels": [n, n],
            },
        )
    )
    g.connect(("src", "out"), ("lens", "in"))
    g.connect(("lens", "out"), ("cam", "in"), transport)
    return Executor(g, reg, ev, available=WITH_SPECTRA)


def test_w77_w78_detector_behind_a_spectrum_link() -> None:
    fieldtracing_focus = pytest.importorskip("fieldtracing.focus")
    ex = focus_graph("S")
    im: Image = ex.evaluate("cam", "image")
    key = "lens.out -> cam.in"
    assert ex.plan is not None and ex.plan.representations[key] == "S"
    assert any(key in line and "H -> S by the PFT" in line for line in ex.context.log)
    # since R30b the detector takes S: the image is the sum of the waves (W79), no S -> H
    assert not any("S -> H" in line for line in ex.context.log)
    with_h: Image = focus_graph("H").evaluate("cam", "image")
    assert np.abs(im.values - with_h.values).max() < 1e-3 * with_h.values.max()  # 7.8e-4
    # W13: against Richards–Wolf (Bessel integrals)
    n, half = 41, 2.5e-3
    u = (np.arange(n) + 0.5) / n * 2 * half - half
    uu, vv = np.meshgrid(u, u)
    plane = Plane.make([0, 0, 3.0], [0, 0, -1])
    pts = uu.reshape(-1, 1) * plane.u[:2] + vv.reshape(-1, 1) * plane.v[:2]
    k_source = 2 * np.pi / (0.6328e-3)  # the source's default wavelength
    ref = np.sum(np.abs(fieldtracing_focus.richards_wolf(pts, 2.0, k_source, 0.5)) ** 2, axis=1)
    ref = ref.reshape(n, n)
    assert im.values.max() == pytest.approx(ref.max(), rel=2e-3)
    assert np.abs(im.values - ref).max() / ref.max() < 5e-3


# ------------------------------------------------------------------------------ W78 focus


def focal_field(w: PlaneWaves, points: FloatArray, plane: Plane | None = None) -> ComplexArray:
    """Σ field_j e^(i k_j·r) at points (global; local in `plane` if given)."""
    rot = np.eye(3) if plane is None else np.stack([plane.u, plane.v, plane.normal])
    origin = np.zeros(3) if plane is None else plane.origin
    kz = np.real(w.port.kz(w.kx, w.ky))
    ez = w.ez()
    out = []
    for p in points:
        q = rot @ (p - origin)
        ph = np.exp(1j * (w.kx * q[0] + w.ky * q[1] + kz * (q[2] - w.port.z)))
        out.append([np.sum(w.field[:, 0] * ph), np.sum(w.field[:, 1] * ph), np.sum(ez * ph)])
    return np.asarray(out, dtype=np.complex128)


def test_w78_focus_from_the_spectrum_against_richards_wolf() -> None:
    from rtt_fusion.ports import spectrum_power

    na, f = 0.8, 2.0
    sb = sp.rays_to_spectrum(lens_rays(na, f))
    assert sp.spectrum_power(sb) == pytest.approx(np.pi * (na * f) ** 2, rel=1e-3)
    rw = richards_wolf(Port(-0.002, 1, 1.0, LAM), na, (0, 0, 0), (1.0, 0.0), 40, 80)
    pts = np.array([[0, 0, 0], [0.3, 0, 0], [0, 0.3, 0], [0, 0, 0.5], [0.2, -0.1, -0.3]]) * 1e-3
    ours = focal_field(sb.waves, pts + np.array([0.0, 0.0, f]), sb.plane) / np.sqrt(
        spectrum_power(sb.waves)
    )
    theirs = focal_field(rw, pts) / np.sqrt(spectrum_power(rw))
    scale = theirs[0, 0] / ours[0, 0]
    assert abs(scale) == pytest.approx(1.0, abs=1e-3)  # measured 6e-5
    # as in W29: the Debye picture has −i between pupil and focus, couplings.focus leaves it out
    assert np.angle(scale) == pytest.approx(np.pi / 2, abs=1e-6)
    assert np.abs(ours * scale - theirs).max() < 2e-3 * abs(theirs[0, 0])  # measured 9.3e-4


# ------------------------------------------------------------------------------ conversion block


def test_convert_block_to_spectrum_and_back() -> None:
    reg = builtin_registry()
    ev = {**optics.register(reg), **conversion.register(reg)}
    g = Graph()
    g.add(Node("src", "optics.collimated_source", params={"radius": 0.25, "rings": 12}))
    g.add(
        Node(
            "lens",
            "optics.aplanatic_lens",
            params={"focal_length": 2.0, "aperture_radius": 0.3, "point": [0, 0, 1]},
        )
    )
    g.add(Node("toS", "fusion.convert", params={"to": "S"}))
    g.add(Node("toR", "fusion.convert", params={"to": "R"}))
    g.connect(("src", "out"), ("lens", "in"))
    g.connect(("lens", "out"), ("toS", "in"))
    g.connect(("toS", "out"), ("toR", "in"), "S")
    ex = Executor(g, reg, ev, available=WITH_SPECTRA)
    spectra = ex.evaluate("toS", "out")
    assert spectra and all(isinstance(s, sp.SpectrumBundle) for s in spectra)
    rays = ex.evaluate("toR", "out")
    lens = ex.evaluate("lens", "out")
    ref = field_to_rays(rays_to_field(lens[0]))
    assert isinstance(rays[0], RayBundle)
    assert np.abs(np.asarray(rays[0].direction) - np.asarray(ref.direction)).max() < 1e-12
    assert any("toS: src: H -> S" in line for line in ex.context.log)
    assert reg.get("fusion.convert") is not None
    to = next(p for p in reg.get("fusion.convert").params if p.name == "to")  # type: ignore[union-attr]
    assert to.default == "H"  # plan S6


# ------------------------------------------------------------------------------ planner


def registry(src: str, det: str) -> NodeTypeRegistry:
    reg = NodeTypeRegistry()
    reg.add(NodeType("t.src", Category.SOURCE, (light_out("out", src),)))
    sockets = (light_in("in", det), data_out("image", "result"))
    reg.add(NodeType("t.det", Category.DETECTOR, sockets))
    return reg


def chain() -> Graph:
    g = Graph()
    g.add(Node("src", "t.src"))
    g.add(Node("det", "t.det"))
    g.connect(("src", "out"), ("det", "in"))
    return g


@pytest.mark.parametrize(
    ("natives", "expected"),
    [
        ({"src": "RS", "det": "RS"}, "S"),  # both take S natively
        ({"src": "RS", "det": "R"}, None),  # one end would convert: S not chosen by "auto"
    ],
)
def test_auto_takes_s_only_between_native_blocks(
    natives: dict[str, str], expected: str | None
) -> None:
    g = chain()
    req = Requirements(available=WITH_SPECTRA, native=lambda n: frozenset(natives[n]))
    p = plan(g, registry("S", "S"), req)
    if expected is None:
        assert not p.ok and "S only between blocks that take it natively" in str(
            p.diagnostics[0].message
        )
    else:
        assert p.ok and p.representation(g.links[0]) == expected
    # a fixed S is the user's choice, natively or not
    g.set_transport(g.links[0], "S")
    assert plan(g, registry("S", "S"), req).representation(g.links[0]) == "S"
    # where R is admissible, "auto" stays R between blocks that take S if nothing says the light
    # is a spectrum (R30f: a source of rays; spectra are kept only behind S)
    g.set_transport(g.links[0], "auto")
    p = plan(g, registry("RS", "RS"), Requirements(available=WITH_SPECTRA, native=req.native))
    assert p.representation(g.links[0]) == "R"


def test_existing_plans_are_unchanged() -> None:
    from rtt_fusion.graph.workflow import available, native_of

    assert available() == WITH_SPECTRA
    g = graph_from_raytatouille(DATA / "linnik_structure.rtt.json")
    reg = builtin_registry()
    for native in (None, native_of(g)):
        p = plan(g, reg, Requirements(available=WITH_SPECTRA, native=native))
        assert p.ok and set(p.representations.values()) == {"R"} and not p.conversions
    # the analytic focus: "auto" is R with spectra available, as before
    ex = focus_graph("auto")
    ex.check()
    assert ex.plan is not None and set(ex.plan.representations.values()) == {"R"}


def test_transport_s_needs_spectra() -> None:
    ex = focus_graph("S")
    ex.available = WITH_FIELDS
    with pytest.raises(ExecutionError, match="exec.transport_unsupported"):
        ex.check()
    ex.available = WITH_SPECTRA
    ex.check()
    assert ex.plan is not None
    assert ex.plan.representations[link_key(ex.graph.links[1])] == "S"


def test_spectrum_link_into_a_ray_block() -> None:
    # S on a link into a block that takes rays only (a thin lens): S -> H -> R at its port
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(Node("src", "optics.collimated_source", params={"radius": 0.3, "rings": 12}))
    g.add(Node("lens", "optics.thin_lens", params={"point": [0, 0, 1], "focal_length": 50.0}))
    g.connect(("src", "out"), ("lens", "in"), "S")
    via_s: Any = Executor(g, reg, ev, available=WITH_SPECTRA).evaluate("lens", "out")
    g.set_transport(g.links[0], "H")
    via_h: Any = Executor(g, reg, ev, available=WITH_SPECTRA).evaluate("lens", "out")
    assert isinstance(via_s[0], RayBundle) and isinstance(via_h[0], RayBundle)
    for name in ("position", "direction", "opl"):
        a, b = np.asarray(getattr(via_s[0], name)), np.asarray(getattr(via_h[0], name))
        assert np.abs(a - b).max() < 1e-6 * max(1.0, float(np.abs(b).max())), name


# ------------------------------------------------------------------------------ used on links


def test_used_representation_per_link() -> None:
    # R30e: what a link carried in the run (`Context.used`), for the GUI's "planned R, used S"
    ex = focus_graph("S", n=9)
    ex.evaluate("cam", "image")
    used = ex.context.used
    assert "src.out -> lens.in" not in used  # plain rays on an R link: nothing to show
    s_link = used["lens.out -> cam.in"]
    assert s_link["representation"] == "S" and s_link["chain"] == "R→H→S"
    assert s_link["reason"] == "lens gave rays; spectra on the link; cam takes them as they are"
    # a spectrum from fusion.convert reaches the detector as it is; since R30f the link on
    # "auto" is planned S too (before: planned R, used S)
    ex = focus_graph("auto", n=9)
    ev = {**ex.evaluators, **conversion.register(ex.registry)}
    g = ex.graph
    g.add(Node("toS", "fusion.convert", params={"to": "S"}))
    g.disconnect(g.links[1])
    g.connect(("lens", "out"), ("toS", "in"))
    g.connect(("toS", "out"), ("cam", "in"))
    ex = Executor(g, ex.registry, ev, available=WITH_SPECTRA)
    ex.evaluate("cam", "image")
    assert ex.plan is not None and ex.plan.representations["toS.out -> cam.in"] == "S"
    assert ex.context.used["toS.out -> cam.in"] == {
        "representation": "S",
        "reason": "toS gave spectra; cam takes them as they are",
        "chain": "S",
    }
    assert "lens.out -> toS.in" not in ex.context.used


def test_used_representation_into_a_ray_block() -> None:
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(Node("src", "optics.collimated_source", params={"radius": 0.3, "rings": 6}))
    g.add(Node("lens", "optics.thin_lens", params={"point": [0, 0, 1], "focal_length": 50.0}))
    g.connect(("src", "out"), ("lens", "in"), "S")
    ex = Executor(g, reg, ev, available=WITH_SPECTRA)
    ex.evaluate("lens", "out")
    used = ex.context.used["src.out -> lens.in"]
    assert used["representation"] == "S" and used["chain"] == "R→H→S→H→R"
    assert used["reason"].endswith("converted to rays at the port of lens")


def test_with_used_carries_the_chain() -> None:
    from rtt_fusion.graph.workflow import with_used

    plan_dict = {"links": {"a.out -> b.in": {"representation": "R", "reason": "cheapest"}}}
    used = {"a.out -> b.in": {"representation": "S", "reason": "why", "chain": "S"}}
    shown = with_used(plan_dict, used)["links"]["a.out -> b.in"]
    assert (shown["used"], shown["used_reason"], shown["used_chain"]) == ("S", "why", "S")
    assert "used" not in plan_dict["links"]["a.out -> b.in"]  # the plan itself is unchanged
