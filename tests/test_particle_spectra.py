"""Particle with plane-wave spectra S (R30d, docs/plan-s-spektren.md): `fusion.particle` takes the
incident angular spectrum as it is and gives the backscattered far field as a spectrum.

References, independent of the code under test:
- W84 W29–W31 with S on the link into the particle. (a) The incident spectrum CliffordBem gets
  from an S link (R -> H -> S at the lens, moved onto the particle's port in the spectrum)
  against the Richards–Wolf spectrum of `couplings.focus` (Gauss–Legendre, no rays) as in W29:
  focal field after power normalisation, the constant +i of the Debye picture, power = the rays'.
  (b) A plane wave as a discrete spectrum on the sphere: the backscattered power in the half space
  against Mie (tests/mie_reference.py) with the Richardson extrapolation of W31, and equal to the
  ray path of W31 (the same plane wave) to the solver's rounding.
- W85 K3 without rays between the focus and the particle. (a) The far field of a dipole (analytic,
  E∞ = k²(r̂ × p) × r̂) as the particle gives it: S -> H on a plane at distance D equals the
  spherical wave E∞ e^(ikr)/r of the far-field expansion (stationary phase is exact for it) and
  differs from the exact dipole field by its near-field term, 1/(kD) relative on the axis; the
  power of the spectrum is the analytic one; rays from it carry the power per direction.
  (b) Pilot K3 in the graph with S on both links (source → lens → [S] → particle → [S] →
  detector): the power fraction backscattered into the NA against CliffordBem alone
  (`BeamField.focused`, far field, no adapter) as in W30, and the detector image in the far
  field against |E∞|²/r² of that same direct computation.
- Planner: the plans with "auto" around a particle stay as they were (R); a spectrum leaving the
  particle reaches a detector on an "auto" link as a spectrum.
- A continuous spectrum without sampling (behind a condensed group) into a block that needs rays:
  `exec.spectrum_support` (the open point of R30c), not an error deep in the conversion.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion import portdata
from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.convert import far_field_to_spectrum, hemisphere
from rtt_fusion.couplings.focus import richards_wolf
from rtt_fusion.graph import Category, Graph, Node, NodeType, builtin_registry, optics, particle
from rtt_fusion.graph.execute import WITH_SPECTRA, ExecutionError, Executor, ResultCache
from rtt_fusion.graph.light import CoherenceKey, Plane, RayBundle
from rtt_fusion.graph.model import light_out
from rtt_fusion.graph.planner import link_key
from rtt_fusion.ports import PlaneWaves, Port, spectrum_power

pytest.importorskip("fieldtracing")
from rtt_fusion.graph import spectrum as sp  # noqa: E402
from rtt_fusion.graph.field import FieldError, field_to_rays  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from mie_reference import S12  # noqa: E402
from test_particle import focal_field, focused  # noqa: E402

LAM, NA, F, R, N_P = 0.6, 0.8, 2.0, 0.25, 1.5
K = 2 * np.pi / (LAM * 1e-3)  # 1/mm
KEY = CoherenceKey("w")


def emitting(bundles: list[Any]) -> tuple[Graph, Any, dict[str, Any]]:
    """A graph with a source block "src" that emits the given spectra (S only)."""
    reg = builtin_registry()
    ev: dict[str, Any] = dict(optics.register(reg))
    ev.update(particle.register(reg))
    reg.add(NodeType("t.waves", Category.SOURCE, (light_out("out", "S"),)))

    def emit(ctx: Any, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
        return list(bundles)

    emit.representations = frozenset("S")  # type: ignore[attr-defined]
    ev["t.waves"] = emit
    g = Graph()
    g.add(Node("src", "t.waves"))
    return g, reg, ev


def power_in_na(sb: sp.SpectrumBundle, na: float) -> float:
    """Power (units of `field.field_power`) of the waves of a continuous spectrum within the NA
    about its axis, as `spectrum.spectrum_power` counts it."""
    w = sb.waves
    assert w.solid_angle is not None
    inside = np.hypot(w.kx, w.ky) <= na * sb.refractive_index * sb.k0
    part = PlaneWaves(
        w.port, w.kx[inside], w.ky[inside], w.field[inside], solid_angle=w.solid_angle[inside]
    )
    return float(spectrum_power(part)) / sb.k0


# ------------------------------------------------------------------------------ W84 (a)


def test_w84_incident_spectrum_of_an_s_link_against_richards_wolf() -> None:
    """W29 with S: what CliffordBem gets from an S link (no Debye per ray)."""
    from test_spectrum_links import lens_rays  # no ray of the rim lost: the axis is z

    rays = lens_rays(NA, F)
    sb = sp.rays_to_spectrum(rays)  # what `Executor._convert` makes on a link planned S
    notes: list[str] = []
    waves, rot = particle.incident_waves(sb, np.array([0.0, 0.0, F]), -0.002, 1.0, notes)
    assert np.allclose(rot[2], [0.0, 0.0, 1.0], atol=1e-9) and not notes
    power = particle.ray_powers(rays, np.asarray(rays.alive)).sum()
    assert spectrum_power(waves) == pytest.approx(K * power, rel=1e-3)  # measured 0 (rounding)
    rw = richards_wolf(Port(-0.002, 1, 1.0, LAM), NA, (0, 0, 0), (1.0, 0.0), 40, 80)
    pts = np.array([[0, 0, 0], [0.3, 0, 0], [0, 0.3, 0], [0, 0, 0.5], [0.2, -0.1, -0.3]]) * 1e-3
    ours = focal_field(waves, pts) / np.sqrt(spectrum_power(waves))
    theirs = focal_field(rw, pts) / np.sqrt(spectrum_power(rw))
    scale = theirs[0, 0] / ours[0, 0]
    assert abs(scale) == pytest.approx(1.0, abs=1e-3)  # measured 6.5e-5
    assert np.angle(scale) == pytest.approx(np.pi / 2, abs=1e-6)  # the −i of the Debye picture
    assert np.abs(ours * scale - theirs).max() < 2e-3 * abs(theirs[0, 0])  # measured 9.3e-4
    # the same focus as the Debye way of the ray path (W29)
    debye, _ = particle.rays_to_spectrum(rays, np.array([0.0, 0.0, F]), rot[2], -0.002)
    ref = focal_field(debye, pts) / np.sqrt(spectrum_power(debye))
    assert np.abs(ours - ref).max() < 1e-9 * abs(ref[0, 0])  # measured 5e-13


def test_w84_spectrum_in_another_medium_is_refused() -> None:
    plane = Plane.make([0, 0, -0.1], [0, 0, 1])
    e = np.ones((1, 2), dtype=complex)
    waves = PlaneWaves(Port(0.0, 1, 1.33 + 0j, LAM), np.zeros(1), np.zeros(1), e)
    with pytest.raises(ExecutionError, match="exec.particle"):
        particle.incident_waves(sp.SpectrumBundle(KEY, plane, waves), np.zeros(3), -0.002, 1.0, [])


def test_w84_spectrum_is_moved_onto_the_port_exactly() -> None:
    """A discrete wave tilted by 0.3 rad, 0.1 mm before the particle and off the axis: on the
    port its phase is that of the plane wave there, e^(ik·(P − o)) – the propagation is exact."""
    o = np.array([0.02, -0.01, -0.1])
    plane = Plane.make(o, [0, 0, 1])
    kpar = K * np.sin(0.3) * np.array([np.cos(0.4), np.sin(0.4)])
    e = np.array([[0.6 + 0.2j, -0.3j]])
    waves = PlaneWaves(Port(0.0, 1, 1.0 + 0j, LAM), kpar[:1].copy(), kpar[1:].copy(), e)
    sb = sp.SpectrumBundle(KEY, plane, waves)
    centre = np.array([0.0, 0.0, 0.0])
    got, rot = particle.incident_waves(sb, centre, -0.002, 1.0, [])
    kvec = np.array([kpar[0], kpar[1], np.sqrt(K * K - kpar @ kpar)]) @ rot  # global
    port_origin = centre - 0.002 * rot[2]
    expected = e[0] * np.exp(1j * kvec @ (port_origin - o))
    assert np.abs(got.field[0] - expected).max() < 1e-12
    assert got.port.z == -0.002 and got.solid_angle is None


# ------------------------------------------------------------------------------ W84 (b)


def _plane_wave_backscatter(level: int, store: ResultCache, spectrum: bool) -> float:
    """Power backscattered into the half space by the sphere for a unit plane wave along +z."""
    reg_params = {"point": [0, 0, 0], "radius_um": R, "index": N_P, "refinement": level}
    if spectrum:
        waves = PlaneWaves(
            Port(0.0, 1, 1.0 + 0j, LAM), np.zeros(1), np.zeros(1), np.array([[1.0 + 0j, 0.0]])
        )
        sb = sp.SpectrumBundle(KEY, Plane.make([0, 0, -0.1], [0, 0, 1]), waves)
        g, reg, ev = emitting([sb])
        g.add(Node("p", "fusion.particle", params=reg_params))
        g.connect(("src", "out"), ("p", "in"))
        out = Executor(g, reg, ev, available=WITH_SPECTRA, operators=store).evaluate("p", "out")
        assert isinstance(out[0], sp.SpectrumBundle)
        return sp.spectrum_power(out[0])
    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(particle.register(reg))
    g = Graph()
    src = {"radius": 0.01, "rings": 4, "center": [0, 0, -0.1], "wavelengths_um": [LAM]}
    g.add(Node("src", "optics.collimated_source", params=src))
    g.add(Node("p", "fusion.particle", params=reg_params))
    g.connect(("src", "out"), ("p", "in"))
    back: RayBundle = Executor(g, reg, ev, operators=store).evaluate("p", "out")[0]
    # the power of the rays with the far-field quadrature they come from (dσ = cos θ dΩ):
    # `ray_powers` would count it with the Delaunay areas of the directions (1 % off here)
    _, w = hemisphere(40, 80, False)
    cos = np.abs(np.asarray(back.direction)[:, 2])
    return float(np.sum(np.abs(np.asarray(back.amplitude)) ** 2 * w * cos) * back.irradiance)


def test_w84_plane_wave_spectrum_backscatter_against_mie() -> None:
    """W31 with S: the same Mie check, and the same number as the ray path."""
    cb = pytest.importorskip("cliffordbem")
    store = ResultCache()
    powers = {}
    for level in (6, 8):
        powers[len(np.asarray(cb.make_icosphere(level).triangles))] = _plane_wave_backscatter(
            level, store, True
        )
    (n1, p1), (n2, p2) = sorted(powers.items())
    extrapolated = (n2 * p2 - n1 * p1) / (n2 - n1)
    x = 2 * np.pi * R / LAM
    t = np.linspace(0, np.pi / 2, 8001)
    s1, s2 = S12(N_P, x, np.pi - t)
    mie = np.pi / K**2 * np.trapezoid((np.abs(s1) ** 2 + np.abs(s2) ** 2) * np.sin(t), t)
    assert p2 == pytest.approx(mie, rel=5e-2)
    assert abs(p2 / mie - 1) < abs(p1 / mie - 1)
    assert extrapolated == pytest.approx(mie, rel=1.5e-2)
    # the ray path of W31 sends the same plane wave (the field of the central ray) and gets the
    # same far field back
    rays = _plane_wave_backscatter(6, store, False)
    assert p1 == pytest.approx(rays, rel=1e-9)
    assert extrapolated == pytest.approx(mie, rel=5e-3)  # measured 3.3e-3 (rays: 7e-3, W31)


# ------------------------------------------------------------------------------ W85 (a)

DIPOLE = np.array([1.0, 0.3, 0.0])  # linear (rays carry linear polarisation)


def dipole_far(d: FloatArray, p: ComplexArray = DIPOLE) -> ComplexArray:
    """E∞ of a dipole p (Gaussian units, 1/(4πε) = 1): k²(r̂ × p) × r̂."""
    return np.asarray(K**2 * np.cross(np.cross(d, p), d), dtype=np.complex128)


def dipole_exact(r: FloatArray, p: ComplexArray = DIPOLE) -> ComplexArray:
    """The full field of the dipole: far term and near terms (1/r³ − ik/r²)(3r̂(r̂·p) − p)."""
    dist = np.linalg.norm(r, axis=1)
    rh = r / dist[:, None]
    near = (3 * rh * (rh @ p)[:, None] - p) * (1 / dist**3 - 1j * K / dist**2)[:, None]
    return np.asarray((dipole_far(rh, p) / dist[:, None] + near) * np.exp(1j * K * dist)[:, None])


def dipole_spectrum(
    centre: FloatArray, distance: float, p: ComplexArray = DIPOLE, n_theta: int = 40
) -> sp.SpectrumBundle:
    """The dipole's far field as the particle gives it: CliffordBem's front port (2 µm before the
    centre, direction −1, local frame = global axes) -> `particle.scattered_spectrum`."""
    dirs, w = hemisphere(n_theta, 2 * n_theta, False)
    waves = far_field_to_spectrum(
        Port(-0.002, -1, 1.0 + 0j, LAM), np.zeros(3), dirs, w, dipole_far(dirs, p)
    )
    incident = sp.SpectrumBundle(
        KEY, Plane.make(centre - [0, 0, 1], [0, 0, 1]), waves, ("src",), ("b", 3)
    )
    return particle.scattered_spectrum(
        waves, centre, np.eye(3), distance, incident, ("src", "p[back]")
    )


def field_on_plane(fb: Any) -> tuple[FloatArray, ComplexArray]:
    """Global sample points and global field of a gridless field bundle."""
    f = fb.field
    rot = np.stack([fb.plane.u, fb.plane.v, fb.plane.normal])
    pts = np.asarray(f.points)
    glob = fb.plane.origin + pts[:, :1] * fb.plane.u + pts[:, 1:] * fb.plane.v
    e = (np.asarray(f.amplitude) * np.exp(1j * np.asarray(f.phase))[:, None]) @ rot
    return np.asarray(glob), np.asarray(e)


@pytest.mark.parametrize("distance", [0.1, 1.0])
def test_w85_far_field_spectrum_is_the_spherical_wave(distance: float) -> None:
    p = np.array([1.0, 0.3j, 0.0])  # elliptical: the field keeps it
    centre = np.array([0.1, -0.2, 2.0])
    sb = dipole_spectrum(centre, distance, p)
    assert np.allclose(sb.plane.normal, [0, 0, -1]) and np.allclose(sb.plane.v, [0, -1, 0])
    assert sb.plane.origin @ [0, 0, 1] == pytest.approx(centre[2] - distance, abs=1e-15)
    assert sb.path == ("src", "p[back]") and sb.branch == ("b", 3)
    # power: n ∫|E∞|² dΩ over the half space = (4π/3) k⁴ |p|² (Gauss–Legendre is exact here)
    assert sp.spectrum_power(sb) == pytest.approx(4 * np.pi / 3 * K**4 * 1.09, rel=1e-12)
    fb = sp.spectrum_to_field(sb)
    glob, e = field_on_plane(fb)
    rel = glob - centre
    dist = np.linalg.norm(rel, axis=1)
    inside = np.hypot(rel[:, 0], rel[:, 1]) / dist <= NA
    spherical = dipole_far(rel / dist[:, None], p) * (np.exp(1j * K * dist) / dist)[:, None]
    scale = np.abs(spherical[inside]).max()
    assert np.abs(e[inside] - spherical[inside]).max() < 1e-9 * scale  # measured 3e-12
    error = np.abs(e[inside] - dipole_exact(rel[inside], p)).max() / scale
    assert error / (1 / (K * distance)) == pytest.approx(1.0, abs=2e-2)  # the near-field term
    # every sample lies where the ray of its direction from the centre crosses the plane
    d = rel / dist[:, None]
    w = sb.waves
    k_glob = np.column_stack([w.kx, -np.asarray(w.ky), -np.real(w.port.kz(w.kx, w.ky))]) / K
    assert np.abs(d - k_glob).max() < 1e-9


def test_w85_rays_from_the_far_field_spectrum_carry_its_power() -> None:
    """S -> H -> R of the particle's spectrum: rays start on the plane, run away from the
    centre and carry the power per direction (the Delaunay quadrature of the plane, 1e−2)."""
    centre = np.array([0.0, 0.0, 2.0])
    sb = dipole_spectrum(centre, 1.0)
    rays = field_to_rays(sp.spectrum_to_field(sb))
    alive = np.asarray(rays.alive)
    d = np.asarray(rays.direction)[alive]
    pos = np.asarray(rays.position)[alive]
    back_to_centre = pos - d * np.sum((pos - centre) * d, axis=1)[:, None]  # closest approach
    assert np.abs(back_to_centre - centre).max() < 1e-9
    power = particle.ray_powers(rays, alive)
    inside = np.hypot(d[:, 0], d[:, 1]) <= NA
    assert power[inside].sum() == pytest.approx(power_in_na(sb, NA), rel=1e-2)  # measured 4e-3


def test_w85_particle_spectrum_through_port_data(tmp_path: Path) -> None:
    sb = dipole_spectrum(np.array([0.0, 0.0, 2.0]), 1.0)
    back = portdata.read(portdata.write(sb, tmp_path), tmp_path)
    assert isinstance(back, sp.SpectrumBundle)
    for name in ("points", "phase", "gradient", "hessian", "amplitude"):
        assert np.array_equal(getattr(back.support, name), getattr(sb.support, name))
    a, b = field_on_plane(sp.spectrum_to_field(sb)), field_on_plane(sp.spectrum_to_field(back))
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


# ------------------------------------------------------------------------------ W85 (b)


@pytest.fixture(scope="module")
def k3() -> dict[str, Any]:
    """Pilot K3 in the graph with S on the link into the particle and "auto" behind it (into a
    detector in the far field), the ray path of W30, and CliffordBem alone."""
    cb = pytest.importorskip("cliffordbem")
    from rtt_fusion.adapters import cliffordbem as cba

    level = 4
    store = ResultCache()
    out: dict[str, Any] = {}
    det = {"point": [0, 0, F - 1.0], "normal": [0, 0, -1], "half_size": [0.6, 0.6]}
    for transport in ("S", "R"):
        g, ex = focused()
        # no ray of the rim lost by rounding: the axis is z, the spectrum's plane parallel to the
        # detector (a tilted one is resampled through rays, which keep only the linear part of
        # the elliptically polarised backscatter off the axes – an open point of H -> R)
        g.nodes["lens"].params["aperture_radius"] = NA * F * (1 + 1e-6)
        g.add(
            Node(
                "p",
                "fusion.particle",
                params={"point": [0, 0, F], "radius_um": R, "index": N_P, "refinement": level},
            )
        )
        g.add(Node("cam", "optics.detector", params={**det, "pixels": [9, 9]}))
        g.connect(("lens", "out"), ("p", "in"), transport)
        g.connect(("p", "out"), ("cam", "in"))
        ex = Executor(g, ex.registry, ex.evaluators, available=WITH_SPECTRA, operators=store)
        rays = ex.evaluate("lens", "out")[0]
        out[transport] = {
            "incident": float(particle.ray_powers(rays, np.asarray(rays.alive)).sum()),
            "back": ex.evaluate("p", "out"),
            "executor": ex,
        }
        if transport == "S":
            out["image"] = ex.evaluate("cam", "image")
    m0 = cb.make_icosphere(level)
    mesh = cb.TriangleMesh(np.asarray(m0.points) * R, m0.triangles)
    k_um = 2 * np.pi / LAM
    handle = cba.prepare(cba.Structure(mesh, N_P**2, 1.0, (0, 0, 0.002), 0.0, 0.004), [LAM])
    beam = cb.BeamField.focused(
        cb.Medium(eps=1.0), k_um, np.zeros(3), NA, 1e6, np.array([1, 0, 0], complex), 40, 80
    )
    b = np.asarray(beam.project(mesh, cb.Medium(eps=1.0)))
    r = cba._problem(handle, LAM).solve_rhs(b, cb.SolveOptions(tol=1e-8))
    h_scat = np.asarray(r.h) - b

    def far(dirs: FloatArray) -> ComplexArray:
        e = np.asarray(cb.far_field_E(mesh, h_scat, k_um, 1.0, dirs))
        return np.asarray(e - np.sum(e * dirs, axis=1, keepdims=True) * dirs)

    dirs, w = hemisphere(40, 80, False)
    e = far(dirs)
    inside = np.sqrt(1 - dirs[:, 2] ** 2) <= NA
    out["direct"] = np.sum(w[inside] * np.sum(np.abs(e[inside]) ** 2, axis=1)) / (2 * beam.power)
    out["far"], out["beam_power"] = far, float(beam.power)
    return out


def test_w85_k3_in_the_graph_without_rays(k3: dict[str, Any]) -> None:
    s = k3["S"]
    back = s["back"]
    assert len(back) == 1 and isinstance(back[0], sp.SpectrumBundle)
    sb = back[0]
    log = "\n".join(s["executor"].context.log)
    assert "incident waves (spectrum)" in log and "S -> R" not in log
    plan = s["executor"].plan
    assert plan.representations["lens.out -> p.in"] == "S"
    assert plan.representations["p.out -> cam.in"] == "S"  # "auto" behind a spectrum (R30f)
    assert [c.node for c in plan.conversions] == ["lens"]  # R -> S at the lens, no S -> R
    assert sb.plane.normal @ [0, 0, 1] == pytest.approx(-1.0, abs=1e-12)  # the axis of the lens
    # the scattered spectrum's plane: 100 wavelengths in front of the particle, facing back
    in_front = (np.array([0.0, 0.0, F]) - sb.plane.origin) @ -np.asarray(sb.plane.normal)
    assert in_front == pytest.approx(particle.OUT_WAVELENGTHS * LAM * 1e-3, rel=1e-12)
    fraction = power_in_na(sb, NA) / s["incident"]
    assert fraction == pytest.approx(k3["direct"], rel=5e-3)  # as W30; measured 4.3e-4
    assert 0.03 < fraction < 0.05
    rays = k3["R"]["back"][0]
    alive = np.asarray(rays.alive)
    d = np.asarray(rays.direction)[alive]
    p = particle.ray_powers(rays, alive)
    ray_fraction = p[np.hypot(d[:, 0], d[:, 1]) <= NA].sum() / k3["R"]["incident"]
    # the ray path counts with the Delaunay areas of the directions (W30: 2.2e-3 off)
    assert fraction == pytest.approx(ray_fraction, rel=5e-3)  # measured 1.8e-3


def test_w85_k3_detector_in_the_far_field(k3: dict[str, Any]) -> None:
    """The image of the backscattered spectrum 1 mm in front of the particle against the far
    field of CliffordBem alone, |E∞|²/r² per incident power."""
    image = k3["image"]
    log = "\n".join(k3["S"]["executor"].context.log)
    assert "S -> H" in log  # far from the particle the sum would alias: through H, exactly
    values = np.asarray(image.values) / k3["S"]["incident"]
    ny, nx = values.shape
    u = (np.arange(nx) + 0.5) / nx * 1.2 - 0.6
    v = (np.arange(ny) + 0.5) / ny * 1.2 - 0.6
    plane = Plane.make([0, 0, F - 1.0], [0, 0, -1])
    uu, vv = np.meshgrid(u, v)
    pts = plane.origin + uu.reshape(-1, 1) * plane.u + vv.reshape(-1, 1) * plane.v
    rel = (pts - [0, 0, F]) * 1e3  # µm, the frame of CliffordBem (centre at the origin)
    dist = np.linalg.norm(rel, axis=1)
    e = k3["far"](rel / dist[:, None])  # E∞ in CliffordBem's units (r in µm)
    per_um2 = np.sum(np.abs(e) ** 2, axis=1) / dist**2 / (2 * k3["beam_power"])
    ref = (per_um2 * 1e6).reshape(ny, nx)  # per mm², as the graph's powers
    error = np.abs(values - ref).max() / ref.max()
    assert error < 1e-2  # measured 5.5e-3 (interpolation of the field between its samples)


# ------------------------------------------------------------------------------ planner


def test_auto_plans_around_the_particle_are_unchanged() -> None:
    """Sockets RHS, evaluator R and S: "auto" stays R on every link (S would be admissible only
    between native blocks, and costs more), with no conversions – as before R30d."""
    g, ex = focused()
    g.add(Node("p", "fusion.particle", params={"point": [0, 0, F]}))
    g.add(Node("cam", "optics.detector", params={"point": [0, 0, 1.0], "normal": [0, 0, -1]}))
    g.connect(("lens", "out"), ("p", "in"))
    g.connect(("p", "out"), ("cam", "in"))
    ex = Executor(g, ex.registry, ex.evaluators, available=WITH_SPECTRA)
    ex.check()
    p = ex.plan
    assert p is not None and p.ok and set(p.representations.values()) == {"R"}
    assert not p.conversions
    assert p.reasons[link_key(g.links[1])].endswith("among RH")  # the lens does not take S
    assert p.reasons["p.out -> cam.in"].endswith("among RHS")  # S admitted (native), R cheaper
    assert ex.native("p") == frozenset("RS")


# ------------------------------------------------------------------------------ R30c open point


def test_spectrum_without_sampling_gives_a_clear_error() -> None:
    """A continuous spectrum without sampling (as behind a condensed group) into a block that
    takes rays: `exec.spectrum_support` naming the block; a detector still takes it."""
    w = PlaneWaves(
        Port(0.0, 1, 1.0 + 0j, LAM),
        np.array([0.0, 100.0, 0.0]),
        np.array([0.0, 0.0, 100.0]),
        np.ones((3, 2), dtype=complex),
        solid_angle=np.full(3, 1e-6),
    )
    sb = sp.SpectrumBundle(KEY, Plane.make([0, 0, 0], [0, 0, 1]), w, ("src", "cond"))
    g, reg, ev = emitting([sb])
    g.add(Node("lens", "optics.thin_lens", params={"point": [0, 0, 1], "focal_length": 5.0}))
    g.connect(("src", "out"), ("lens", "in"), "S")
    with pytest.raises(ExecutionError, match="from cond carries no sampling") as e:
        Executor(g, reg, ev, available=WITH_SPECTRA).evaluate("lens", "out")
    assert e.value.code == "exec.spectrum_support"
    g2, reg2, ev2 = emitting([sb])
    g2.add(Node("cam", "optics.detector", params={"point": [0, 0, 1], "pixels": [4, 4]}))
    g2.connect(("src", "out"), ("cam", "in"))
    image = Executor(g2, reg2, ev2, available=WITH_SPECTRA).evaluate("cam", "image")
    assert float(np.asarray(image.values).max()) > 0
    # discrete waves keep their own message (no solid angles)
    with pytest.raises(FieldError, match="discrete plane waves"):
        sp.spectrum_to_field(replace(sb, waves=replace(w, solid_angle=None)))
