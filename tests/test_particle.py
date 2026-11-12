"""Particle as the sample (R19c): CliffordBem in the working graph.

References, independent of the code under test:
- W29 converging rays of an aplanatic lens (NA 0.8) as an angular spectrum: the focal field equals
  the Richards–Wolf spectrum of `couplings.focus` (Gauss–Legendre quadrature, no rays) up to one
  complex constant with modulus 1 after power normalisation; `spectrum_power` = k·Σ P_rays and
  Σ P_rays = the source power (energy).
- W30 pilot K3 in the graph: source → aplanatic lens → particle (sphere R = 0.25 µm, n = 1.5); the
  backscattered power within the NA divided by the incident power equals the same quantity
  computed with CliffordBem alone (its own Richards–Wolf beam `BeamField.focused`, far field,
  no adapter, no rays) on the same mesh.
- W31 plane wave on the sphere: the power scattered into the backward hemisphere against Mie
  (tests/mie_reference.py); flat triangles make the BEM converge like 1/N_triangles, so two
  meshes are extrapolated (Richardson) before the comparison.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.couplings.focus import richards_wolf
from rtt_fusion.graph import Graph, Node, builtin_registry, optics, particle
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.light import RayBundle
from rtt_fusion.graph.particle import ray_powers, rays_to_spectrum
from rtt_fusion.ports import PlaneWaves, Port, spectrum_power

sys.path.insert(0, str(Path(__file__).parent))
from mie_reference import S12

LAM, NA, F, R, N_P = 0.6, 0.8, 2.0, 0.25, 1.5
K = 2 * np.pi / (LAM * 1e-3)  # 1/mm


def focused(rings: int = 30) -> tuple[Graph, Executor]:
    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(particle.register(reg))
    g = Graph()
    a = NA * F  # aplanatic: h = f sin θ
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": a, "rings": rings, "center": [0, 0, -1], "wavelengths_um": [LAM]},
        )
    )
    g.add(Node("lens", "optics.aplanatic_lens", params={"focal_length": F, "aperture_radius": a}))
    g.connect(("src", "out"), ("lens", "in"))
    return g, Executor(g, reg, ev)


def focal_field(w: PlaneWaves, points: FloatArray) -> ComplexArray:
    kz = np.real(w.port.kz(w.kx, w.ky))
    ez = w.ez()
    out = []
    for p in points:
        ph = np.exp(1j * (w.kx * p[0] + w.ky * p[1] + kz * (p[2] - w.port.z)))
        out.append([np.sum(w.field[:, 0] * ph), np.sum(w.field[:, 1] * ph), np.sum(ez * ph)])
    return np.asarray(out, dtype=np.complex128)


def test_w29_rays_to_spectrum_equals_richards_wolf() -> None:
    _, ex = focused()
    rays: RayBundle = ex.evaluate("lens", "out")[0]
    waves, _ = rays_to_spectrum(rays, np.array([0.0, 0.0, F]), np.array([0.0, 0.0, 1.0]), -0.002)
    power = ray_powers(rays, np.asarray(rays.alive)).sum()
    assert power == pytest.approx(np.pi * (NA * F) ** 2, rel=1e-3)  # source power, unit irradiance
    assert spectrum_power(waves) == pytest.approx(K * power, rel=1e-12)
    rw = richards_wolf(Port(-0.002, 1, 1.0, LAM), NA, (0, 0, 0), (1.0, 0.0), 40, 80)
    pts = np.array([[0, 0, 0], [0.3, 0, 0], [0, 0.3, 0], [0, 0, 0.5], [0.2, -0.1, -0.3]]) * 1e-3
    ours = focal_field(waves, pts) / np.sqrt(spectrum_power(waves))
    theirs = focal_field(rw, pts) / np.sqrt(spectrum_power(rw))
    scale = theirs[0, 0] / ours[0, 0]
    assert abs(scale) == pytest.approx(1.0, abs=1e-3)
    # the Debye picture has −i between the field on the sphere and the focus; couplings.focus
    # leaves that constant out (its docstring), so the scale is +i
    assert np.angle(scale) == pytest.approx(np.pi / 2, abs=1e-6)
    assert np.abs(ours * scale - theirs).max() < 2e-3 * abs(theirs[0, 0])


def test_w30_pilot_k3_in_the_graph() -> None:
    cb = pytest.importorskip("cliffordbem")
    from rtt_fusion.adapters import cliffordbem as cba
    from rtt_fusion.convert import hemisphere

    level = 4
    g, ex = focused()
    g.add(
        Node(
            "p",
            "fusion.particle",
            params={"point": [0, 0, F], "radius_um": R, "index": N_P, "refinement": level},
        )
    )
    g.connect(("lens", "out"), ("p", "in"))
    incident = ray_powers(
        ex.evaluate("lens", "out")[0], np.asarray(ex.evaluate("lens", "out")[0].alive)
    ).sum()
    back: RayBundle = ex.evaluate("p", "out")[0]
    alive = np.asarray(back.alive)
    p = ray_powers(back, alive)
    d = np.asarray(back.direction)[alive]
    assert np.all(d[:, 2] < 1e-3)  # back towards the lens (the last quadrature node is grazing)
    graph_fraction = p[np.hypot(d[:, 0], d[:, 1]) <= NA].sum() / incident
    # CliffordBem alone, as in pilot K3, on the same mesh
    m0 = cb.make_icosphere(level)
    mesh = cb.TriangleMesh(np.asarray(m0.points) * R, m0.triangles)
    k_um = 2 * np.pi / LAM
    handle = cba.prepare(cba.Structure(mesh, N_P**2, 1.0, (0, 0, 0.002), 0.0, 0.004), [LAM])
    beam = cb.BeamField.focused(
        cb.Medium(eps=1.0), k_um, np.zeros(3), NA, 1e6, np.array([1, 0, 0], complex), 40, 80
    )
    b = np.asarray(beam.project(mesh, cb.Medium(eps=1.0)))
    r = cba._problem(handle, LAM).solve_rhs(b, cb.SolveOptions(tol=1e-8))
    dirs, w = hemisphere(40, 80, False)
    e = np.asarray(cb.far_field_E(mesh, np.asarray(r.h) - b, k_um, 1.0, dirs))
    e = e - np.sum(e * dirs, axis=1, keepdims=True) * dirs
    inside = np.sqrt(1 - dirs[:, 2] ** 2) <= NA
    direct = np.sum(w[inside] * np.sum(np.abs(e[inside]) ** 2, axis=1)) / (2 * beam.power)
    assert graph_fraction == pytest.approx(direct, rel=5e-3)
    assert 0.03 < graph_fraction < 0.05  # about 4 % return into the objective


def test_w31_plane_wave_backscatter_against_mie() -> None:
    cb = pytest.importorskip("cliffordbem")
    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(particle.register(reg))
    powers = {}
    for level in (6, 8):
        g = Graph()
        g.add(
            Node(
                "src",
                "optics.collimated_source",
                params={
                    "radius": 0.01,
                    "rings": 4,
                    "center": [0, 0, -0.1],
                    "wavelengths_um": [LAM],
                },
            )
        )
        g.add(
            Node(
                "p",
                "fusion.particle",
                params={"point": [0, 0, 0], "radius_um": R, "index": N_P, "refinement": level},
            )
        )
        g.connect(("src", "out"), ("p", "in"))
        back: RayBundle = Executor(g, reg, ev).evaluate("p", "out")[0]
        n = len(np.asarray(cb.make_icosphere(level).triangles))
        powers[n] = ray_powers(back, np.asarray(back.alive)).sum()
    (n1, p1), (n2, p2) = sorted(powers.items())
    extrapolated = (n2 * p2 - n1 * p1) / (n2 - n1)  # flat-triangle error ∝ 1/N
    x = 2 * np.pi * R / LAM
    t = np.linspace(0, np.pi / 2, 8001)
    s1, s2 = S12(N_P, x, np.pi - t)
    mie = np.pi / K**2 * np.trapezoid((np.abs(s1) ** 2 + np.abs(s2) ** 2) * np.sin(t), t)
    assert p2 == pytest.approx(mie, rel=5e-2)  # finer mesh within the flat-triangle error
    assert abs(p2 / mie - 1) < abs(p1 / mie - 1)  # and converging
    assert extrapolated == pytest.approx(mie, rel=1.5e-2)
