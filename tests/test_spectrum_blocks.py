"""Gratings, layer stacks and condensed groups with plane-wave spectra S (R30c,
docs/plan-s-spektren.md).

References, independent of the code under test:
- W81 two gratings in sequence. (a) A plane wave (discrete S) through two lamellar gratings 2 µm
  apart in glass, the chain of the two blocks against FourierGrating with both gratings and the
  spacer in one structure: equal up to the multiple reflections between the gratings that a
  chain of blocks leaves out – they scale with the square of the index contrast, which is checked
  at two contrasts. (b) A Gaussian beam (continuous S, a grid spectrum) against the H chain of
  the same blocks (the field operator of F1, a different code path: a grid per block and
  interpolation back to the samples): with the second grating in the plane of the first equal to
  the sampling of the H path; with a gap of z the H chain's geometric step between the gratings
  leaves out the diffraction of the beam, a phase of order z/z_R, which S keeps (exact
  propagation in the spectrum).
- W82 energy: a lossless grating sends the incident power into its orders, Σ R_n + Σ T_n = 1, for
  every wave of a continuous spectrum (per wave k_z|E|²/Δ²k with Δ²k kept); a lossless stack
  R + T = 1.
- W83 the stack against the Airy formula of one layer (Fresnel coefficients, p̂ = ŝ × k̂) per
  wave: complex r_s, t_s, r_p, t_p for discrete waves of any azimuth and the waves of a
  continuous spectrum, from both sides, and dΩ behind the stack (Δ²k kept); the condensed plate
  (R29, three inner branches) against the first three terms of the Airy sum of a plate.
- Port data: spectra (PFT, grid, gridded, discrete) bit for bit through files, S -> H after the
  round trip equal; the grating's S path across a process boundary equal to the local one.
- Planner and diagnostics: existing "auto" plans (Linnik with a grating sample) unchanged; S on
  the links between blocks that take it natively, chosen by "auto" behind a spectrum (R30f);
  S -> R of orders that overlap reported.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion import portdata
from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.graph import (
    Category,
    Graph,
    Node,
    NodeType,
    Requirements,
    builtin_registry,
    grating,
    optics,
    plan,
    stack,
)
from rtt_fusion.graph.execute import WITH_SPECTRA, Context, Executor
from rtt_fusion.graph.light import CoherenceKey, Plane
from rtt_fusion.graph.model import light_out
from rtt_fusion.graph.planner import link_key
from rtt_fusion.ports import PlaneWaves, Port

fg = pytest.importorskip("fouriergrating")
pytest.importorskip("fieldtracing")

from rtt_fusion.adapters import fouriergrating as fga  # noqa: E402
from rtt_fusion.graph import spectrum as sp  # noqa: E402
from rtt_fusion.graph.field import FieldBundle, rays_to_field  # noqa: E402

DATA = Path(__file__).parent / "data"
LAM = 0.55
K0 = 2e3 * np.pi / LAM  # 1/mm
KEY = CoherenceKey("s")
DOWN = Plane.make([0.0, 0.0, 0.0], [0.0, 0.0, -1.0])  # light along −z onto blocks at z = 0


def glob(sb: sp.SpectrumBundle) -> tuple[FloatArray, ComplexArray]:
    """Global wave vectors (N, 3) and fields (N, 3) of a spectrum."""
    w = sb.waves
    rot = np.stack([sb.plane.u, sb.plane.v, sb.plane.normal])
    kz = np.real(w.port.kz(w.kx, w.ky))
    k = np.column_stack([w.kx, w.ky, kz]) @ rot
    e = np.column_stack([np.asarray(w.field), w.ez()]) @ rot
    return np.asarray(k), np.asarray(e)


def wave_set(spectra: list[sp.SpectrumBundle], digits: int = 6) -> dict[Any, ComplexArray]:
    """The waves of several spectra by global k (rounded), fields of equal k added."""
    out: dict[Any, ComplexArray] = {}
    for sb in spectra:
        k, e = glob(sb)
        for j in range(len(k)):
            key: tuple[float, ...] = tuple(float(x) for x in np.round(k[j], digits))
            out[key] = np.asarray(out.get(key, 0) + e[j])
    return out


def discrete(plane: Plane, n: float, kpar: FloatArray, e: ComplexArray) -> sp.SpectrumBundle:
    waves = PlaneWaves(Port(0.0, 1, complex(n), LAM), kpar[:, 0].copy(), kpar[:, 1].copy(), e)
    return sp.SpectrumBundle(KEY, plane, waves)


def gaussian_spectrum(
    plane: Plane, n: float, centre: FloatArray, width: float, count: int = 11
) -> sp.SpectrumBundle:
    """A continuous spectrum on a regular k∥ grid (count² waves within ±3 widths of `centre`):
    Gaussian amplitude, a polarisation that turns across the spectrum, dΩ = Δ²k/(k k_z)."""
    k = n * K0
    a = np.linspace(-3.0, 3.0, count) * width
    ku, kv = np.meshgrid(centre[0] + a, centre[1] + a)
    ku, kv = ku.ravel(), kv.ravel()
    amp = np.exp(-((ku - centre[0]) ** 2 + (kv - centre[1]) ** 2) / width**2)
    angle = 0.4 + 0.3 * (ku - centre[0]) / width
    e = amp[:, None] * np.column_stack([np.cos(angle), np.sin(angle) * np.exp(0.7j)])
    dk2 = (a[1] - a[0]) ** 2
    kz = np.sqrt(k * k - ku * ku - kv * kv)
    waves = PlaneWaves(
        Port(0.0, 1, complex(n), LAM), ku, kv, e.astype(complex), solid_angle=dk2 / (k * kz)
    )
    return sp.SpectrumBundle(KEY, plane, waves)


def lamellar(t: float, edge: float, materials: tuple[str, str]) -> dict[str, Any]:
    return {
        "type": "lamellar1d",
        "thickness_um": t,
        "breakpoints_um": [0.0, edge, 1.0],
        "materials": list(materials),
    }


# ------------------------------------------------------------------------------ W81 (a)


def _two_gratings(contrast: float) -> float:
    """Largest deviation of the chain from the joint structure, relative to the largest order."""
    n, gap = 1.5, 2.0  # glass, µm
    base = {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 12},
        "period_um": {"x": 1.0},
        "materials": {"g": {"n": n}, "h": {"n": n + contrast}},
        "superstrate": "g",
        "substrate": "g",
    }
    g1, g2 = lamellar(1.5, 0.4, ("h", "g")), lamellar(1.0, 0.3, ("h", "g"))
    # FourierGrating lists the layers from the substrate upwards: g2, spacer, g1 (light from top)
    joint = dict(
        base, layers=[g2, {"type": "homogeneous", "thickness_um": gap, "material": "g"}, g1]
    )
    a = Node("g1", "fg.grating", params={"config": dict(base, layers=[g1]), "side": "transmitted"})
    # the blocks are thin: the transmitted orders of g1 refer to its bottom, so g2 is `gap` below
    b = Node(
        "g2",
        "fg.grating",
        params={
            "config": dict(base, layers=[g2]),
            "side": "transmitted",
            "point": [0, 0, -gap * 1e-3],
        },
    )
    theta, phi = np.radians(3.0), np.radians(25.0)
    k_in = (
        n
        * K0
        * np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), -np.cos(theta)])
    )
    kpar = np.array([[k_in @ DOWN.u, k_in @ DOWN.v]])
    e_in = np.array([[0.6 + 0.2j, -0.3 + 0.5j]])
    ctx = Context(Graph())
    first = grating.grating(ctx, a, "out", {"in": [discrete(DOWN, n, kpar, e_in)]})
    chain = wave_set(grating.grating(ctx, b, "out", {"in": first}))
    assert all(not sb.continuous and sb.support is None for sb in first)
    # the joint structure, the incident wave in its frame (x = grating vector, z along the light)
    x, z = np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, -1.0])
    frame = np.stack([x, np.cross(z, x), z])
    handle = fga.prepare(fga.Structure(joint), (LAM,))
    rot = np.stack([DOWN.u, DOWN.v, DOWN.normal])
    kl = np.array([kpar[0, 0], kpar[0, 1], -k_in[2]])
    e_glob = np.array([*e_in[0], -(kpar[0] @ e_in[0]) / kl[2]]) @ rot
    kf, ef = frame @ k_in, frame @ e_glob
    inc = PlaneWaves(fga.front_port(handle, LAM), kf[:1], kf[1:2], ef[None, :2])
    ref = fga.apply_many(handle, [inc])[0]["back"]
    kz = np.real(ref.port.kz(ref.kx, ref.ky))
    k_ref = np.column_stack([ref.kx, ref.ky, kz]) @ frame
    e_ref = np.column_stack([np.asarray(ref.field), ref.ez()]) @ frame
    assert len(chain) == len(k_ref) == 5  # orders −2 … 2 propagate in the glass
    scale = float(np.abs(e_ref).max())
    worst = 0.0
    for kk, ee in zip(k_ref, e_ref, strict=True):
        got = chain[tuple(np.round(kk, 6))]
        worst = max(worst, float(np.abs(got - ee).max()) / scale)
    return worst


def test_w81_two_gratings_as_a_chain_equal_the_joint_structure() -> None:
    """W81 (a): measured 5.9e−5 at Δn = 0.03 and 5.7e−6 at Δn = 0.01 – the multiple reflections
    between the gratings (∝ Δn²), which a chain of blocks does not have."""
    strong, weak = _two_gratings(0.03), _two_gratings(0.01)
    assert strong < 1e-4 and weak < 1.5e-5
    assert strong / weak > 5.0  # second order in the contrast: the chain itself adds nothing


# ------------------------------------------------------------------------------ W81 (b)


def _phase_grating(edge: float, nx: int = 6) -> dict[str, Any]:
    return {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": nx},
        "period_um": {"x": 10.0},
        "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}},
        "superstrate": "air",
        "substrate": "air",
        "layers": [
            {
                "type": "lamellar1d",
                "thickness_um": 0.4,
                "breakpoints_um": [0, edge, 10.0],
                "materials": ["glass", "air"],
            }
        ],
    }


def _beam() -> FieldBundle:
    from test_grating_field import _graph_beam

    return rays_to_field(_graph_beam(5.0, radius=0.1, rings=10, waist=0.035, lam=LAM))


def _field_at(fb: FieldBundle) -> tuple[ComplexArray, FloatArray, FloatArray]:
    """Global field and global points of the samples of a field, and their in-plane points."""
    f: Any = fb.field
    rot = np.stack([fb.plane.u, fb.plane.v, fb.plane.normal])
    p = np.asarray(f.points)
    e = (np.asarray(f.amplitude) * np.exp(1j * np.asarray(f.phase))[:, None]) @ rot
    return e, p @ rot[:2] + fb.plane.origin, p


def test_w81_beam_through_two_gratings_equals_the_h_chain() -> None:
    """W81 (b): a Gaussian beam (w₀ = 35 µm, 5° incidence) through two phase gratings (Λ = 10
    µm, orders −1, 0, 1 each), S chain against H chain, per pair of orders on the H field's
    samples (the S field there by the sum of its waves). Second grating in the plane of the first:
    measured 5.8e−4 (the H path's sampling: 1.8e−3 with 6 rings instead of 10). 20 µm apart:
    the H chain's geometric step leaves out the beam's diffraction, measured 2.7e−3 against
    z/z_R = 2.9e−3."""
    fb = _beam()
    params = {"side": "transmitted", "orders": [-1, 0, 1], "field_transfer": "interpolated"}
    g1 = Node("g1", "fg.grating", params={"config": _phase_grating(5.0), **params})
    ctx_h, ctx_s = Context(Graph()), Context(Graph())
    h1 = grating.grating(ctx_h, g1, "out", {"in": [fb]})
    s1 = grating.grating(ctx_s, g1, "out", {"in": [sp.field_to_spectrum(fb)]})
    assert all(isinstance(b, sp.SpectrumBundle) and b.continuous for b in s1)
    assert any("verified to" in line for line in ctx_s.log)  # the interpolated transfer
    z_r = K0 * 0.035**2 / 2
    for gap, limit in ((0.0, 1e-3), (0.02, None)):
        g2 = Node(
            "g2",
            "fg.grating",
            params={"config": _phase_grating(3.0), **params, "point": [0, 0, -gap]},
        )
        h2 = grating.grating(ctx_h, g2, "out", {"in": h1})
        s2 = grating.grating(ctx_s, g2, "out", {"in": s1})
        assert [b.path for b in h2] == [b.path for b in s2] and len(s2) == 9
        worst = 0.0
        for a, b in zip(h2, s2, strict=True):
            assert isinstance(a, FieldBundle) and isinstance(b, sp.SpectrumBundle)
            assert b.support is not None  # S -> H works behind the second grating
            e_h, points, local = _field_at(a)
            core = np.hypot(*(local - local.mean(axis=0)).T) < 0.07  # inside the beam
            e_s = sp.wave_sum(b, points)
            worst = max(worst, float(np.abs(e_h - e_s)[core].max() / np.abs(e_h).max()))
        if limit is not None:
            assert worst < limit, (gap, worst)
        else:
            assert 0.5 < worst / (gap / z_r) < 1.5, (worst, gap / z_r)


# ------------------------------------------------------------------------------ W82


LAMELLAR = {
    "dimension": "1D",
    "solver": "rcwa",
    "truncation": {"Nx": 8},
    "period_um": {"x": 0.7},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}},
    "superstrate": "air",
    "substrate": "glass",
    "layers": [
        {
            "type": "lamellar1d",
            "thickness_um": 0.25,
            "breakpoints_um": [0, 0.35, 0.7],
            "materials": ["glass", "air"],
        }
    ],
}


def _incident(n: float = 1.0) -> sp.SpectrumBundle:
    """A beam of about 20 µm waist at 10° in the x–z plane (and 4° out of it), on a plane 1 µm
    above the blocks at z = 0."""
    plane = Plane.make([0.0, 0.0, 1e-3], [0.0, 0.0, -1.0])
    k = n * K0 * np.array([np.sin(np.radians(10.0)), np.sin(np.radians(4.0)), 0.0])
    return gaussian_spectrum(plane, n, np.array([k @ plane.u, k @ plane.v]), 100.0)


def test_w82_orders_of_both_sides_carry_the_incident_power() -> None:
    """W82: Σ over the reflected and transmitted orders = incident power, measured 4e−15 (the
    spectrum's 121 waves, FourierGrating's energy balance per wave)."""
    sb = _incident()
    ctx = Context(Graph())
    total = 0.0
    orders = 0
    for side in ("reflected", "transmitted"):
        node = Node("gr", "fg.grating", params={"config": LAMELLAR, "side": side})
        out = grating.grating(ctx, node, "out", {"in": [sb]})
        assert all(b.continuous and b.waves.orders is not None for b in out)
        total += sum(sp.spectrum_power(b) for b in out)
        orders += len(out)
    assert orders >= 4  # reflected 0, −1; transmitted 0, −1 (and more)
    assert total / sp.spectrum_power(sb) == pytest.approx(1.0, abs=1e-10)


BRAGG = [(2.4, LAM / (4 * 2.4)), (1.46, LAM / (4 * 1.46))] * 2


def _stack_config(n_in: float, layers: list[tuple[float, float]], n_out: float) -> dict[str, Any]:
    """`layers` from the front (incident side) on; FourierGrating lists them from the back."""
    materials: dict[str, Any] = {"m_in": {"n": n_in}, "m_out": {"n": n_out}}
    cfg_layers = []
    for k, (n, d) in reversed(list(enumerate(layers))):
        materials[f"l{k}"] = {"n": n}
        cfg_layers.append({"type": "homogeneous", "thickness_um": d, "material": f"l{k}"})
    return {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 0},
        "period_um": {"x": 1.0},
        "materials": materials,
        "superstrate": "m_in",
        "substrate": "m_out",
        "layers": cfg_layers,
    }


def test_w82_stack_reflects_and_transmits_the_incident_power() -> None:
    """W82 for a lossless Bragg stack: R + T = 1 over a continuous spectrum, measured 7e−16."""
    node = Node("st", "fg.stack", params={"config": _stack_config(1.0, BRAGG, 1.5)})
    sb = _incident()
    ctx = Context(Graph())
    r = stack.stack(ctx, node, "front_out", {"front_in": [sb]})[0]
    t = stack.stack(ctx, node, "back_out", {"front_in": [sb]})[0]
    ratio = (sp.spectrum_power(r) + sp.spectrum_power(t)) / sp.spectrum_power(sb)
    assert ratio == pytest.approx(1.0, abs=1e-12)
    assert 0.5 < sp.spectrum_power(r) / sp.spectrum_power(sb) < 1.0  # a mirror


# ------------------------------------------------------------------------------ W83


def airy(
    ns: tuple[float, float, float], d: float, sin0: float, pol: str
) -> tuple[complex, complex]:
    """r and t of one layer (index ns[1], thickness d µm) between ns[0] and ns[2] at n₀ sin θ₀ =
    ns[0]·sin0: Fresnel coefficients (s; p with p̂ = ŝ × k̂ for every wave, r_p = (n₂c₁ − n₁c₂)/
    (n₂c₁ + n₁c₂)), Airy sum; the transmitted phase refers to the exit face."""
    s0 = ns[0] * sin0
    c = [np.sqrt(1 - (s0 / n) ** 2 + 0j) for n in ns]

    def rt(i: int, j: int) -> tuple[complex, complex]:
        ni, nj, ci, cj = ns[i], ns[j], c[i], c[j]
        if pol == "s":
            return (ni * ci - nj * cj) / (ni * ci + nj * cj), 2 * ni * ci / (ni * ci + nj * cj)
        return (nj * ci - ni * cj) / (nj * ci + ni * cj), 2 * ni * ci / (nj * ci + ni * cj)

    r01, t01 = rt(0, 1)
    r12, t12 = rt(1, 2)
    beta = 2 * np.pi / LAM * ns[1] * d * c[1]
    den = 1 + r01 * r12 * np.exp(2j * beta)
    return complex((r01 + r12 * np.exp(2j * beta)) / den), complex(
        t01 * t12 * np.exp(1j * beta) / den
    )


def expected(
    k_in: FloatArray,
    e_in: ComplexArray,
    k_out: FloatArray,
    normal: FloatArray,
    cs: complex,
    cp: complex,
) -> ComplexArray:
    """c_s (ŝ·E) ŝ + c_p (p̂_in·E) p̂_out with ŝ ∝ k × n̂, p̂ = ŝ × k̂."""
    s = np.cross(k_in, normal)
    s = (
        s / np.linalg.norm(s)
        if np.linalg.norm(s) > 1e-9 * np.linalg.norm(k_in)
        else np.array([0.0, 1.0, 0.0])
    )
    p_in = np.cross(s, k_in / np.linalg.norm(k_in))
    p_out = np.cross(s, k_out / np.linalg.norm(k_out))
    return np.asarray(cs * (s @ e_in) * s + cp * (p_in @ e_in) * p_out)


@pytest.mark.parametrize("from_back", [False, True])
def test_w83_stack_per_wave_against_airy(from_back: bool) -> None:
    """W83: a film (n = 2.0, 0.3 µm) between air and glass, every wave of a discrete set (θ up to
    50°, any azimuth) and of a continuous spectrum against the Airy formula, both sides of the
    stack; measured 1.9e−13 at most; dΩ behind the stack Δ²k/(k k_z)."""
    ns = (1.0, 2.0, 1.5)
    node = Node("st", "fg.stack", params={"config": _stack_config(ns[0], [(ns[1], 0.3)], ns[2])})
    n_in = ns[2] if from_back else ns[0]
    order = (ns[2], ns[1], ns[0]) if from_back else ns
    normal = np.array([0.0, 0.0, 1.0])  # the stack's normal: towards the front (air)
    plane = Plane.make([0, 0, 0], normal if from_back else -normal)
    theta = np.radians([0.0, 10.0, 30.0, 50.0])
    phi = np.radians([0.0, 20.0, 70.0, -120.0])
    kpar = n_in * K0 * np.column_stack([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi)])
    rng = np.random.default_rng(1)
    e = rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
    inputs = [
        discrete(plane, n_in, kpar, e),
        gaussian_spectrum(plane, n_in, np.array([0.3, -0.2]) * n_in * K0, 0.05 * K0, 7),
    ]
    side = "back" if from_back else "front"
    for sb in inputs:
        k_in, e_in = glob(sb)
        for socket, reflect in (
            (f"{side}_out", True),
            (f"{'front' if from_back else 'back'}_out", False),
        ):
            ctx = Context(Graph())
            out = stack.stack(ctx, node, socket, {f"{side}_in": [sb]})[0]
            k_out, e_out = glob(out)
            sin_in = np.linalg.norm(k_in[:, :2], axis=1) / (n_in * K0)
            # beyond the critical angle (from the glass) nothing is transmitted
            passes = np.ones(len(k_in), dtype=bool) if reflect else n_in * sin_in < order[2]
            assert len(k_out) == int(passes.sum()) >= 3
            worst = 0.0
            for j, i in enumerate(np.flatnonzero(passes)):
                r_s, t_s = airy(order, 0.3, float(sin_in[i]), "s")
                r_p, t_p = airy(order, 0.3, float(sin_in[i]), "p")
                cs, cp = (r_s, r_p) if reflect else (t_s, t_p)
                ref = expected(k_in[i], e_in[i], k_out[j], normal, cs, cp)
                worst = max(worst, float(np.abs(e_out[j] - ref).max() / np.abs(ref).max()))
                assert np.allclose(k_out[j, :2], k_in[i, :2], rtol=0, atol=1e-9 * K0)
            assert worst < 1e-11, (sb.continuous, socket, worst)
            if sb.continuous:
                w_in, w_out = sb.waves, out.waves
                assert w_in.solid_angle is not None and w_out.solid_angle is not None
                k1 = out.refractive_index * K0
                kz1 = np.real(w_out.port.kz(w_out.kx, w_out.ky))
                kz0 = np.real(w_in.port.kz(w_in.kx, w_in.ky))
                dk2 = (np.asarray(w_in.solid_angle) * n_in * K0 * kz0)[passes]
                np.testing.assert_allclose(
                    np.asarray(w_out.solid_angle), dk2 / (k1 * kz1), rtol=1e-12
                )


def test_w83_condensed_plate_per_wave_against_airy() -> None:
    """W83 for `fusion.condensed`: the plate of R29 (n = 1.5, three inner branches) condensed into
    its plane-wave table; plane waves at directions on the table's grid (where the table is exact)
    against the first three terms t₀₁t₁₂e^(iβ)·Σ (r₁₀r₁₂e^(2iβ))^m of the Airy sum, s and p,
    phase referred to the exit face; measured 2e−12."""
    from test_condense import D_MM, LAM_UM, N_PLATE, _registry, plate_model

    from rtt_fusion.graph import replace_block, unroll

    reg, ev = _registry()
    reg.add(NodeType("t.waves", Category.SOURCE, (light_out("out", "S"),)))
    k0 = 2e3 * np.pi / LAM_UM
    sigma = np.array([[0.1, 0.0], [0.0, -0.2], [0.1, 0.1], [0.0, 0.0]])
    plane = Plane.make([0.0, 0.0, -0.5], [0.0, 0.0, 1.0])
    e = np.array([[1.0, 0.5j], [0.3, 1.0], [0.7 - 0.2j, -0.4], [1.0, 0.0]], dtype=complex)
    waves = PlaneWaves(Port(0.0, 1, 1.0 + 0j, LAM_UM), sigma[:, 0] * k0, sigma[:, 1] * k0, e)
    sb = sp.SpectrumBundle(KEY, plane, waves)

    def emit(ctx: Any, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
        return [sb]

    emit.representations = frozenset("S")  # type: ignore[attr-defined]
    ev["t.waves"] = emit
    g = unroll(plate_model({"sigma_max": 0.2, "count": 5}), reg)
    replace_block(g, reg, "src", Node("src", "t.waves"))
    ex = Executor(g, reg, ev, available=WITH_SPECTRA)
    out = ex.evaluate("plate", "out_out")
    assert len(out) == 3 and all(isinstance(b, sp.SpectrumBundle) for b in out)
    assert ex.plan is not None and ex.plan.representations["src.out -> plate.in_in"] == "S"
    got = wave_set(out, 4)
    k_in, e_in = glob(sp.propagate(sb, Plane.make([0, 0, 0], [0, 0, 1.0])))
    normal = np.array([0.0, 0.0, 1.0])
    ns = (1.0, N_PLATE, 1.0)
    for j in range(len(k_in)):
        sin0 = float(np.linalg.norm(sigma[j]))
        terms = []
        for pol in ("s", "p"):
            c = [np.sqrt(1 - (sin0 / n) ** 2 + 0j) for n in ns]

            def fresnel(i: int, m: int, pol: str = pol, c: Any = c) -> tuple[complex, complex]:
                ni, nj, ci, cj = ns[i], ns[m], c[i], c[m]
                if pol == "s":
                    return (ni * ci - nj * cj) / (ni * ci + nj * cj), 2 * ni * ci / (
                        ni * ci + nj * cj
                    )
                return (nj * ci - ni * cj) / (nj * ci + ni * cj), 2 * ni * ci / (nj * ci + ni * cj)

            _, t01 = fresnel(0, 1)
            r12, t12 = fresnel(1, 2)
            beta = 2 * np.pi / LAM_UM * N_PLATE * D_MM * 1e3 * c[1]
            q = r12 * r12 * np.exp(2j * beta)  # r₁₀ = r₁₂ (glass to air on both faces)
            terms.append(t01 * t12 * np.exp(1j * beta) * (1 + q + q * q))
        ref = expected(k_in[j], e_in[j], k_in[j], normal, terms[0], terms[1])
        e_got = got[tuple(np.round(k_in[j], 4))]
        assert float(np.abs(e_got - ref).max() / np.abs(ref).max()) < 1e-9


# ------------------------------------------------------------------------------ port data


def _spectra() -> list[sp.SpectrumBundle]:
    from fieldtracing.field import GriddedField
    from test_spectrum_links import gaussian

    pft = sp.field_to_spectrum(gaussian(100 * 0.13), method="pft")
    grid = sp.field_to_spectrum(gaussian(0.0, samples=31))
    n = 32
    x = (np.arange(n) - n // 2) * 0.002
    xx, yy = np.meshgrid(x, x)
    values = np.zeros((n, n, 3), dtype=complex)
    values[..., 0] = np.exp(-(xx**2 + yy**2) / 0.01**2)
    gridded = FieldBundle(
        KEY,
        Plane.make([0, 0, 0], [0, 0, 1]),
        GriddedField(0.002, values, 0.6, 1.0, np.zeros((2, 2))),
    )
    kpar = np.array([[100.0, 0.0], [-50.0, 20.0]])
    plain = discrete(DOWN, 1.0, kpar, np.ones((2, 2), dtype=complex))
    return [pft, grid, sp.field_to_spectrum(gridded), plain]


def test_spectra_through_port_data(tmp_path: Path) -> None:
    from test_portdata import same

    for i, sb in enumerate(_spectra()):
        handle = json.loads(json.dumps(portdata.write(sb, tmp_path / str(i)), allow_nan=False))
        back = portdata.read(handle, tmp_path / str(i))
        assert same(back, sb), i
        if sb.support is not None:
            a, b = sp.spectrum_to_field(sb).field, sp.spectrum_to_field(back).field
            for name in ("amplitude", "values"):
                if hasattr(a, name):
                    assert np.array_equal(
                        np.asarray(getattr(a, name)), np.asarray(getattr(b, name))
                    )
    kinds = [type(sb.support).__name__ for sb in _spectra()]
    assert kinds == ["GridlessField", "GridSupport", "GridSupport", "NoneType"]


def test_grating_spectrum_in_another_process(tmp_path: Path) -> None:
    """GF6 for S: the grating's spectrum path behind a JSON channel (the FourierGrating plugin of
    Raytatouille Fusion) gives every order as the local one, with its sampling."""
    from rtt_fusion.graph.remote import BlockProvider, register_remote

    node = Node("gr", "fg.grating", params={"config": LAMELLAR, "side": "transmitted"})
    sb = sp.field_to_spectrum(_beam_small())
    local = grating.grating(Context(Graph()), node, "out", {"in": [sb]})
    provider = BlockProvider([grating.TYPE], dict(grating.EVALUATORS))
    assert provider.catalogue()[0]["representations"] == ["H", "R", "S"]

    def call(name: str, request: dict[str, Any]) -> dict[str, Any]:
        wire = json.loads(json.dumps(request, allow_nan=False))
        answer: dict[str, Any] = json.loads(json.dumps(provider.handle(wire), allow_nan=False))
        return answer

    reg = builtin_registry()
    catalogue = [{"provider": "fg", "type": e} for e in provider.catalogue()]
    remote = register_remote(reg, catalogue, call, tmp_path / "portdata")["fg.grating"]
    again = remote(Context(Graph()), node, "out", {"in": [sb]})
    assert [b.path for b in again] == [b.path for b in local] and len(local) >= 2
    for a, b in zip(again, local, strict=True):
        assert isinstance(a, sp.SpectrumBundle) and isinstance(a.support, sp.GridSupport)
        np.testing.assert_array_equal(np.asarray(a.waves.field), np.asarray(b.waves.field))
        fa: Any = sp.spectrum_to_field(a).field
        fb: Any = sp.spectrum_to_field(b).field
        np.testing.assert_array_equal(np.asarray(fa.amplitude), np.asarray(fb.amplitude))


def _beam_small() -> FieldBundle:
    from test_grating_field import _graph_beam

    return rays_to_field(_graph_beam(10.0, radius=0.02, rings=4, waist=0.008, lam=LAM))


# ------------------------------------------------------------------------------ planner


def test_existing_auto_plans_with_a_grating_sample_are_unchanged() -> None:
    """The Linnik with a grating sample (W22/W23, GF5): with gratings and stacks taking S, "auto"
    stays R everywhere (and H where the user fixed it, with the same conversions as before)."""
    from rtt_fusion.graph import backend_rtt, graph_from_raytatouille, replace_block
    from rtt_fusion.graph.workflow import native_of

    d = json.loads((DATA / "linnik_mini.rtt.json").read_text("utf-8"))
    reg = builtin_registry()
    ev = dict(backend_rtt.EVALUATORS)
    ev.update(grating.register(reg))
    for fields in (False, True):
        g = graph_from_raytatouille(d)
        replace_block(g, reg, "sample", Node("sample", "fg.grating"))
        if fields:
            for lk in list(g.links):
                if "sample" in (lk.source[0], lk.target[0]):
                    g.set_transport(lk, "H")
        for native in (None, native_of(g, ev)):
            p = plan(g, reg, Requirements(available=WITH_SPECTRA, native=native))
            assert p.ok
            expected = {
                link_key(lk): "H" if fields and "sample" in (lk.source[0], lk.target[0]) else "R"
                for lk in g.links
                if link_key(lk) in p.representations
            }
            assert p.representations == expected
            chains = [(c.node, c.chain) for c in p.conversions]
            assert chains == ([("objective", "RH"), ("objective", "HR")] if fields else [])


def _chain_graph() -> tuple[Executor, Graph]:
    """A source that emits a spectrum (S only) -> grating (reflected orders, upwards) -> stack
    (through its back) -> detector."""
    reg = builtin_registry()
    ev: dict[str, Any] = dict(optics.register(reg))
    ev.update(grating.register(reg))
    ev.update(stack.register(reg))
    reg.add(NodeType("t.waves", Category.SOURCE, (light_out("out", "S"),)))
    sb = _incident()

    def emit(ctx: Any, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
        return [sb]

    emit.representations = frozenset("S")  # type: ignore[attr-defined]
    ev["t.waves"] = emit
    g = Graph()
    g.add(Node("src", "t.waves"))
    g.add(Node("gr", "fg.grating", params={"config": LAMELLAR, "side": "reflected"}))
    film = _stack_config(1.5, [(2.0, 0.3)], 1.0)  # the reflected orders enter from its back (air)
    g.add(Node("st", "fg.stack", params={"config": film, "point": [0, 0, 0.05]}))
    det = {
        "point": [0, 0, 0.1],
        "normal": [0, 0, -1],
        "half_size": [0.05, 0.05],
        "pixels": [16, 16],
    }
    g.add(Node("cam", "optics.detector", params=det))
    g.connect(("src", "out"), ("gr", "in"))
    g.connect(("gr", "out"), ("st", "back_in"))
    g.connect(("st", "front_out"), ("cam", "in"))
    return Executor(g, reg, ev, available=WITH_SPECTRA), g


def test_spectra_flow_between_native_blocks() -> None:
    """Decision 3: between blocks that take S natively "auto" admits S, and since R30f chooses it
    outright behind a spectrum (before: R by cost, the spectrum kept anyway) – the orders of the
    grating reach the stack and the detector as spectra, without a conversion, and the plan says
    so."""
    ex, g = _chain_graph()
    image = ex.evaluate("cam", "image")
    p = ex.plan
    assert p is not None
    assert p.representations["src.out -> gr.in"] == "S"
    assert p.representations["gr.out -> st.back_in"] == "S"
    assert p.representations["st.front_out -> cam.in"] == "S"
    assert p.reasons["gr.out -> st.back_in"].startswith("spectra from gr; S kept")
    assert not p.conversions
    log = "\n".join(ex.context.log)
    assert "S -> H" not in log and "S -> R" not in log
    assert "image from the sum of" in log
    assert float(np.asarray(image.values).max()) > 0
    g.set_transport(g.links[1], "S")  # fixed S: the same light (and the same plan)
    ex2 = Executor(g, ex.registry, ex.evaluators, available=WITH_SPECTRA)
    assert np.allclose(ex2.evaluate("cam", "image").values, image.values, rtol=0, atol=0)


def test_orders_that_overlap_are_reported_as_rays() -> None:
    """S -> R behind a grating (plan S3, `exec.spectrum_rays`): orders whose spectra overlap – a
    narrow beam on a coarse grating – are reported; separated orders are not."""
    from test_spectrum_links import gaussian

    sb = sp.field_to_spectrum(gaussian(0.0, w0=0.005, samples=15))
    for period, overlap in ((10.0, True), (1.0, False)):
        cfg = copy.deepcopy(_phase_grating(period / 2, 4))
        cfg["period_um"] = {"x": period}
        cfg["layers"][0]["breakpoints_um"] = [0, period / 2, period]
        params = {"config": cfg, "side": "transmitted", "field_transfer": "interpolated"}
        node = Node("gr", "fg.grating", params={**params, "point": [0, 0, 0.001]})
        out = grating.grating(Context(Graph()), node, "out", {"in": [sb]})
        notes = sp.overlapping_orders(out)
        assert len(out) >= 3
        assert bool(notes) is overlap, (period, notes)
    # in the executor: a spectrum link into a block that takes rays only (S -> H -> R)
    reg = builtin_registry()
    ev: dict[str, Any] = dict(optics.register(reg))
    ev.update(grating.register(reg))
    reg.add(NodeType("t.waves", Category.SOURCE, (light_out("out", "S"),)))

    def emit(ctx: Any, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
        return [sb]

    emit.representations = frozenset("S")  # type: ignore[attr-defined]
    ev["t.waves"] = emit
    g = Graph()
    g.add(Node("src", "t.waves"))
    cfg = _phase_grating(5.0, 4)
    params = {"config": cfg, "side": "transmitted", "field_transfer": "interpolated"}
    g.add(Node("gr", "fg.grating", params={**params, "point": [0, 0, 0.001]}))
    g.add(Node("lens", "optics.thin_lens", params={"point": [0, 0, 0.002], "focal_length": 5.0}))
    g.connect(("src", "out"), ("gr", "in"))
    g.connect(("gr", "out"), ("lens", "in"))
    ex = Executor(g, reg, ev, available=WITH_SPECTRA)
    rays = ex.evaluate("lens", "out")
    assert rays
    assert any("exec.spectrum_rays: the orders" in line for line in ex.context.log)
