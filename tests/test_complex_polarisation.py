"""Complex polarisation of rays (stage R31, docs/plan-r31-polarisation.md).

References, independent of the code under test:
- W86 (a) the far field of a dipole with elliptical p = (1, 0.3i, 0) (E∞ = k²(r̂ × p) × r̂) as
  the particle gives it, on a detector tilted by 1.5 mrad against the spectrum's plane (the field
  is resampled H -> R -> H): the spherical wave E∞ e^(ikr)/r at the pixels, as for the untilted
  detector. The sphere of R30d (CliffordBem) is W86 (b) in `test_particle_spectra.py`.
- W87 Fresnel rhomb (Born/Wolf 1.5.4): n = 1.51, two total internal reflections at 54.6°, each a
  phase difference δ between s and p with tan(δ/2) = cos θ sqrt(sin²θ − 1/n²)/sin²θ (45.0°);
  light linear at 45° to the plane of incidence leaves circular. The field per ray against the
  boundary conditions of Maxwell's equations solved directly (tangential E and H continuous,
  E ⟂ k), not the Fresnel formulas of the code.
- W88 an elliptically polarised Gaussian beam H -> R -> H on its plane: the Stokes parameters
  per sample and integrated are those of the field before (the conversion propagates nothing).
- W89 linear polarisation as before R31 (the split of the field into a scalar and a real
  direction of R18–R30, written out here) and real polarisation data (port data from before R31)
  read as complex with imaginary part 0.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion import portdata
from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.graph import Graph, Node, builtin_registry, optics
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.light import (
    CoherenceKey,
    Plane,
    RayBundle,
    physical_split,
    split_along,
    split_field,
)

LAM = 0.6
K = 2 * np.pi / (LAM * 1e-3)
N_GLASS, THETA = 1.51, np.radians(54.6)


def stokes(e: ComplexArray, u: FloatArray, v: FloatArray) -> FloatArray:
    """(S0, S1, S2, S3) of field vectors e (N, 3) in the transverse axes u, v (S3 > 0: the field
    turns from u to v, i.e. e ∝ u + i v)."""
    eu, ev = e @ u, e @ v
    return np.stack(
        [
            np.abs(eu) ** 2 + np.abs(ev) ** 2,
            np.abs(eu) ** 2 - np.abs(ev) ** 2,
            2 * np.real(np.conj(eu) * ev),
            2 * np.imag(np.conj(eu) * ev),
        ],
        axis=-1,
    )


# ------------------------------------------------------------------------------ split


def _split_before_r31(f: ComplexArray) -> tuple[ComplexArray, FloatArray]:
    """The split of R18–R30 (`field.field_to_rays`, `grating.split_prt`): the phase of the
    largest component taken out, the real part normalised, the amplitude along it."""
    big = np.argmax(np.abs(f), axis=1)
    alpha = np.angle(f[np.arange(len(f)), big])
    e = np.real(f * np.exp(-1j * alpha)[:, None])
    e = e / np.linalg.norm(e, axis=1)[:, None]
    return np.sum(f * e, axis=1), e


def test_split_field_keeps_every_polarisation() -> None:
    rng = np.random.default_rng(31)
    f = rng.normal(size=(500, 3)) + 1j * rng.normal(size=(500, 3))
    a, e = split_field(f)
    assert np.abs(a[:, None] * e - f).max() < 1e-14 * np.abs(f).max()
    assert np.abs(np.linalg.norm(e, axis=1) - 1).max() < 1e-14
    assert np.allclose(np.abs(a), np.linalg.norm(f, axis=1), rtol=1e-14)
    # circular light and zero fields
    a, e = split_field(np.array([[1.0, 1j, 0.0], [0.0, 0.0, 0.0]]))
    assert abs(a[0]) == pytest.approx(np.sqrt(2), rel=1e-15) and a[1] == 0
    assert np.allclose(e[0] * a[0], [1, 1j, 0]) and np.all(e[1] == 0)


def test_w89_linear_polarisation_splits_as_before() -> None:
    """Linear polarisation (any phase, also two equally large components): ê real and the
    amplitude of R18–R30 to rounding."""
    rng = np.random.default_rng(89)
    d = rng.normal(size=(400, 3))
    d[:50, 2] = 0.0
    d[:50, 1] = -d[:50, 0]  # two equally large components
    phase = np.exp(1j * rng.uniform(-np.pi, np.pi, 400))
    f = d * (rng.uniform(0.1, 2.0, 400) * phase)[:, None]
    a, e = split_field(f)
    a0, e0 = _split_before_r31(f)
    assert np.abs(e.imag).max() < 1e-15
    assert np.abs(e.real - e0).max() < 1e-15
    assert np.abs(a - a0).max() < 1e-15 * np.abs(a0).max()


def test_w89_interface_factor_for_s_and_p_light_as_before() -> None:
    """`physical_split` (interfaces, condensed groups): for s or p light the scalar factor and
    the carried polarisation of R27 (factor = carried†·P·ê/|carried|², prt -> G·prt)."""
    rng = np.random.default_rng(27)
    n = 40
    k = np.tile([0.0, 0.0, 1.0], (n, 1))
    ang = rng.uniform(0, np.pi, n)
    s = np.stack([np.cos(ang), np.sin(ang), 0 * ang], axis=1)  # s ⟂ k
    p_in = np.cross(s, k)
    p_out = np.cross(s, k)  # transmitted at normal incidence: the identity transversally
    geometric = np.einsum("ni,nj->nij", s, s) + np.einsum("ni,nj->nij", p_out, p_in)
    c_s = rng.normal(size=n) + 1j * rng.normal(size=n)
    c_p = rng.normal(size=n) + 1j * rng.normal(size=n)
    physical = c_s[:, None, None] * np.einsum("ni,nj->nij", s, s) + c_p[:, None, None] * np.einsum(
        "ni,nj->nij", p_out, p_in
    )
    for pol, coeff in ((s, c_s), (p_in, c_p)):
        b = _bundle(np.zeros((n, 3)), k, np.ones(n, complex))
        b = b.evolve(prt=np.einsum("ni,j->nij", pol, [1.0, 0.0, 0.0]).astype(complex))
        c, prt = physical_split(physical, geometric, b)
        assert np.abs(c - coeff).max() < 1e-14 * np.abs(coeff).max()
        carried = np.einsum("nij,nj->ni", geometric, b.field_direction())
        assert np.abs(np.einsum("nij,j->ni", prt, [1.0, 0, 0]) - carried).max() < 1e-14
    # mixed light: the full field P·ê, the factor along the carried polarisation
    b = _bundle(np.zeros((n, 3)), k, np.ones(n, complex))
    b = b.evolve(prt=np.einsum("ni,j->nij", s + 1j * p_in, [1.0, 0, 0]).astype(complex))
    c, prt = physical_split(physical, geometric, b)
    full = np.einsum("nij,nj->ni", physical, b.field_direction())
    assert np.abs(c[:, None] * np.einsum("nij,j->ni", prt, [1.0, 0, 0]) - full).max() < 1e-13
    carried = np.einsum("nij,nj->ni", geometric, b.field_direction())
    assert np.allclose(split_along(full, carried), c)


def _bundle(pos: FloatArray, d: FloatArray, amp: ComplexArray) -> RayBundle:
    n = len(pos)
    return RayBundle(
        key=CoherenceKey("t"),
        wavelength_um=LAM,
        position=pos,
        direction=d,
        opl=np.zeros(n),
        amplitude=amp,
        source_xy=pos[:, :2].copy(),
        alive=np.ones(n, dtype=bool),
        jones=np.array([1.0, 0.0, 0.0]),
    )


def test_w89_real_polarisation_data_reads_as_complex(tmp_path: Path) -> None:
    """Port data written before R31 could carry real prt/jones: read as complex128, the same
    field."""
    n = 7
    rng = np.random.default_rng(1)
    rot = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    b = _bundle(rng.normal(size=(n, 3)), np.tile([0.0, 0, 1], (n, 1)), np.ones(n, complex))
    assert b.jones is not None and b.jones.dtype == np.complex128  # real given, complex kept
    b = b.evolve(prt=np.tile(rot, (n, 1, 1)))
    assert b.prt is not None and b.prt.dtype == np.complex128
    assert np.array_equal(b.prt.imag, np.zeros((n, 3, 3)))
    handle = portdata.write(b, tmp_path)
    # as before R31: the arrays of prt and jones real on disk
    fields = handle["root"]["fields"]
    for name in ("prt", "jones"):
        file = tmp_path / fields[name]["file"]
        np.save(file, np.load(file).real, allow_pickle=False)
        assert np.load(file).dtype == np.float64
    back = portdata.read(handle, tmp_path)
    assert isinstance(back, RayBundle)
    assert back.prt is not None and back.prt.dtype == np.complex128
    assert back.jones is not None and back.jones.dtype == np.complex128
    assert np.array_equal(back.field_vectors(), b.field_vectors())


# ------------------------------------------------------------------------------ W86 (a)


def test_w86_elliptical_far_field_on_a_tilted_detector() -> None:
    """The far field of an elliptical dipole on a detector tilted by 1.5 mrad (H -> R -> H)
    equals the spherical wave, as on the untilted detector; before R31 the resampling kept the
    linear part only (up to 13 % less at the diagonals for the sphere of R30d)."""
    pytest.importorskip("fieldtracing")
    sys.path.insert(0, str(Path(__file__).parent))
    from test_particle_spectra import dipole_far, dipole_spectrum

    from rtt_fusion.graph import spectrum as sp
    from rtt_fusion.graph.field import field_image

    p = np.array([1.0, 0.3j, 0.0])
    centre = np.array([0.0, 0.0, 2.0])
    fb = sp.spectrum_to_field(dipole_spectrum(centre, 1.0, p))
    half, n = 0.4, 15
    errors = {}
    for tilt in (0.0, 1.5e-3):
        normal = np.array([0.0, np.sin(tilt), -np.cos(tilt)])
        plane = Plane.make([0.0, 0.0, 1.0], normal, (0.0, 1.0, 0.0))
        log: list[str] = []
        image = field_image([fb], plane, half, half, n, n, log)
        assert ("resampled" in "\n".join(log)) == (tilt > 0)
        u = (np.arange(n) + 0.5) / n * 2 * half - half
        uu, vv = np.meshgrid(u, u)
        pts = plane.origin + uu.reshape(-1, 1) * plane.u + vv.reshape(-1, 1) * plane.v
        rel = pts - centre
        dist = np.linalg.norm(rel, axis=1)
        ref = np.sum(np.abs(dipole_far(rel / dist[:, None], p)) ** 2, axis=1) / dist**2
        ref = ref.reshape(n, n)
        errors[tilt] = float(np.abs(image.values - ref).max() / ref.max())
        ratio = image.values / ref
        assert ratio.min() > 0.99  # before R31 tilted: 0.906 (0.918 at the diagonal)
    # the sampling error of the field (as W85 (b)): measured 3.9e-3 untilted, 5.0e-3 tilted
    # (before R31 tilted: 8.7e-2; a linear dipole p = (1, 0.3, 0) tilted: 4.9e-3)
    assert errors[0.0] < 1e-2 and errors[1.5e-3] < 1e-2
    assert errors[1.5e-3] == pytest.approx(errors[0.0], abs=2e-3)


# ------------------------------------------------------------------------------ W87


def maxwell_reflection(
    e_in: ComplexArray, k_in: FloatArray, normal: FloatArray, n1: float, n2: float
) -> ComplexArray:
    """The reflected field of a plane wave (field e_in, unit direction k_in in medium n1) at a
    plane interface to n2 (`normal` into medium 2), from the boundary conditions solved directly:
    tangential E and tangential H ∝ k × E continuous, E ⟂ k for both new waves (the
    transmitted one evanescent under total internal reflection, Im k_z > 0)."""
    nn = normal / np.linalg.norm(normal)
    kt_par = n1 * (k_in - (k_in @ nn) * nn)  # tangential wave vector / k0, conserved
    k_r = n1 * k_in - 2 * n1 * (k_in @ nn) * nn
    kz2 = np.sqrt(complex(n2**2 - kt_par @ kt_par))
    kz2 = kz2 if kz2.imag > 0 or (kz2.imag == 0 and kz2.real >= 0) else -kz2
    k_t = kt_par.astype(complex) + kz2 * nn
    k_i = n1 * k_in
    t1 = np.cross(nn, [1.0, 0.0, 0.0] if abs(nn[0]) < 0.9 else [0.0, 1.0, 0.0])
    t1 = t1 / np.linalg.norm(t1)
    t2 = np.cross(nn, t1)

    def cross_matrix(k: Any) -> Any:
        return np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]], dtype=complex)

    rows, rhs = [], []
    for t in (t1, t2):  # tangential E: E_i + E_r = E_t
        rows.append(np.concatenate([t, -t]))
        rhs.append(-(t @ e_in))
    for t in (t1, t2):  # tangential H ∝ k × E
        rows.append(np.concatenate([t @ cross_matrix(k_r), -(t @ cross_matrix(k_t))]))
        rhs.append(-(t @ cross_matrix(k_i) @ e_in))
    rows.append(np.concatenate([k_r, np.zeros(3)]))  # transversality
    rhs.append(0.0)
    rows.append(np.concatenate([np.zeros(3), k_t]))
    rhs.append(0.0)
    x = np.linalg.solve(np.array(rows, dtype=complex), np.array(rhs, dtype=complex))
    return np.asarray(x[:3])


def rhomb(pol: Any) -> tuple[Graph, Any, Any]:
    """Source in glass → TIR at face 1 → TIR at the parallel face 2 (a Fresnel rhomb)."""
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "center": [0.0, 0.0, -3.0],
                "radius": 0.3,
                "rings": 3,
                "wavelengths_um": [LAM],
                "polarisation": pol,
            },
        )
    )
    n1 = [np.sin(THETA), 0.0, -np.cos(THETA)]  # towards the glass, the light comes from there
    k1 = [np.sin(2 * THETA), 0.0, -np.cos(2 * THETA)]
    for name, point, normal in (
        ("f1", [0.0, 0.0, 0.0], n1),
        ("f2", list(6.0 * np.asarray(k1)), [-x for x in n1]),
    ):
        params = {"point": point, "normal": normal, "n_front": N_GLASS, "n_back": 1.0}
        g.add(Node(name, "optics.interface", params=params))
    g.connect(("src", "out"), ("f1", "front_in"))
    g.connect(("f1", "front_out"), ("f2", "front_in"))
    return g, reg, ev


def test_w87_fresnel_rhomb_makes_circular_light() -> None:
    sin2 = np.sin(THETA) ** 2
    delta = 2 * np.arctan(np.cos(THETA) * np.sqrt(sin2 - 1 / N_GLASS**2) / sin2)
    assert np.degrees(delta) == pytest.approx(45.0, abs=0.05)  # Born/Wolf 1.5.4
    s_axis, p_axis = np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0])
    pol = (s_axis + p_axis) / np.sqrt(2)
    g, reg, ev = rhomb([float(x) for x in pol])
    out: RayBundle = Executor(g, reg, ev).evaluate("f2", "front_out")[0]
    assert np.all(out.alive)
    k = np.asarray(out.direction)
    assert np.abs(k - [0.0, 0.0, 1.0]).max() < 1e-12  # parallel faces: the direction back
    # the reference: Maxwell's boundary conditions at both faces
    n1 = np.array([np.sin(THETA), 0.0, -np.cos(THETA)])
    k1 = np.array([np.sin(2 * THETA), 0.0, -np.cos(2 * THETA)])
    e1 = maxwell_reflection(pol.astype(complex), np.array([0.0, 0, 1]), -n1, N_GLASS, 1.0)
    e2 = maxwell_reflection(e1, k1, n1, N_GLASS, 1.0)
    got = out.field_vectors()
    assert np.abs(got - e2).max() < 1e-12  # measured
    st = stokes(got, p_axis, s_axis)
    assert np.abs(st[:, 0] - 1).max() < 1e-12  # lossless
    # circular: S3/S0 = ±sin(2δ), the phase difference of two reflections
    assert np.abs(np.abs(st[:, 3]) - abs(np.sin(2 * delta))).max() < 1e-12
    assert np.abs(st[:, 3]).min() > 0.99999
    # the same light as a field (R -> H) keeps it
    from rtt_fusion.graph.field import rays_to_field

    pytest.importorskip("fieldtracing")
    fb = rays_to_field(out)
    rot = np.stack([fb.plane.u, fb.plane.v, fb.plane.normal])
    st_h = stokes(np.asarray(fb.field.amplitude) @ rot, p_axis, s_axis)
    assert np.abs(st_h[:, 3] / st_h[:, 0] - st[0, 3]).max() < 1e-12


def test_w87_half_a_rhomb_is_a_quarter_of_the_way() -> None:
    """One total internal reflection: the field against the boundary conditions, the phase
    difference δ (45°) between s and p."""
    g, reg, ev = rhomb([np.sqrt(0.5), np.sqrt(0.5), 0.0])
    out: RayBundle = Executor(g, reg, ev).evaluate("f1", "front_out")[0]
    n1 = np.array([np.sin(THETA), 0.0, -np.cos(THETA)])
    e1 = maxwell_reflection(
        np.array([np.sqrt(0.5), np.sqrt(0.5), 0.0], complex), np.array([0.0, 0, 1]), -n1, N_GLASS, 1
    )
    assert np.abs(out.field_vectors() - e1).max() < 1e-12
    eu = out.field_vectors() @ np.cross([0.0, 1, 0], np.asarray(out.direction)[0])
    ev_ = out.field_vectors() @ np.array([0.0, 1, 0])
    sin2 = np.sin(THETA) ** 2
    delta = 2 * np.arctan(np.cos(THETA) * np.sqrt(sin2 - 1 / N_GLASS**2) / sin2)
    assert (
        np.abs(np.abs(np.angle(eu / ev_)) - delta).max() < 1e-12
        or np.abs(np.abs(np.angle(-eu / ev_)) - delta).max() < 1e-12
    )


# ------------------------------------------------------------------------------ W88


def test_w88_elliptical_gaussian_round_trip_keeps_stokes() -> None:
    pytest.importorskip("fieldtracing")
    from rtt_fusion.graph.field import field_power, field_to_rays, rays_to_field

    chi, psi = 0.35, 0.6  # ellipticity angle and orientation of the ellipse
    jones = np.array(
        [
            np.cos(psi) * np.cos(chi) - 1j * np.sin(psi) * np.sin(chi),
            np.sin(psi) * np.cos(chi) + 1j * np.cos(psi) * np.sin(chi),
            0.0,
        ]
    )
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    params = {"radius": 0.25, "rings": 20, "gaussian_waist": 0.1, "polarisation": list(jones)}
    g.add(Node("src", "optics.collimated_source", params=params))
    bundle = Executor(g, reg, ev).evaluate("src", "out")[0]
    plane = Plane.make([0, 0, 1.0], [0, 0, 1])
    fb = rays_to_field(bundle, plane)
    back = rays_to_field(field_to_rays(fb), plane)
    rot = np.stack([plane.u, plane.v, plane.normal])
    st0 = stokes(np.asarray(fb.field.amplitude) @ rot, plane.u, plane.v)
    st1 = stokes(np.asarray(back.field.amplitude) @ rot, plane.u, plane.v)
    pts0, pts1 = np.asarray(fb.field.points), np.asarray(back.field.points)
    assert len(pts1) == len(pts0) and np.abs(pts1 - pts0).max() < 1e-15
    assert np.abs(st1 - st0).max() < 1e-12 * st0[:, 0].max()  # measured
    # the state: S3/S0 = sin 2χ, S2/S1 = tan 2ψ – a beam of one polarisation
    expected = np.array([np.cos(2 * chi) * np.cos(2 * psi), np.cos(2 * chi) * np.sin(2 * psi)])
    assert np.abs(st1[:, 1:3] / st1[:, :1] - expected).max() < 1e-12
    assert np.abs(st1[:, 3] / st1[:, 0] - np.sin(2 * chi)).max() < 1e-12
    weight = st0[:, 0] / st0[:, 0].sum()
    total0, total1 = weight @ st0, weight @ st1
    assert np.abs(total1 - total0).max() < 1e-12 * total0[0]
    assert field_power(back) == pytest.approx(field_power(fb), rel=1e-12)
