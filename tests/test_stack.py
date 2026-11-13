"""Layer stack `fg.stack` (ADR grating as field operator and layer stacks, stage F4a, ST0–ST2).

References independent of FourierGrating: the transfer (characteristic) matrix of thin-film optics
(Macleod, *Thin-Film Optical Filters*, ch. 2) written out here for s and p, Fresnel's equations
for the bare interface, and RTT-Fusion's own `optics.interface` for the full field.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import Graph, Node, builtin_registry, optics, stack
from rtt_fusion.graph.execute import Context, ExecutionError, Executor
from rtt_fusion.graph.light import CoherenceKey, RayBundle

pytest.importorskip("fouriergrating")


# ------------------------------------------------------------------ reference: transfer matrix
def tmm(
    n_in: complex,
    layers: Sequence[tuple[complex, float]],
    n_out: complex,
    theta: float,
    lam: float,
    pol: str,
) -> tuple[float, float]:
    """Power reflectance and transmittance of a stack (Macleod's characteristic matrix; indices
    N = n − ik in his convention, i.e. the conjugate of FourierGrating's n + ik)."""

    def conj(n: complex) -> complex:
        return complex(n).conjugate()

    n0 = conj(n_in)
    s0 = n0 * np.sin(theta)  # invariant n sin θ

    def cos_of(n: complex) -> complex:
        c = complex(np.sqrt(1 - (s0 / n) ** 2 + 0j))  # Macleod: the root with Re > 0
        return -c if c.real < 0 else c

    def eta(n: complex) -> complex:
        return n * cos_of(n) if pol == "s" else n / cos_of(n)

    m = np.eye(2, dtype=complex)
    for n, d in layers:
        nn = conj(n)
        delta = 2 * np.pi * nn * d * cos_of(nn) / lam
        e = eta(nn)
        m = m @ np.array(
            [[np.cos(delta), 1j * np.sin(delta) / e], [1j * e * np.sin(delta), np.cos(delta)]]
        )
    e0, es = eta(n0), eta(conj(n_out))
    b, c = m @ np.array([1.0, es])
    r = (e0 * b - c) / (e0 * b + c)
    t_pow = 4 * e0.real * es.real / abs(e0 * b + c) ** 2
    return float(abs(r) ** 2), float(t_pow)


# ------------------------------------------------------------------ the block
def _config(
    tmp_path: Path, n_in: Any, layers: list[tuple[Any, float]], n_out: Any, name: str
) -> str:
    """`layers` from the incident side (front) on, as the transfer matrix takes them; written the
    way FourierGrating lists them, from the substrate upwards."""
    materials = {"m_in": {"n": n_in}, "m_out": {"n": n_out}}
    cfg_layers = []
    for k, (n, d) in reversed(list(enumerate(layers))):
        materials[f"l{k}"] = {"n": n}
        cfg_layers.append({"type": "homogeneous", "thickness_um": d, "material": f"l{k}"})
    cfg = {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 0},
        "period_um": {"x": 1.0},
        "materials": materials,
        "superstrate": "m_in",
        "substrate": "m_out",
        "layers": cfg_layers,
    }
    path = tmp_path / f"{name}.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    return str(path)


def _rays(
    theta: float, pol: str, lam: float, from_back: bool = False, index: float = 1.0
) -> RayBundle:
    """Two rays in the x-z plane towards the plane z = 0 (normal +z) from the front (z > 0) or the
    back; s = y, p in the plane of incidence."""
    sz = 1.0 if from_back else -1.0
    d = np.array([np.sin(theta), 0.0, sz * np.cos(theta)])
    jones = (
        np.array([0.0, 1.0, 0.0])
        if pol == "s"
        else np.array([np.cos(theta), 0.0, -sz * np.sin(theta)])
    )
    start = -d * 1.0 + np.array([[0.0, 0.0, 0.0], [0.01, 0.02, 0.0]])
    return RayBundle(
        key=CoherenceKey("t"),
        wavelength_um=lam,
        position=start,
        direction=np.tile(d, (2, 1)),
        opl=np.zeros(2),
        amplitude=np.ones(2, dtype=complex),
        source_xy=np.zeros((2, 2)),
        alive=np.ones(2, dtype=bool),
        jones=jones.astype(complex),
        index=index,
        source_index=index,
    )


def _power(b: RayBundle | None) -> float:
    if b is None:
        return 0.0
    f = np.asarray(b.amplitude)[:, None] * b.field_direction()
    p = np.where(np.asarray(b.alive), np.sum(np.abs(f) ** 2, axis=1), 0.0)
    assert abs(p[0] - p[1]) < 1e-12  # both rays alike (laterally invariant)
    return float(p[0])


def _run(path: str, rays: RayBundle, from_back: bool, out: str) -> RayBundle | None:
    node = Node("st", stack.TYPE_ID, params={"fouriergrating": path})
    ctx = Context(Graph())
    result = stack.stack(ctx, node, out, {"back_in" if from_back else "front_in": [rays]})
    return result[0] if result else None


def _rt(
    path: str, theta: float, pol: str, lam: float, from_back: bool = False, index: float = 1.0
) -> tuple[float, float]:
    rays = _rays(theta, pol, lam, from_back, index)
    near, far = ("back_out", "front_out") if from_back else ("front_out", "back_out")
    return _power(_run(path, rays, from_back, near)), _power(_run(path, rays, from_back, far))


@pytest.mark.parametrize("pol", ["s", "p"])
@pytest.mark.parametrize("theta_deg", [0.0, 30.0, 60.0, 80.0])
def test_st0_bare_interface_is_fresnel(tmp_path: Path, pol: str, theta_deg: float) -> None:
    path = _config(tmp_path, 1.0, [], 1.5, "bare")
    theta = np.radians(theta_deg)
    r, t = _rt(path, theta, pol, 0.6)
    r_ref, t_ref = tmm(1.0, [], 1.5, theta, 0.6, pol)
    # 1e-11: at grazing incidence (80°) FourierGrating rounds to a few 1e-12
    assert r == pytest.approx(r_ref, abs=1e-11) and t == pytest.approx(t_ref, abs=1e-11)
    assert r + t == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("pol", ["s", "p"])
def test_st0_full_field_equals_optics_interface(tmp_path: Path, pol: str) -> None:
    """The field leaving the bare stack equals RTT-Fusion's Fresnel interface, direction, index and
    complex field per ray (reflected and transmitted)."""
    path = _config(tmp_path, 1.0, [], 1.5, "bare")
    theta = np.radians(40.0)
    rays = _rays(theta, pol, 0.6)
    plane = optics._plane(Node("i", "optics.interface", params={}))
    for out, reflect in (("front_out", True), ("back_out", False)):
        mine = _run(path, rays, False, out)
        ref = optics._fresnel(rays, plane, 1.0, 1.5, reflect, "i")
        assert mine is not None
        np.testing.assert_allclose(mine.direction, ref.direction, atol=1e-12)
        assert mine.index == ref.index
        f_mine = np.asarray(mine.amplitude)[:, None] * mine.field_direction()
        f_ref = np.asarray(ref.amplitude)[:, None] * ref.field_direction()
        np.testing.assert_allclose(f_mine, f_ref, atol=1e-12)


def test_st1_quarter_wave_antireflection(tmp_path: Path) -> None:
    lam, n_sub = 0.55, 1.52
    n_c = float(np.sqrt(n_sub))
    path = _config(tmp_path, 1.0, [(n_c, lam / (4 * n_c))], n_sub, "ar")
    for pol in ("s", "p"):
        r, t = _rt(path, 0.0, pol, lam)
        assert r < 1e-12 and t == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("pol", ["s", "p"])
def test_st2_bragg_mirror_against_the_transfer_matrix(tmp_path: Path, pol: str) -> None:
    lam0, n_h, n_l, n_sub = 0.55, 2.4, 1.46, 1.52
    layers = [(n_h, lam0 / (4 * n_h)), (n_l, lam0 / (4 * n_l))] * 5 + [(n_h, lam0 / (4 * n_h))]
    path = _config(tmp_path, 1.0, layers, n_sub, "bragg")
    worst = 0.0
    for lam in (0.50, 0.55, 0.60):
        for theta_deg in (0.0, 25.0, 50.0):
            theta = np.radians(theta_deg)
            r, t = _rt(path, theta, pol, lam)
            r_ref, t_ref = tmm(1.0, layers, n_sub, theta, lam, pol)
            worst = max(worst, abs(r - r_ref), abs(t - t_ref))
            assert r + t == pytest.approx(1.0, abs=1e-10)
    assert worst < 1e-10, worst
    assert _rt(path, 0.0, pol, lam0)[0] > 0.99  # content guard: it is a mirror


def test_st2_absorbing_layer_and_light_from_the_back(tmp_path: Path) -> None:
    layers = [(1.38, 0.1), ([2.0, 0.3], 0.05)]  # MgF2-like and an absorbing layer (n + ik)
    path = _config(tmp_path, 1.0, layers, 1.52, "absorbing")
    theta = np.radians(20.0)
    for pol in ("s", "p"):
        r, t = _rt(path, theta, pol, 0.6)
        r_ref, t_ref = tmm(1.0, [(1.38, 0.1), (2.0 + 0.3j, 0.05)], 1.52, theta, 0.6, pol)
        assert r == pytest.approx(r_ref, abs=1e-10) and t == pytest.approx(t_ref, abs=1e-10)
        assert r + t < 0.999  # absorbed
        # from the substrate: the stack turned round, incident medium glass
        theta_b = np.radians(15.0)
        rb, tb = _rt(path, theta_b, pol, 0.6, from_back=True, index=1.52)
        rb_ref, tb_ref = tmm(1.52, [(2.0 + 0.3j, 0.05), (1.38, 0.1)], 1.0, theta_b, 0.6, pol)
        assert rb == pytest.approx(rb_ref, abs=1e-10) and tb == pytest.approx(tb_ref, abs=1e-10)


def test_st2_total_internal_reflection(tmp_path: Path) -> None:
    path = _config(tmp_path, 1.0, [], 1.5, "bare")
    theta = np.radians(60.0)  # beyond the critical angle 41.8° from glass
    r, t = _rt(path, theta, "s", 0.6, from_back=True, index=1.5)
    assert r == pytest.approx(1.0, abs=1e-12) and t == 0.0


def test_a_grating_layer_is_refused(tmp_path: Path) -> None:
    cfg = {
        "materials": {"a": {"n": 1.0}},
        "superstrate": "a",
        "substrate": "a",
        "layers": [
            {
                "type": "lamellar1d",
                "thickness_um": 0.1,
                "breakpoints_um": [0, 0.5, 1],
                "materials": ["a", "a"],
            }
        ],
    }
    path = tmp_path / "grating.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    with pytest.raises(ExecutionError, match="homogeneous layers only"):
        _run(str(path), _rays(0.0, "s", 0.6), False, "front_out")


def test_the_stack_in_a_graph(tmp_path: Path) -> None:
    """A collimated beam through an AR-coated window onto a detector: in the registry, with the
    executor, the light arrives with the coating's transmittance."""
    lam, n_sub = 0.55, 1.52
    n_c = float(np.sqrt(n_sub))
    path = _config(tmp_path, 1.0, [(n_c, lam / (4 * n_c))], n_sub, "ar")
    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(stack.register(reg))
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": 0.5, "rings": 4, "center": [0, 0, -1], "wavelengths_um": [lam]},
        )
    )
    g.add(Node("st", stack.TYPE_ID, params={"fouriergrating": path, "normal": [0, 0, -1]}))
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, 1.0],
                "normal": [0, 0, -1],
                "half_size": [1, 1],
                "pixels": [8, 8],
            },
        )
    )
    g.connect(("src", "out"), ("st", "front_in"))
    g.connect(("st", "back_out"), ("cam", "in"))
    ex = Executor(g, reg, ev)
    image = ex.evaluate("cam", "image")
    assert float(np.asarray(image.values).sum()) > 0
    out: list[RayBundle] = ex.evaluate("st", "back_out")
    assert _power(out[0]) == pytest.approx(1.0, abs=1e-12)  # AR: everything through
    assert out[0].index == pytest.approx(n_sub)
