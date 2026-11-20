"""The emitter response table M of a planar stack for BindOptiX (pilot K5, stage K5d).

References (W97), independent of the table's assembly:
- (a) quadratic form: for a tilted dipole μ the rates computed directly (`dipole_layers.rates`,
  `band_power`) equal μᵀAμ of the tensors built from an x and a z dipole only – the stack is
  invariant about z, so there are no cross terms and a_xx = a_yy (also with a gold film);
- (b) DP4: at a near-perfect mirror a_tot is the rate of the image dipole (full near field written
  out, as W94) for both orientations;
- (c) DP6: BindOptiX accepts the table (symmetric, positive semidefinite, A_tot − A_rad and
  A_rad − Σ A_det positive semidefinite), and its quantum yield and collection equal
  a_rad/(a_tot + (1 − q₀)/q₀) and band/(up + down) of the direct rates (DP5);
- (d) the chain `emitter.response_map` → `bindoptix.signal` gives the rate written out from these
  tensors; wrong heights or regions are diagnosed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from rtt_fusion.adapters import bindoptix as bxa
from rtt_fusion.couplings import dipole_layers as dl
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import ExecutionError, Executor

LAM = 0.67
GOLD = dl.Stack(1.515, ((0.16 + 3.8j, 50e-6),), 1.333, -50e-6)
GLASS = dl.Stack(1.515, (), 1.333, 0.0)
REGIONS = [
    ("camera", "substrate", 0.0, 1.45),
    ("saf", "substrate", 1.45, 1.515),  # beyond the camera's NA: bands must not overlap
    ("top", "solution", 0.0, 1.0),
]


@pytest.mark.parametrize("stack", [GLASS, GOLD], ids=["glass", "gold"])
def test_w97_tensors_are_the_quadratic_form(stack: dl.Stack) -> None:
    heights = [0.005, 0.03]  # µm
    tot, rad, det, n_host = bxa.response_tensors(stack, LAM, heights, REGIONS)
    assert n_host == pytest.approx(1.333)
    mu = np.array([0.48, -0.36, 0.8])  # unit, tilted
    for i, h in enumerate(heights):
        r = dl.rates(stack, LAM, h * 1e-3, mu)
        assert mu @ tot[i] @ mu == pytest.approx(r.total, rel=1e-9)
        assert mu @ rad[i] @ mu == pytest.approx(r.up + r.down, rel=1e-9)
        for d, (_, side, lo, hi) in enumerate(REGIONS):
            side_ = "bottom" if side == "substrate" else "top"
            band = dl.band_power(stack, LAM, h * 1e-3, mu, side_, lo, hi)
            assert mu @ det[i, d] @ mu == pytest.approx(band, rel=1e-9, abs=1e-15)


def test_w97_mirror_is_the_image_dipole() -> None:
    k0 = 2 * np.pi / (LAM * 1e-3)
    mirror = dl.Stack(1e5j, (), 1.0, 0.0)
    heights = [0.05, 0.2]
    tot, _, _, _ = bxa.response_tensors(mirror, LAM, heights, [("top", "solution", 0.0, 1.0)])
    for i, h in enumerate(heights):
        d = h * 1e-3
        for axis, sign in ((0, -1.0), (2, 1.0)):
            p = np.zeros(3)
            p[axis] = 1.0
            img = p * np.array([-1.0, -1.0, 1.0])
            r, nh = 2 * d, np.array([0.0, 0.0, 1.0])
            e = (
                (
                    k0**2 * np.cross(np.cross(nh, img), nh) / r
                    + (3 * nh * (nh @ img) - img) * (1 / r**3 - 1j * k0 / r**2)
                )
                * np.exp(1j * k0 * r)
                / (4 * np.pi)
            )
            want = 1 + 6 * np.pi / k0**3 * np.imag(p @ e)
            assert tot[i, axis, axis] == pytest.approx(want, abs=1e-4), (h, axis, sign)


def test_w97_bindoptix_accepts_the_table_and_agrees() -> None:
    pytest.importorskip("bindoptix")
    heights = [0.005, 0.01, 0.03]
    table = bxa.response_table(GOLD, LAM, heights, REGIONS)  # BindOptiX checks DP6 here
    assert table.detectors == ("camera", "saf", "top")
    mu = np.array([[0.48, -0.36, 0.8]])
    q0 = 0.6
    for h in heights:
        pt = np.array([[0.0, 0.0, h]])
        r = dl.rates(GOLD, LAM, h * 1e-3, mu[0])
        assert table.quantum_yield(pt, mu, q0)[0] == pytest.approx(
            (r.up + r.down) / (r.total + (1 - q0) / q0), rel=1e-9
        )
        band = dl.band_power(GOLD, LAM, h * 1e-3, mu[0], "bottom", 0.0, 1.45)
        assert table.collection(pt, mu)[0, 0] == pytest.approx(band / (r.up + r.down), rel=1e-9)


def test_w97_response_map_feeds_the_signal(tmp_path: Path) -> None:
    pytest.importorskip("bindoptix")
    import bindoptix.binding
    import bindoptix.signal

    from rtt_fusion.adapters import biooptomd as bmd

    reg, ev = workflow.registry()
    field = bindoptix.signal.uniform_field_map(2.0e7, 0.633, 1.333)
    (tmp_path / "f.json").write_text(field.to_json(), encoding="utf-8")
    occ = bindoptix.binding.Occupancy(np.array([0.0]), (0,), ("free",), np.ones((1, 1, 1)))
    (tmp_path / "b.json").write_text(occ.to_json(), encoding="utf-8")
    series = np.array([[[0.0, 0.0, 0.0, 0.6, 0.0, 0.8, 10.0]]])  # one tilted dipole, 10 nm
    bmd.write_layer(tmp_path / "run", dipoles=bmd.DipoleSeries({"free": series}))
    g = Graph()
    g.add(Node("F", "data.field_map", params={"file": str(tmp_path / "f.json")}))
    g.add(Node("B", "data.occupancy", params={"file": str(tmp_path / "b.json")}))
    g.add(Node("Q", "biooptomd.layer", params={"directory": str(tmp_path / "run")}))
    g.add(
        Node(
            "M",
            "emitter.response_map",
            params={
                "wavelength_um": LAM,
                "n_substrate": 1.515,
                "n_solution": 1.333,
                "layers": [[50.0, 0.16, 3.8]],
                "heights_nm": [5.0, 10.0, 20.0],
                "regions": "camera=substrate:0:1.45",
            },
        )
    )
    dye = {"epsilon_per_molar_cm": 80000.0, "quantum_yield": 0.6, "excitation_um": 0.633}
    g.add(Node("S", "bindoptix.signal", params={**dye, "emission_um": LAM, "power_mw": 1.0}))
    for a, b in (("F", "excitation"), ("M", "response"), ("Q", "dipoles"), ("B", "occupancy")):
        g.connect(
            (a, {"F": "field_map", "M": "response_map", "Q": "dipoles", "B": "occupancy"}[a]),
            ("S", b),
        )
    res = Executor(g, reg, ev).evaluate("S", "rates")
    mu = [0.6, 0.0, 0.8]
    r = dl.rates(GOLD, LAM, 1e-5, mu)
    band = dl.band_power(GOLD, LAM, 1e-5, mu, "bottom", 0.0, 1.45)
    sigma = np.log(10) * 80000.0 * 1000.0 / 6.02214076e23 * 1e-4  # m²
    r_exc = sigma * 2.0e7 * 3.0 * 0.6**2 / (6.62607015e-34 * 299792458.0 / 0.633e-6)
    want = r_exc * (r.up + r.down) / (r.total + (1 - 0.6) / 0.6) * band / (r.up + r.down)
    assert res.rates[0, 0, 0] == pytest.approx(want, rel=1e-6)  # BindOptiX's constants
    assert res.labels == ("camera",)


def test_w97_wrong_parameters_are_diagnosed() -> None:
    pytest.importorskip("bindoptix")
    reg, ev = workflow.registry()
    for params, words in (
        ({"heights_nm": [-2.0, 5.0]}, "heights_nm"),
        ({"regions": "camera=side:0:1"}, "substrate|solution"),
        ({"regions": "a=substrate:0:1, a=solution:0:1"}, "twice"),
        ({"regions": "a=substrate:0:1.45, b=substrate:1.3:1.5"}, "overlap"),
    ):
        g = Graph()
        g.add(Node("M", "emitter.response_map", params=params))
        with pytest.raises(ExecutionError) as e:
            Executor(g, reg, ev).evaluate("M", "response_map")
        assert e.value.code == "exec.param" and words in str(e.value), params
