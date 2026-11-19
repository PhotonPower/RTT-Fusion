"""BindOptiX as a data node (D13 step 3, RTT-Fusion#106): `bindoptix.signal` with F, M, Q and B
from files and from `biooptomd.layer`.

References, independent of the code under test (W93):
- (a) the detected rate per receptor written out from the signal model of BindOptiX (ADR draft
  "Datenports", section 3; Lakowicz for σ): R = P · σ · ½ n c ε₀ · 3|μ̂·E|² / (h c/λ_ex) · q · η_d,
  σ = ln 10 · ε · 1000 cm³/L / N_A, with a uniform field (½ n c ε₀ |E|² = I per mW) along x and a
  table A_tot = A_rad = 1 (q = q₀), A_det = η_d · 1; per state the mean over the samples of Q;
- (b) DP7: the occupancy mixes the states linearly, R(t) = Σ_s f_s(t) R_s; only "free" → R_free;
- (c) DP8: Q goes in unconverted – with a field linear in z the rate follows the heights d of the
  tracker above the surface (sensor frame, solution z > 0), not mirrored ones;
- (d) the regions of the table are bound by the parameter, errors for unknown or doubled ones.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion._types import FloatArray
from rtt_fusion.adapters import biooptomd as bmd
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import ExecutionError, Executor

bx = pytest.importorskip("bindoptix")
import bindoptix.binding  # noqa: E402
import bindoptix.ports  # noqa: E402
import bindoptix.signal  # noqa: E402

H, C, NA, EPS0 = 6.62607015e-34, 299792458.0, 6.02214076e23, 8.8541878128e-12
EPS_DYE, Q0, LAM_EX, LAM_EM = 80000.0, 0.8, 0.633, 0.67
N_HOST, I_PER_MW, POWER = 1.333, 2.5e7, 3.0  # W/m² per mW, mW
ETA = (0.3, 0.2)  # collection of the two regions


def _series(seed: int, frames: int = 30) -> FloatArray:
    rng = np.random.default_rng(seed)
    r = rng.uniform(-5.0, 5.0, (2, frames, 3))
    mu = rng.normal(size=(2, frames, 3))
    mu /= np.linalg.norm(mu, axis=-1, keepdims=True)
    d = rng.uniform(2.0, 18.0, (2, frames, 1))  # nm above the surface
    return np.concatenate([r, mu, d], axis=-1)


def _point_grid() -> Any:
    return bindoptix.ports.RegularGrid((np.zeros(1), np.zeros(1), np.zeros(1)))


def _table(regions: tuple[str, ...] = ("substrate", "solution")) -> Any:
    eye = np.eye(3).reshape(1, 1, 1, 3, 3)
    det = np.stack([e * eye for e in ETA[: len(regions)]], axis=3)
    return bindoptix.ports.EmitterResponseTable(
        _point_grid(), LAM_EM, N_HOST, regions, eye, eye.copy(), det, {"solver": "test"}
    )


def _occupancy() -> Any:
    t = np.array([0.0, 10.0, 20.0])
    f = np.zeros((3, 2, 2))
    f[:, 0, 0] = 1.0  # spot 0: always free
    f[:, 1, 1] = [0.0, 0.4, 0.9]  # spot 1: binds over time
    f[:, 1, 0] = 1.0 - f[:, 1, 1]
    return bindoptix.binding.Occupancy(t, (0, 1), ("free", "bound"), f)


def _write(
    tmp: Path, field: Any, table: Any, occupancy: Any, q: bmd.DipoleSeries
) -> dict[str, str]:
    paths = {"f": tmp / "f.json", "m": tmp / "m.json", "b": tmp / "b.json"}
    paths["f"].write_text(field.to_json(), encoding="utf-8")
    paths["m"].write_text(table.to_json(), encoding="utf-8")
    paths["b"].write_text(occupancy.to_json(), encoding="utf-8")
    bmd.write_layer(tmp / "run", dipoles=q)
    return {k: str(v) for k, v in paths.items()}


def _graph(paths: dict[str, str], tmp: Path, **signal: Any) -> Executor:
    reg, ev = workflow.registry()
    g = Graph()
    g.add(Node("F", "data.field_map", params={"file": paths["f"]}))
    g.add(Node("M", "data.response_map", params={"file": paths["m"]}))
    g.add(Node("B", "data.occupancy", params={"file": paths["b"]}))
    g.add(Node("Q", "biooptomd.layer", params={"directory": str(tmp / "run")}))
    params = {
        "epsilon_per_molar_cm": EPS_DYE,
        "quantum_yield": Q0,
        "excitation_um": LAM_EX,
        "emission_um": LAM_EM,
        "power_mw": POWER,
        "regions": "camera=substrate, ring=solution",
        **signal,
    }
    g.add(Node("sig", "bindoptix.signal", params=params))
    g.connect(("F", "field_map"), ("sig", "excitation"))
    g.connect(("M", "response_map"), ("sig", "response"))
    g.connect(("Q", "dipoles"), ("sig", "dipoles"))
    g.connect(("B", "occupancy"), ("sig", "occupancy"))
    return Executor(g, reg, ev)


def _rate_per_state(series: FloatArray, e_of_z: Any) -> FloatArray:
    """(a): mean over the samples of P σ I_eff/(hν) q₀ η_d, per region."""
    a = series.reshape(-1, 7)
    sigma = math.log(10.0) * EPS_DYE * 1000.0 / NA * 1e-4  # m²
    e = e_of_z(a[:, 6] * 1e-3)  # (N, 3) V/m per mW at the height d (µm)
    i_eff = 0.5 * N_HOST * C * EPS0 * 3.0 * np.abs(np.einsum("ni,ni->n", a[:, 3:6], e)) ** 2
    r = POWER * sigma * i_eff / (H * C / (LAM_EX * 1e-6)) * Q0
    return np.asarray([r.mean() * eta for eta in ETA])


def _uniform_e(z_um: FloatArray) -> FloatArray:
    amp = math.sqrt(2.0 * I_PER_MW / (N_HOST * C * EPS0))
    return np.tile([amp, 0.0, 0.0], (len(z_um), 1))


def test_w93_rates_follow_the_signal_model_and_the_occupancy(tmp_path: Path) -> None:
    free, bound = _series(1), _series(2)
    q = bmd.DipoleSeries({"free": free, "bound": bound}, dt_ps=5.0)
    field = bindoptix.signal.uniform_field_map(I_PER_MW, LAM_EX, N_HOST)
    paths = _write(tmp_path, field, _table(), _occupancy(), q)
    res = _graph(paths, tmp_path).evaluate("sig", "rates")
    assert res.labels == ("camera", "ring") and res.spots == (0, 1)
    r_free, r_bound = _rate_per_state(free, _uniform_e), _rate_per_state(bound, _uniform_e)
    f = _occupancy().fractions
    want = f[..., 0:1] * r_free + f[..., 1:2] * r_bound  # (b) DP7
    assert res.rates.shape == (3, 2, 2)
    assert np.abs(res.rates - want).max() < 1e-12 * want.max()
    assert np.abs(res.rates[:, 0, :] - r_free).max() < 1e-12 * r_free.max()  # only free


def test_w93_q_enters_in_the_sensor_frame(tmp_path: Path) -> None:
    """(c) DP8: E_x = a + b z on z ∈ [0, 0.02] µm; the rate follows the tracker heights d."""
    a, b = 1.0e4, 4.0e5  # V/m, V/m per µm
    z = np.array([0.0, 0.02])
    e = np.zeros((1, 1, 2, 3), dtype=complex)
    e[0, 0, :, 0] = a + b * z
    grid = bindoptix.ports.RegularGrid((np.zeros(1), np.zeros(1), z))
    field = bindoptix.ports.FieldMap(grid, LAM_EX, e, source={"solver": "test"})
    s = _series(3)
    q = bmd.DipoleSeries({"free": s, "bound": s})
    paths = _write(tmp_path, field, _table(), _occupancy(), q)
    res = _graph(paths, tmp_path).evaluate("sig", "rates")

    def linear(z_um: FloatArray) -> FloatArray:
        out = np.zeros((len(z_um), 3))
        out[:, 0] = a + b * z_um
        return out

    want = _rate_per_state(s, linear)
    assert np.abs(res.rates[0, 0] - want).max() < 1e-12 * want.max()
    mirrored = _rate_per_state(s * [1, 1, 1, 1, 1, -1, -1], linear)  # z → −z would differ
    assert np.abs(mirrored - want).max() > 0.1 * want.max()


def test_w93_regions_are_bound_by_the_parameter(tmp_path: Path) -> None:
    q = bmd.DipoleSeries({"free": _series(4), "bound": _series(5)})
    field = bindoptix.signal.uniform_field_map(I_PER_MW, LAM_EX, N_HOST)
    paths = _write(tmp_path, field, _table(), _occupancy(), q)
    one = _graph(paths, tmp_path, regions="ring=solution").evaluate("sig", "rates")
    both = _graph(paths, tmp_path).evaluate("sig", "rates")
    assert one.labels == ("ring",) and np.array_equal(one.rates[..., 0], both.rates[..., 1])
    for text, words in (
        ("", "bind them"),
        ("camera=nowhere", "not in the response table"),
        ("a=substrate, a=solution", "given twice"),
        ("a=substrate, b=substrate", "bound twice"),
    ):
        with pytest.raises(ExecutionError) as e:
            _graph(paths, tmp_path, regions=text).evaluate("sig", "rates")
        assert e.value.code == "exec.data" and words in str(e.value), text
    # a table with one region needs no binding
    single = _table(("camera",))
    (tmp_path / "m1.json").write_text(single.to_json(), encoding="utf-8")
    paths1 = {**paths, "m": str(tmp_path / "m1.json")}
    res = _graph(paths1, tmp_path, regions="").evaluate("sig", "rates")
    assert res.labels == ("camera",)


def test_w93_wrong_inputs_are_diagnosed(tmp_path: Path) -> None:
    q = bmd.DipoleSeries({"free": _series(6)})  # no "bound" ensemble
    field = bindoptix.signal.uniform_field_map(I_PER_MW, LAM_EX, N_HOST)
    paths = _write(tmp_path, field, _table(), _occupancy(), q)
    with pytest.raises(ExecutionError) as e:
        _graph(paths, tmp_path).evaluate("sig", "rates")
    assert e.value.code == "exec.data" and "bound" in str(e.value)
    swapped = {**paths, "f": paths["m"]}  # a response table where a field map belongs
    with pytest.raises(ExecutionError) as e:
        _graph(swapped, tmp_path).evaluate("sig", "rates")
    assert e.value.code == "exec.data" and "EmitterResponseTable" in str(e.value)


def test_w93_rates_plot_over_time(tmp_path: Path) -> None:
    from rtt_fusion.graph.display import plot_view

    q = bmd.DipoleSeries({"free": _series(7), "bound": _series(8)})
    field = bindoptix.signal.uniform_field_map(I_PER_MW, LAM_EX, N_HOST)
    paths = _write(tmp_path, field, _table(), _occupancy(), q)
    res = _graph(paths, tmp_path).evaluate("sig", "rates")
    view = plot_view(None, Node("p", "display.plot"), "view", {"in": res})  # type: ignore[arg-type]
    assert len(view.curves) == 4 and view.x_label == "time (s)"
    assert view.curves[3].y == [float(x) for x in res.rates[:, 1, 1]]
