"""The data block `biooptomd.layer` (D13 step 2, RTT-Fusion#106): a BioOptoMD run directory as
data ports N (index profile) and Q (dipole series per receptor state).

References, independent of the code under test:
- W92 (a) the files written as BioOptoMD's formats come back unchanged: the profile as the CSV of
  `export_refractive_index` read with `read_profile_csv`, the tracker series bit for bit with
  replicas and frames (the mixing by occupancy needs them, D13);
- W92 (b) Q stays in BioOptoMD's sensor frame; the conversion into the RTT frame is the mirror
  rule written out here (BioOptoMD docs/UNITS_AND_CONVENTIONS.md: nm, surface z = 0, solution
  z > 0): "substrate" x, y, z = R_x, R_y, z_s + d, μ̂ as is; "solution" z = z_s − d, μ̂_z → −μ̂_z;
- W92 (c) a rewritten run is read again by the same executor and cache (files in the key);
- W92 (d) missing or wrong files are `exec.data`, a missing directory `exec.param`.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from rtt_fusion._types import FloatArray
from rtt_fusion.adapters import biooptomd as bmd
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import ExecutionError, Executor, ResultCache

DZ, N_W, LAM = 0.05, 1.333, 0.633


def _profile() -> bmd.Profile:
    z = (np.arange(200) + 0.5) * DZ
    dn = 0.12 * np.exp(-z / 3.0)
    return bmd.Profile(z, np.vstack([N_W + dn, N_W + 0.97 * dn]), DZ, N_W, LAM, "synthetic")


def _series(seed: int, replicas: int = 3, frames: int = 50) -> FloatArray:
    rng = np.random.default_rng(seed)
    r = rng.uniform(-20.0, 20.0, (replicas, frames, 3))
    mu = rng.normal(size=(replicas, frames, 3))
    mu /= np.linalg.norm(mu, axis=-1, keepdims=True)
    d = rng.uniform(1.0, 15.0, (replicas, frames, 1))
    return np.concatenate([r, mu, d], axis=-1)


def _run(directory: str | Path, cache: ResultCache | None = None, **params: object) -> Executor:
    reg, ev = workflow.registry()
    g = Graph()
    g.add(Node("bmd", "biooptomd.layer", params={"directory": str(directory), **params}))
    return Executor(g, reg, ev, cache=cache)


def test_w92_run_directory_comes_back_unchanged(tmp_path: Path) -> None:
    free, bound = _series(1), _series(2, frames=40)
    q = bmd.DipoleSeries({"free": free, "bound": bound}, dt_ps=10.0, source="synthetic")
    bmd.write_layer(tmp_path, _profile(), q, source="BioOptoMD synthetic")
    ex = _run(tmp_path)
    got_q = ex.evaluate("bmd", "dipoles")
    assert isinstance(got_q, bmd.DipoleSeries)
    assert got_q.states == ["free", "bound"]
    assert np.array_equal(got_q.series["free"], free)
    assert np.array_equal(got_q.series["bound"], bound)
    assert got_q.series["bound"].shape == (3, 40, 7) and got_q.dt_ps == 10.0
    profile = ex.evaluate("bmd", "profile")
    ref = bmd.read_profile_csv(tmp_path / "n_profile.csv", DZ, N_W, LAM)
    assert np.array_equal(profile.n, ref.n) and np.array_equal(profile.z_nm, ref.z_nm)
    assert np.abs(profile.n - _profile().n).max() < 1e-15
    # the manifest itself is accepted as the parameter
    same = _run(tmp_path / bmd.MANIFEST).evaluate("bmd", "dipoles")
    assert np.array_equal(same.series["free"], free)


@pytest.mark.parametrize("light_from", ["substrate", "solution"])
def test_w92_q_converts_with_the_mirror_rule(tmp_path: Path, light_from: str) -> None:
    s = _series(3)
    bmd.write_layer(tmp_path, dipoles=bmd.DipoleSeries({"bound": s}))
    q = _run(tmp_path).evaluate("bmd", "dipoles")
    got = q.in_rtt_frame("bound", surface_z_mm=0.2, light_from=light_from)
    a = s.reshape(-1, 7)
    sign = 1.0 if light_from == "substrate" else -1.0
    want_pos = np.column_stack([a[:, 0] * 1e-6, a[:, 1] * 1e-6, 0.2 + sign * a[:, 6] * 1e-6])
    want_mu = a[:, 3:6] * [1.0, 1.0, sign]
    assert np.abs(got.positions_mm - want_pos).max() < 1e-18
    assert np.abs(got.orientations - want_mu).max() < 1e-15


def test_w92_a_rewritten_run_is_read_again(tmp_path: Path) -> None:
    cache = ResultCache()
    bmd.write_layer(tmp_path, dipoles=bmd.DipoleSeries({"free": _series(4)}))
    first = _run(tmp_path, cache).evaluate("bmd", "dipoles")
    again = _run(tmp_path, cache).evaluate("bmd", "dipoles")
    assert again is first  # from the cache
    bmd.write_layer(tmp_path, dipoles=bmd.DipoleSeries({"free": _series(5)}))
    npy = tmp_path / "dipoles_0.npy"
    st = npy.stat()
    os.utime(npy, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))  # coarse file clocks
    new = _run(tmp_path, cache).evaluate("bmd", "dipoles")
    assert np.array_equal(new.series["free"], _series(5))


def test_w92_wrong_runs_are_diagnosed(tmp_path: Path) -> None:
    with pytest.raises(ExecutionError) as e:
        _run("").evaluate("bmd", "dipoles")
    assert e.value.code == "exec.param"
    with pytest.raises(ExecutionError) as e:
        _run(tmp_path / "nowhere").evaluate("bmd", "dipoles")
    assert e.value.code == "exec.data" and "manifest" in str(e.value)
    bmd.write_layer(tmp_path, profile=_profile())
    with pytest.raises(ExecutionError) as e:
        _run(tmp_path).evaluate("bmd", "dipoles")
    assert e.value.code == "exec.data" and "no dipole series" in str(e.value)
    np.save(tmp_path / "bad.npy", np.zeros((2, 5, 6)))
    (tmp_path / bmd.MANIFEST).write_text(
        '{"format": "biooptomd-layer", "version": 1, "dipoles": {"states": {"x": "bad.npy"}}}'
    )
    with pytest.raises(ExecutionError) as e:
        _run(tmp_path).evaluate("bmd", "dipoles")
    assert e.value.code == "exec.data" and "(n_replicas, n_frames, 7)" in str(e.value)
    (tmp_path / bmd.MANIFEST).write_text('{"format": "other", "version": 1}')
    with pytest.raises(ExecutionError) as e:
        _run(tmp_path).evaluate("bmd", "profile")
    assert e.value.code == "exec.data" and "format" in str(e.value)
