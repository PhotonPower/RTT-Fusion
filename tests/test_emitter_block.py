"""The emitter block `emitter.dipole_layer` (pilot K5, stage K5b).

References (W95), independent of the block's own assembly:
- (a) the power of each spectrum is `power_w` times the share `rates(...).collected` of
  `couplings.dipole_layers` within the block's NA (checked there against the free dipole, the
  image dipole and energy conservation, W94);
- (b) every wave is the angular spectrum of `dipole_layers.angular_spectrum` in the sensor frame,
  mirrored (x, y, z) → (u, v, −n) for the substrate side (the rule of `in_rtt_frame`, DP8), times
  its solid angle, with the phase of the move onto the port plane;
- (c) Q: one coherence group per emitter; two emitters that differ only laterally by Δ have
  spectra that differ by e^(−i k∥·Δ) (Fourier shift theorem); the power adds over the emitters;
- (d) N: an index profile equal to the solvent changes nothing.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from rtt_fusion.adapters import biooptomd as bmd
from rtt_fusion.couplings import dipole_layers as dl
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import Executor

pytest.importorskip("fieldtracing")
from rtt_fusion.graph.spectrum import spectrum_power  # noqa: E402

LAM, NS, NW = 0.6, 1.52, 1.333
BASE = {"wavelength_um": LAM, "n_substrate": NS, "n_solution": NW, "n_theta": 48, "na": 1.45}


def _run(params: dict[str, Any], q: bmd.DipoleSeries | None = None, tmp: Any = None) -> list[Any]:
    reg, ev = workflow.registry()
    g = Graph()
    g.add(Node("e", "emitter.dipole_layer", params={**BASE, **params}))
    if q is not None:
        bmd.write_layer(tmp, dipoles=q)
        g.add(Node("q", "biooptomd.layer", params={"directory": str(tmp)}))
        g.connect(("q", "dipoles"), ("e", "dipoles"))
    out: list[Any] = Executor(g, reg, ev, available=workflow.available()).evaluate("e", "out")
    return out


@pytest.mark.parametrize("dipole", [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.6, 0.0, 0.8]])
def test_w95_power_is_the_collected_share(dipole: list[float]) -> None:
    for layers in ([], [[50.0, 0.197, 3.44]]):
        (sb,) = _run({"dipole": dipole, "height_nm": 15.0, "layers": layers, "power_w": 2.0})
        stack = dl.Stack(
            NS,
            tuple((complex(r[1], r[2]), r[0] * 1e-6) for r in layers),
            NW,
            -sum(r[0] for r in layers) * 1e-6,
        )
        want = 2.0 * dl.rates(stack, LAM, 15e-6, dipole, na=1.45).collected
        assert spectrum_power(sb) == pytest.approx(want, rel=2e-4), layers


def test_w95_waves_are_the_mirrored_angular_spectrum() -> None:
    dip = np.array([0.6, 0.0, 0.8])
    (sb,) = _run({"dipole": list(dip), "height_nm": 15.0, "distance": 0.0, "power_w": 1.0})
    w = sb.waves
    k = NS * 2 * np.pi / (LAM * 1e-3)
    kz = np.sqrt(k**2 - w.kx**2 - w.ky**2)
    theta = np.arcsin(np.hypot(w.kx, w.ky) / k)
    phi = np.arctan2(w.ky, w.kx)
    stack = dl.Stack(NS, (), NW, 0.0)
    a = dl.angular_spectrum(stack, LAM, 15e-6, dip, theta, phi, "bottom")
    scale = np.sqrt(1.0 / (2 * dl.free_power(stack, LAM, 15e-6)))
    want = a[:, :2] * scale * np.asarray(w.solid_angle)[:, None]  # u, v = x, y (mirror)
    assert np.abs(np.asarray(w.field) - want).max() < 1e-9 * np.abs(want).max()
    # the longitudinal part follows from the mirror (A_n = −A_z)
    ez = -(w.kx * w.field[:, 0] + w.ky * w.field[:, 1]) / kz
    want_z = -a[:, 2] * scale * np.asarray(w.solid_angle)
    assert np.abs(ez - want_z).max() < 1e-9 * np.abs(want_z).max()


def test_w95_emitters_are_incoherent_and_shifted(tmp_path: Any) -> None:
    mu = np.array([0.0, 0.6, 0.8])
    rows = np.array(
        [[[10.0, 20.0, 0.0, *mu, 8.0], [14.0, 17.0, 0.0, *mu, 8.0], [12.0, 18.5, 0.0, *mu, 12.0]]]
    )  # (1 replica, 3 frames, 7), nm
    out = _run({"samples": 3, "distance": 0.0}, bmd.DipoleSeries({"bound": rows}), tmp_path)
    assert len({sb.key for sb in out}) == 3
    a, b = out[0].waves, out[1].waves
    d = (rows[0, 1, :2] - rows[0, 0, :2]) * 1e-6  # mm, (u, v) = (x, y)
    shift = np.exp(-1j * (a.kx * d[0] + a.ky * d[1]))
    assert np.abs(np.asarray(b.field) - np.asarray(a.field) * shift[:, None]).max() < 1e-9 * (
        np.abs(np.asarray(a.field)).max()
    )
    stack = dl.Stack(NS, (), NW, 0.0)
    want = sum(dl.rates(stack, LAM, r[6] * 1e-6, r[3:6], na=1.45).collected for r in rows[0])
    assert sum(spectrum_power(sb) for sb in out) == pytest.approx(want, rel=2e-4)


def test_w95_solvent_profile_changes_nothing(tmp_path: Any) -> None:
    z = (np.arange(100) + 0.5) * 0.1
    flat = bmd.Profile(z, np.full((1, 100), NW), 0.1, NW, LAM)
    reg, ev = workflow.registry()
    bmd.write_layer(tmp_path, profile=flat)
    g = Graph()
    g.add(Node("e", "emitter.dipole_layer", params={**BASE, "height_nm": 6.0}))
    g.add(Node("q", "biooptomd.layer", params={"directory": str(tmp_path)}))
    g.connect(("q", "profile"), ("e", "profile"))
    (with_n,) = Executor(g, reg, ev, available=workflow.available()).evaluate("e", "out")
    (without,) = _run({"height_nm": 6.0})
    assert np.abs(np.asarray(with_n.waves.field) - np.asarray(without.waves.field)).max() < 1e-12
