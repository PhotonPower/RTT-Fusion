"""Pilot K1 (R3): FourierGrating efficiencies as a table in Raytatouille's sequential tracer.

Bench: lamellar transmission grating (period 0.7 µm, n = 2.1 / air, 0.25 µm deep) in vacuum,
λ = 0.6 µm, paths for the orders −1, 0, +1 (tests/data/pilot_k1_grating.rtt.json, derived from
Raytatouille tests/reference/m4/grating_transmission.rtt.json). Checks:

1. direction of every order traced by Raytatouille = direction of the same order from
   FourierGrating (grating equation in both), at normal and oblique incidence – this fixes the
   mapping of the order numbers (Raytatouille m = FourierGrating n);
2. path transmission of every order in Raytatouille = FourierGrating efficiency (s, p,
   unpolarised): no double counting of interface losses;
3. energy: Σ transmitted (Raytatouille) + Σ reflected (FourierGrating) = 1 for the lossless grating.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np
import pytest

try:
    import fouriergrating  # noqa: F401
    import raytatouille as rt
except ImportError:
    if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
        raise
    pytest.skip("fouriergrating or raytatouille not installed", allow_module_level=True)

from rtt_fusion.adapters import fouriergrating as fga
from rtt_fusion.couplings.grating_table import order_efficiencies, with_efficiencies
from rtt_fusion.ports import PlaneWaves

BENCH = Path(__file__).parent / "data" / "pilot_k1_grating.rtt.json"
LAM = 0.6
K0 = 2 * np.pi / (LAM * 1e-3)
CONFIG = {
    "dimension": "1D",
    "solver": "rcwa",
    "truncation": {"Nx": 15, "Nz": 0},
    "period_um": {"x": 0.7, "z": 0.0},
    "materials": {"air": {"n": 1.0}, "hi": {"n": 2.1}},
    "superstrate": "air",
    "substrate": "air",
    "layers": [
        {
            "type": "lamellar1d",
            "thickness_um": 0.25,
            "breakpoints_um": [0.0, 0.3, 0.7],
            "materials": ["hi", "air"],
        }
    ],
}


def _system() -> Any:
    return rt.System.from_json(BENCH.read_text("utf-8"))


def _rays(theta: float, n: int = 5) -> Any:
    rays = rt.trace.RayBatch(n)
    for k in range(n):
        rays.pos_x[k] = -2.0 + 4.0 * k / (n - 1)
        rays.pos_z[k] = -5.0
        rays.dir_x[k] = np.sin(theta)
        rays.dir_z[k] = np.cos(theta)
    return rays


@pytest.mark.parametrize("theta", [0.0, 0.08])
def test_order_directions_agree(theta: float) -> None:
    handle = fga.prepare(fga.Structure(CONFIG), [LAM])
    cs = rt.compile(_system(), rt.MaterialLibrary())
    kx_in = K0 * np.sin(theta)
    eff = order_efficiencies(handle, LAM, kx_in, 0.0, "s")["transmitted"]
    port = fga.front_port(handle, LAM)
    fg_out = fga.apply(
        handle,
        PlaneWaves(
            port, np.array([kx_in]), np.array([0.0]), np.array([[0.0, 1.0]], dtype=np.complex128)
        ),
    )["back"]
    assert fg_out.orders is not None
    fg_kx = {int(o[0]): float(fg_out.kx[j]) for j, o in enumerate(fg_out.orders)}
    for m in (-1, 0, 1):
        rays = _rays(theta)
        rt.trace.trace(cs, rays, path=f"order {m:+d}" if m else "order 0")
        if m in eff:
            assert rt.trace.RayStatus(rays.status[0]) == rt.trace.RayStatus.ALIVE
            # FourierGrating order n = m: k_x = k_x,in + m·2π/Λ
            sx = (kx_in + m * 2 * np.pi / 0.7e-3) / K0
            assert abs(rays.dir_x[0] - sx) < 1e-12
            # the order FourierGrating calls n = m leaves in the same direction
            assert abs(rays.dir_x[0] - fg_kx[m] / K0) < 1e-9
        else:
            assert rt.trace.RayStatus(rays.status[0]) != rt.trace.RayStatus.ALIVE


@pytest.mark.parametrize("polarisation", ["s", "p", "unpolarised"])
@pytest.mark.parametrize("theta", [0.0, 0.08])
def test_path_transmission_equals_efficiency(polarisation: str, theta: float) -> None:
    handle = fga.prepare(fga.Structure(CONFIG), [LAM])
    eff = order_efficiencies(handle, LAM, K0 * np.sin(theta), 0.0, polarisation)
    table = {m: eta for m, eta in eff["transmitted"].items() if m in (-1, 0, 1)}
    cs = rt.compile(with_efficiencies(_system(), "G", table), rt.MaterialLibrary())
    for m in (-1, 0, 1):
        name = f"order {m:+d}" if m else "order 0"
        t = rt.analysis.path_transmission(cs, name, start=_rays(theta))
        assert abs(t.mean - table.get(m, 0.0)) < 1e-12
    if theta == 0.0:  # symmetric grating at normal incidence: η(−1) = η(+1)
        assert abs(table[-1] - table[1]) < 1e-9
    else:
        assert abs(table[-1] - table[1]) > 1e-3  # oblique: the order mapping matters


@pytest.mark.parametrize("theta", [0.0, 0.08])
def test_energy_over_all_orders(theta: float) -> None:
    handle = fga.prepare(fga.Structure(CONFIG), [LAM])
    for pol in ("s", "p"):
        eff = order_efficiencies(handle, LAM, K0 * np.sin(theta), 0.0, pol)
        total = sum(eff["transmitted"].values()) + sum(eff["reflected"].values())
        assert abs(total - 1.0) < 1e-9
        assert set(eff["transmitted"]) <= {-1, 0, 1, -2}  # only propagating orders
