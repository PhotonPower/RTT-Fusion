"""Grating block with hp-FEM in the loop (R35): `hpfem.grating` builds its Jones table from an
hp-FEM job.

References, independent of the block:
- W99 a flat interface (no scatterer, n = 1 → 1.5): the reflected and transmitted order 0 carry
  the Fresnel powers |r_s|², |r_p|² and 1 − |r_s|², 1 − |r_p|² (analytic);
- W100 a dielectric bar (lossless): every order's power equals the efficiency of hp-FEM's own
  `grating.solve` at the same incidence (the solver's result without the adapter, the table and
  the PRT chain), and the orders of both sides carry all the incident power (energy);
- W101 a table around a bundle of directions: its interpolation equals the direct Jones matrices
  at a direction between its nodes (`jones_table.direct_jones`), and a bundle wider than
  `max_nodes` stops before solving.
Needs hp-FEM (skipped without it, except with RTT_FUSION_REQUIRE_MODULES=1); small cells, a few
seconds.
"""

from __future__ import annotations

import copy
import os
from typing import Any

import numpy as np
import pytest

try:
    import hpfem  # noqa: F401
except ImportError:
    if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
        raise
    pytest.skip("hpfem not installed", allow_module_level=True)

from rtt_fusion.adapters import hpfem as hpa
from rtt_fusion.couplings import jones_table as jt
from rtt_fusion.graph import Graph, Node, builtin_registry, hpfem_grating, optics
from rtt_fusion.graph.execute import Context, ExecutionError, Executor

LAMBDA_UM = 0.5
N_SUB = 1.5

#: the cell of the hp-FEM plugin's example (nm): air over glass, optionally a bar of n = 1.5
FLAT: dict[str, Any] = {
    "version": 1,
    "problem": "grating",
    "unit": 1e-9,
    "model": {
        "period": 400,
        "y_bottom": -800,
        "y_top": 950,
        "background_tag": 1,
        "slabs": [{"tag": 2, "y0": -800, "y1": 0}],
        "shapes": [],
    },
    "mesh": {"structured": {"nx": 16, "rows": [[-800, 0, 24], [0, 148, 6], [148, 950, 24]]}},
    "materials": {"1": {"n": [1.0, 0.0]}, "2": {"n": [N_SUB, 0.0]}, "3": {"eps": [2.25, 0.0]}},
    "stack": {"incidence": {"n": [1.0, 0.0]}, "substrate": {"n": [N_SUB, 0.0]}, "layers": []},
    "incidence": {"polarisation": "p", "theta_deg": 50, "wavelength": 400},
    "solver": {"order": 3, "pml": {"top": 419, "bottom": 419}, "orders_max": 2},
}
BAR = copy.deepcopy(FLAT)
BAR["model"]["shapes"] = [
    {"kind": "rectangle", "tag": 3, "params": {"x": -100, "y": 0, "width": 200, "height": 148}}
]


def run_block(
    job: dict[str, Any], side: str, alpha_deg: float, polarisation: list[float]
) -> tuple[dict[int, float], float, list[str]]:
    """Order powers behind the block for a collimated beam at `alpha_deg` in the x-z plane, the
    beam's σ_x and the log."""
    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(hpfem_grating.register(reg))
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": 0.01,
                "rings": 3,
                "center": [0, 0, -0.1],
                "wavelengths_um": [LAMBDA_UM],
                "field_angles_deg": [[alpha_deg, 0.0]],
                "polarisation": polarisation,
            },
        )
    )
    g.add(Node("gr", "hpfem.grating", params={"job_data": job, "side": side, "orders": []}))
    g.connect(("src", "out"), ("gr", "in"))
    ex = Executor(g, reg, ev)
    sigma_x = float(np.asarray(ex.evaluate("src", "out")[0].direction)[0, 0])
    powers: dict[int, float] = {}
    for b in ex.evaluate("gr", "out"):
        m = int(b.path[-1].split("[")[1].rstrip("]"))
        alive = np.asarray(b.alive)
        field = np.asarray(b.amplitude)[alive][:, None] * np.einsum(
            "nij,j->ni", np.asarray(b.prt)[alive], np.asarray(b.jones)
        )
        powers[m] = float(np.mean(np.sum(np.abs(field) ** 2, axis=1)))
    return powers, sigma_x, ex.context.log


def fresnel(theta: float, n1: float, n2: float) -> tuple[complex, complex]:
    """Amplitude reflection coefficients r_s, r_p (Born & Wolf 1.5.2, eq. 21a)."""
    ci = np.cos(theta)
    ct = np.sqrt(1 - (n1 / n2 * np.sin(theta)) ** 2)
    rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct)
    rp = (n2 * ci - n1 * ct) / (n2 * ci + n1 * ct)
    return complex(rs), complex(rp)


TE = [0.0, 1.0, 0.0]  # along the lines (the job's z): s
TM = [1.0, 0.0, 0.0]  # in the plane of incidence: p


@pytest.mark.parametrize(("pol", "which"), [(TE, 0), (TM, 1)])
def test_w99_flat_interface_is_fresnel(pol: list[float], which: int) -> None:
    alpha = 30.0
    r = abs(fresnel(np.radians(alpha), 1.0, N_SUB)[which]) ** 2
    reflected, _, log = run_block(FLAT, "reflected", alpha, pol)
    transmitted, _, _ = run_block(FLAT, "transmitted", alpha, pol)
    assert abs(reflected[0] - r) < 1e-6
    assert abs(transmitted[0] - (1 - r)) < 1e-6
    assert all(p < 1e-8 for m, p in {**reflected, **transmitted}.items() if m != 0)
    assert any("Jones matrices from hp-FEM at one incidence" in line for line in log)


@pytest.mark.parametrize(("pol", "name"), [(TE, "s"), (TM, "p")])
def test_w100_bar_equals_hpfem_and_conserves_energy(pol: list[float], name: str) -> None:
    from hpfem import grating, run, units

    alpha = 20.0
    reflected, sigma_x, _ = run_block(BAR, "reflected", alpha, pol)
    transmitted, _, _ = run_block(BAR, "transmitted", alpha, pol)
    # hp-FEM directly, without adapter, table and PRT: efficiencies of the same cell
    unit = BAR["unit"]
    omega = units.angular_frequency(wavelength=LAMBDA_UM * 1e-6)
    cell = run._unit_cell(BAR["model"], unit)
    mats = {int(t): run._material(m, "m") for t, m in BAR["materials"].items()}
    stack = run._stack(BAR["stack"], omega, unit)
    mesh = run._build_mesh(BAR, cell, mats, omega, unit, ".")
    phi = 0.0 if sigma_x > 0 else np.pi
    res = grating.solve(
        mesh, mats, stack, name, np.radians(alpha), phi, omega, 3,
        pml={"top": 419e-9, "bottom": 419e-9}, orders_max=2,
    )  # fmt: skip
    for mine, theirs in ((reflected, res.R_orders), (transmitted, res.T_orders)):
        expected = {o.m: o.efficiency for o in theirs if o.propagating}
        assert {m for m in mine if mine[m] > 1e-12} <= set(expected)
        for m, eff in expected.items():
            if m in (-1, 0, 1):
                assert abs(mine.get(m, 0.0) - eff) < 1e-6 * max(1.0, eff), (m, mine, expected)
    total = sum(reflected.values()) + sum(transmitted.values())
    assert abs(total - 1.0) < 1e-3  # lossless; the FEM balance residual (here ≤ 7e-5)


def _ctx() -> Context:
    return Context(Graph(), shared={}, operators={})


def test_w101_table_interpolates_direct_jones() -> None:
    job = copy.deepcopy(BAR)
    job["solver"]["order"] = 2
    node = Node(
        "gr",
        "hpfem.grating",
        params={"job_data": job, "side": "reflected", "sigma_step": 0.02, "sigma_margin": 0.01},
    )
    # σ_x 0.04 … 0.08, far from the Rayleigh anomalies at ±0.25 (orders ±1 grazing in air and glass)
    sx = np.linspace(0.04, 0.08, 5)
    d = np.stack([sx, np.zeros_like(sx), np.sqrt(1 - sx**2)], axis=-1)
    ctx = _ctx()
    table = hpfem_grating.hpfem_table(ctx, node, LAMBDA_UM, d)
    assert len(table.sigma_x) == 4 and len(table.sigma_y) == 4
    assert any("Jones table from hp-FEM, 4×4 nodes" in line for line in ctx.log)
    handle = hpa.prepare(hpa.structure_from_job(hpfem_grating.job_of(node), LAMBDA_UM), [0.5])
    direct = jt.direct_jones(handle, LAMBDA_UM, 0.065, 0.003, "reflected")
    assert set(direct) == {0}  # Λ = 0.4 µm < λ / (1 + σ): only the specular order propagates in air
    for m in direct:
        got = jt.interpolate(table, m, LAMBDA_UM, 0.065, 0.003)
        assert np.max(np.abs(got - direct[m])) < 2e-4 * np.max(np.abs(direct[m])) + 1e-6
    # the same request again comes from the operator store, without solving
    again = hpfem_grating.hpfem_table(ctx, node, LAMBDA_UM, d)
    assert again is table


def test_w101_too_wide_a_bundle_stops_before_solving() -> None:
    node = Node(
        "gr", "hpfem.grating", params={"job_data": FLAT, "sigma_step": 0.01, "max_nodes": 16}
    )
    sx = np.linspace(0.0, 0.3, 7)
    d = np.stack([sx, np.zeros_like(sx), np.sqrt(1 - sx**2)], axis=-1)
    with pytest.raises(ExecutionError, match="more than 16"):
        hpfem_grating.hpfem_table(_ctx(), node, LAMBDA_UM, d)


def test_structure_ports_on_the_bar_and_refusals() -> None:
    s = hpa.structure_from_job(BAR, LAMBDA_UM)
    assert s.y_front == pytest.approx(148e-9) and s.y_back == pytest.approx(0.0)
    assert s.n_cover == 1.0 and s.n_substrate == N_SUB
    assert dict(s.options)["pml"] == {"top": pytest.approx(419e-9), "bottom": pytest.approx(419e-9)}
    h = hpa.prepare(s, [LAMBDA_UM])
    assert h.period_um == pytest.approx(0.4) and h.front_index == 1.0
    lossy = copy.deepcopy(FLAT)
    lossy["stack"]["incidence"] = {"n": [1.0, 0.1]}
    with pytest.raises(ValueError, match="lossless"):
        hpa.structure_from_job(lossy, LAMBDA_UM)
    with pytest.raises(ValueError, match="grating"):
        hpa.structure_from_job({**FLAT, "problem": "cavity"}, LAMBDA_UM)
    node = Node("gr", "hpfem.grating", params={})
    with pytest.raises(ExecutionError, match="needs a job"):
        hpfem_grating.job_of(node)
