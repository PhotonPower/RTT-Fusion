"""Grating block with FourierGrating in the loop (R19d): the block builds its Jones table itself.

Grating: the lamellar 1D grating of reference/jones_table_lamellar (n = 2.1 / air, Λ = 0.7 µm,
RCWA), λ = 0.6 µm. References, independent of the block:
- W26 the order powers of the block equal FourierGrating's direct solution at the same incidence
  (`jones_table.direct_jones`, 2 solves per point – no table, no interpolation).
- W27 energy: the grating is lossless, so the reflected and transmitted orders carry all the
  incident power, Σ R_m + Σ T_m = 1.
Needs the solver (skipped without it); 1D, a few seconds.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import Graph, Node, builtin_registry, grating, optics
from rtt_fusion.graph.execute import Context, Executor
from rtt_fusion.results import load_json

pytest.importorskip("fouriergrating")

from rtt_fusion.adapters import fouriergrating as fga  # noqa: E402
from rtt_fusion.couplings import jones_table as jt  # noqa: E402

REFERENCE = Path(__file__).parent.parent / "reference" / "jones_table_lamellar.result.json"
LAMBDA_UM = 0.6


@pytest.fixture(scope="module")
def config(tmp_path_factory: pytest.TempPathFactory) -> tuple[dict[str, object], str]:
    _, d = load_json(REFERENCE.read_text(encoding="utf-8"))
    cfg: dict[str, object] = json.loads(d["parameters"]["fouriergrating_config"])
    path = tmp_path_factory.mktemp("fg") / "lamellar.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    return cfg, str(path)


def order_powers(path: str, side: str, alpha_deg: float) -> tuple[dict[int, float], list[str]]:
    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(grating.register(reg))
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": 0.01,
                "rings": 6,
                "center": [0, 0, -0.1],
                "wavelengths_um": [LAMBDA_UM],
                "field_angles_deg": [[alpha_deg, 0.0]],
                "polarisation": [0, 1, 0],  # TE: along the grooves
            },
        )
    )
    g.add(
        Node(
            "gr",
            "fg.grating",
            params={"fouriergrating": path, "side": side, "orders": [-1, 0, 1]},
        )
    )
    g.connect(("src", "out"), ("gr", "in"))
    ex = Executor(g, reg, ev)
    out = ex.evaluate("gr", "out")
    powers: dict[int, float] = {}
    for b in out:
        m = int(b.path[-1].split("[")[1].rstrip("]"))
        alive = np.asarray(b.alive)
        field = np.asarray(b.amplitude)[alive][:, None] * np.einsum(
            "nij,j->ni", np.asarray(b.prt)[alive], np.asarray(b.jones)
        )
        powers[m] = float(np.mean(np.sum(np.abs(field) ** 2, axis=1)))
    return powers, ex.context.log


# 7.87°: σ = 0.137, 0.006 below the Rayleigh anomaly of order +1 (σ = 1 − λ/Λ = 0.143), closer
# than a table's margin – a plane wave is computed at its own σ (no table), so nothing is dropped
@pytest.mark.parametrize("alpha_deg", [0.0, 3.0, 7.874])
def test_w26_block_equals_direct_fouriergrating(
    config: tuple[dict[str, object], str], alpha_deg: float
) -> None:
    cfg, path = config
    handle = fga.prepare(fga.Structure(cfg), (LAMBDA_UM,))
    sx = float(np.sin(np.radians(alpha_deg)))
    for side in ("reflected", "transmitted"):
        powers, log = order_powers(path, side, alpha_deg)
        assert any("from FourierGrating at one incidence" in line for line in log), log
        assert not any("outside the table" in line for line in log), log
        direct = jt.direct_jones(handle, LAMBDA_UM, sx, 0.0, side)
        for m, jones in direct.items():
            p = jt.prt_matrix(
                jones, (sx, 0.0), (sx + m * LAMBDA_UM / 0.7, 0.0), 1.0, 1.0, side == "reflected"
            )
            expected = float(np.sum(np.abs(p @ np.array([0.0, 1.0, 0.0])) ** 2))
            assert powers[m] == pytest.approx(expected, rel=1e-4, abs=1e-6)


@pytest.mark.parametrize("alpha_deg", [0.0, 3.0])
def test_w27_lossless_grating_conserves_energy(
    config: tuple[dict[str, object], str], alpha_deg: float
) -> None:
    _, path = config
    reflected, _ = order_powers(path, "reflected", alpha_deg)
    transmitted, _ = order_powers(path, "transmitted", alpha_deg)
    assert sum(reflected.values()) + sum(transmitted.values()) == pytest.approx(1.0, abs=1e-5)


def _fanned(path: str, step: float, spread: float, ctx: Context | None = None) -> tuple[int, Any]:
    """The grating with a bundle whose directions fan out by ±spread in σ_x around 2°: a table
    is built around them; (its σ_x nodes, the outgoing bundles)."""
    reg = builtin_registry()
    ev = optics.register(reg)
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
                "field_angles_deg": [[2.0, 0.0]],
                "polarisation": [0, 1, 0],
            },
        )
    )
    b = Executor(g, reg, ev).evaluate("src", "out")[0]
    d = np.asarray(b.direction, dtype=float).copy()
    d[:, 0] += np.linspace(-spread, spread, len(d))
    d[:, 2] = np.sqrt(1.0 - d[:, 0] ** 2 - d[:, 1] ** 2)
    params = {"fouriergrating": path, "side": "reflected", "orders": [-1, 0, 1], "sigma_step": step}
    ctx = ctx if ctx is not None else Context(Graph())
    out = grating.grating(
        ctx, Node("gr", "fg.grating", params=params), "out", {"in": [b.evolve(direction=d)]}
    )
    line = next(x for x in ctx.log if "Jones table from FourierGrating" in x)
    return int(line.split(", ")[1].split()[0].split("×")[1]), out  # 4×nx×ny nodes


def test_sigma_step_sets_the_table_grid(config: tuple[dict[str, object], str]) -> None:
    """The ray path's table resolution is a block parameter: a finer σ step gives more nodes and
    the same orders (the interpolation is already exact to 1e-4 at the default step)."""
    _, path = config
    coarse, a = _fanned(path, 0.01, 0.005)
    fine, b = _fanned(path, 0.004, 0.005)
    assert fine > coarse
    for x, y in zip(a, b, strict=True):
        assert np.allclose(np.asarray(x.amplitude), np.asarray(y.amplitude), rtol=1e-4, atol=1e-8)


def test_tables_outlive_a_request(config: tuple[dict[str, object], str]) -> None:
    """A built table is kept in the operator store, which outlives a provider's request (the bild
    and the pupil branch of a microscope, a second run): the second request builds nothing."""
    _, path = config
    store: dict[str, Any] = {}
    _fanned(path, 0.01, 0.005, Context(Graph(), operators=store))
    assert any(k.startswith("jones|") for k in store)
    again = Context(Graph(), operators=store)
    with pytest.raises(StopIteration):  # no "Jones table from FourierGrating" this time
        _fanned(path, 0.01, 0.005, again)


def test_point_table_is_the_direct_solution(config: tuple[dict[str, object], str]) -> None:
    cfg, _ = config
    handle = fga.prepare(fga.Structure(cfg), (LAMBDA_UM,))
    t = jt.point_table(handle, LAMBDA_UM, 0.137, 0.02, (-1, 0, 1), "reflected")
    direct = jt.direct_jones(handle, LAMBDA_UM, 0.137, 0.02, "reflected")
    for m in (-1, 0, 1):
        assert np.allclose(jt.interpolate(t, m, LAMBDA_UM, 0.137, 0.02), direct[m], rtol=1e-12)
    assert not jt.near_rayleigh_anomaly(t, LAMBDA_UM, 0.137, 0.02)


def _source_bundle(alpha_deg: float = 0.0) -> Any:
    reg = builtin_registry()
    ev = optics.register(reg)
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
                "polarisation": [0, 1, 0],
            },
        )
    )
    return Executor(g, reg, ev).evaluate("src", "out")[0]


def _air_to_glass(glass: Any) -> dict[str, Any]:
    """A shallow lamellar grating, air above, glass below (Λ = 1.7 µm: no order grazes at normal
    incidence in air or glass)."""
    return {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 8, "Nz": 0},
        "period_um": {"x": 1.7, "z": 0.0},
        "materials": {"air": {"n": 1.0}, "glass": {"n": glass}},
        "superstrate": "air",
        "substrate": "glass",
        "layers": [
            {
                "type": "lamellar1d",
                "thickness_um": 0.1,
                "breakpoints_um": [0.0, 0.85, 1.7],
                "materials": ["glass", "air"],
            }
        ],
    }


def test_orders_carry_the_index_of_their_medium() -> None:
    """RTT-Fusion#121: a transmitted order runs in the substrate – its bundle carries the
    substrate's index (as fg.stack), a reflected one the incident medium's."""
    b = _source_bundle()
    for side, n in (("transmitted", 1.5), ("reflected", 1.0)):
        node = Node(
            "gr",
            "fg.grating",
            params={"config": _air_to_glass(1.5), "side": side, "orders": [-1, 0, 1]},
        )
        ctx = Context(Graph())
        out = grating.grating(ctx, node, "out", {"in": [b]})
        assert {bo.path[-1] for bo in out} == {"gr[-1]", "gr[0]", "gr[1]"}
        for bo in out:
            assert float(bo.index) == n
            m = int(bo.path[-1].split("[")[1].rstrip("]"))
            d = np.asarray(bo.direction)[np.asarray(bo.alive)]
            # the grating equation in the order's medium: n·sin θ_m = m λ/Λ
            assert np.allclose(n * d[:, 0], m * LAMBDA_UM / 1.7, atol=1e-9)
        assert not any("rays arrive in index" in line for line in ctx.log)


def test_rays_in_another_index_than_the_table_expects_are_logged() -> None:
    b = _source_bundle().evolve(index=1.33)
    node = Node("gr", "fg.grating", params={"config": _air_to_glass(1.5), "side": "transmitted"})
    ctx = Context(Graph())
    grating.grating(ctx, node, "out", {"in": [b]})
    assert any("rays arrive in index 1.33, the grating expects 1.0" in line for line in ctx.log)


def test_a_bundle_without_live_rays_passes_nothing() -> None:
    """RTT-Fusion#122: the blocked group behind an ideal polariser – no table, no crash."""
    b = _source_bundle()
    dead = b.evolve(alive=np.zeros(b.size, dtype=bool))
    node = Node("gr", "fg.grating", params={"config": _air_to_glass(1.5), "side": "transmitted"})
    ctx = Context(Graph())
    assert grating.grating(ctx, node, "out", {"in": [dead]}) == []
    assert not any("Jones" in line for line in ctx.log)
    # beside a live bundle: only the live one is diffracted
    out = grating.grating(ctx, node, "out", {"in": [dead, b]})
    assert out and all(np.any(bo.alive) for bo in out)


def test_a_dispersive_substrate_is_evaluated_at_the_wavelength() -> None:
    """RTT-Fusion#123: {"n": {"sellmeier": …}} at a port, as FourierGrating computes it."""
    fg = pytest.importorskip("fouriergrating")
    if not hasattr(fg, "material_eps"):
        cfg = _air_to_glass({"sellmeier": {"A": 1.0, "B": [1.25], "C_um2": [0.01]}})
        with pytest.raises(ValueError, match="dispersive"):
            fga.prepare(fga.Structure(cfg), (LAMBDA_UM,))
        return
    sellmeier = {
        "A": 1.0,
        "B": [0.6961663, 0.4079426, 0.8974794],
        "C_um2": [0.00467914826, 0.0135120631, 97.9340025],
    }  # fused silica
    cfg = _air_to_glass({"sellmeier": sellmeier})
    lam2 = LAMBDA_UM**2
    n = float(
        np.sqrt(
            1.0
            + sum(
                bk * lam2 / (lam2 - ck)
                for bk, ck in zip(sellmeier["B"], sellmeier["C_um2"], strict=True)
            )
        )
    )
    handle = fga.prepare(fga.Structure(cfg), (LAMBDA_UM,))
    assert handle.back_index.real == pytest.approx(n, abs=1e-9)
    assert handle.back_index.imag == 0.0
    b = _source_bundle()
    node = Node("gr", "fg.grating", params={"config": cfg, "side": "transmitted", "orders": [0]})
    out = grating.grating(Context(Graph()), node, "out", {"in": [b]})
    assert len(out) == 1 and float(out[0].index) == pytest.approx(n, abs=1e-9)
