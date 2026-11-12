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

import numpy as np
import pytest

from rtt_fusion.graph import Graph, Node, builtin_registry, grating, optics
from rtt_fusion.graph.execute import Executor
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
            "fusion.grating",
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


@pytest.mark.parametrize("alpha_deg", [0.0, 3.0])
def test_w26_block_equals_direct_fouriergrating(
    config: tuple[dict[str, object], str], alpha_deg: float
) -> None:
    cfg, path = config
    handle = fga.prepare(fga.Structure(cfg), (LAMBDA_UM,))
    sx = float(np.sin(np.radians(alpha_deg)))
    for side in ("reflected", "transmitted"):
        powers, log = order_powers(path, side, alpha_deg)
        assert any("Jones table from FourierGrating" in line for line in log)
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
