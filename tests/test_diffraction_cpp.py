"""The C++ table (stage R12) in the grating block: `jones_table.prt_many` gives per ray what the
Python loop gives – the same rays dropped (outside the table, a Rayleigh anomaly over the
stencil, a node not propagating) and the same PRT matrices – and the grating block's output is
the same with and without it."""

from __future__ import annotations

import numpy as np
import pytest

from rtt_fusion.couplings import jones_table as jt

pytest.importorskip("rtt_fusion._diffraction")


def _table(reflected: bool) -> jt.JonesTable:
    rng = np.random.default_rng(3)
    orders = (-1, 0, 1)
    lam = np.linspace(0.55, 0.65, 5)
    sx = np.linspace(-0.4, 0.4, 17)
    sy = np.linspace(-0.3, 0.3, 13)
    jones = rng.normal(size=(3, 5, 17, 13, 2, 2)) + 1j * rng.normal(size=(3, 5, 17, 13, 2, 2))
    period = 1.2
    for o, m in enumerate(orders):  # NaN where the order does not propagate
        for a, la in enumerate(lam):
            rho = np.hypot(sx[:, None] + m * la / period, sy[None, :])
            jones[o, a][rho >= 1.0] = np.nan
    return jt.JonesTable(
        orders, "reflected" if reflected else "transmitted", lam, sx, sy, period, 1.0, 1.0, jones
    )


@pytest.mark.parametrize("reflected", [False, True])
@pytest.mark.parametrize("order", [-1, 0, 1])
def test_prt_many_equals_the_python_loop(reflected: bool, order: int) -> None:
    table = _table(reflected)
    rng = np.random.default_rng(order + 7)
    sigma = np.column_stack([rng.uniform(-0.45, 0.45, 400), rng.uniform(-0.35, 0.35, 400)])
    alive = rng.uniform(size=400) > 0.1
    lam = 0.6
    got = jt.prt_many(table, order, lam, sigma, alive)
    assert got is not None
    prt, ok = got
    for k in range(len(sigma)):
        sx, sy = sigma[k]
        inside = bool(
            table.sigma_x.min() - 1e-12 <= sx <= table.sigma_x.max() + 1e-12
            and table.sigma_y.min() - 1e-12 <= sy <= table.sigma_y.max() + 1e-12
        )
        want = alive[k] and inside and not jt.near_rayleigh_anomaly(table, lam, sx, sy)
        if want:
            jm = jt.interpolate(table, order, lam, sx, sy)
            want = bool(np.all(np.isfinite(jm)))
        assert ok[k] == want, k
        if want:
            p = jt.prt_matrix(
                jm, (sx, sy), (sx + order * lam / table.period_um, sy), 1.0, 1.0, reflected
            )
            assert np.allclose(prt[k], p, rtol=1e-12, atol=1e-12), k
    assert 0 < ok.sum() < alive.sum()  # some rays dropped near the anomalies, most kept


def test_without_the_extension_the_python_loop_is_used(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(jt, "USE_CPP", False)
    assert jt.prt_many(_table(False), 0, 0.6, np.zeros((1, 2)), np.ones(1, bool)) is None
