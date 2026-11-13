"""Grating as the sample of a Raytatouille system (R19b): the path is split at the sample.

Test system `tests/data/linnik_mini.rtt.json`: one Linnik arm with ideal lenses (Köhler
illumination, objective f = 2 mm, tube lens f = 20 mm, magnification M = 10, sample in the front
focal plane: a 4f image of the sample on the camera). The sample mirror is replaced by
`fg.grating` with tables built here from closed forms:

- W22 a "mirror table" (order 0, J = −I): the split computation (Raytatouille up to the sample,
  the block, Raytatouille on) must give the image of the whole path with the mirror – for rays
  and for fields.
- W23 a thin sinusoidal phase grating, Raman–Nath: e^{iφ₀ sin Kx} = Σ J_m(φ₀) e^{imKx}; keeping
  the orders −1, 0, 1 the image is I/I_mirror = J₀² + 4 J₁² sin²(K′x) with K′ = 2π/(MΛ): mean
  J₀² + 2 J₁², extremes J₀² and J₀² + 4 J₁², fringes of period MΛ/2 and no component at 1/(MΛ).
  Λ = 10 µm keeps the orders paraxial (sin θ = 0.06), where the ideal lenses image without
  aberration (they are collinear maps, not aplanatic).

The tables are power-normalised for the full field (stage R6): the component in the plane of
incidence is scaled by cos θ_m / sqrt(cos² θ_m + σ_m²) to take out the E_z that the lift adds.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pytest
from scipy.special import jv

from rtt_fusion.graph import (
    Graph,
    Node,
    builtin_registry,
    graph_from_raytatouille,
    grating,
    replace_block,
)
from rtt_fusion.graph.execute import WITH_FIELDS, Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.results import FIELD_REFERENCE, to_json

rt = pytest.importorskip("raytatouille")

SYSTEM = Path(__file__).parent / "data" / "linnik_mini.rtt.json"
PHI0, PERIOD_UM, MAGNIFICATION = 1.0, 10.0, 10.0


def write_table(path: Path, amplitudes: dict[int, float], period_um: float) -> str:
    lam = np.linspace(0.59, 0.61, 4)
    sx = np.linspace(-0.03, 0.03, 7)
    sy = np.linspace(-0.03, 0.03, 7)
    orders = sorted(amplitudes)
    jones = np.zeros((len(orders), len(lam), len(sx), len(sy), 2, 2), dtype=complex)
    for o, m in enumerate(orders):
        for a, wl in enumerate(lam):
            for b, x in enumerate(sx):
                for c, y in enumerate(sy):
                    smx = x + m * wl / period_um
                    cz = float(np.sqrt(max(1.0 - float(smx) ** 2 - float(y) ** 2, 0.0)))
                    jones[o, a, b, c] = amplitudes[m] * np.diag([cz / np.hypot(cz, smx), 1.0])
    data = {
        "parameters": {
            "orders": orders,
            "side": "reflected",
            "period_um": period_um,
            "n_in": 1.0,
            "n_out": 1.0,
        },
        "inputs": {"wavelengths_um": lam, "sigma_x": sx, "sigma_y": sy, "jones": jones},
    }
    path.write_text(to_json(FIELD_REFERENCE, data), encoding="utf-8")
    return str(path)


def image(
    system_json: dict[str, object], table: str | None, transport: str = "auto", pixels: int = 512
) -> Image:
    from rtt_fusion.graph import backend_rtt

    system = rt.System.from_json(json.dumps(system_json))
    g: Graph = graph_from_raytatouille(system_json)
    reg = builtin_registry()
    ev = dict(backend_rtt.EVALUATORS)
    ev.update(grating.register(reg))
    if table is not None:
        replace_block(g, reg, "sample", Node("sample", "fg.grating", params={"table": table}))
    for link in g.links_to("camera"):
        g.set_transport(link, transport)
    ex = Executor(
        g,
        reg,
        ev,
        {"system": system, "rings": 15, "pixels": [pixels, 8]},
        available=WITH_FIELDS,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result: Image = ex.evaluate("camera", "image")
    return result


def load_system() -> dict[str, object]:
    d: dict[str, object] = json.loads(SYSTEM.read_text(encoding="utf-8"))
    return d


def test_w22_mirror_table_equals_mirror(tmp_path: Path) -> None:
    table = write_table(tmp_path / "mirror.json", {0: -1.0}, PERIOD_UM)
    d = load_system()
    whole = image(d, None)
    split = image(d, table)
    assert len(split.groups) == 2  # two input polarisations of the unpolarised source
    assert np.abs(split.values - whole.values).max() < 1e-9 * whole.values.max()


def test_w22_mirror_table_equals_mirror_fields(tmp_path: Path) -> None:
    pytest.importorskip("fieldtracing")
    table = write_table(tmp_path / "mirror.json", {0: -1.0}, PERIOD_UM)
    d = load_system()
    d["aperture"] = {"type": "epd", "value": 0.2}  # a small beam keeps the FFT small
    cam = d["root"]["children"][-1]["surfaces"][0]["aperture"]  # type: ignore[index]
    cam["half_width_x"] = cam["half_width_y"] = 0.25
    whole = image(d, None, "H", pixels=64)
    split = image(d, table, "H", pixels=64)
    assert np.abs(split.values - whole.values).max() < 1e-3 * whole.values.max()


def test_w23_phase_grating_imaged_by_the_linnik(tmp_path: Path) -> None:
    amplitudes = {m: -float(jv(m, PHI0)) for m in (-1, 0, 1)}
    table = write_table(tmp_path / "phase.json", amplitudes, PERIOD_UM)
    d = load_system()
    mirror = image(d, None)
    grated = image(d, table)
    row = grated.values[4] / mirror.values[4].max()
    x = np.linspace(grated.extent[0], grated.extent[1], len(row))
    inside = np.abs(x) < 0.8  # well inside the image of the illuminated field (radius 1 mm)
    j0, j1 = float(jv(0, PHI0)), float(jv(1, PHI0))
    assert row[inside].mean() == pytest.approx(j0**2 + 2 * j1**2, abs=5e-3)
    assert row[inside].max() == pytest.approx(j0**2 + 4 * j1**2, abs=5e-3)
    assert row[inside].min() == pytest.approx(j0**2, abs=5e-3)
    width = grated.extent[1] - grated.extent[0]
    n = 1 << 15
    spectrum = np.abs(np.fft.rfft((row - row[inside].mean()) * inside * np.hanning(len(row)), n=n))
    freq = np.fft.rfftfreq(n, width / len(row))  # 1/mm
    k = 1.0 / (MAGNIFICATION * PERIOD_UM * 1e-3)
    near = np.abs(freq - 2 * k) < 0.3 * k
    peak = freq[np.argmax(np.where(near, spectrum, 0.0))]
    assert peak == pytest.approx(2 * k, rel=1e-2)
    first = spectrum[np.abs(freq - k) < 0.3 * k].max()
    assert first < 0.02 * spectrum[near].max()
