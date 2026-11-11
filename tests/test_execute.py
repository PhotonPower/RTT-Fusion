"""R15: execution without planner, coherence and the coherent image; pilot W9.

Expected values are analytic (two-beam interference of plane waves: fringe period λ / sin α on a
detector perpendicular to one beam, I = I₁ + I₂ + 2√(I₁I₂) cos Δφ, energy) or follow from the
definitions (coherence groups, cache, mute); never from the code under test. The Raytatouille
test uses the Michelson of Raytatouille's references (tests/data/michelson.rtt.json) and is
skipped without the module.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.optimize import curve_fit

from rtt_fusion._types import FloatArray
from rtt_fusion.graph import Graph, Node, builtin_registry, graph_from_raytatouille, optics
from rtt_fusion.graph.execute import ExecutionError, Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.model import light_in, light_out

WAVELENGTH_MM = 0.6328e-3
R = 2.0  # beam radius, mm
I_ARM = 0.25  # |r t|² of a 50/50 splitter passed twice


def michelson(
    theta: float = 0.0,
    delta_mm: float = 0.0,
    wavelengths: tuple[float, ...] = (0.6328,),
    fields: tuple[tuple[float, float], ...] = ((0.0, 0.0),),
) -> Executor:
    """Michelson with plane elements: arm A along +z, arm B along +x, mirror B tilted by θ about z
    (fringes along y on the camera), arm B longer by δ."""
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "center": [0, 0, -10],
                "direction": [0, 0, 1],
                "radius": R,
                "rings": 30,
                "wavelengths_um": list(wavelengths),
                "field_angles_deg": [list(f) for f in fields],
            },
        )
    )
    g.add(Node("bs", "optics.beam_splitter", params={"point": [0, 0, 0], "normal": [-1, 0, 1]}))
    g.add(Node("mA", "optics.plane_mirror", params={"point": [0, 0, 20.0], "normal": [0, 0, -1]}))
    g.add(
        Node(
            "mB",
            "optics.plane_mirror",
            params={"point": [20.0 + delta_mm, 0, 0], "normal": [-np.cos(theta), np.sin(theta), 0]},
        )
    )
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [-30, 0, 0],
                "normal": [1, 0, 0],
                "half_size": [2.5, 2.5],
                "pixels": [256, 256],
            },
        )
    )
    g.connect(("src", "out"), ("bs", "in"))
    g.connect(("bs", "t"), ("mA", "in"))
    g.connect(("mA", "out"), ("bs", "t_in"))
    g.connect(("bs", "r"), ("mB", "in"))
    g.connect(("mB", "out"), ("bs", "r_in"))
    g.connect(("bs", "merged"), ("cam", "in"))
    return Executor(g, reg, ev)


def image(ex: Executor) -> Image:
    result: Image = ex.evaluate("cam", "image")
    return result


def fringe_period(
    values: FloatArray, extent: tuple[float, float, float, float], radius: float
) -> float:
    """Period of the fringes along the axis with the larger variation, by a cosine fit through
    the centre (inside `radius`)."""
    ny, nx = values.shape
    u0, u1, v0, v1 = extent
    along_v = values[:, nx // 2].std() > values[ny // 2, :].std()
    line = values[:, nx // 2] if along_v else values[ny // 2, :]
    a, b = (v0, v1) if along_v else (u0, u1)
    n = len(line)
    x = a + (np.arange(n) + 0.5) * (b - a) / n
    inside = np.abs(x) < radius
    f = np.abs(np.fft.rfft((line - line[inside].mean()) * inside))
    f0 = (np.argmax(f[1:]) + 1) / (b - a)

    def model(x: FloatArray, c: float, d: float, fr: float, ph: float) -> FloatArray:
        return np.asarray(c + d * np.cos(2 * np.pi * fr * x + ph))

    p: Any = curve_fit(
        model, x[inside], line[inside], p0=[line[inside].mean(), line[inside].std() * 1.4, f0, 0.0]
    )[0]
    return float(1.0 / abs(p[2]))


# ---------------------------------------------------------------------------- pilot W9 (analytic)


def test_w9_fringe_period_and_contrast() -> None:
    theta = 1e-3
    im = image(michelson(theta=theta))
    period = fringe_period(im.values, im.extent, 1.5)
    assert period == pytest.approx(WAVELENGTH_MM / np.sin(2 * theta), rel=1e-6)
    centre = im.values[64:192, 64:192]
    contrast = (centre.max() - centre.min()) / (centre.max() + centre.min())
    assert contrast == pytest.approx(1.0, abs=2e-3)


@pytest.mark.parametrize(
    ("delta", "expected"),
    [(0.0, 4 * I_ARM), (WAVELENGTH_MM / 8, 2 * I_ARM), (WAVELENGTH_MM / 4, 0.0)],
)
def test_w9_path_difference(delta: float, expected: float) -> None:
    # Δφ = 2 k₀ δ: I = 2 I_arm (1 + cos Δφ)
    im = image(michelson(delta_mm=delta))
    assert im.values[100:156, 100:156].mean() == pytest.approx(expected, abs=1e-3)


def test_w9_energy() -> None:
    # quarter-wave difference: the interference term vanishes, ∫I = 2 · I_arm · π R²
    im = image(michelson(delta_mm=WAVELENGTH_MM / 8))
    assert im.power == pytest.approx(2 * I_ARM * np.pi * R * R, rel=5e-3)
    assert im.ray_power == pytest.approx(2 * I_ARM * np.pi * R * R, rel=5e-3)


# ---------------------------------------------------------------------------- coherence


def test_different_wavelengths_and_fields_add_incoherently() -> None:
    both = image(michelson(theta=1e-3, wavelengths=(0.6328, 0.5)))
    a = image(michelson(theta=1e-3, wavelengths=(0.6328,)))
    b = image(michelson(theta=1e-3, wavelengths=(0.5,)))
    assert len(both.groups) == 2
    np.testing.assert_allclose(both.values, a.values + b.values, atol=1e-12)
    two_fields = image(michelson(fields=((0.0, 0.0), (0.05, 0.0))))
    one = image(michelson(fields=((0.0, 0.0),)))
    other = image(michelson(fields=((0.05, 0.0),)))
    assert len(two_fields.groups) == 2
    np.testing.assert_allclose(two_fields.values, one.values + other.values, atol=1e-12)


# ---------------------------------------------------------------------------- execution


def test_cache_recomputes_only_downstream() -> None:
    ex = michelson()
    image(ex)
    first = dict(ex.calls)
    image(ex)
    assert ex.calls == first  # nothing recomputed
    ex.graph.nodes["mB"].params["point"] = [20.0 + WAVELENGTH_MM / 4, 0, 0]
    image(ex)
    assert ex.calls[("src", "out")] == first[("src", "out")]  # upstream: cached
    assert ex.calls[("mA", "out")] == first[("mA", "out")]  # the other arm: cached
    assert ex.calls[("mB", "out")] == first[("mB", "out")] + 1
    assert ex.calls[("cam", "image")] == first[("cam", "image")] + 1


def test_execution_follows_the_plan() -> None:
    # R17: the executor plans before running; with "auto" every light link is R (as in R15)
    ex = michelson()
    image(ex)
    assert ex.plan is not None and ex.plan.ok
    assert set(ex.plan.representations.values()) == {"R"}
    assert len(ex.plan.representations) == len(ex.graph.links)
    assert not ex.plan.conversions
    assert len([p for p in ex.plan.paths if p.output == ("cam", "image")]) == 2  # two arms


def test_muted_mirror_passes_light_through() -> None:
    # muted reference mirror: its light runs on and misses the camera – only arm A remains
    ex = michelson()
    ex.graph.set_muted("mB")
    assert image(ex).values[100:156, 100:156].mean() == pytest.approx(I_ARM, abs=1e-3)


def test_execution_errors() -> None:
    ex = michelson()
    ex.graph.set_transport(ex.graph.links[0], "H")  # the optics sockets take R only
    with pytest.raises(ExecutionError, match="graph.transport_unsupported"):
        image(ex)
    # sockets that take H, but R15 runs R only
    g = Graph()
    g.add(Node("a", "rtt.section", sockets=(light_out("o", "RH"),)))
    g.add(
        Node(
            "b",
            "rtt.section",
            sockets=(light_in("i", "RH"), light_out("o", "RH")),
            through=(("o", ("i",)),),
        )
    )
    g.connect(("a", "o"), ("b", "i"), "H")
    ex2 = Executor(g, builtin_registry(), {"rtt.section": lambda ctx, node, socket, inputs: []})
    with pytest.raises(ExecutionError, match="exec.transport_unsupported"):
        ex2.evaluate("b", "o")
    ex = michelson()
    del ex.evaluators["optics.plane_mirror"]
    with pytest.raises(ExecutionError, match="exec.no_evaluator"):
        image(ex)
    ex = michelson()
    ex.graph.disconnect(ex.graph.links[0])
    ex.graph.connect(("cam", "image"), ("bs", "in"))  # data into light
    with pytest.raises(ExecutionError, match="graph.kind_mismatch"):
        image(ex)


# ---------------------------------------------------------------------------- Raytatouille (W9)


def test_w9_michelson_through_raytatouille() -> None:
    rt = pytest.importorskip("raytatouille")
    from rtt_fusion.graph import backend_rtt

    d = json.loads(
        (Path(__file__).parent / "data" / "michelson.rtt.json").read_text(encoding="utf-8")
    )
    theta = 1e-3
    test_mirror = next(c for c in d["root"]["children"] if c["name"] == "test mirror")
    test_mirror["pose"]["rotation_deg"][0] = -90.0 + np.degrees(theta)
    g = graph_from_raytatouille(d)
    assert {n.type_id for n in g.nodes.values()} >= {"rtt.splitter", "rtt.mirror", "rtt.detector"}
    system = rt.System.from_json(json.dumps(d))
    ex = Executor(
        g,
        builtin_registry(),
        backend_rtt.EVALUATORS,
        {"system": system, "rings": 30, "pixels": [320, 240]},
    )
    im: Image = ex.evaluate("camera", "image")
    assert len(im.groups) == 1  # one field, one wavelength: both arms in one coherence group
    period = fringe_period(im.values, im.extent, 2.0)
    assert period == pytest.approx(WAVELENGTH_MM / np.sin(2 * theta), rel=1e-2)
    ny, nx = im.values.shape
    centre = im.values[ny // 2 - 40 : ny // 2 + 40, nx // 2 - 40 : nx // 2 + 40]
    assert (centre.max() - centre.min()) / (centre.max() + centre.min()) > 0.99
