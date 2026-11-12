"""Display blocks and progress (R21).

References: definitions and closed forms –
- W32 `display.image` shows exactly the camera image; `display.plot` of a sweep shows exactly the
  sweep's measure over its values.
- W33 `display.beam` attached to a collimated source (radius r, irradiance I₀): I₀ in the centre
  and the power π r² I₀ – as rays; as a field (Gaussian beam at its waist) I₀ at the peak.
- W34 progress: nondecreasing in [0, 1], one report per sweep value with the parameter in the
  text; the running blocks are reported.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

from rtt_fusion.graph import Graph, Node, builtin_registry, display, meta, optics
from rtt_fusion.graph.execute import WITH_FIELDS, Executor
from rtt_fusion.graph.imaging import Image

sys.path.insert(0, str(Path(__file__).parent))
from test_execute import WAVELENGTH_MM, michelson


def with_displays() -> Executor:
    ex = michelson()
    ex.evaluators.update(display.register(ex.registry))
    ex.evaluators.update(meta.register(ex.registry))
    g = ex.graph
    g.add(Node("show", "display.image", name="Camera"))
    g.connect(("cam", "image"), ("show", "in"))
    g.add(
        Node(
            "sweep",
            "meta.sweep",
            params={
                "parameter": "mB.point[0]",
                "collect": "cam.image",
                "measure": "centre",
                "values": [20.0 + k * WAVELENGTH_MM / 8 for k in range(5)],
            },
        )
    )
    g.add(Node("curve", "display.plot"))
    g.connect(("sweep", "results"), ("curve", "in"))
    return ex


def test_w32_image_and_plot_displays() -> None:
    ex = with_displays()
    view = ex.evaluate("show", "view")
    camera: Image = ex.evaluate("cam", "image")
    assert view.kind == "image" and view.title == "Camera"
    assert view.image is not None and np.array_equal(view.image.values, camera.values)
    plot = ex.evaluate("curve", "view")
    sweep = ex.evaluate("sweep", "results")
    assert plot.kind == "curves" and plot.x_label == "mB.point[0]"
    assert plot.curves[0].x == sweep.values and plot.curves[0].y == sweep.measures
    profile = display.plot_view(ex.context, Node("p", "display.plot"), "view", {"in": camera})
    assert len(profile.curves) == 2 and len(profile.curves[0].x) == camera.values.shape[1]


@pytest.mark.parametrize("transport", ["R", "H"])
def test_w33_beam_probe(transport: str) -> None:
    if transport == "H":
        pytest.importorskip("fieldtracing")
    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(display.register(reg))
    g = Graph()
    params: dict[str, object] = {"radius": 0.3, "rings": 20, "irradiance": 2.0}
    if transport == "H":
        params.update(gaussian_waist=0.12)
    g.add(Node("src", "optics.collimated_source", params=params))
    g.add(Node("probe", "display.beam", params={"pixels": 64}))
    g.connect(("src", "out"), ("probe", "in"), transport)
    view = Executor(g, reg, ev, available=WITH_FIELDS).evaluate("probe", "view")
    assert view.kind == "image" and view.image is not None
    v = view.image.values
    centre = v[31:33, 31:33].mean()
    u0, u1, _, _ = view.image.extent
    pixel = (u1 - u0) / 64
    r2 = 2 * (pixel / 2) ** 2  # the four central pixel centres
    expected = 2.0 * (np.exp(-2 * r2 / 0.12**2) if transport == "H" else 1.0)
    assert centre == pytest.approx(expected, rel=1e-2)  # I₀ (the Gaussian at its waist)
    if transport == "R":
        assert view.image.power == pytest.approx(np.pi * 0.3**2 * 2.0, rel=3e-2)


def test_w34_progress_and_running_blocks() -> None:
    ex = with_displays()
    parts: list[float] = []
    texts: list[str] = []
    nodes: list[str] = []

    def record(p: float, t: str) -> None:
        parts.append(p)
        texts.append(t)

    ex.context.progress = record
    ex.context.on_node = nodes.append
    ex.evaluate("curve", "view")
    assert parts and all(0.0 <= p <= 1.0 for p in parts)
    assert all(b >= a - 1e-12 for a, b in itertools.pairwise(parts))
    assert sum("mB.point[0] =" in t for t in texts) == 5
    assert "mB" in nodes and "cam" in nodes and nodes[-1] == "curve"
