"""Grating block in the working graph (R19): a periodic sample with diffraction orders as channels.

The interaction comes from the stored Jones table of stage R6 (reference/jones_table_lamellar,
FourierGrating RCWA of a lamellar grating, Λ = 0.7 µm, n = 2.1 / air, orders −1, 0, 1,
transmitted); no solver is needed. References, independent of the table values:
- W19 grating equation: sin θ_m = sin α + m·λ/Λ in the plane of incidence; at normal incidence
  the symmetric lamellar profile gives equal power in the orders +1 and −1.
- W20 interference of the orders behind the grating: the intensity is periodic with the grating
  period – components at 1/Λ (orders 0 with ±1) and 2/Λ (−1 with +1), nothing in between.
- W21 working graph: mirror -> grating with `replace_node` keeps the links and the graph runs;
  a field link (H) into the grating gives the same orders as rays.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from rtt_fusion._types import FloatArray
from rtt_fusion.graph import Graph, Node, builtin_registry, grating, optics
from rtt_fusion.graph.execute import WITH_FIELDS, Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.light import Plane, RayBundle
from rtt_fusion.graph.topology import Model

TABLE = str(Path(__file__).parent.parent / "reference" / "jones_table_lamellar.result.json")
LAMBDA_UM, PERIOD_UM = 0.6, 0.7


def setup(alpha_deg: float = 0.0, transport: str = "auto") -> tuple[Graph, Executor]:
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
                "rings": 12,
                "center": [0, 0, -0.1],
                "wavelengths_um": [LAMBDA_UM],
                "field_angles_deg": [[alpha_deg, 0.0]],
            },
        )
    )
    g.add(Node("gr", "fg.grating", params={"table": TABLE}))
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, 0.002],
                "normal": [0, 0, -1],
                "half_size": [0.0035, 0.0035],
                "pixels": [400, 4],
            },
        )
    )
    g.connect(("src", "out"), ("gr", "in"), transport)
    g.connect(("gr", "out"), ("cam", "in"))
    return g, Executor(g, reg, ev, available=WITH_FIELDS)


def orders(ex: Executor) -> dict[int, RayBundle]:
    out: list[RayBundle] = ex.evaluate("gr", "out")
    return {int(b.path[-1].split("[")[1].rstrip("]")): b for b in out}


def mean_direction(b: RayBundle) -> FloatArray:
    alive = np.asarray(b.alive)
    return np.asarray(np.asarray(b.direction)[alive].mean(axis=0), dtype=np.float64)


@pytest.mark.parametrize("alpha_deg", [0.0, 3.0, -2.0])
def test_w19_grating_equation(alpha_deg: float) -> None:
    _, ex = setup(alpha_deg)
    found = orders(ex)
    assert sorted(found) == [-1, 0, 1]
    # the field angle tilts the beam about the source frame's v axis, in the plane of its u axis
    u = Plane.make([0, 0, 0], [0, 0, 1]).u
    d0 = np.tan(np.radians(alpha_deg)) * u + np.array([0.0, 0.0, 1.0])
    s_in = float((d0 / np.linalg.norm(d0))[0])
    for m, b in found.items():
        d = mean_direction(b)
        assert d[0] == pytest.approx(s_in + m * LAMBDA_UM / PERIOD_UM, abs=1e-9)
        assert d[1] == pytest.approx(0.0, abs=1e-12)
        assert np.ptp(np.asarray(b.direction)[np.asarray(b.alive)], axis=0).max() < 1e-12
    if alpha_deg == 0.0:
        power = {m: np.mean(np.abs(np.asarray(b.amplitude)) ** 2) for m, b in found.items()}
        assert power[1] == pytest.approx(power[-1], rel=1e-9)
        assert 0 < power[0] + power[1] + power[-1] < 1  # the rest is reflected


def test_w20_order_interference_has_the_grating_period() -> None:
    _, ex = setup()
    image: Image = ex.evaluate("cam", "image")
    row = image.values[2]
    width = image.extent[1] - image.extent[0]
    n = 1 << 14
    spectrum = np.abs(np.fft.rfft((row - row.mean()) * np.hanning(len(row)), n=n))
    freq = np.fft.rfftfreq(n, width / len(row))  # 1/mm
    k1, k2 = 1e3 / PERIOD_UM, 2e3 / PERIOD_UM
    peak1 = freq[np.argmax(np.where(np.abs(freq - k1) < 0.3 * k1, spectrum, 0))]
    peak2 = freq[np.argmax(np.where(np.abs(freq - k2) < 0.3 * k1, spectrum, 0))]
    assert peak1 == pytest.approx(k1, rel=1e-2)
    assert peak2 == pytest.approx(k2, rel=1e-2)
    # nothing else: away from the main lobes of the two peaks (10 periods in the window, Hann
    # window: lobes of ±0.3/Λ) the spectrum is below 5 % of the peaks
    away = (freq > 0.3 * k1) & (freq < 3.0 * k1)
    away &= (np.abs(freq - k1) > 0.3 * k1) & (np.abs(freq - k2) > 0.3 * k1)
    peaks = max(
        spectrum[np.abs(freq - k1) < 0.1 * k1].max(), spectrum[np.abs(freq - k2) < 0.1 * k1].max()
    )
    assert spectrum[away].max() < 0.05 * peaks


def test_w21_mirror_replaced_by_grating() -> None:
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
                "rings": 8,
                "center": [0, 0, -0.1],
                "wavelengths_um": [LAMBDA_UM],
            },
        )
    )
    g.add(Node("sample", "optics.plane_mirror", params={"point": [0, 0, 0], "normal": [0, 0, -1]}))
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, 0.002],
                "normal": [0, 0, -1],
                "half_size": [0.003, 0.003],
                "pixels": [32, 32],
            },
        )
    )
    g.connect(("src", "out"), ("sample", "in"))
    g.connect(("sample", "out"), ("cam", "in"))
    links = [(lk.source, lk.target) for lk in g.links]
    g.replace_node("sample", Node("sample", "fg.grating", params={"table": TABLE}))
    assert [(lk.source, lk.target) for lk in g.links] == links
    image: Image = Executor(g, reg, ev).evaluate("cam", "image")
    assert image.values.max() > 0
    # a field link into the grating: converted to rays at its port, same orders
    _, rays = setup()
    _, fields = setup(transport="H")
    a, b = orders(rays), orders(fields)
    assert sorted(a) == sorted(b)
    for m in a:
        assert np.allclose(mean_direction(a[m]), mean_direction(b[m]), atol=1e-9)


def test_w72_old_grating_type_loads_as_fg_grating() -> None:
    """CAD ADR 0010 W4: `fusion.grating` (files before the rename) loads as `fg.grating` – graph,
    model and the variants of a group – and is saved with the new id."""
    from rtt_fusion.graph import io
    from rtt_fusion.graph.workflow import registry as workflow_registry

    old = Node("gr", "fusion.grating", params={"table": TABLE})
    graph_dict = io.to_dict(Graph())
    graph_dict["nodes"] = [{**io._node_to_dict(old)}]
    graph = io.from_dict(graph_dict)
    assert graph.nodes["gr"].type_id == "fg.grating"
    assert '"fusion.grating"' not in io.dumps(graph) and '"fg.grating"' in io.dumps(graph)

    reg, _ = workflow_registry()
    assert reg.get("fg.grating") is not None and reg.get("fusion.grating") is None
    model_dict = io.model_to_dict(io.model_from_dict(io.model_to_dict(_model_with(old)), reg))
    assert [n["type"] for n in model_dict["nodes"]] == ["fg.grating"]


def _model_with(node: Node) -> Model:  # a model holding one node, with its type id as given
    m = Model()
    m.add(node)
    return m
