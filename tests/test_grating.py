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
    pytest.importorskip("fieldtracing")  # H needs Fieldtracing (requirements-modules.txt)
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


# ------------------------------------------------- the back port (#120): gratings in transmission


def _transmission_chain(back: bool) -> Graph:
    """Source → G1 → G2 → camera, both gratings in transmission: through their back port (two
    neighbours each, #120) or, as before, through `out` with side = transmitted."""
    g = Graph()
    g.add(Node("src", "optics.collimated_source", params={"radius": 0.01, "rings": 2}))
    for name, z in (("G1", 1.0), ("G2", 2.0)):
        g.add(
            Node(
                name,
                "fg.grating",
                params={
                    "table": TABLE,
                    "side": "reflected" if back else "transmitted",
                    "point": [0, 0, z],
                    "normal": [0, 0, -1],  # towards the light (front)
                },
            )
        )
    g.add(Node("cam", "optics.detector", params={"point": [0, 0, 3.0], "normal": [0, 0, -1]}))
    out = "back_out" if back else "out"
    g.connect(("src", "out"), ("G1", "in"))
    g.connect(("G1", out), ("G2", "in"))
    g.connect(("G2", out), ("cam", "in"))
    return g


def test_transmission_gratings_are_steps_of_a_branch() -> None:
    """Version 2 (the model CAD writes): a grating with a front and a back port is a step between
    two neighbours; with its single front port it is not (a port has one connection)."""
    from rtt_fusion.graph.topology import model_from_graph, unroll

    reg = builtin_registry()
    optics.register(reg)
    grating.register(reg)
    m = model_from_graph(_transmission_chain(back=True), reg)
    assert len(m.connections) == 3
    steps = [(s.node, s.port_in, s.port_out) for s in m.branches[0].steps]
    assert ("G1", "front", "back") in steps and ("G2", "front", "back") in steps
    u = unroll(m, reg)
    assert {(lk.source, lk.target) for lk in u.links} >= {
        (("G1", "back_out"), ("G2", "in")),
        (("G2", "back_out"), ("cam", "in")),
    }
    # the version 1 form (out → in on "transmitted") is read through the back port as well
    assert len(model_from_graph(_transmission_chain(back=False), reg).connections) == 3


def test_the_back_port_needs_transmitted_orders(monkeypatch: pytest.MonkeyPatch) -> None:
    """A table of reflected orders cannot feed the back port: a clear message, not reflected
    light leaving at the back."""
    from dataclasses import replace

    from rtt_fusion.graph import grating as module
    from rtt_fusion.graph.execute import ExecutionError

    held = module._table
    monkeypatch.setattr(
        module, "_table", lambda ctx, node: replace(held(ctx, node), side="reflected")
    )
    _, ex = setup()
    assert ex.evaluate("gr", "out")  # the front port takes it
    with pytest.raises(ExecutionError, match="back port needs a table of the transmitted"):
        ex.evaluate("gr", "back_out")


def test_port_frames_face_both_ways() -> None:
    frames = grating.port_frames(
        Node("gr", "fg.grating", params={"point": [0, 0, 1], "normal": [0, 0, 2]})
    )
    assert frames["front"].normal == (0.0, 0.0, 1.0) and frames["back"].normal == (0.0, 0.0, -1.0)
    assert frames["front"].origin == frames["back"].origin == (0.0, 0.0, 1.0)


def test_a_version_1_chain_of_transmission_gratings_becomes_a_model() -> None:
    """A version 1 graph (as scripts write it) chains transmission gratings out → in; as a model
    the light goes on through the back port (#120), the branch steps front → back, and the model
    unrolls to the same computation (the same blocks and links up to the socket name)."""
    from rtt_fusion.graph.topology import model_from_graph, unroll

    reg = builtin_registry()
    optics.register(reg)
    grating.register(reg)
    g = _transmission_chain(back=False)
    m = model_from_graph(g, reg)
    steps = [(s.node, s.port_in, s.port_out) for s in m.branches[0].steps]
    assert ("G1", "front", "back") in steps and ("G2", "front", "back") in steps
    u = unroll(m, reg)
    assert {(lk.source, lk.target) for lk in u.links} >= {
        (("G1", "back_out"), ("G2", "in")),
        (("G2", "back_out"), ("cam", "in")),
    }
    # a grating on "reflected" that sends light on to another block stays refused
    for name in ("G1", "G2"):
        g.nodes[name].params["side"] = "reflected"
    with pytest.raises(ValueError, match="a port has one connection"):
        model_from_graph(g, reg)


def test_a_model_document_of_gratings_reads_back() -> None:
    """The workflow document of a fresh model (no ports stored on the nodes yet, as a script or
    the conversion of an old graph makes it) keeps the gratings' declared ports front and back: the
    file reads back and unrolls (the ports were derived from the sockets as in / out before)."""
    import json

    from rtt_fusion.graph import io, workflow
    from rtt_fusion.graph.topology import model_from_graph, unroll

    reg, _ = workflow.registry()
    m = model_from_graph(_transmission_chain(back=False), reg)
    doc = workflow.model_document(m, reg)
    assert not [d for d in doc["diagnostics"] if d["severity"] == "error"], doc["diagnostics"]
    g1 = next(n for n in doc["model"]["nodes"] if n["id"] == "G1")
    assert [p["name"] for p in g1["ports"]] == ["front", "back"]
    back = io.model_from_dict(json.loads(json.dumps(doc["model"])), reg)
    links = {(lk.source, lk.target) for lk in unroll(back, reg).links}
    assert (("G1", "back_out"), ("G2", "in")) in links


def test_replace_block_in_a_model_moves_ports_by_the_branches() -> None:
    """`replace_block` on a model (version 2): a thin lens between source and camera replaced by a
    transmission grating – the port the light enters through becomes `front`, the one it leaves
    through `back`; connections and branch steps follow, the model unrolls through `back_out`."""
    from rtt_fusion.graph import Model, replace_block, unroll, validate_model

    reg = builtin_registry()
    optics.register(reg)
    grating.register(reg)
    m = Model()
    m.add(Node("src", "optics.collimated_source", params={"radius": 0.01, "rings": 2}))
    m.add(Node("lens", "optics.thin_lens", params={"point": [0, 0, 1.0]}, position=(200.0, 0.0)))
    m.add(Node("cam", "optics.detector", params={"point": [0, 0, 3.0], "normal": [0, 0, -1]}))
    m.connect(("src", "out"), ("lens", "in"))
    m.connect(("lens", "out"), ("cam", "in"))
    m.branch("main", ("src", None, "out"), ("lens", "in", "out"), ("cam", "in", None))
    placed = replace_block(
        m,
        reg,
        "lens",
        Node(
            "lens",
            "fg.grating",
            params={
                "table": TABLE,
                "side": "transmitted",
                "point": [0, 0, 1.0],
                "normal": [0, 0, -1],
            },
        ),
    )
    assert placed.type_id == "fg.grating" and placed.position == (200.0, 0.0)
    assert {(c.a, c.b) for c in m.connections} == {
        (("src", "out"), ("lens", "front")),
        (("lens", "back"), ("cam", "in")),
    }
    assert ("lens", "front", "back") in [
        (s.node, s.port_in, s.port_out) for s in m.branches[0].steps
    ]
    assert validate_model(m, reg) == []
    links = {(lk.source, lk.target) for lk in unroll(m, reg).links}
    assert (("lens", "back_out"), ("cam", "in")) in links


def test_replace_block_in_a_model_mirror_to_grating_keeps_the_system(tmp_path: Path) -> None:
    """The Linnik's sample mirror (one port, entered and left through it) replaced by a grating in
    the model: its port becomes `front`, the block keeps the mirror's Raytatouille data (`rtt`,
    `runs` keyed by the new port), and the image equals the mirror's (W22 in the model)."""
    import json as _json

    from rtt_fusion.graph import backend_rtt, model_from_raytatouille, replace_block, unroll
    from rtt_fusion.graph.execute import Executor as _Executor

    rt = pytest.importorskip("raytatouille")
    from test_grating_rtt import PERIOD_UM, image, load_system, write_table

    d = load_system()
    table = write_table(tmp_path / "mirror.json", {0: -1.0}, PERIOD_UM)
    reg = builtin_registry()
    ev = dict(backend_rtt.EVALUATORS)
    ev.update(grating.register(reg))
    m = model_from_raytatouille(d)
    sample = next(n for n in m.nodes.values() if "sample" in n.id.lower())
    old_ports = {c.a[1] for c in m.connections if c.a[0] == sample.id} | {
        c.b[1] for c in m.connections if c.b[0] == sample.id
    }
    assert len(old_ports) == 1
    placed = replace_block(
        m, reg, sample.id, Node(sample.id, "fg.grating", params={"table": table})
    )
    assert "rtt" in placed.params and all("front" in k for k in placed.params.get("runs", {}))
    assert all(
        (s.port_in, s.port_out) == ("front", "front")
        for b in m.branches
        for s in b.steps
        if s.node == sample.id
    )
    g = unroll(m, reg)
    ex = _Executor(
        g, reg, ev, {"system": rt.System.from_json(_json.dumps(d)), "rings": 15, "pixels": [512, 8]}
    )
    camera = next(n for n in g.nodes if n.startswith("camera"))
    split = ex.evaluate(camera, "image")
    whole = image(d, None)
    assert np.abs(split.values - whole.values).max() < 1e-9 * whole.values.max()
