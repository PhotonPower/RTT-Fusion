"""R24: ports, connections and branches in the model; the unrolled socket form (ADR ports, poses,
branches and groups; acceptance PA1, PA2, PA10 and the unrolling of section 8).

Expected values are independent of the code under test: the version 1 graphs written by the
importer before R24 (tests/data/import_v1, generated from main at b0e3799), the two-beam
interference I = 2 I_arm (1 + cos 2k₀δ) of W5, and – for the unrolling – the branches
themselves: a light model whose light records its way returns, at the detector, exactly the
list of ways the modeller wrote (decision P6), with multiple visits and junctions.
"""

from __future__ import annotations

import glob
import json
import random
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import (
    Graph,
    Model,
    ModelError,
    Node,
    NodeType,
    Step,
    builtin_registry,
    display,
    graph_from_raytatouille,
    io,
    meta,
    model_from_dict,
    model_from_graph,
    model_from_raytatouille,
    model_to_dict,
    optics,
    unroll,
    validate_model,
)
from rtt_fusion.graph.execute import Context, Executor
from rtt_fusion.graph.model import Category, PortSpec, data_out, light_in, light_out
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.graph.validate import errors, node_sockets, node_through, validate

DATA = Path(__file__).parent / "data"
WAVELENGTH_MM = 0.6328e-3
I_ARM = 0.25


def _registry() -> tuple[NodeTypeRegistry, dict[str, Any]]:
    reg = builtin_registry()
    ev = dict(optics.register(reg))
    ev.update(meta.register(reg))
    ev.update(display.register(reg))
    return reg, ev


def _effective(g: Graph, reg: NodeTypeRegistry) -> tuple[Any, ...]:
    """What the computation sees: nodes with type, parameters and sockets, their dependencies
    between linked sockets, and the links."""
    linked = {lk.source for lk in g.links} | {lk.target for lk in g.links}
    nodes = {}
    for n in g.nodes.values():
        deps = {
            (o, tuple(sorted(i for i in ins if (n.id, i) in linked)))
            for o, ins in node_through(n, reg).items()
            if (n.id, o) in linked
        }
        nodes[n.id] = (
            n.type_id,
            n.name,
            json.dumps(n.params, sort_keys=True),
            n.muted,
            tuple(s.to_dict()["name"] for s in node_sockets(n, reg)),
            frozenset(deps),
        )
    links = frozenset((lk.source, lk.target, lk.transport) for lk in g.links)
    return nodes, links, g.meta.get("paths")


# ---------------------------------------------------------------------------- the Michelson


def _michelson_nodes(delta: float, theta: float = 0.0) -> list[Node]:
    return [
        Node(
            "src",
            "optics.collimated_source",
            params={
                "center": [0, 0, -10],
                "direction": [0, 0, 1],
                "radius": 2.0,
                "rings": 30,
                "wavelengths_um": [0.6328],
            },
        ),
        Node("bs", "optics.beam_splitter", params={"point": [0, 0, 0], "normal": [-1, 0, 1]}),
        Node("mA", "optics.plane_mirror", params={"point": [0, 0, 20.0], "normal": [0, 0, -1]}),
        Node(
            "mB",
            "optics.plane_mirror",
            params={"point": [20.0 + delta, 0, 0], "normal": [-np.cos(theta), np.sin(theta), 0]},
        ),
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [-30, 0, 0],
                "normal": [1, 0, 0],
                "half_size": [2.5, 2.5],
                "pixels": [256, 256],
            },
        ),
    ]


def michelson_v1(delta: float = 0.0, theta: float = 0.0) -> Graph:
    """The Michelson of W4/W5 in version 1: six directed links."""
    g = Graph()
    for n in _michelson_nodes(delta, theta):
        g.add(n)
    g.connect(("src", "out"), ("bs", "in"))
    g.connect(("bs", "t"), ("mA", "in"))
    g.connect(("mA", "out"), ("bs", "t_in"))
    g.connect(("bs", "r"), ("mB", "in"))
    g.connect(("mB", "out"), ("bs", "r_in"))
    g.connect(("bs", "merged"), ("cam", "in"))
    return g


def michelson_v2(delta: float = 0.0, arms: str = "AB") -> Model:
    """The same Michelson as a model: one connection per arm, the arms as branches."""
    m = Model()
    for n in _michelson_nodes(delta):
        m.add(n)
    m.connect(("src", "out"), ("bs", "in"))
    m.connect(("bs", "t"), ("mA", "front"))
    m.connect(("bs", "r"), ("mB", "front"))
    m.connect(("bs", "merged"), ("cam", "in"))
    for arm, port, mirror in (("A", "t", "mA"), ("B", "r", "mB")):
        if arm in arms:
            m.branch(
                f"arm {arm}",
                ("src", None, "out"),
                ("bs", "in", port),
                (mirror, "front", "front"),
                ("bs", port, "merged"),
                ("cam", "in", None),
            )
    return m


def _centre(g: Graph) -> float:
    reg, ev = _registry()
    im = Executor(g, reg, ev).evaluate("cam", "image")
    return float(im.values[100:156, 100:156].mean())


# ---------------------------------------------------------------------------- W39 (PA1)


@pytest.mark.parametrize(
    "path", sorted(glob.glob(str(DATA / "*.rtt.json"))), ids=lambda p: Path(p).name
)
def test_w39_import_unchanged(path: str) -> None:
    """The importer now builds the model and unrolls it: the socket graph is the one of the
    importer before R24."""
    reg, _ = _registry()
    before = io.load(DATA / "import_v1" / Path(path).name.replace(".rtt.json", ".graph.json"))
    after = graph_from_raytatouille(Path(path).read_text(encoding="utf-8"))
    assert _effective(after, reg) == _effective(before, reg)


@pytest.mark.parametrize(
    "path", sorted(glob.glob(str(DATA / "import_v1" / "*.json"))), ids=lambda p: Path(p).name
)
def test_w39_version_1_round_trip(path: str) -> None:
    """Version 1 -> model -> file version 2 -> model -> socket form: the same graph."""
    reg, _ = _registry()
    g = io.load(path)
    m = model_from_graph(g, reg)
    text = json.dumps(model_to_dict(m))
    again = model_from_dict(json.loads(text), reg)
    assert model_to_dict(again) == model_to_dict(m)
    assert _effective(unroll(again, reg), reg) == _effective(g, reg)


def test_w39_michelson_with_observer_and_drive() -> None:
    reg, _ = _registry()
    g = michelson_v1(delta=WAVELENGTH_MM / 8)
    g.add(Node("probe", "display.beam"))
    g.connect(("bs", "t"), ("probe", "in"))  # an observer: a second link at a light output
    g.add(Node("sweep", "meta.sweep", params={"values": [0.0, 1e-4]}))
    g.connect(("sweep", "drive"), ("mB", "drive"))
    g.connect(("cam", "image"), ("sweep", "observe"))
    m = model_from_graph(g, reg)
    assert len(m.connections) == 4 and [b.name for b in m.branches] == ["branch 1", "branch 2"]
    assert {lk.target for lk in m.links} == {("probe", "in"), ("mB", "drive"), ("sweep", "observe")}
    u = unroll(m, reg)
    assert _effective(u, reg)[:2] == _effective(g, reg)[:2]
    assert not errors(validate(u, reg))
    assert _centre(u) == _centre(g)  # bit for bit: the same computation


# ---------------------------------------------------------------------------- W40 (PA2)


@pytest.mark.parametrize(
    ("delta", "expected"),
    [(0.0, 4 * I_ARM), (WAVELENGTH_MM / 8, 2 * I_ARM), (WAVELENGTH_MM / 4, 0.0)],
)
def test_w40_michelson_from_ports_and_branches(delta: float, expected: float) -> None:
    # W5: I = 2 I_arm (1 + cos 2k₀δ), with one connection per arm instead of two links
    m = michelson_v2(delta)
    assert len(m.connections) == 4
    assert _centre(unroll(m, _registry()[0])) == pytest.approx(expected, abs=1e-3)


def test_w40_exactly_the_branches() -> None:
    # only arm A in the model: one beam, |r t|² = 1/4, no interference, whatever δ is
    for delta in (0.0, WAVELENGTH_MM / 4):
        g = unroll(michelson_v2(delta, arms="A"), _registry()[0])
        assert not g.links_to("mB")  # the other arm is not computed at all
        assert _centre(g) == pytest.approx(I_ARM, abs=1e-3)


def test_w40_linnik_two_connections_per_arm() -> None:
    m = model_from_raytatouille((DATA / "linnik_structure.rtt.json").read_text(encoding="utf-8"))
    v1 = io.load(DATA / "import_v1" / "linnik_structure.graph.json")
    for objective in ("Objective_A", "Objective_B"):
        assert len([c for c in m.connections if objective in (c.a[0], c.b[0])]) == 2
        assert len([lk for lk in v1.links if objective in (lk.source[0], lk.target[0])]) == 4
    assert [b.name for b in m.branches] == ["measurement", "reference", "measurement_pupil"]
    measurement = m.branches[0].steps
    assert [s.node for s in measurement].count("Objective_A") == 2  # there and back
    assert m.meta["raytatouille_paths"] == ["measurement", "reference", "measurement_pupil"]


# ---------------------------------------------------------------------------- W41 (PA10)


def _codes(m: Model) -> set[str]:
    return {d.code for d in validate_model(m, _registry()[0]) if d.severity == "error"}


def test_w41_no_branch() -> None:
    m = michelson_v2(arms="")
    assert _codes(m) == {"graph.no_branch"}
    with pytest.raises(ModelError) as e:
        unroll(m, _registry()[0])
    assert {d.code for d in e.value.diagnostics} == {"graph.no_branch"}
    m = michelson_v2()
    m.add(Node("cam2", "optics.detector"))  # a detector no branch reaches
    assert _codes(m) == {"graph.no_branch"}


def test_w41_impossible_branches() -> None:
    m = michelson_v2(arms="A")
    m.branch("through", ("src", None, "out"), ("bs", "in", "merged"), ("cam", "in", None))
    assert _codes(m) == {"graph.branch_impossible"}  # the splitter gives nothing in -> merged
    m = michelson_v2()
    m.connections = [c for c in m.connections if c.b != ("mB", "front")]
    assert _codes(m) == {"graph.branch_impossible"}  # arm B over a port without connection
    m = michelson_v2(arms="A")
    m.branch("short", ("src", None, "out"), ("bs", "in", "t"))
    assert _codes(m) == {"graph.branch_impossible"}  # does not end at a detector
    m = michelson_v2()
    m.branches.append(m.branches[0])
    assert _codes(m) == {"graph.branch_duplicate"}
    m = michelson_v2()
    m.connect(("mA", "front"), ("mB", "front"))
    assert "graph.port_links" in _codes(m)


def test_w41_nothing_computed_on_error() -> None:
    m = michelson_v2(arms="")
    with pytest.raises(ModelError):
        unroll(m, _registry()[0])


# ---------------------------------------------------------------------------- W42 (unrolling)

_TWO = PortSpec("front", "front_in", "front_out", "front"), PortSpec("back", "back_in", "back_out")


def _way_registry() -> tuple[NodeTypeRegistry, dict[str, Any]]:
    """Blocks whose light is the record of its way: a source gives the empty way, a block appends
    its step, a detector returns the ways that arrive."""
    reg = builtin_registry()
    reg.add(NodeType("way.source", Category.SOURCE, (light_out("out", "R"),)))
    reg.add(
        NodeType(
            "way.block",
            Category.DOMAIN,
            (
                light_in("front_in", "R"),
                light_out("front_out", "R"),
                light_in("back_in", "R"),
                light_out("back_out", "R"),
            ),
            through=(
                ("front_out", ("front_in", "back_in")),
                ("back_out", ("front_in", "back_in")),
            ),
            ports=_TWO,
            scatter=tuple((i, o, "full", ()) for i in ("front", "back") for o in ("front", "back")),
        )
    )
    reg.add(
        NodeType(
            "way.detector",
            Category.DETECTOR,
            (light_in("in", "R"), data_out("image", "result")),
            through=(("image", ("in",)),),
        )
    )

    def port(socket: str) -> str:
        return socket.rsplit("_", 1)[0]

    def source(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
        return [(Step(node.id, None, "out"),)]

    def block(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
        base = node.id.split("@")[0]
        return [
            (*w, Step(base, port(i), port(socket))) for i, ws in inputs.items() for w in ws or []
        ]

    def detector(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
        return sorted(((*w, Step(node.id, "in", None)) for w in inputs["in"]), key=repr)

    return reg, {"way.source": source, "way.block": block, "way.detector": detector}


def _chain(n: int) -> Model:
    """source — b1 — b2 — … — bn — detector, every block reflecting and transmitting."""
    m = Model()
    m.add(Node("s", "way.source"))
    for k in range(1, n + 1):
        m.add(Node(f"b{k}", "way.block"))
    m.add(Node("d", "way.detector"))
    m.connect(("s", "out"), ("b1", "front"))
    for k in range(1, n):
        m.connect((f"b{k}", "back"), (f"b{k + 1}", "front"))
    m.connect((f"b{n}", "back"), ("d", "in"))
    return m


def _random_way(rng: random.Random, n: int, bounces: int) -> tuple[Step, ...]:
    """A way through the chain with up to `bounces` reflections."""
    steps = [Step("s", None, "out")]
    k, forward, left = 1, True, bounces
    while True:
        reflect = left > 0 and rng.random() < 0.3
        port_in = "front" if forward else "back"
        if reflect:
            forward, left = not forward, left - 1
        port_out = "back" if forward else "front"
        steps.append(Step(f"b{k}", port_in, port_out))
        k += 1 if forward else -1
        if k > n:
            steps.append(Step("d", "in", None))
            return tuple(steps)
        if k < 1:  # back at the source side: turn around in b1 instead
            raise ValueError


def _ways(rng: random.Random, n: int, count: int) -> list[tuple[Step, ...]]:
    out: list[tuple[Step, ...]] = []
    for _ in range(200):  # a short chain has fewer different ways than asked for
        if len(out) == count:
            break
        try:
            w = _random_way(rng, n, bounces=rng.choice((0, 2, 2, 4)))
        except ValueError:
            continue
        if w not in out:
            out.append(w)
    return out


@pytest.mark.parametrize("seed", range(40))
def test_w42_detector_sees_exactly_the_branches(seed: int) -> None:
    rng = random.Random(seed)
    reg, ev = _way_registry()
    n = rng.choice((1, 2, 3, 4))
    m = _chain(n)
    ways = _ways(rng, n, rng.choice((1, 2, 3, 5, 8)))
    for k, w in enumerate(ways):
        m.branch(f"w{k}", *w)
    g = unroll(m, reg)
    assert not errors(validate(g, reg))  # a DAG, one link per input
    seen = Executor(g, reg, ev).evaluate("d", "image")
    assert seen == sorted(ways, key=repr)  # each branch once, nothing else


def test_w42_junction() -> None:
    """Ways that meet again after different histories need a junction block: b2 is passed
    straight by w0 and w1, but only w0's history also reflects back (w2)."""
    reg, ev = _way_registry()
    m = _chain(3)
    s, d = Step("s", None, "out"), Step("d", "in", None)
    straight = [
        Step("b1", "front", "back"),
        Step("b2", "front", "back"),
        Step("b3", "front", "back"),
    ]
    bounce = [
        Step("b1", "front", "back"),
        Step("b2", "front", "front"),
        Step("b1", "back", "back"),
    ]
    w0 = (s, *straight, d)
    w1 = (s, *bounce, Step("b2", "front", "back"), Step("b3", "front", "back"), d)
    w2 = (s, *bounce, Step("b2", "front", "front"), Step("b1", "back", "back"), *straight[1:], d)
    for k, w in enumerate((w0, w1, w2)):
        m.branch(f"w{k}", *w)
    g = unroll(m, reg)
    assert not errors(validate(g, reg))
    assert Executor(g, reg, ev).evaluate("d", "image") == sorted((w0, w1, w2), key=repr)
    assert [n.id for n in g.nodes.values() if "@" in n.id]  # multiple visits: copies


def test_w42_copies_keep_data_links() -> None:
    reg, _ = _registry()
    m = michelson_v2()
    m.add(Node("sweep", "meta.sweep", params={"values": [0.0]}))
    from rtt_fusion.graph.model import Link

    m.links.append(Link(("sweep", "drive"), ("mB", "drive")))
    g = unroll(m, reg)
    assert g.links_to("mB", "drive")
    assert g.meta["paths"] == {
        "arm A": ["bs", "mA", "bs", "cam"],
        "arm B": ["bs", "mB", "bs", "cam"],
    }
