"""R14: working graph of the node workflow (docs/konzept-node-workflow.md).

Expected values come from the concept (section 4: the Linnik graph) and from the definitions of
the checks, not from the code under test. No solver modules needed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from rtt_fusion.graph import (
    Graph,
    GraphError,
    Node,
    builtin_registry,
    dumps,
    errors,
    find_cycle,
    graph_from_raytatouille,
    loads,
    node_type_from_describe,
    validate,
)
from rtt_fusion.graph.model import Direction, Link, data_in, data_out, light_in, light_out

LINNIK = Path(__file__).parent / "data" / "linnik_structure.rtt.json"


def codes(graph: Graph) -> list[str]:
    return [d.code for d in validate(graph, builtin_registry())]


# ---------------------------------------------------------------------------- import (concept 4)


def test_linnik_import_gives_the_concept_graph() -> None:
    g = graph_from_raytatouille(LINNIK)
    types = {n.name: n.type_id for n in g.nodes.values()}
    assert types == {
        "Source": "rtt.source",
        "Illumination stop": "rtt.section",
        "Beam splitter": "rtt.splitter",
        "Objective A": "rtt.section",  # stop and lenses of the assembly in one block
        "Sample": "rtt.mirror",  # exchangeable: a block of its own inside the assembly
        "Objective B": "rtt.section",
        "Reference": "rtt.mirror",
        "Pupil camera": "rtt.detector",
        "Tube and camera: Tube lens … Tube lens": "rtt.section",
        "Camera": "rtt.detector",
    }
    assert errors(validate(g, builtin_registry())) == []
    assert find_cycle(g, builtin_registry()) == []


def test_linnik_objective_is_passed_both_ways() -> None:
    g = graph_from_raytatouille(LINNIK)
    a = g.nodes["Objective_A"]
    names = {s.name: s for s in a.sockets}
    assert set(names) == {"from_Beam_splitter", "to_Sample", "from_Sample", "to_Beam_splitter"}
    # forward from the splitter (front in, back out), backward from the sample
    assert names["from_Beam_splitter"].side == "front"
    assert names["to_Sample"].side == "back"
    assert names["from_Sample"].side == "back"
    assert names["to_Beam_splitter"].side == "front"
    assert dict(a.through) == {
        "to_Sample": ("from_Beam_splitter",),
        "to_Beam_splitter": ("from_Sample",),
    }


def test_linnik_splitter_splits_and_combines() -> None:
    g = graph_from_raytatouille(LINNIK)
    through = dict(g.nodes["Beam_splitter"].through)
    tube = next(n.id for n in g.nodes.values() if n.name.startswith("Tube and camera"))
    # the light of both arms is combined towards the tube lens (interference at the camera)
    assert set(through[f"to_{tube}"]) == {"from_Objective_A", "from_Objective_B"}
    assert through["to_Objective_A"] == ("from_Illumination_stop",)
    assert through["to_Objective_B"] == ("from_Illumination_stop",)
    # one link per connection, shared by the paths
    assert len(g.links) == 13
    assert g.meta["paths"]["measurement"][-1] == "Camera"
    assert g.meta["paths"]["measurement_pupil"][-1] == "Pupil_camera"


def test_auto_path_is_a_chain() -> None:
    d = {
        "root": {
            "children": [
                {"type": "stop", "name": "STO", "surfaces": [{"id": "STO"}]},
                {"type": "lens", "name": "L1", "surfaces": [{"id": "L1.S1"}, {"id": "L1.S2"}]},
                {"type": "detector", "name": "IMG", "surfaces": [{"id": "IMG"}]},
            ]
        },
        "paths": [{"name": "main", "events": "auto"}],
    }
    g = graph_from_raytatouille(d)
    assert [n.type_id for n in g.nodes.values()] == [
        "rtt.source",
        "rtt.section",
        "rtt.section",
        "rtt.detector",
    ]
    assert len(g.links) == 3
    assert errors(validate(g, builtin_registry())) == []


def test_unknown_surface_in_a_path_is_reported() -> None:
    d = {
        "root": {"children": [{"type": "lens", "name": "L", "surfaces": [{"id": "S1"}]}]},
        "paths": [{"name": "p", "events": [{"surface": "S9"}]}],
    }
    with pytest.raises(ValueError, match="S9"):
        graph_from_raytatouille(d)


# ---------------------------------------------------------------------------- checks


def two_blocks(out_reps: str, in_reps: str, transport: str = "auto") -> Graph:
    g = Graph()
    g.add(Node("a", "rtt.section", sockets=(light_out("o", out_reps),)))
    g.add(Node("b", "rtt.section", sockets=(light_in("i", in_reps),)))
    g.connect(("a", "o"), ("b", "i"), transport)
    return g


def test_representations_must_meet() -> None:
    assert codes(two_blocks("G", "GS")) == []
    assert codes(two_blocks("G", "S")) == [
        "graph.no_common_representation"
    ]  # e.g. BWPM -> FourierGrating
    assert codes(two_blocks("RH", "RHS", "S")) == ["graph.transport_unsupported"]
    assert codes(two_blocks("RH", "RHS", "H")) == []


def test_kinds_data_types_direction_and_inputs() -> None:
    g = Graph()
    g.add(Node("a", "rtt.section", sockets=(light_out("o"), data_out("n", "N"), light_in("i"))))
    g.add(Node("b", "rtt.section", sockets=(light_in("i"), data_in("q", "Q"), data_in("n", "N"))))
    g.connect(("a", "o"), ("b", "n"))  # light into data
    g.connect(("a", "n"), ("b", "q"))  # N into Q
    g.connect(("b", "i"), ("a", "i"))  # input to input
    g.connect(("a", "x"), ("b", "i"))  # unknown socket
    assert sorted(codes(g)) == sorted(
        [
            "graph.kind_mismatch",
            "graph.data_type_mismatch",
            "graph.link_direction",
            "graph.socket_unknown",
        ]
    )
    g2 = Graph()
    for nid in ("a", "c"):
        g2.add(Node(nid, "rtt.section", sockets=(light_out("o"),)))
    g2.add(Node("b", "rtt.section", sockets=(light_in("i"),)))
    g2.connect(("a", "o"), ("b", "i"))
    g2.connect(("c", "o"), ("b", "i"))
    assert codes(g2) == ["graph.input_multiple"]


def test_feedback_is_a_cycle_but_a_return_path_is_not() -> None:
    reg = builtin_registry()
    # two mirrors facing each other: a resonator, i.e. feedback
    g = Graph()
    g.add(Node("m1", "rtt.mirror"))
    g.add(Node("m2", "rtt.mirror"))
    g.connect(("m1", "out"), ("m2", "in"))
    g.connect(("m2", "out"), ("m1", "in"))
    assert "graph.cycle" in codes(g)
    assert find_cycle(g, reg)
    # an objective passed forth and back via a mirror: no cycle on the socket level
    g = Graph()
    g.add(Node("obj", "rtt.section"))
    g.add(Node("m", "rtt.mirror"))
    g.connect(("obj", "back_out"), ("m", "in"))
    g.connect(("m", "out"), ("obj", "back_in"))
    assert find_cycle(g, reg) == []


def test_unknown_type_and_parameter() -> None:
    g = Graph()
    g.add(Node("x", "nobody.knows"))
    g.add(Node("c", "fusion.convert", params={"to": "S", "colour": "red"}))
    assert sorted(codes(g)) == ["graph.node_type_unknown", "graph.param_unknown"]


# ---------------------------------------------------------------------------- editing


def test_insert_replace_mute() -> None:
    g = graph_from_raytatouille(LINNIK)
    # exchange the sample mirror for a grating block (concept 1: Linnik microscope); links stay
    grating = node_type_from_describe(
        {
            "name": "fouriergrating",
            "representation": "S",
            "ports": [{"name": "front"}, {"name": "back"}],
        }
    )
    before = len(g.links)
    g.replace_node(
        "Sample",
        Node("Grating", grating.type_id),
        {"from_Objective_A": "front_in", "to_Objective_A": "front_out"},
    )
    assert len(g.links) == before
    assert g.links_to("Grating", "front_in") and g.links_from("Grating", "front_out")
    # insert a conversion on a link (drop on the link)
    link = g.links_to("Grating", "front_in")[0]  # objective -> grating
    assert link.source == ("Objective_A", "to_Sample")  # the objective keeps its socket
    g.insert_on_link(link, Node("conv", "fusion.convert", params={"to": "S"}), "in", "out")
    assert g.links_to("conv", "in") and g.links_from("conv", "out")
    assert len(g.links) == before + 1
    g.set_muted("conv")
    assert g.nodes["conv"].muted
    with pytest.raises(GraphError):
        g.connect(("conv", "out"), ("Grating", "front_in"))  # the same link twice
    with pytest.raises(GraphError):
        g.remove("nobody")


def test_transport_is_kept_and_checked() -> None:
    g = two_blocks("RHS", "RHS")
    link = g.links[0]
    fixed = g.set_transport(link, "H")
    assert fixed.transport == "H" and g.links == [fixed]
    with pytest.raises(ValueError):
        Link(("a", "o"), ("b", "i"), "X")


# ---------------------------------------------------------------------------- adapter types (A1/A2)


def test_adapter_without_socket_description() -> None:
    t = node_type_from_describe({"name": "bwpm", "representation": "G", "method": "BWPM"})
    assert [(s.name, s.direction, "".join(sorted(s.accepts))) for s in t.sockets] == [
        ("front_in", Direction.IN, "G"),
        ("front_out", Direction.OUT, "G"),
        ("back_out", Direction.OUT, "G"),
    ]
    back = node_type_from_describe(
        {"name": "x", "representation": "S", "ports": [{"name": "back", "accepts_incidence": True}]}
    )
    assert back.socket("back_in") is not None


def test_adapter_with_socket_and_parameter_description() -> None:
    t = node_type_from_describe(
        {
            "name": "cliffordbem",
            "sockets": [
                light_in("in", "S").to_dict(),
                light_out("scattered", "S").to_dict(),
                data_in("emitters", "Q").to_dict(),
            ],
            "through": {"scattered": ["in", "emitters"]},
            "parameters": [
                {"name": "radius", "type": "number", "default": 0.05, "unit": "µm", "min": 0}
            ],
        }
    )
    assert t.socket("emitters") is not None
    assert t.params[0].unit == "µm" and t.params[0].minimum == 0
    assert dict(t.through) == {"scattered": ("in", "emitters")}


# ---------------------------------------------------------------------------- file format


def test_round_trip_and_version() -> None:
    g = graph_from_raytatouille(LINNIK)
    g.set_transport(g.links[0], "R")
    g.set_muted("Reference")
    text = dumps(g)
    back = loads(text)
    assert dumps(back) == text
    assert validate(back, builtin_registry()) == validate(g, builtin_registry())
    d = json.loads(text)
    assert d["format"] == "rtt-fusion-graph" and d["version"] == 1
    d["version"] = 2
    with pytest.raises(ValueError, match="version 2"):
        loads(json.dumps(d))
    with pytest.raises(ValueError, match="not a"):
        loads(json.dumps({"format": "something"}))
