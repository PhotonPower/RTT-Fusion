"""Planner stages 1 and 2 (R17): decision tables and the R15 equivalence.

The expected choices follow from the rules of docs/fusionsplan.md 7.2 and table 6.3 together with
the provisional cost ranking of the planner (R 1, P 1.5, H 3, S 4, G 6, T 8; operators 0.5–2); each
case states the arithmetic. The Linnik counts come from the structure of concept section 4 (two
arms meet at the camera).
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import pytest

from rtt_fusion.graph import (
    Category,
    Graph,
    Link,
    Node,
    NodeType,
    NodeTypeRegistry,
    Requirements,
    builtin_registry,
    graph_from_raytatouille,
    plan,
)
from rtt_fusion.graph.model import data_out, light_in, light_out
from rtt_fusion.graph.planner import CHAINS, link_key

DATA = Path(__file__).parent / "data"
ALL = frozenset("RHGSTP")


def registry(src: str = "RHGSTP", dom: str = "RHGSTP", det: str = "RHGSTP") -> NodeTypeRegistry:
    reg = NodeTypeRegistry()
    reg.add(NodeType("t.src", Category.SOURCE, (light_out("out", src),)))
    reg.add(NodeType("t.dom", Category.DOMAIN, (light_in("in", dom), light_out("out", dom))))
    reg.add(
        NodeType(
            "t.split",
            Category.SPLITTER,
            (light_in("in", dom), light_out("a", dom), light_out("b", dom)),
            through=(("a", ("in",)), ("b", ("in",))),
        )
    )
    reg.add(
        NodeType(
            "t.merge",
            Category.SPLITTER,
            (light_in("a", dom), light_in("b", dom), light_out("out", dom)),
            through=(("out", ("a", "b")),),
        )
    )
    reg.add(
        NodeType("t.det", Category.DETECTOR, (light_in("in", det), data_out("image", "result")))
    )
    return reg


def chain(t1: str = "auto", t2: str = "auto") -> Graph:
    """src -> dom -> det with transports t1, t2."""
    g = Graph()
    g.add(Node("src", "t.src"))
    g.add(Node("dom", "t.dom"))
    g.add(Node("det", "t.det"))
    g.connect(("src", "out"), ("dom", "in"), t1)
    g.connect(("dom", "out"), ("det", "in"), t2)
    return g


def reps(p: object, g: Graph) -> list[str]:
    from rtt_fusion.graph.planner import Plan

    assert isinstance(p, Plan)
    return [p.representation(lk) for lk in g.links]


def codes(p: object) -> list[str]:
    return [d.code for d in getattr(p, "diagnostics")]  # noqa: B009


# ------------------------------------------------------------------------------ decision tables


def test_auto_prefers_rays() -> None:
    # everything accepted and available: R is cheapest on every link, no conversion
    g = chain()
    p = plan(g, registry())
    assert p.ok and reps(p, g) == ["R", "R"] and not p.conversions
    assert p.cost == 2.0
    assert "lowest cost" in p.reasons[link_key(g.links[0])]


def test_field_only_detector_converts_upstream() -> None:
    # detector takes only H: link 2 = H (only admissible); link 1: R (1 + R->H 1 = 2) beats
    # H (3 + 0 = 3), so the domain converts rays to a field
    g = chain()
    p = plan(g, registry(det="H"))
    assert p.ok and reps(p, g) == ["R", "H"]
    assert [(c.node, c.chain) for c in p.conversions] == [("dom", "RH")]
    assert [op.name for op in p.conversions[0].operators] == ["rays_to_field"]
    assert p.cost == 1.0 + 3.0 + 1.0
    assert p.reasons[link_key(g.links[1])] == "only H admissible"


def test_fixed_transport_is_kept() -> None:
    # fixed H on link 1 although R would be cheaper; link 2 stays R: conversion H -> R
    g = chain("H")
    p = plan(g, registry())
    assert p.ok and reps(p, g) == ["H", "R"]
    assert [c.chain for c in p.conversions] == ["HR"]
    assert p.reasons[link_key(g.links[0])] == "fixed by the user (H)"


def test_fixed_transport_not_available() -> None:
    # stage R15/R17 executes only R: a fixed H cannot run yet
    p = plan(chain("H"), registry(), Requirements(available=frozenset("R")))
    assert not p.ok and codes(p) == ["plan.representation_unavailable"]


def test_auto_without_available_representation() -> None:
    # detector takes only H, only R available: link 2 has no admissible representation
    p = plan(chain(), registry(det="H"), Requirements(available=frozenset("R")))
    assert not p.ok and "plan.no_representation" in codes(p)


def test_conversion_chain() -> None:
    # G fixed, then R fixed: G->H->R (1 + 1 = 2) beats G->S->R (2 + 1 = 3)
    g = chain("G", "R")
    p = plan(g, registry())
    assert p.ok and [c.chain for c in p.conversions] == ["GHR"]


def test_no_operator_after_phase_loss() -> None:
    # P has no phase: P -> R ends the chain, a coherent field H cannot follow
    p = plan(chain("P", "H"), registry())
    assert not p.ok and codes(p) == ["plan.no_operator"]
    assert CHAINS[("P", "R")][0] == 0.5 and CHAINS[("P", "H")][0] == float("inf")


def test_validity_hook() -> None:
    # the hook rejects rays on link 2 (e.g. a caustic at the port): H instead, converted in dom
    def validity(link: Link, rep: str) -> str | None:
        return "caustic at the port" if link.target[0] == "det" and rep == "R" else None

    g = chain()
    p = plan(g, registry(), Requirements(validity=validity))
    assert p.ok and reps(p, g) == ["R", "H"] and [c.chain for c in p.conversions] == ["RH"]
    # with only R available the reasons are reported
    p2 = plan(g, registry(), Requirements(available=frozenset("R"), validity=validity))
    assert not p2.ok
    assert "caustic at the port" in p2.diagnostics[0].message


def test_shared_input_sums_the_demands() -> None:
    # a splitter input feeds two outputs: one consumer needs H, the other takes anything.
    # input R: 1 + (R->H 1) + 0 = 2; input H: 3 + 0 + (H->R 1) = 4 -> R, one conversion
    reg = registry()
    reg.add(
        NodeType("t.hdet", Category.DETECTOR, (light_in("in", "H"), data_out("image", "result")))
    )
    g = Graph()
    for n in (
        Node("src", "t.src"),
        Node("sp", "t.split"),
        Node("d1", "t.hdet"),
        Node("d2", "t.det"),
    ):
        g.add(n)
    g.connect(("src", "out"), ("sp", "in"))
    g.connect(("sp", "a"), ("d1", "in"))
    g.connect(("sp", "b"), ("d2", "in"))
    p = plan(g, reg)
    assert p.ok and reps(p, g) == ["R", "H", "R"]
    assert [(c.node, c.out_socket, c.chain) for c in p.conversions] == [("sp", "a", "RH")]
    assert len(p.paths) == 2


def test_incoherent_merge_warns() -> None:
    # two arms meet; one arrives as P (no phase): intensities, not amplitudes -> warning
    g = Graph()
    for n in (
        Node("src", "t.src"),
        Node("sp", "t.split"),
        Node("m", "t.merge"),
        Node("det", "t.det"),
    ):
        g.add(n)
    g.connect(("src", "out"), ("sp", "in"))
    g.connect(("sp", "a"), ("m", "a"), "P")
    g.connect(("sp", "b"), ("m", "b"))
    g.connect(("m", "out"), ("det", "in"))
    p = plan(g, registry())
    assert p.ok
    assert codes(p) == ["plan.incoherent_merge"]
    assert p.diagnostics[0].location == "det.image"
    # both arms coherent: no warning
    g.set_transport(g.links[1], "auto")
    assert codes(plan(g, registry())) == []


# ------------------------------------------------------------------------------ Linnik, R15


def linnik() -> Graph:
    return graph_from_raytatouille(DATA / "linnik_structure.rtt.json")


def test_linnik_topology() -> None:
    p = plan(linnik(), builtin_registry(), Requirements(available=frozenset("R")))
    assert p.ok
    to_camera = [path for path in p.paths if path.output == ("Camera", "image")]
    assert len(to_camera) == 2  # sample arm and reference arm
    arms = {
        next(lk.target[0] for lk in path.links if lk.target[0] in ("Sample", "Reference"))
        for path in to_camera
    }
    assert arms == {"Sample", "Reference"}
    assert all(path.source == ("source", "to_Illumination_stop") for path in p.paths)


def test_linnik_auto_is_rays() -> None:
    # acceptance R17: with "auto" the plan is R on all 13 links, as executed in R15
    g = linnik()
    for available in (frozenset("R"), ALL):
        p = plan(g, builtin_registry(), Requirements(available=available))
        assert p.ok and set(reps(p, g)) == {"R"} and len(p.representations) == 13
        assert not p.conversions


def test_path_limit() -> None:
    p = plan(linnik(), builtin_registry(), Requirements(max_paths=1))
    assert len(p.paths) == 1 and "plan.paths_truncated" in codes(p)


def test_plan_as_json() -> None:
    g = chain()
    p = plan(g, registry(det="H"))
    d = json.loads(json.dumps(p.to_dict()))
    assert d["links"][link_key(g.links[1])]["representation"] == "H"
    assert d["conversions"][0]["operators"] == ["rays_to_field"]
    assert d["paths"] == [[link_key(g.links[0]), link_key(g.links[1])]]


@pytest.mark.parametrize("pair", sorted(CHAINS))
def test_chains_follow_table_6_3(pair: tuple[str, str]) -> None:
    # every chain is a sequence of operators of the table, joined end to start, with phase loss
    # only at the end
    cost, ops = CHAINS[pair]
    if not ops:
        return
    assert ops[0].source == pair[0] and ops[-1].target == pair[1]
    assert all(a.target == b.source for a, b in itertools.pairwise(ops))
    assert all(op.coherent for op in ops[:-1])
    assert cost == sum(op.cost for op in ops)
