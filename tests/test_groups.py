"""R28: groups and variants (ADR ports, poses, branches, section 7; acceptance PA5 open/closed,
PA6).

Expected values are independent of the code under test: the same set-up without the group (PA5,
bit for bit), the mirror tilt of W47 (fringe period λ / sin 2α) through the group's pose input and
through an exported parameter, and matrix optics for the variants (PA6): a single thin lens f and
a pair f₁, f₂ at distance d with f = f₁f₂ / (f₁ + f₂ − d) bend a collimated ray at height h to
the slope −h/f.
"""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import (
    Group,
    Model,
    Node,
    Variant,
    builtin_registry,
    meta,
    model_from_dict,
    model_to_dict,
    optics,
    set_variant,
    unroll,
    validate_model,
)
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.model import Link
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.graph.topology import Branch, Connection, Step

sys.path.insert(0, str(Path(__file__).parent))
from test_execute import fringe_period
from test_topology import michelson_v2

WAVELENGTH_MM = 0.6328e-3


def _registry() -> tuple[NodeTypeRegistry, dict[str, Any]]:
    reg = builtin_registry()
    ev = dict(optics.register(reg))
    ev.update(meta.register(reg))
    return reg, ev


def _image(m: Model, node: str = "cam", socket: str = "image") -> Any:
    reg, ev = _registry()
    return Executor(unroll(m, reg), reg, ev).evaluate(node, socket)


def grouped(
    group_pose: dict[str, Any] | None = None,
    mirror: dict[str, Any] | None = None,
    params: dict[str, str] | None = None,
    nested: bool = False,
) -> Model:
    """W5's Michelson with arm B's mirror inside the group "armB" (port "face")."""
    m = michelson_v2()
    flat_mirror = m.nodes.pop("mB")
    inner = Model()
    if nested:
        holder = Model()
        holder.add(Node("M", flat_mirror.type_id, params=mirror or dict(flat_mirror.params)))
        holder.branch("bounce", ("M", "front", "front"))
        Group({"only": Variant(holder, {"face": ("M", "front")})}, "only").add_to(
            inner, Node("holder", "group")
        )
        inner.branch("bounce", ("holder", "face", "face", "bounce"))
        ports = {"face": ("holder", "face")}
    else:
        inner.add(Node("M", flat_mirror.type_id, params=mirror or dict(flat_mirror.params)))
        inner.branch("bounce", ("M", "front", "front"))
        ports = {"face": ("M", "front")}
    Group({"plain": Variant(inner, ports, params or {})}, "plain").add_to(
        m, Node("armB", "group", params={"pose": group_pose} if group_pose else {})
    )
    m.connections = [
        Connection(c.a, ("armB", "face")) if c.b == ("mB", "front") else c for c in m.connections
    ]
    m.branches = [
        Branch(
            b.name,
            tuple(Step("armB", "face", "face", "bounce") if s.node == "mB" else s for s in b.steps),
        )
        for b in m.branches
    ]
    return m


# ---------------------------------------------------------------------------- W53 (PA5)


def test_w53_group_closed_equals_open() -> None:
    flat: Image = _image(michelson_v2(WAVELENGTH_MM / 8))
    m = michelson_v2(WAVELENGTH_MM / 8)
    closed = grouped(mirror=dict(m.nodes["mB"].params))
    np.testing.assert_array_equal(_image(closed).values, flat.values)
    nested = grouped(mirror=dict(m.nodes["mB"].params), nested=True)
    np.testing.assert_array_equal(_image(nested).values, flat.values)
    reg, _ = _registry()
    assert "armB/holder/M" in unroll(nested, reg).nodes  # path addresses


def test_w53_blocks_lie_relative_to_their_group() -> None:
    flat: Image = _image(michelson_v2(WAVELENGTH_MM / 8))
    # the mirror at the group's origin, the group where the flat mirror was
    m = grouped(
        group_pose={"position": [20.0 + WAVELENGTH_MM / 8, 0.0, 0.0], "rotation": [0, 0, 0]},
        mirror={"point": [0, 0, 0], "normal": [-1, 0, 0]},
    )
    assert np.abs(_image(m).values - flat.values).max() < 1e-9 * flat.values.max()


def _tilted_period(m: Model) -> float:
    result = _image(m, "tilt", "results")
    im: Image = result.results[0]
    return fringe_period(im.values, im.extent, 1.5)


def test_w53_group_pose_input_and_exported_parameter() -> None:
    alpha = 1e-3
    expected = WAVELENGTH_MM / np.sin(2 * alpha)
    # a sweep drives the group's pose input (rz about the group's origin = the mirror)
    m = grouped(
        group_pose={"position": [20.0, 0, 0], "rotation": [0, 0, 0]},
        mirror={"point": [0, 0, 0], "normal": [-1, 0, 0]},
    )
    m.add(Node("tilt", "meta.sweep", params={"values": [alpha], "component": "rz"}))
    m.links += [
        Link(("tilt", "drive"), ("armB", "pose")),
        Link(("cam", "image"), ("tilt", "observe")),
    ]
    assert _tilted_period(m) == pytest.approx(expected, rel=1e-6)
    # an exported parameter "tilt" of the group, swept by its address
    m = grouped(
        mirror={"point": [20.0, 0, 0], "normal": [-1, 0, 0]}, params={"tilt": "M.pose_offset[5]"}
    )
    m.add(Node("tilt", "meta.sweep", params={"values": [alpha], "parameter": "armB.tilt"}))
    m.links.append(Link(("cam", "image"), ("tilt", "observe")))
    assert _tilted_period(m) == pytest.approx(expected, rel=1e-6)


def test_w53_file_round_trip() -> None:
    m = grouped(params={"tilt": "M.pose_offset[5]"}, nested=True)
    reg, _ = _registry()
    again = model_from_dict(json.loads(json.dumps(model_to_dict(m))), reg)
    assert model_to_dict(again) == model_to_dict(m)
    np.testing.assert_array_equal(_image(again).values, _image(m).values)


def test_w53_group_errors() -> None:
    reg, _ = _registry()
    m = grouped()
    m.branches[1] = Branch(
        m.branches[1].name,
        tuple(replace_channel(s, "missing") for s in m.branches[1].steps),
    )
    assert {d.code for d in validate_model(m, reg)} == {"graph.group_branch"}
    m = grouped()
    other = copy.deepcopy(m.groups["armB"].variants["plain"])
    other.ports = {"side": other.ports["face"]}
    m.groups["armB"].variants["other"] = other
    assert "graph.group_ports" in {d.code for d in validate_model(m, reg)}


def replace_channel(s: Step, channel: str) -> Step:
    return Step(s.node, s.port_in, s.port_out, channel) if s.node == "armB" else s


# ---------------------------------------------------------------------------- W54 (PA6)

F_SINGLE, F1, F2, D = 50.0, 100.0, 80.0, 20.0  # f₁f₂ / (f₁ + f₂ − d) = 50


def objective() -> Model:
    m = Model()
    m.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"center": [0, 0, -5], "radius": 2.0, "rings": 5, "polarisation": [0, 1, 0]},
        )
    )
    single = Model()
    single.add(Node("L", "optics.thin_lens", params={"point": [0, 0, 0], "focal_length": F_SINGLE}))
    single.branch("through", ("L", "in", "out"))
    pair = Model()
    pair.add(Node("L1", "optics.thin_lens", params={"point": [0, 0, 0], "focal_length": F1}))
    pair.add(Node("L2", "optics.thin_lens", params={"point": [0, 0, D], "focal_length": F2}))
    pair.connect(("L1", "out"), ("L2", "in"))
    pair.branch("through", ("L1", "in", "out"), ("L2", "in", "out"))
    Group(
        {
            "ideal": Variant(single, {"in": ("L", "in"), "out": ("L", "out")}),
            "pair": Variant(pair, {"in": ("L1", "in"), "out": ("L2", "out")}),
        },
        "ideal",
    ).add_to(m, Node("obj", "group"))
    m.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, 300.0],
                "normal": [0, 0, -1],
                "half_size": [5, 5],
                "pixels": [8, 8],
            },
        )
    )
    m.connect(("src", "out"), ("obj", "in"))
    m.connect(("obj", "out"), ("cam", "in"))
    m.branch("main", ("src", None, "out"), ("obj", "in", "out", "through"), ("cam", "in", None))
    return m


def _focal_length(m: Model) -> Any:
    reg, ev = _registry()
    g = unroll(m, reg)
    feeding = g.links_to("cam")[0].source
    b = Executor(g, reg, ev).evaluate(*feeding)[0]
    start = Executor(g, reg, ev).evaluate("src", "out")[0]
    h = np.asarray(start.position)[:, 1]
    d = np.asarray(b.direction)
    slope = d[:, 1] / d[:, 2]
    use = np.abs(h) > 0.1
    return -h[use] / slope[use]


def test_w54_variants_keep_ports_and_focal_length() -> None:
    m = objective()
    connections = list(m.connections)
    f_ideal = _focal_length(m)
    set_variant(m, "obj", "pair")
    assert m.connections == connections  # the outside does not change
    assert set(m.groups["obj"].variants["pair"].ports) == set(
        m.groups["obj"].variants["ideal"].ports
    )
    f_pair = _focal_length(m)
    np.testing.assert_allclose(f_ideal, F_SINGLE, rtol=1e-12)
    np.testing.assert_allclose(f_pair, F1 * F2 / (F1 + F2 - D), rtol=1e-12)
    with pytest.raises(ValueError):
        set_variant(m, "obj", "unknown")


# ---------------------------------------------------------------------------- W60


def _grouped(m: Model, ids: set[str], gid: str) -> Model:
    """A group of blocks as the node editor makes it: connections across become exported
    ports "<block>.<port>", each stretch of a branch inside an inner branch the step names."""
    inner = Model()
    for i in ids:
        inner.add(m.nodes.pop(i))
    ports: dict[str, tuple[str, str]] = {}
    outer = []
    for c in m.connections:
        a, b = c.a[0] in ids, c.b[0] in ids
        if a and b:
            inner.connections.append(c)
            continue
        if a or b:
            end = c.a if a else c.b
            name = f"{end[0]}.{end[1]}"
            ports[name] = end
            c = (
                Connection((gid, name), c.b, c.transport)
                if a
                else Connection(c.a, (gid, name), c.transport)
            )
        outer.append(c)
    m.connections = outer
    branches = []
    for br in m.branches:
        steps: list[Step] = []
        k = 0
        while k < len(br.steps):
            if br.steps[k].node not in ids:
                steps.append(br.steps[k])
                k += 1
                continue
            run = []
            while k < len(br.steps) and br.steps[k].node in ids:
                run.append(br.steps[k])
                k += 1
            entry = next(p for p, e in ports.items() if e == (run[0].node, run[0].port_in))
            exit_ = next(p for p, e in ports.items() if e == (run[-1].node, run[-1].port_out))
            inner_name = next((b.name for b in inner.branches if list(b.steps) == run), None)
            if inner_name is None:
                inner_name = br.name
                inner.branches.append(Branch(inner_name, tuple(run)))
            steps.append(Step(gid, entry, exit_, inner_name))
        branches.append(Branch(br.name, tuple(steps)))
    m.branches = branches
    Group({"main": Variant(inner, ports)}, "main").add_to(m, Node(gid, "group"))
    return m


@pytest.mark.parametrize("kept", [True, False], ids=["blocks hold their runs", "import list"])
def test_w60_raytatouille_blocks_in_a_group_compute(kept: bool) -> None:
    rt = pytest.importorskip("raytatouille")
    from rtt_fusion.graph import backend_rtt, model_from_raytatouille

    d = json.loads(
        (Path(__file__).parent / "data" / "linnik_mini.rtt.json").read_text(encoding="utf-8")
    )
    reg, _ = _registry()

    def image(m: Model) -> Image:
        if (
            not kept
        ):  # a model from before S1: no data in the blocks, the import's list and the file
            for n in m.nodes.values():
                for key in ("rtt", "runs", "rtt_system", "rtt_tree", "rtt_paths"):
                    n.params.pop(key, None)
        g = unroll(m, reg)
        ex = Executor(
            g,
            reg,
            backend_rtt.EVALUATORS,
            {"system": rt.System.from_json(json.dumps(d)), "rings": 12, "pixels": [32, 32]},
        )
        result: Image = ex.evaluate("camera", "image")
        return result

    plain = image(model_from_raytatouille(d))
    m = _grouped(model_from_raytatouille(d), {"objective", "sample"}, "arm")
    assert not validate_model(m, reg)
    np.testing.assert_allclose(image(m).values, plain.values, atol=1e-12 * plain.values.max())


# ---------------------------------------------------------------------------- W61


def _sections(d: dict[str, Any]) -> list[str]:
    from rtt_fusion.graph import model_from_raytatouille

    m = model_from_raytatouille(d)
    return [
        n.id
        for n in m.nodes.values()
        if n.type_id == "rtt.section"
        and sum(1 for e in n.params.get("rtt", {}).values() if e.get("surfaces")) > 1
    ]


def _linnik_two_element_objective() -> dict[str, Any]:
    """linnik_mini with the objective as an assembly of two elements: a weak front lens (a
    window, f = 1 m) and the objective's lens – a section of several elements that computes."""
    d: dict[str, Any] = json.loads(
        (Path(__file__).parent / "data" / "linnik_mini.rtt.json").read_text(encoding="utf-8")
    )
    lens = d["root"]["children"][2]
    front = copy.deepcopy(lens)
    front["name"] = "front lens"
    front["pose"] = {"position": [0.0, 0.0, -0.5]}
    front["surfaces"][0]["id"] = "WIN"
    front["surfaces"][0]["interaction"] = {"type": "ideal_lens", "focal_length": 1000.0}
    lens = {**lens, "name": "objective lens", "pose": {"position": [0.0, 0.0, 0.0]}}
    d["root"]["children"][2] = {
        "type": "assembly",
        "name": "objective",
        "pose": {"position": [0.0, 0.0, 12.0]},
        "children": [front, lens],
    }
    events = d["paths"][0]["events"]
    k = next(i for i, e in enumerate(events) if e["surface"] == "OBJ")
    events.insert(k, {"surface": "WIN", "kind": "transmit"})
    k = max(i for i, e in enumerate(events) if e["surface"] == "OBJ")
    events.insert(k + 1, {"surface": "WIN", "kind": "transmit"})
    return d


@pytest.mark.parametrize("name", ["feature_tour", "linnik"])
def test_w61_section_to_group_keeps_the_system(name: str) -> None:
    from rtt_fusion.graph import model_from_raytatouille, section_to_group
    from rtt_fusion.graph.importer import system_from_model

    if name == "feature_tour":  # Raytatouille's reference system (MIT)
        path = Path(__file__).parent / "data" / "raytatouille" / "feature_tour.rtt.json"
        d = json.loads(path.read_text(encoding="utf-8"))
    else:
        d = _linnik_two_element_objective()
    targets = _sections(d)
    assert targets
    reg, _ = _registry()
    m = model_from_raytatouille(d)
    before = system_from_model(m)
    for t in targets:
        section_to_group(m, t)
        assert m.nodes[t].type_id == "group"
        assert len(m.groups[t].variants[m.groups[t].active].model.nodes) > 1
    assert not [x for x in validate_model(m, reg) if x.severity == "error"]
    again = model_from_dict(json.loads(json.dumps(model_to_dict(m))), reg)  # through the file
    assert system_from_model(again) == before  # the same system, exactly (groups resolved)


def test_w61_section_to_group_keeps_the_image() -> None:
    rt = pytest.importorskip("raytatouille")
    from rtt_fusion.graph import backend_rtt, model_from_raytatouille, section_to_group

    d = _linnik_two_element_objective()
    reg, _ = _registry()

    def image(m: Model) -> Any:
        ex = Executor(
            unroll(m, reg),
            reg,
            backend_rtt.EVALUATORS,
            {"system": rt.System.from_json(json.dumps(d)), "rings": 12, "pixels": [32, 32]},
        )
        return ex.evaluate("camera", "image").values

    plain = image(model_from_raytatouille(d))
    assert plain.max() > 0
    grouped = image(section_to_group(model_from_raytatouille(d), "objective"))
    np.testing.assert_allclose(grouped, plain, atol=1e-12 * plain.max())


def test_w61_refusals() -> None:
    from rtt_fusion.graph import model_from_raytatouille, section_to_group

    m = model_from_raytatouille(
        json.loads(
            (Path(__file__).parent / "data" / "linnik_mini.rtt.json").read_text(encoding="utf-8")
        )
    )
    with pytest.raises(ValueError, match="one element"):
        section_to_group(m, "objective")
    with pytest.raises(ValueError, match="no Raytatouille section"):
        section_to_group(m, "camera")
