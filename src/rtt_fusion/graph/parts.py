"""New Raytatouille parts from the library (ADR ports, poses, branches, Nachtrag 1, stage S4).

A part is a Raytatouille element made in the workflow, not imported: its element (type, pose,
surfaces with shape, aperture, interaction, glass) is the node's parameter `rtt_new`, editable
like the data of an imported block. Its ports are fixed per kind; the events of a traversal follow
from the ports – leaving by the side the light came in is a reflection, by the other side a
transmission. Before unrolling (`prepare_parts`) each part gets its place in the system (a child
of the root after the imported ones), surface ids unique in the system (`<block>.<surface>`), its
`surfaces` and its traversals (`runs`); from then on it is a Raytatouille block like an imported
one: the generated system holds it, a branch through it is a path of that system.
"""

from __future__ import annotations

import json
from dataclasses import replace
from typing import Any

from rtt_fusion.graph.model import (
    Category,
    Node,
    NodeType,
    ParamSpec,
    PortSpec,
    SocketSpec,
    data_out,
    light_in,
    light_out,
)

_REPS = "RH"
_APERTURE = {"type": "circular", "radius": 12.5}
_POSE = {"position": [0.0, 0.0, 0.0], "rotation_deg": [0.0, 0.0, 0.0]}


def _port_sockets(names: tuple[str, ...]) -> tuple[tuple[SocketSpec, ...], tuple[PortSpec, ...]]:
    sockets: tuple[SocketSpec, ...] = ()
    ports: tuple[PortSpec, ...] = ()
    for name in names:
        side = "back" if name.startswith("back") else "front"
        sockets += (light_in(f"{name}_in", _REPS, side), light_out(f"{name}_out", _REPS, side))
        ports += (PortSpec(name, f"{name}_in", f"{name}_out", side),)
    return sockets, ports


#: kind -> (ports, element template, category)
PARTS: dict[str, tuple[tuple[str, ...], dict[str, Any], Category]] = {
    "mirror": (
        ("front",),
        {
            "type": "mirror",
            "name": "Mirror",
            "pose": _POSE,
            "surfaces": [
                {"id": "M", "aperture": _APERTURE, "interaction": {"type": "ideal_mirror"}}
            ],
        },
        Category.DOMAIN,
    ),
    "splitter": (
        ("front_a", "front_b", "back_a", "back_b"),
        {
            "type": "thin_element",
            "name": "Beam splitter",
            "pose": {"position": [0.0, 0.0, 0.0], "rotation_deg": [45.0, 0.0, 0.0]},
            "surfaces": [
                {"id": "BS", "aperture": _APERTURE, "interaction": {"type": "ideal_beam_splitter"}}
            ],
        },
        Category.SPLITTER,
    ),
    "lens": (
        ("front", "back"),
        {
            "type": "lens",
            "name": "Lens",
            "pose": _POSE,
            "material": "CONST:1.5168",
            "surfaces": [
                {
                    "id": "S1",
                    "shape": {"base": {"type": "conic", "radius": 51.68}},
                    "aperture": _APERTURE,
                },
                {
                    "id": "S2",
                    "pose": {"position": [0.0, 0.0, 4.0]},
                    "shape": {"base": {"type": "conic", "radius": -51.68}},
                    "aperture": _APERTURE,
                },
            ],
        },
        Category.DOMAIN,
    ),
    "detector": (
        ("in",),
        {
            "type": "detector",
            "name": "Detector",
            "pose": _POSE,
            "surfaces": [
                {
                    "id": "IMG",
                    "aperture": {"type": "rectangular", "half_width_x": 5.0, "half_width_y": 5.0},
                }
            ],
        },
        Category.DETECTOR,
    ),
}


def part_types() -> list[NodeType]:
    """The block types of the parts ("rtt.part.<kind>")."""
    out = []
    for kind, (ports, template, category) in PARTS.items():
        if kind == "detector":
            sockets: tuple[SocketSpec, ...] = (
                light_in("in_in", _REPS, "front"),
                data_out("image", "result", "image"),
            )
            port_specs: tuple[PortSpec, ...] = (PortSpec("in", "in_in", "", "front"),)
            through: tuple[tuple[str, tuple[str, ...]], ...] = (("image", ("in_in",)),)
        else:
            sockets, port_specs = _port_sockets(ports)
            through = ()
        out.append(
            NodeType(
                f"rtt.part.{kind}",
                category,
                sockets,
                (
                    ParamSpec(
                        "rtt_new",
                        "data",
                        template,
                        description="the Raytatouille element (pose, surfaces, glass)",
                    ),
                ),
                through,
                ports=port_specs,
                description=f"Raytatouille {template['name'].lower()} (a new element)",
            )
        )
    return out


def part_events(kind: str, surfaces: list[str], port_in: str, port_out: str) -> list[Any] | None:
    """The events of a traversal of a part, None if there is none."""
    if kind == "detector":
        return [{"surface": surfaces[0], "kind": "transmit"}] if port_out == "" else None
    if kind == "mirror":
        return [{"surface": surfaces[0], "kind": "reflect"}] if port_in == port_out else None
    if kind == "lens":
        if {port_in, port_out} != {"front", "back"}:
            return None
        order = surfaces if port_in == "front" else list(reversed(surfaces))
        return [{"surface": s} for s in order]  # refracted (Raytatouille's default for a lens)
    if kind == "splitter":
        side_in, line_in = port_in.split("_")
        side_out, line_out = port_out.split("_") if port_out else ("", "")
        if side_in == side_out and line_in != line_out:
            return [{"surface": surfaces[0], "kind": "reflect"}]
        if side_in != side_out and line_in == line_out:
            return [{"surface": surfaces[0], "kind": "transmit"}]
    return None


def is_part(node: Node) -> bool:
    return node.type_id.startswith("rtt.part.") and "rtt_new" in node.params


def prepare_parts(nodes: dict[str, Node]) -> dict[str, Node]:
    """The nodes with every part placed in the system: a location after the root's children,
    unique surface ids, `rtt`, `surfaces`, `location`, `element_type` and `runs` – a part is then a
    Raytatouille block like an imported one. Other nodes unchanged (the same objects)."""
    parts = [n for n in nodes.values() if is_part(n)]
    if not parts:
        return nodes
    taken = 0
    for n in nodes.values():
        for location in list(n.params.get("rtt", {})) + list(n.params.get("rtt_tree", {})):
            pieces = location.split("/")
            if len(pieces) == 4 and pieces[2] == "children":  # /root/children/i
                taken = max(taken, int(pieces[3]) + 1)
    out = dict(nodes)
    for k, node in enumerate(parts):
        kind = node.type_id.rsplit(".", 1)[-1]
        element: dict[str, Any] = json.loads(json.dumps(node.params["rtt_new"]))
        for s in element.get("surfaces", []):
            s["id"] = f"{node.id}.{s['id']}"
        surfaces = [str(s["id"]) for s in element.get("surfaces", [])]
        location = f"/root/children/{taken + k}"
        ports = PARTS[kind][0]
        runs: dict[str, Any] = {}
        for p_in in ports:
            for p_out in (*ports, ""):
                events = part_events(kind, surfaces, p_in, p_out)
                if events is not None:
                    runs[f"{p_in}>{p_out}|"] = events
        params = {
            **node.params,
            "rtt": {location: element},
            "surfaces": surfaces,
            "location": location,
            "element_type": str(element.get("type", kind)),
            "runs": runs,
        }
        out[node.id] = replace(node, params=params)
    return out


__all__ = ["PARTS", "is_part", "part_events", "part_types", "prepare_parts"]
