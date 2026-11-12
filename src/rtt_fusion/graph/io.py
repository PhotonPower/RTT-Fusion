"""File format `rtt-fusion-graph` (concept W7 and L12): versioned JSON, parameters are data.

    {"format": "rtt-fusion-graph", "version": 1, "meta": {…},
     "nodes": [{"id", "type", "name", "params", "position": [x, y], "muted",
                "sockets": [...], "through": {"out": ["in", …]}}],  # dynamic types only
     "links": [{"from": [node, socket], "to": [node, socket], "transport": "auto"}]}

Version 2 (ADR ports, poses, branches; R24) is the model: topology and branches.

    {"format": "rtt-fusion-graph", "version": 2, "meta": {…},
     "nodes": [... as in version 1, plus "ports": [{"name", "in", "out", "side"}] ...],
     "connections": [{"ports": [[node, port], [node, port]],
                      "transport": "auto" | "R" | {"forward": "R", "back": "H"}}],
     "links": [... data and observer links as in version 1, "branch": name for an observer ...],
     "branches": [{"name": "...", "steps": [[node, port in | null, port out | null, channel?]]}]}

A group (R28) is a node of type "group" with "variant": the active one and "variants": {name:
{nodes, connections, links, branches (the inner branches), "ports": {port: [node, port]},
"params": {name: "node.param[i]"}}}; a step through it names the inner branch as its channel.

`to_dict`/`from_dict` read and write the socket form (version 1); `model_to_dict` writes
version 2, `model_from_dict` reads both (version 1 via `topology.model_from_graph`). Readers reject
other formats and newer versions with a message.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from rtt_fusion.graph.model import Graph, Link, Node, PortSpec, SocketSpec
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.graph.topology import (
    Group,
    Model,
    Variant,
    branch_from_dict,
    branch_to_dict,
    model_from_graph,
)

FORMAT = "rtt-fusion-graph"
VERSION = 1  # socket form
MODEL_VERSION = 2


def _node_to_dict(n: Node) -> dict[str, Any]:
    d: dict[str, Any] = {"id": n.id, "type": n.type_id}
    if n.name:
        d["name"] = n.name
    if n.params:
        d["params"] = n.params
    if n.position != (0.0, 0.0):
        d["position"] = [n.position[0], n.position[1]]
    if n.muted:
        d["muted"] = True
    if n.sockets:
        d["sockets"] = [s.to_dict() for s in n.sockets]
    if n.through:
        d["through"] = {o: list(ins) for o, ins in n.through}
    if n.ports:
        d["ports"] = [p.to_dict() for p in n.ports]
    return d


def _node_from_dict(nd: dict[str, Any]) -> Node:
    pos = nd.get("position", [0.0, 0.0])
    return Node(
        id=str(nd["id"]),
        type_id=str(nd["type"]),
        name=str(nd.get("name", "")),
        params=dict(nd.get("params", {})),
        position=(float(pos[0]), float(pos[1])),
        muted=bool(nd.get("muted", False)),
        sockets=tuple(SocketSpec.from_dict(s) for s in nd.get("sockets", [])),
        through=tuple(
            (str(o), tuple(str(i) for i in ins)) for o, ins in nd.get("through", {}).items()
        ),
        ports=tuple(PortSpec.from_dict(p) for p in nd.get("ports", [])),
    )


def _link_to_dict(lk: Link) -> dict[str, Any]:
    e: dict[str, Any] = {"from": list(lk.source), "to": list(lk.target)}
    if lk.transport != "auto":
        e["transport"] = lk.transport
    return e


def _link_from_dict(ld: dict[str, Any]) -> Link:
    return Link(
        (str(ld["from"][0]), str(ld["from"][1])),
        (str(ld["to"][0]), str(ld["to"][1])),
        str(ld.get("transport", "auto")),
    )


def _check_header(d: dict[str, Any], newest: int) -> int:
    if d.get("format") != FORMAT:
        raise ValueError(f"not a {FORMAT} document (format {d.get('format')!r})")
    version = int(d.get("version", 0))
    if version < 1 or version > newest:
        raise ValueError(f"{FORMAT} version {version} is not supported (this reader: {newest})")
    return version


def _group_to_dict(group: Group) -> dict[str, Any]:
    variants = {}
    for name, v in group.variants.items():
        body = model_to_dict(v.model)
        del body["format"], body["version"]
        body["ports"] = {p: list(end) for p, end in v.ports.items()}
        if v.params:
            body["params"] = dict(v.params)
        variants[name] = body
    return {"variant": group.active, "variants": variants}


def _group_from_dict(d: dict[str, Any], registry: NodeTypeRegistry) -> Group:
    variants = {}
    for name, body in d["variants"].items():
        inner = model_from_dict({"format": FORMAT, "version": MODEL_VERSION, **body}, registry)
        ports = {str(p): (str(e[0]), str(e[1])) for p, e in body.get("ports", {}).items()}
        variants[str(name)] = Variant(
            inner, ports, {str(k): str(v) for k, v in body.get("params", {}).items()}
        )
    return Group(variants, str(d["variant"]))


def model_to_dict(model: Model) -> dict[str, Any]:
    """Version 2: topology and branches; a group node carries its variants (R28)."""
    out: dict[str, Any] = {"format": FORMAT, "version": MODEL_VERSION}
    if model.meta:
        out["meta"] = model.meta
    out["nodes"] = []
    for n in model.nodes.values():
        nd = _node_to_dict(n)
        if n.id in model.groups:
            nd.update(_group_to_dict(model.groups[n.id]))
        out["nodes"].append(nd)
    connections = []
    for c in model.connections:
        e: dict[str, Any] = {"ports": [list(c.a), list(c.b)]}
        if c.transport_back is not None and c.transport_back != c.transport:
            e["transport"] = {"forward": c.transport, "back": c.transport_back}
        elif c.transport != "auto":
            e["transport"] = c.transport
        connections.append(e)
    out["connections"] = connections
    links = []
    for lk in model.links:
        e = _link_to_dict(lk)
        if lk.target in model.observe_branch:
            e["branch"] = model.observe_branch[lk.target]
        links.append(e)
    out["links"] = links
    out["branches"] = [branch_to_dict(b) for b in model.branches]
    return out


def model_from_dict(d: dict[str, Any], registry: NodeTypeRegistry) -> Model:
    """Version 2, or version 1 converted (`topology.model_from_graph`)."""
    if _check_header(d, MODEL_VERSION) == 1:
        return model_from_graph(from_dict(d), registry)
    m = Model(meta=dict(d.get("meta", {})))
    for nd in d.get("nodes", []):
        m.add(_node_from_dict(nd))
        if "variants" in nd:
            m.groups[str(nd["id"])] = _group_from_dict(nd, registry)
    for cd in d.get("connections", []):
        a, b = cd["ports"]
        t = cd.get("transport", "auto")
        forward, back = (t["forward"], t["back"]) if isinstance(t, dict) else (t, None)
        m.connect((str(a[0]), str(a[1])), (str(b[0]), str(b[1])), str(forward), back)
    for ld in d.get("links", []):
        lk = _link_from_dict(ld)
        m.links.append(lk)
        if "branch" in ld:
            m.observe_branch[lk.target] = str(ld["branch"])
    m.branches = [branch_from_dict(b) for b in d.get("branches", [])]
    return m


def to_dict(graph: Graph) -> dict[str, Any]:
    """The socket form (version 1)."""
    out: dict[str, Any] = {"format": FORMAT, "version": VERSION}
    if graph.meta:
        out["meta"] = graph.meta
    out["nodes"] = [_node_to_dict(n) for n in graph.nodes.values()]
    out["links"] = [_link_to_dict(lk) for lk in graph.links]
    return out


def from_dict(d: dict[str, Any]) -> Graph:
    """The socket form (version 1)."""
    if d.get("format") == FORMAT and int(d.get("version", 0)) == MODEL_VERSION:
        raise ValueError(
            f"{FORMAT} version 2 is a model: read it with model_from_dict and topology.unroll"
        )
    _check_header(d, VERSION)
    g = Graph(meta=dict(d.get("meta", {})))
    for nd in d.get("nodes", []):
        g.add(_node_from_dict(nd))
    g.links = [_link_from_dict(ld) for ld in d.get("links", [])]
    return g


def dumps(graph: Graph) -> str:
    return json.dumps(to_dict(graph), indent=2, ensure_ascii=False)


def loads(text: str) -> Graph:
    return from_dict(json.loads(text))


def save(graph: Graph, path: str | Path) -> None:
    Path(path).write_text(dumps(graph) + "\n", encoding="utf-8")


def load(path: str | Path) -> Graph:
    return loads(Path(path).read_text(encoding="utf-8"))
