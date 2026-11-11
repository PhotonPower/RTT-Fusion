"""File format `rtt-fusion-graph` (concept W7 and L12): versioned JSON, parameters are data.

    {"format": "rtt-fusion-graph", "version": 1, "meta": {…},
     "nodes": [{"id", "type", "name", "params", "position": [x, y], "muted",
                "sockets": [...], "through": {"out": ["in", …]}}],  # dynamic types only
     "links": [{"from": [node, socket], "to": [node, socket], "transport": "auto"}]}

Readers reject other formats and newer versions with a message.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from rtt_fusion.graph.model import Graph, Link, Node, SocketSpec

FORMAT = "rtt-fusion-graph"
VERSION = 1


def to_dict(graph: Graph) -> dict[str, Any]:
    nodes = []
    for n in graph.nodes.values():
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
        nodes.append(d)
    links = []
    for lk in graph.links:
        e: dict[str, Any] = {"from": list(lk.source), "to": list(lk.target)}
        if lk.transport != "auto":
            e["transport"] = lk.transport
        links.append(e)
    out: dict[str, Any] = {"format": FORMAT, "version": VERSION}
    if graph.meta:
        out["meta"] = graph.meta
    out["nodes"] = nodes
    out["links"] = links
    return out


def from_dict(d: dict[str, Any]) -> Graph:
    if d.get("format") != FORMAT:
        raise ValueError(f"not a {FORMAT} document (format {d.get('format')!r})")
    version = int(d.get("version", 0))
    if version < 1 or version > VERSION:
        raise ValueError(f"{FORMAT} version {version} is not supported (this reader: {VERSION})")
    g = Graph(meta=dict(d.get("meta", {})))
    for nd in d.get("nodes", []):
        pos = nd.get("position", [0.0, 0.0])
        g.add(
            Node(
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
            )
        )
    for ld in d.get("links", []):
        g.links.append(
            Link(
                (str(ld["from"][0]), str(ld["from"][1])),
                (str(ld["to"][0]), str(ld["to"][1])),
                str(ld.get("transport", "auto")),
            )
        )
    return g


def dumps(graph: Graph) -> str:
    return json.dumps(to_dict(graph), indent=2, ensure_ascii=False)


def loads(text: str) -> Graph:
    return from_dict(json.loads(text))


def save(graph: Graph, path: str | Path) -> None:
    Path(path).write_text(dumps(graph) + "\n", encoding="utf-8")


def load(path: str | Path) -> Graph:
    return loads(Path(path).read_text(encoding="utf-8"))
