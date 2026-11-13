"""The workflow document a front end shows and edits (CAD ADR 0010, stage W1).

What the node editor of Raytatouille Fusion receives from the workflow engine – and from a plugin
that makes a workflow, e.g. by importing a system – is one JSON-like document: the model (every
block with its sockets and ports, so the editor needs no type knowledge), the block library, the
diagnostics of the model and of the unrolled graph and the plan, per link of the unrolled graph and
per connection of the model. This module builds it; it knows no kernel: the editable data fields of
a kernel's blocks come from a callback of the plugin that provides them (ADR "Bidirektionale
Ports, Posen, Äste und Gruppen", Nachtrag 1: the plugin describes what its blocks hold).
"""

from __future__ import annotations

import re
from collections.abc import Callable
from importlib import import_module
from typing import Any

from rtt_fusion.graph import io
from rtt_fusion.graph.execute import EXECUTABLE
from rtt_fusion.graph.model import Graph, Node
from rtt_fusion.graph.planner import Requirements, plan
from rtt_fusion.graph.registry import NodeTypeRegistry, builtin_registry
from rtt_fusion.graph.topology import Model, node_ports, unroll, validate_model
from rtt_fusion.graph.validate import node_sockets, validate

#: modules with block types beyond the built-in ones: grating (R19), particle (R20), meta blocks
#: and displays (R21); each registers what its optional packages allow
OPTIONAL_MODULES = (
    "rtt_fusion.graph.grating",
    "rtt_fusion.graph.particle",
    "rtt_fusion.graph.meta",
    "rtt_fusion.graph.display",
)

#: the editable data fields of one node (a dict as in the model's file) – the providing plugin's
#: description; empty for a node without such data
Describe = Callable[[dict[str, Any]], list[dict[str, Any]]]


def registry() -> tuple[NodeTypeRegistry, dict[str, Any]]:
    """The block types of the workflow (built-in and optional modules) and the evaluators of the
    optional modules."""
    reg = builtin_registry()
    extra: dict[str, Any] = {}
    for module in OPTIONAL_MODULES:
        try:
            extra.update(import_module(module).register(reg))
        except ImportError:
            pass
    return reg, extra


def available() -> frozenset[str]:
    """Representations the workflow can run: R, and H (fields, R18b) with fieldtracing."""
    execute = import_module("rtt_fusion.graph.execute")
    try:
        import_module("fieldtracing")
    except ImportError:
        return EXECUTABLE
    return frozenset(getattr(execute, "WITH_FIELDS", EXECUTABLE))


def catalogue(reg: NodeTypeRegistry) -> list[dict[str, Any]]:
    """The block library for the editor: type, category, description, sockets, parameters,
    through-paths and the light ports (one dot per port, R24)."""
    out = []
    for t in reg.types():
        out.append(
            {
                "type": t.type_id,
                "category": t.category.value,
                "description": t.description,
                "sockets": [s.to_dict() for s in t.sockets],
                "params": [p.to_dict() for p in t.params],
                "through": {o: list(ins) for o, ins in t.through},
                "ports": [p.to_dict() for p in node_ports(Node("x", t.type_id), reg)],
            }
        )
    return out


def checked(graph: Graph, reg: NodeTypeRegistry) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Diagnostics of the graph check and of the planner (with what can run), and the plan."""
    diagnostics = [dict(d.__dict__) for d in validate(graph, reg)]
    p = plan(graph, reg, Requirements(available=available()))
    diagnostics += [dict(d.__dict__) for d in p.diagnostics]
    return diagnostics, p.to_dict()


def model_id(unrolled: str) -> str:
    """The block of the model behind a block of the unrolled graph: a copy "x@2" is "x", a block
    "g/x" inside a group is the group "g"."""
    return unrolled.split("@")[0].split("/")[0]


def connection_plan(
    graph: Graph, model: Model, reg: NodeTypeRegistry, p: dict[str, Any]
) -> dict[str, dict[str, str]]:
    """The planned representation per connection of the model and direction ("forward" a -> b,
    "back")."""
    out: dict[str, dict[str, str]] = {}
    ports = {nid: node_ports(n, reg) for nid, n in model.nodes.items()}
    for lk in graph.links:
        key = f"{lk.source[0]}.{lk.source[1]} -> {lk.target[0]}.{lk.target[1]}"
        if key not in p.get("links", {}):
            continue
        a, b = model_id(lk.source[0]), model_id(lk.target[0])
        pa = next((q.name for q in ports.get(a, ()) if q.out_socket == lk.source[1]), None)
        pb = next((q.name for q in ports.get(b, ()) if q.in_socket == lk.target[1]), None)
        for i, c in enumerate(model.connections):
            if {c.a, c.b} == {(a, pa), (b, pb)}:
                side = "forward" if c.a == (a, pa) else "back"
                out.setdefault(str(i), {})[side] = p["links"][key]["representation"]
    return out


def model_document(
    model: Model,
    reg: NodeTypeRegistry,
    describe: Describe | None = None,
    keep_fields: dict[str, list[dict[str, Any]]] | None = None,
    providers: dict[str, dict[str, str]] | None = None,
    recorded: dict[str, dict[str, str]] | None = None,
) -> dict[str, Any]:
    """The workflow document of a model: {"model", "types", "diagnostics", "plan", "providers"}.
    The fields of a node come from `describe` (the providing plugin), else from `keep_fields` (by
    node id: what the document held before – the engine keeps what the providers described).

    `providers`: block type → {"plugin", "version"} of the types other plugins provide now
    (CAD ADR 0010, D7); `recorded`: the same from the document before. The document records the
    provider of every such type its model uses; a block whose type nobody provides but whose
    provider is recorded gets `graph.provider_missing` (the plugin is missing or too old) instead
    of `graph.node_type_unknown` – it stays in the model with its data."""
    for node in model.nodes.values():  # with the inputs every block has (drive, pose)
        node.sockets = tuple(node_sockets(node, reg))
        if not node.ports:
            node.ports = node_ports(node, reg)
    diagnostics = [dict(d.__dict__) for d in validate_model(model, reg)]
    p: dict[str, Any] = {"links": {}, "conversions": [], "connections": {}}
    if not any(d["severity"] == "error" for d in diagnostics):
        graph = unroll(model, reg)
        more, p = checked(graph, reg)
        diagnostics += more
        p["connections"] = connection_plan(graph, model, reg, p)
    for d in diagnostics:  # copies and blocks inside groups: the block of the model
        d["location"] = re.sub(r"@\d+", "", str(d["location"]))
    now, before = providers or {}, recorded or {}
    used: set[str] = set()

    def types_of(m: Model) -> None:
        for n in m.nodes.values():
            used.add(n.type_id)
        for g in m.groups.values():
            for v in g.variants.values():
                types_of(v.model)

    types_of(model)
    record = {t: dict(now.get(t) or before[t]) for t in sorted(used) if t in now or t in before}
    for d in diagnostics:
        if d["code"] != "graph.node_type_unknown":
            continue
        t = next((x for x in record if f"{x!r}" in str(d["message"])), None)
        if t is not None and t not in now:
            r = record[t]
            d["code"] = "graph.provider_missing"
            d["message"] = (
                f"block type {t!r} comes from the plugin {r.get('plugin', '?')} "
                f"(version {r.get('version', '?')}), which is not available – the block keeps its "
                "data; install or start the plugin"
            )
    data = io.model_to_dict(model)
    kept = keep_fields or {}

    def fields(nodes: list[dict[str, Any]]) -> None:  # also inside groups (every variant)
        for nd in nodes:
            got = describe(nd) if describe is not None else kept.get(str(nd.get("id")), [])
            if got:
                nd["fields"] = got
            for body in (nd.get("variants") or {}).values():
                fields(body.get("nodes", []))

    fields(data["nodes"])
    types = catalogue(reg)
    for entry in types:  # the editor groups the block library by plugin
        if entry["type"] in now:
            entry["provider"] = now[entry["type"]].get("plugin", "")
    return {
        "model": data,
        "types": types,
        "diagnostics": diagnostics,
        "plan": p,
        "providers": record,
    }


def kept_fields(document: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """The fields a document holds, by node id (also inside groups)."""
    out: dict[str, list[dict[str, Any]]] = {}

    def walk(nodes: list[dict[str, Any]]) -> None:
        for nd in nodes:
            if nd.get("fields"):
                out[str(nd.get("id"))] = list(nd["fields"])
            for body in (nd.get("variants") or {}).values():
                walk(body.get("nodes", []))

    walk((document.get("model") or {}).get("nodes", []))
    return out


def block_states(raw: dict[str, bool], model_ids: set[str] | None) -> dict[str, bool]:
    """Up to date per block of the unrolled graph → per block of the model (all its copies and the
    blocks of a group must be up to date); unchanged without a model (an old socket graph)."""
    if model_ids is None:
        return dict(raw)
    merged: dict[str, bool] = {}
    for k, v in raw.items():
        m = model_id(k)
        if m in model_ids:
            merged[m] = merged.get(m, True) and v
    return merged


__all__ = [
    "OPTIONAL_MODULES",
    "Describe",
    "available",
    "block_states",
    "catalogue",
    "checked",
    "connection_plan",
    "kept_fields",
    "model_document",
    "model_id",
    "registry",
]
