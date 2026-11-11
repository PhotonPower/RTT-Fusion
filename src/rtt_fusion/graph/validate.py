"""Checks of a working graph (concept sections 3 and 4) with diagnostic codes.

Codes (severity): graph.node_type_unknown, graph.socket_unknown, graph.link_direction,
graph.kind_mismatch, graph.data_type_mismatch, graph.no_common_representation,
graph.transport_unsupported, graph.input_multiple (errors); graph.cycle (error);
graph.param_unknown, graph.detector_unconnected (warnings).

Cycles are looked for on the level of sockets: links plus each node's internal dependencies
(`through`), so an objective passed in both directions is no cycle (concept section 4); a real
feedback (multiple reflection, resonator) is one and belongs into an iteration meta block.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from rtt_fusion.graph.model import Category, Direction, Graph, Node, SocketKind, SocketSpec
from rtt_fusion.graph.registry import NodeTypeRegistry


@dataclass(frozen=True)
class Diagnostic:
    code: str
    severity: str  # "error" or "warning"
    message: str
    location: str  # node id, "node.socket" or "a.s -> b.t"


def node_sockets(node: Node, registry: NodeTypeRegistry) -> tuple[SocketSpec, ...]:
    """The sockets of a node: its own for types with dynamic sockets, else the type's."""
    t = registry.get(node.type_id)
    if node.sockets:
        return node.sockets
    return t.sockets if t is not None else ()


def node_through(node: Node, registry: NodeTypeRegistry) -> dict[str, tuple[str, ...]]:
    """Internal dependencies output -> inputs; default: every output on every light input."""
    t = registry.get(node.type_id)
    pairs = node.through or (t.through if t is not None else ())
    if pairs:
        return dict(pairs)
    sockets = node_sockets(node, registry)
    ins = tuple(s.name for s in sockets if s.direction is Direction.IN)
    return {s.name: ins for s in sockets if s.direction is Direction.OUT}


def validate(graph: Graph, registry: NodeTypeRegistry) -> list[Diagnostic]:
    out: list[Diagnostic] = []

    def err(code: str, msg: str, loc: str) -> None:
        out.append(Diagnostic(code, "error", msg, loc))

    def warn(code: str, msg: str, loc: str) -> None:
        out.append(Diagnostic(code, "warning", msg, loc))

    sockets: dict[str, dict[str, SocketSpec]] = {}
    for node in graph.nodes.values():
        t = registry.get(node.type_id)
        if t is None:
            err("graph.node_type_unknown", f"unknown block type {node.type_id!r}", node.id)
            sockets[node.id] = {s.name: s for s in node.sockets}
            continue
        sockets[node.id] = {s.name: s for s in node_sockets(node, registry)}
        known = {p.name for p in t.params}
        if known and not t.dynamic_sockets:
            for name in node.params:
                if name not in known:
                    warn(
                        "graph.param_unknown",
                        f"parameter {name!r} is not one of {t.type_id}",
                        node.id,
                    )

    for lk in graph.links:
        loc = f"{lk.source[0]}.{lk.source[1]} -> {lk.target[0]}.{lk.target[1]}"
        a = sockets.get(lk.source[0], {}).get(lk.source[1])
        b = sockets.get(lk.target[0], {}).get(lk.target[1])
        if a is None or b is None:
            missing = lk.source if a is None else lk.target
            err("graph.socket_unknown", f"no socket {missing[1]!r} on {missing[0]!r}", loc)
            continue
        if a.direction is not Direction.OUT or b.direction is not Direction.IN:
            err("graph.link_direction", "a link runs from an output to an input", loc)
            continue
        if a.kind is not b.kind:
            err(
                "graph.kind_mismatch",
                f"{a.kind.value} socket linked to a {b.kind.value} socket",
                loc,
            )
            continue
        if a.kind is SocketKind.DATA:
            if a.accepts != b.accepts:
                err(
                    "graph.data_type_mismatch",
                    f"data type {next(iter(a.accepts))} does not fit {next(iter(b.accepts))}",
                    loc,
                )
            continue
        common = a.accepts & b.accepts
        if not common:
            err(
                "graph.no_common_representation",
                f"no common representation ({''.join(sorted(a.accepts))} / "
                f"{''.join(sorted(b.accepts))});"
                " insert a conversion",
                loc,
            )
        elif lk.transport != "auto" and lk.transport not in common:
            err(
                "graph.transport_unsupported",
                f"transport {lk.transport} not accepted by both sockets "
                f"({''.join(sorted(common))})",
                loc,
            )

    for (node_id, socket), n in Counter(lk.target for lk in graph.links).items():
        if n > 1:
            err(
                "graph.input_multiple",
                f"{n} links into one input; use a combining block",
                f"{node_id}.{socket}",
            )

    for node in graph.nodes.values():
        t = registry.get(node.type_id)
        if t is not None and t.category is Category.DETECTOR and not graph.links_to(node.id):
            warn("graph.detector_unconnected", "detector without light", node.id)

    cycle = find_cycle(graph, registry)
    if cycle:
        err(
            "graph.cycle",
            "feedback loop: " + " -> ".join(cycle) + " (use an iteration block)",
            cycle[0],
        )
    return out


def find_cycle(graph: Graph, registry: NodeTypeRegistry) -> list[str]:
    """A cycle on the socket level as a list of "node.socket", or [] if the graph is acyclic."""
    edges: dict[str, list[str]] = {}
    for lk in graph.links:
        edges.setdefault(f"{lk.source[0]}.{lk.source[1]}", []).append(
            f"{lk.target[0]}.{lk.target[1]}"
        )
    for node in graph.nodes.values():
        for out_name, ins in node_through(node, registry).items():
            for i in ins:
                edges.setdefault(f"{node.id}.{i}", []).append(f"{node.id}.{out_name}")
    state: dict[str, int] = {}  # 1: on the stack, 2: done
    stack: list[str] = []

    def visit(v: str) -> list[str]:
        state[v] = 1
        stack.append(v)
        for w in edges.get(v, ()):
            if state.get(w) == 1:
                return [*stack[stack.index(w) :], w]
            if w not in state:
                found = visit(w)
                if found:
                    return found
        stack.pop()
        state[v] = 2
        return []

    for v in list(edges):
        if v not in state:
            found = visit(v)
            if found:
                return found
    return []


def errors(diagnostics: list[Diagnostic]) -> list[Diagnostic]:
    return [d for d in diagnostics if d.severity == "error"]
