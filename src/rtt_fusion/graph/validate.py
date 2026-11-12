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

from rtt_fusion.graph.model import (
    Category,
    Direction,
    Graph,
    GraphError,
    Node,
    SocketKind,
    SocketSpec,
)
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
        extra = t.extra_sockets if t is not None else ()
        names = {s.name for s in node.sockets}
        own = tuple(node.sockets) + tuple(s for s in extra if s.name not in names)
    else:
        own = t.sockets if t is not None else ()
    if t is not None and t.category in _POSED and not any(s.name == "pose" for s in own):
        own = (*own, POSE_INPUT)  # every block with a pose has the offset input (R26)
    return own


#: the offset input of every block with a pose (R26): six values, driven by meta blocks
POSE_INPUT = SocketSpec("pose", Direction.IN, SocketKind.DATA, frozenset({"drive"}), "", "pose")
_POSED = (Category.SOURCE, Category.DOMAIN, Category.SPLITTER, Category.DETECTOR)


def node_through(node: Node, registry: NodeTypeRegistry) -> dict[str, tuple[str, ...]]:
    """Internal dependencies output -> inputs; default: every output on every light input."""
    t = registry.get(node.type_id)
    pairs = node.through or (t.through if t is not None else ())
    if pairs:
        return dict(pairs)
    sockets = node_sockets(node, registry)
    ins = tuple(s.name for s in sockets if s.direction is Direction.IN)
    return {s.name: ins for s in sockets if s.direction is Direction.OUT}


def replace_block(graph: Graph, registry: NodeTypeRegistry, node_id: str, new: Node) -> Node:
    """`Graph.replace_node` with the socket map found from the sockets: linked sockets of the old
    block move to the new block's socket of the same name, else to the only new socket of the same
    direction and kind (a mirror's `from_objective`/`to_objective` -> a grating's `in`/`out`).
    Raises GraphError if a linked socket has no unique counterpart."""
    old = graph.node(node_id)
    new_sockets = node_sockets(new, registry)
    names = {s.name for s in new_sockets}
    mapping: dict[str, str] = {}
    linked = {lk.source[1] for lk in graph.links_from(node_id)} | {
        lk.target[1] for lk in graph.links_to(node_id)
    }
    for spec in node_sockets(old, registry):
        if spec.name not in linked or spec.name in names:
            continue
        same = [
            s.name for s in new_sockets if s.direction is spec.direction and s.kind is spec.kind
        ]
        if len(same) != 1:
            raise GraphError(
                f"no unique socket of {new.type_id} for {old.id}.{spec.name} ({len(same)} fit)"
            )
        mapping[spec.name] = same[0]
    if len(set(mapping.values())) != len(mapping):
        raise GraphError(f"two sockets of {old.id} would share one socket of {new.type_id}")
    return graph.replace_node(node_id, new, mapping)


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
                if (
                    name not in known
                    and name not in ("pose", "pose_offset")
                    and not (name.startswith("pose_by_"))
                ):
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
