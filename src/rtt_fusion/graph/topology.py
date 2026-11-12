"""Ports, connections and branches (ADR ports, poses, branches and groups; stage R24).

The model has three levels (ADR section 1):

- **topology**: nodes with light ports (two channels each, `PortSpec`), undirected connections
  between two ports (`Connection`), data and observer links (`Link`, directed as before);
- **branches**: finite lists of steps (node, port in, port out, channel) from a source to a
  detector (`Branch`). They are given by the modeller only (decision P6): exactly the branches in
  the model are computed; nothing here enumerates, adds or drops a way, there are no limits and no
  defaults, no loops, no series and no feedback;
- the **execution plan** works on the socket graph (`Graph`) that `unroll` derives from the
  branches; planner stage 2, execution, result cache and coherence are unchanged.

Unrolling: the branches are words over steps. Their prefix tree, merged by equal remainders
(the minimal automaton of the finite set of words), accepts exactly the branches. Every
transition is the visit of a node; the visits of a node become as few copies of it (`id`,
`id@2`, …) as the sockets allow – an input socket of a copy takes light of one state only, an
output socket gives light to one state only. Light at a state is the sum of what arrives there;
where it arrives from several copies a junction block (`fusion.sum`) adds it. The socket paths
from the sources to the detectors are then exactly the branches (multiple visits included), and
the socket graph is acyclic because every link leads to a later state.

Channels (diffraction orders, …) are carried with the steps; in this stage one socket carries all
channels of a pair of ports (as in version 1) – which channels a node gives is its parameter, the
Raytatouille backend traces its paths by name.

Connection checks (R26, `connection_diagnostics`): graph.port_facing, graph.port_behind,
graph.port_medium (errors); graph.port_aperture, graph.port_near_field, graph.port_tilt (warnings).

Diagnostic codes (errors): graph.no_branch, graph.branch_impossible, graph.branch_duplicate,
graph.port_unknown, graph.port_links, graph.group_variant, graph.group_ports,
graph.group_params, graph.group_branch (groups, R28), graph.branch_backend (a step through a
Raytatouille block that no path of the system makes), graph.observer_ambiguous.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

from rtt_fusion.graph.model import (
    TRANSPORTS,
    Category,
    Direction,
    Graph,
    Link,
    Node,
    PortFrame,
    PortSpec,
    SocketKind,
    SocketSpec,
    light_in,
    light_out,
)
from rtt_fusion.graph.pose import posed
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.graph.validate import Diagnostic, node_sockets, node_through

#: block type of the junction that adds the light of several node copies at one state
SUM_TYPE = "fusion.sum"


class ModelError(ValueError):
    """A model that cannot be unrolled; `diagnostics` says why."""

    def __init__(self, diagnostics: list[Diagnostic]) -> None:
        self.diagnostics = diagnostics
        super().__init__("; ".join(f"{d.code} at {d.location}: {d.message}" for d in diagnostics))


@dataclass(frozen=True)
class Connection:
    """An undirected connection of two ports; it carries both directions (ADR section 5).

    `transport` holds for the direction a -> b, `transport_back` for b -> a (None: the same)."""

    a: tuple[str, str]  # (node id, port name)
    b: tuple[str, str]
    transport: str = "auto"
    transport_back: str | None = None

    def __post_init__(self) -> None:
        for t in (self.transport, self.transport_back or "auto"):
            if t not in TRANSPORTS:
                raise ValueError(f"transport {t!r} not in {sorted(TRANSPORTS)}")
        if self.a == self.b:
            raise ValueError(f"a connection joins two ports, not {self.a} with itself")

    def other(self, end: tuple[str, str]) -> tuple[str, str]:
        return self.b if end == self.a else self.a

    def transport_from(self, end: tuple[str, str]) -> str:
        """The transport of the direction leaving the port `end`."""
        if end == self.a:
            return self.transport
        return self.transport if self.transport_back is None else self.transport_back


@dataclass(frozen=True)
class Step:
    """A visit of a node: entering through `port_in` (None at a source), leaving through
    `port_out` (None at a detector), with an optional channel (an order, "r"/"t", …)."""

    node: str
    port_in: str | None
    port_out: str | None
    channel: str = ""

    def __str__(self) -> str:
        ch = f" [{self.channel}]" if self.channel else ""
        return f"{self.node}({self.port_in or '·'}→{self.port_out or '·'}){ch}"


@dataclass(frozen=True)
class Branch:
    """A named, finite way from a source to a detector (ADR section 6)."""

    name: str
    steps: tuple[Step, ...]


@dataclass
class Model:
    """Topology plus branches (`rtt-fusion-graph` version 2).

    `links` are directed links of data sockets and observer links (light output of a node to a
    display or analysis block); `observe_branch` names, per observer input (node, socket), the
    branch whose light it reads where the observed node has several copies."""

    nodes: dict[str, Node] = field(default_factory=dict)
    connections: list[Connection] = field(default_factory=list)
    links: list[Link] = field(default_factory=list)
    branches: list[Branch] = field(default_factory=list)
    meta: dict[str, Any] = field(default_factory=dict)
    observe_branch: dict[tuple[str, str], str] = field(default_factory=dict)
    #: groups by the id of their node (type "group"), R28
    groups: dict[str, Group] = field(default_factory=dict)

    def add(self, node: Node) -> Node:
        if node.id in self.nodes:
            raise ValueError(f"node id {node.id!r} exists")
        self.nodes[node.id] = node
        return node

    def connect(
        self,
        a: tuple[str, str],
        b: tuple[str, str],
        transport: str = "auto",
        transport_back: str | None = None,
    ) -> Connection:
        c = Connection(a, b, transport, transport_back)
        self.connections.append(c)
        return c

    def branch(self, name: str, *steps: Step | tuple[Any, ...]) -> Branch:
        b = Branch(name, tuple(s if isinstance(s, Step) else Step(*s) for s in steps))
        self.branches.append(b)
        return b

    def connection_at(self, end: tuple[str, str]) -> Connection | None:
        return next((c for c in self.connections if end in (c.a, c.b)), None)


@dataclass
class Variant:
    """One content of a group (ADR section 7): a sub-model whose `branches` are the inner
    branches (from an exported port in to an exported port out, named), the exported ports
    (name -> inner node and port) and the exported parameters (name -> "node.param" or
    "node.param[i]" inside)."""

    model: Model
    ports: dict[str, tuple[str, str]]
    params: dict[str, str] = field(default_factory=dict)


@dataclass
class Group:
    """A group: variants with the same exported ports and parameters, one of them active.
    Its node (type "group") carries the pose of the group and the values of the exported
    parameters; the blocks inside lie relative to it."""

    variants: dict[str, Variant]
    active: str
    #: condensed (R29): {"sigma_max": …, "count": …} – the group is one block applying plane-wave
    #: tables of its inner branches (`graph.condense`); None: resolved into its blocks
    condense: dict[str, Any] | None = None

    def add_to(self, model: Model, node: Node) -> Node:
        model.add(node)
        model.groups[node.id] = self
        return node


def set_variant(model: Model, group_id: str, name: str) -> None:
    """Switches a group to another variant (replaces `replace_block` for groups)."""
    group = model.groups[group_id]
    if name not in group.variants:
        raise ValueError(f"group {group_id!r} has no variant {name!r}")
    group.active = name


def _group_checks(model: Model) -> list[Diagnostic]:
    out: list[Diagnostic] = []
    for gid, g in model.groups.items():
        if g.active not in g.variants:
            out.append(Diagnostic("graph.group_variant", "error", f"no variant {g.active!r}", gid))
            continue
        names = {n: (set(v.ports), set(v.params)) for n, v in g.variants.items()}
        ports, params = names[g.active]
        for n, (p, q) in names.items():
            if p != ports:
                out.append(
                    Diagnostic(
                        "graph.group_ports",
                        "error",
                        f"variant {n!r} exports the ports {sorted(p)}, "
                        f"{g.active!r} {sorted(ports)}",
                        gid,
                    )
                )
            if q != params:
                out.append(
                    Diagnostic(
                        "graph.group_params",
                        "error",
                        f"variant {n!r} exports the parameters {sorted(q)}",
                        gid,
                    )
                )
    return out


def _flatten(
    model: Model,
) -> tuple[Model, dict[tuple[str, str], tuple[str, str]], dict[tuple[str, str], list[Step]]]:
    """The model with every group replaced by its active variant: inner blocks as
    "<group>/<block>" with `_parent` and `_bind`, inner connections and links, outer
    connections to the inner ports, outer steps through a group replaced by the inner branch
    named in their channel. Returns also the exported ports and the inner branches."""
    flat = Model(meta=dict(model.meta), observe_branch=dict(model.observe_branch))
    ports: dict[tuple[str, str], tuple[str, str]] = {}
    inner_branches: dict[tuple[str, str], list[Step]] = {}
    problems: list[Diagnostic] = []
    for node in model.nodes.values():
        flat.add(replace(node, params=dict(node.params)))
    condensed: set[str] = set()
    for gid, g in model.groups.items():
        variant = g.variants[g.active]
        if g.condense is not None:
            from rtt_fusion.graph.condense import condensed_node

            flat.nodes[gid] = condensed_node(gid, flat.nodes[gid], variant, g.condense)
            condensed.add(gid)
            for b in variant.model.branches:
                inner_branches[(gid, b.name)] = []
            continue
        inner, inner_ports, nested = _flatten(variant.model)
        prefix = f"{gid}/"

        def name(node_id: str, prefix: str = prefix) -> str:
            return prefix + node_id

        bindings: dict[str, dict[str, str]] = {}
        for exported, address in variant.params.items():
            block, _, param = address.partition(".")
            bindings.setdefault(block, {})[param] = exported
        for iid, inode in inner.nodes.items():
            params = dict(inode.params)
            params["_parent"] = name(params["_parent"]) if "_parent" in params else gid
            if iid in bindings:
                params["_bind"] = {**params.get("_bind", {}), **bindings[iid]}
            flat.add(replace(inode, id=name(iid), params=params))
        for c in inner.connections:
            flat.connections.append(replace(c, a=(name(c.a[0]), c.a[1]), b=(name(c.b[0]), c.b[1])))
        for lk in inner.links:
            flat.links.append(
                replace(
                    lk,
                    source=(name(lk.source[0]), lk.source[1]),
                    target=(name(lk.target[0]), lk.target[1]),
                )
            )
        for port, (iid, iport) in variant.ports.items():
            # an exported port of a nested group resolves to the port inside it
            target_node, target_port = inner_ports.get((iid, iport), (iid, iport))
            ports[(gid, port)] = (name(target_node), target_port)
        for b in variant.model.branches:
            steps: list[Step] = []
            for s in b.steps:
                if s.node in variant.model.groups:
                    steps += nested.get((s.node, s.channel), [])
                else:
                    steps.append(s)
            inner_branches[(gid, b.name)] = [replace(s, node=name(s.node)) for s in steps]
    for c in model.connections:
        flat.connections.append(replace(c, a=ports.get(c.a, c.a), b=ports.get(c.b, c.b)))
    flat.links = [*flat.links, *model.links]
    for b in model.branches:
        steps = []
        for k, s in enumerate(b.steps):
            if s.node in condensed:  # one block: the step stays, its channel names the table
                if (s.node, s.channel) not in inner_branches:
                    problems.append(
                        Diagnostic(
                            "graph.group_branch",
                            "error",
                            f"step {k + 1}: group {s.node} has no inner branch {s.channel!r}",
                            f"branch {b.name}",
                        )
                    )
                steps.append(s)
                continue
            if s.node not in model.groups:
                steps.append(s)
                continue
            through = inner_branches.get((s.node, s.channel)) or []
            entry = ports.get((s.node, s.port_in or ""))
            exit_ = ports.get((s.node, s.port_out or ""))
            if (
                not through
                or entry != (through[0].node, through[0].port_in)
                or exit_ != (through[-1].node, through[-1].port_out)
            ):
                problems.append(
                    Diagnostic(
                        "graph.group_branch",
                        "error",
                        f"step {k + 1}: group {s.node} has no inner branch {s.channel!r} from "
                        f"{s.port_in} to {s.port_out}",
                        f"branch {b.name}",
                    )
                )
                continue
            steps += through
        flat.branches.append(Branch(b.name, tuple(steps)))
    if problems:
        raise ModelError(problems)
    return flat, ports, inner_branches


def flatten(model: Model, registry: NodeTypeRegistry | None = None) -> Model:
    """The model with its groups resolved (see `_flatten`); a model without groups as it is."""
    if not model.groups:
        return model
    problems = _group_checks(model)
    if problems:
        raise ModelError(problems)
    return _flatten(model)[0]


# ---------------------------------------------------------------------------- ports and scattering


def _port_key(spec: SocketSpec) -> str:
    n = spec.name
    if spec.direction is Direction.IN:
        if n.startswith("from_"):
            return n[5:]
        if n.endswith("_in") and len(n) > 3:
            return n[:-3]
    else:
        if n.startswith("to_"):
            return n[3:]
        if n.endswith("_out") and len(n) > 4:
            return n[:-4]
    return n


def derive_ports(sockets: Iterable[SocketSpec]) -> tuple[PortSpec, ...]:
    """Ports from the names of light sockets: `X_in`/`X_out`, `from_X`/`to_X` and `X`/`X_in`
    form port X; a socket without a partner is a port with one channel. Types where the names
    do not tell (a mirror's `in`/`out` is one port) declare `NodeType.ports`."""
    ports: dict[str, dict[str, str]] = {}
    order: list[str] = []
    for s in sockets:
        if s.kind is not SocketKind.LIGHT:
            continue
        key = _port_key(s)
        channel = "in" if s.direction is Direction.IN else "out"
        entry = ports.get(key)
        if entry is None or channel in entry:
            if entry is not None:  # a second socket of the same channel: a port of its own
                key = s.name
            entry = ports.setdefault(key, {})
            if key not in order:
                order.append(key)
        entry[channel] = s.name
        if s.side and not entry.get("side"):
            entry["side"] = s.side
    return tuple(
        PortSpec(k, ports[k].get("in", ""), ports[k].get("out", ""), ports[k].get("side", ""))
        for k in order
    )


def node_ports(node: Node, registry: NodeTypeRegistry) -> tuple[PortSpec, ...]:
    """The light ports of a node: its own, its type's, or derived from its light sockets."""
    if node.ports:
        return node.ports
    t = registry.get(node.type_id)
    if t is not None and t.ports and not node.sockets:
        return t.ports
    return derive_ports(node_sockets(node, registry))


def _port_of_socket(ports: Iterable[PortSpec], socket: str) -> PortSpec | None:
    return next((p for p in ports if socket in (p.in_socket, p.out_socket)), None)


def node_scatter(
    node: Node, registry: NodeTypeRegistry
) -> dict[tuple[str, str], tuple[str, tuple[str, ...]]]:
    """Scattering table (port in, port out) -> (kind, channels); pairs missing are "none".

    The type's table if it has one, else "oneway" for every pair of the node's `through`."""
    t = registry.get(node.type_id)
    if t is not None and t.scatter:
        return {(i, o): (k, tuple(ch)) for i, o, k, ch in t.scatter}
    ports = node_ports(node, registry)
    out: dict[tuple[str, str], tuple[str, tuple[str, ...]]] = {}
    if node.type_id == "rtt.section":  # surface reflections inside a section (R27): channel
        for p_in in ports:  # "ghost:<surface>,…" names the reflecting surfaces
            for p_out in ports:
                if p_in.in_socket and p_out.out_socket:
                    out[(p_in.name, p_out.name)] = ("ghost", ())
    for o_socket, ins in node_through(node, registry).items():
        po = _port_of_socket(ports, o_socket)
        if po is None or po.out_socket != o_socket:
            continue
        for i_socket in ins:
            pi = _port_of_socket(ports, i_socket)
            if pi is not None and pi.in_socket == i_socket:
                out[(pi.name, po.name)] = ("oneway", ())
    return out


# ---------------------------------------------------------------------------- checks


def _category(model: Model, registry: NodeTypeRegistry, node_id: str) -> Category | None:
    t = registry.get(model.nodes[node_id].type_id) if node_id in model.nodes else None
    return t.category if t is not None else None


def validate_model(model: Model, registry: NodeTypeRegistry) -> list[Diagnostic]:
    """Checks of the topology and the branches (the socket graph is checked after unrolling by
    `validate.validate`). A model with groups is checked resolved (`flatten`)."""
    if model.groups:
        try:
            return validate_model(flatten(model, registry), registry)
        except ModelError as e:
            return e.diagnostics
    out: list[Diagnostic] = []

    def err(code: str, msg: str, loc: str) -> None:
        out.append(Diagnostic(code, "error", msg, loc))

    ports = {nid: {p.name: p for p in node_ports(n, registry)} for nid, n in model.nodes.items()}

    seen: dict[tuple[str, str], int] = {}
    for c in model.connections:
        loc = f"{c.a[0]}.{c.a[1]} — {c.b[0]}.{c.b[1]}"
        for end in (c.a, c.b):
            if end[0] not in model.nodes or end[1] not in ports[end[0]]:
                err("graph.port_unknown", f"no port {end[1]!r} on {end[0]!r}", loc)
            seen[end] = seen.get(end, 0) + 1
    for end, n in seen.items():
        if n > 1:
            err(
                "graph.port_links",
                f"{n} connections at one port; splitting and joining happen in blocks",
                f"{end[0]}.{end[1]}",
            )

    if not model.branches:
        err("graph.no_branch", "the model has no branch; the modeller gives every branch", "")
    ends = {b.steps[-1].node for b in model.branches if b.steps}
    for nid in model.nodes:
        if _category(model, registry, nid) is Category.DETECTOR and nid not in ends:
            err("graph.no_branch", "no branch ends at this detector", nid)

    names: set[str] = set()
    words: set[tuple[Step, ...]] = set()
    for b in model.branches:
        loc = f"branch {b.name}"
        if b.name in names:
            err("graph.branch_duplicate", f"two branches named {b.name!r}", loc)
        names.add(b.name)
        if b.steps in words:
            err("graph.branch_duplicate", "the same steps as another branch", loc)
        words.add(b.steps)

        def impossible(msg: str, k: int, loc: str = loc) -> None:
            err("graph.branch_impossible", f"step {k + 1}: {msg}", loc)

        if len(b.steps) < 2:
            impossible("a branch runs from a source to a detector", 0)
            continue
        for k, s in enumerate(b.steps):
            if s.node not in model.nodes:
                impossible(f"unknown node {s.node!r}", k)
                continue
            first, last = k == 0, k == len(b.steps) - 1
            cat = _category(model, registry, s.node)
            if first and (s.port_in is not None or cat is not Category.SOURCE):
                impossible("a branch starts at a source (no port in)", k)
            if last and (s.port_out is not None or cat is not Category.DETECTOR):
                impossible("a branch ends at a detector (no port out)", k)
            if (not first and s.port_in is None) or (not last and s.port_out is None):
                impossible(f"{s} needs a port in and a port out", k)
                continue
            p = ports[s.node]
            if s.port_in is not None and (s.port_in not in p or not p[s.port_in].in_socket):
                impossible(f"{s.node} has no port {s.port_in!r} with an in channel", k)
                continue
            if s.port_out is not None and (s.port_out not in p or not p[s.port_out].out_socket):
                impossible(f"{s.node} has no port {s.port_out!r} with an out channel", k)
                continue
            if s.port_in is not None and s.port_out is not None:
                kind, channels = node_scatter(model.nodes[s.node], registry).get(
                    (s.port_in, s.port_out), ("none", ())
                )
                if kind == "none":
                    impossible(f"{s.node} gives nothing from {s.port_in} to {s.port_out}", k)
                elif s.channel and channels and s.channel not in channels:
                    impossible(f"{s.node} has no channel {s.channel!r} ({list(channels)})", k)
            if not last and s.port_out is not None:
                conn = model.connection_at((s.node, s.port_out))
                nxt = b.steps[k + 1]
                if conn is None:
                    impossible(f"port {s.node}.{s.port_out} has no connection", k)
                elif conn.other((s.node, s.port_out)) != (nxt.node, nxt.port_in):
                    impossible(
                        f"{s.node}.{s.port_out} is not connected to {nxt.node}.{nxt.port_in}", k
                    )

    out += connection_diagnostics(model, registry)

    # the Raytatouille backend passes a branch block by block through the traversals the
    # system's paths make (R25); one it does not know (a ghost) follows with R27
    runs = model.meta.get("raytatouille_runs")
    if isinstance(runs, list):
        known = {(str(r[0]), str(r[1]), str(r[2]), str(r[3])) for r in runs}
        for b in model.branches:
            for k, s in enumerate(b.steps):
                block = model.nodes.get(s.node)
                if block is None or not block.type_id.startswith("rtt.") or s.port_in is None:
                    continue
                if "rtt_new" in block.params:  # a new part: its traversals follow from its ports
                    continue
                if s.channel.startswith("ghost:"):
                    surfaces = [x for x in s.channel[6:].split(",") if x]
                    own = {str(x) for x in block.params.get("surfaces", [])}
                    back = s.port_in == s.port_out
                    if not surfaces or not set(surfaces) <= own or (len(surfaces) % 2 == 1) != back:
                        err(
                            "graph.branch_backend",
                            f"step {k + 1}: ghost {s.channel!r} – surfaces of {s.node}, an odd "
                            "number leaves through the port it came in, an even one through the "
                            "other",
                            f"branch {b.name}",
                        )
                    continue
                if (s.node, s.port_in, s.port_out or "", s.channel) not in known:
                    err(
                        "graph.branch_backend",
                        f"step {k + 1}: no path of the system passes {s.node} from "
                        f"{s.port_in} to {s.port_out or 'the detector'}",
                        f"branch {b.name}",
                    )
    return out


def port_frames(model: Model, registry: NodeTypeRegistry) -> dict[tuple[str, str], PortFrame]:
    """Frames of the ports in world coordinates, where the node types give them (R26)."""
    out: dict[tuple[str, str], PortFrame] = {}
    for nid, node in model.nodes.items():
        t = registry.get(node.type_id)
        if t is None or t.port_frames is None:
            continue
        world = posed(node, {p.name: p.default for p in t.params})
        for port, frame in t.port_frames(world).items():
            out[(nid, port)] = frame
    return out


#: tolerances of the connection checks: angle (rad), distance (mm), near field (wavelengths)
TILT_TOL, DISTANCE_TOL, NEAR_FIELD_WAVELENGTHS = 1e-6, 1e-9, 10.0


def connection_diagnostics(model: Model, registry: NodeTypeRegistry) -> list[Diagnostic]:
    """Checks of the relative position of the two ports of each connection (ADR section 5):
    graph.port_facing, graph.port_behind, graph.port_medium (errors); graph.port_aperture,
    graph.port_near_field, graph.port_tilt (warnings, the last for fields H: between tilted
    planes they are resampled geometrically, R18)."""
    frames = port_frames(model, registry)
    wavelengths = [
        float(w)
        for n in model.nodes.values()
        for w in (n.params.get("wavelengths_um") or [])
        if isinstance(w, int | float)
    ]
    lam_mm = max(wavelengths, default=0.0) * 1e-3
    out: list[Diagnostic] = []
    for c in model.connections:
        fa, fb = frames.get(c.a), frames.get(c.b)
        if fa is None or fb is None:
            continue
        loc = f"{c.a[0]}.{c.a[1]} — {c.b[0]}.{c.b[1]}"
        oa, ob = np.asarray(fa.origin), np.asarray(fb.origin)
        na, nb = np.asarray(fa.normal), np.asarray(fb.normal)
        facing = float(-(na @ nb))
        if facing <= 0.0:
            out.append(Diagnostic("graph.port_facing", "error", "the ports do not face", loc))
            continue
        gap = ob - oa
        if float(gap @ na) < -DISTANCE_TOL or float(-gap @ nb) < -DISTANCE_TOL:
            out.append(
                Diagnostic("graph.port_behind", "error", "a port lies behind the other", loc)
            )
            continue
        if fa.index is not None and fb.index is not None and abs(fa.index - fb.index) > 1e-12:
            out.append(
                Diagnostic(
                    "graph.port_medium",
                    "error",
                    f"different media outside the ports ({fa.index} / {fb.index})",
                    loc,
                )
            )
        if fa.radius is not None and fb.radius is not None:
            lateral = gap - (gap @ na) * na
            if float(np.linalg.norm(lateral)) > fa.radius + fb.radius:
                out.append(
                    Diagnostic(
                        "graph.port_aperture", "warning", "the apertures do not overlap", loc
                    )
                )
        distance = float(gap @ na)
        if lam_mm > 0 and DISTANCE_TOL < distance < NEAR_FIELD_WAVELENGTHS * lam_mm:
            out.append(
                Diagnostic(
                    "graph.port_near_field",
                    "warning",
                    f"{distance * 1e3:.3g} µm apart: evanescent fields are not carried",
                    loc,
                )
            )
        tilt = float(np.arccos(np.clip(facing, -1.0, 1.0)))
        h = "H" in (c.transport, c.transport_back or c.transport)
        if h and tilt > TILT_TOL:
            out.append(
                Diagnostic(
                    "graph.port_tilt",
                    "warning",
                    f"fields between planes tilted by {tilt:.3g} rad are resampled geometrically",
                    loc,
                )
            )
    return out


# ---------------------------------------------------------------------------- unrolling


@dataclass
class _Transition:
    src: int
    step: Step
    dst: int
    copy: int = -1


def _automaton(branches: list[Branch]) -> tuple[list[_Transition], dict[str, list[int]]]:
    """Prefix tree of the branches merged by equal remainders. Returns the transitions and, per
    branch, the indices of its transitions in order."""
    children: list[dict[Step, int]] = [{}]
    for b in branches:
        state = 0
        for s in b.steps:
            nxt = children[state].get(s)
            if nxt is None:
                children.append({})
                nxt = len(children) - 1
                children[state][s] = nxt
            state = nxt
    # class of a state = its remainder: the sorted (step, class of child) pairs, bottom up
    cls: dict[int, int] = {}
    signatures: dict[tuple[tuple[str, int], ...], int] = {}

    def classify(q: int) -> int:
        if q in cls:
            return cls[q]
        sig = tuple(sorted((repr(s), classify(c)) for s, c in children[q].items()))
        cls[q] = signatures.setdefault(sig, len(signatures))
        return cls[q]

    stack = [0]
    order: list[int] = []
    while stack:  # iterative post-order, branches may be long
        q = stack.pop()
        order.append(q)
        stack.extend(children[q].values())
    for q in reversed(order):
        classify(q)
    transitions: dict[tuple[int, Step, int], int] = {}
    listed: list[_Transition] = []
    per_branch: dict[str, list[int]] = {}
    for b in branches:
        state = 0
        idx: list[int] = []
        for s in b.steps:
            nxt = children[state][s]
            key = (cls[state], s, cls[nxt])
            if key not in transitions:
                transitions[key] = len(listed)
                listed.append(_Transition(cls[state], s, cls[nxt]))
            idx.append(transitions[key])
            state = nxt
        per_branch[b.name] = idx
    return listed, per_branch


def _copy_id(node_id: str, k: int) -> str:
    return node_id if k == 0 else f"{node_id}@{k + 1}"


def unroll(model: Model, registry: NodeTypeRegistry) -> Graph:
    """The socket graph of the model's branches (see the module docstring). Raises ModelError
    if `validate_model` finds an error."""
    model = flatten(model, registry)
    from rtt_fusion.graph.parts import prepare_parts  # S4: new elements take their place

    model = replace(model, nodes=prepare_parts(model.nodes))
    problems = [d for d in validate_model(model, registry) if d.severity == "error"]
    if problems:
        raise ModelError(problems)
    transitions, per_branch = _automaton(model.branches)
    ports = {nid: {p.name: p for p in node_ports(n, registry)} for nid, n in model.nodes.items()}

    # copies: a copy binds each of its input sockets to one source state and each output socket
    # to one target state
    binds: dict[str, list[tuple[dict[str, int], dict[str, int]]]] = {}
    for t in transitions:
        s = t.step
        copies = binds.setdefault(s.node, [])
        fits = [
            k
            for k, (bound_in, bound_out) in enumerate(copies)
            if (s.port_in is None or bound_in.get(s.port_in, t.src) == t.src)
            and (s.port_out is None or bound_out.get(s.port_out, t.dst) == t.dst)
        ]
        if fits:
            k = fits[0]
        else:
            copies.append(({}, {}))
            k = len(copies) - 1
        if s.port_in is not None:
            copies[k][0][s.port_in] = t.src
        if s.port_out is not None:
            copies[k][1][s.port_out] = t.dst
        t.copy = k

    g = Graph(meta=dict(model.meta))
    for nid, node in model.nodes.items():
        n_copies = max(len(binds.get(nid, ())), 1)
        original = node_through(node, registry)
        light = {s.name for s in node_sockets(node, registry) if s.kind is SocketKind.LIGHT}
        for k in range(n_copies):
            mine = [t for t in transitions if t.step.node == nid and t.copy == k]
            through: dict[str, list[str]] = {}
            if mine:
                for o, ins in original.items():  # data dependencies stay (a mirror's drive)
                    if o not in light:
                        through[o] = list(ins)
                    else:
                        through[o] = [i for i in ins if i not in light]
                for t in mine:
                    st = t.step
                    if st.port_out is None:
                        continue
                    o = ports[nid][st.port_out].out_socket
                    deps = through.setdefault(o, [])
                    if st.port_in is not None:
                        i = ports[nid][st.port_in].in_socket
                        if i not in deps:
                            deps.append(i)
            same = {o: set(i) for o, i in through.items()} == {
                o: set(i) for o, i in original.items()
            }
            sockets = node.sockets
            if node.sockets and mine:  # a ghost may use a channel no path of the import used
                have = {s.name for s in node.sockets}
                accepts = "".join(sorted(next(iter(node.sockets)).accepts))
                for t in mine:
                    for port_name, direction in ((t.step.port_in, "in"), (t.step.port_out, "out")):
                        if port_name is None:
                            continue
                        p = ports[nid][port_name]
                        name = p.in_socket if direction == "in" else p.out_socket
                        if name and name not in have:
                            make = light_in if direction == "in" else light_out
                            sockets = (*sockets, make(name, accepts, p.side, port_name))
                            have.add(name)
            params = dict(node.params)
            if node.type_id == "fusion.condensed":  # which inner branches this copy applies
                channels: dict[str, list[str]] = {}
                for t in mine:
                    key = f"{t.step.port_in}>{t.step.port_out}"
                    if t.step.channel not in channels.setdefault(key, []):
                        channels[key].append(t.step.channel)
                params["_channels"] = channels
            g.add(
                replace(
                    node,
                    id=_copy_id(nid, k),
                    sockets=sockets,
                    params=params,
                    position=(node.position[0], node.position[1] + 120.0 * k),
                    through=node.through
                    if not mine or same
                    else tuple((o, tuple(ins)) for o, ins in through.items()),
                )
            )

    # light links: from the out socket of a copy (bound to state d) to the in sockets bound to d
    feeds: dict[tuple[str, str], list[tuple[str, str]]] = {}  # in socket -> out sockets
    transport: dict[tuple[tuple[str, str], tuple[str, str]], str] = {}
    for t in transitions:
        st = t.step
        if st.port_out is None:
            continue
        out_end = (st.node, st.port_out)
        conn = model.connection_at(out_end)
        assert conn is not None  # checked
        source = (_copy_id(st.node, t.copy), ports[st.node][st.port_out].out_socket)
        for u in transitions:
            if u.src != t.dst or u.step.port_in is None:
                continue
            if (u.step.node, u.step.port_in) != conn.other(out_end):
                continue
            target = (_copy_id(u.step.node, u.copy), ports[u.step.node][u.step.port_in].in_socket)
            if source not in feeds.setdefault(target, []):
                feeds[target].append(source)
                transport[(source, target)] = conn.transport_from(out_end)
    junctions = 0
    for target, sources in feeds.items():
        if len(sources) == 1:
            g.connect(sources[0], target, transport[(sources[0], target)])
            continue
        junctions += 1
        jid = g.unique_id(f"sum_{target[0]}_{target[1]}")
        tr = {transport[(src, target)] for src in sources}
        g.add(
            Node(
                jid,
                SUM_TYPE,
                f"+ {target[0]}.{target[1]}",
                sockets=(
                    *(light_in(f"in_{k + 1}") for k in range(len(sources))),
                    light_out("out"),
                ),
                through=(("out", tuple(f"in_{k + 1}" for k in range(len(sources)))),),
            )
        )
        for k, src in enumerate(sources):
            g.connect(src, (jid, f"in_{k + 1}"), transport[(src, target)])
        g.connect((jid, "out"), target, tr.pop() if len(tr) == 1 else "auto")

    # data and observer links
    for lk in model.links:
        src_node = model.nodes.get(lk.source[0])
        light_src = src_node is not None and any(
            s.name == lk.source[1] and s.kind is SocketKind.LIGHT
            for s in node_sockets(src_node, registry)
        )
        src_copy = lk.source[0]
        if light_src:
            src_copy = _observed_copy(model, lk, transitions, per_branch, ports)
        n_target = max(len(binds.get(lk.target[0], ())), 1)
        for k in range(n_target):
            g.connect(
                (src_copy, lk.source[1]), (_copy_id(lk.target[0], k), lk.target[1]), lk.transport
            )

    # the blocks of each branch after its source (the Raytatouille backend traces these paths)
    paths: dict[str, list[str]] = {}
    for b in model.branches:
        ids = []
        for index in per_branch[b.name]:
            t = transitions[index]
            if _category(model, registry, t.step.node) is not Category.SOURCE:
                ids.append(_copy_id(t.step.node, t.copy))
        paths[b.name] = ids
    g.meta["paths"] = paths
    g.meta["branches"] = [branch_to_dict(b) for b in model.branches]
    if junctions:
        g.meta["junctions"] = junctions
    return g


def _observed_copy(
    model: Model,
    lk: Link,
    transitions: list[_Transition],
    per_branch: dict[str, list[int]],
    ports: dict[str, dict[str, PortSpec]],
) -> str:
    node_id, socket = lk.source
    port = next((p.name for p in ports[node_id].values() if p.out_socket == socket), None)
    used = [t for t in transitions if t.step.node == node_id and t.step.port_out == port]
    wanted = model.observe_branch.get(lk.target)
    if wanted is not None:
        used = [transitions[i] for i in per_branch.get(wanted, ()) if transitions[i] in used]
    copies = sorted({t.copy for t in used})
    if len(copies) > 1:
        raise ModelError(
            [
                Diagnostic(
                    "graph.observer_ambiguous",
                    "error",
                    f"{node_id}.{socket} is passed by several copies; name the branch",
                    f"{lk.target[0]}.{lk.target[1]}",
                )
            ]
        )
    return _copy_id(node_id, copies[0] if copies else 0)


# ---------------------------------------------------------------------------- version 1 -> 2


def _is_topology(registry: NodeTypeRegistry, node: Node) -> bool:
    t = registry.get(node.type_id)
    return t is not None and t.category in (
        Category.SOURCE,
        Category.DOMAIN,
        Category.SPLITTER,
        Category.DETECTOR,
        Category.CONVERSION,
    )


def model_from_graph(graph: Graph, registry: NodeTypeRegistry) -> Model:
    """A version 1 working graph as a model: pairs of sockets become ports, the light links
    between two ports one connection, and the source-to-detector ways of the directed graph –
    which in version 1 *were* the choice of the ways – the branches (named after
    `graph.meta["paths"]` where they match). Light links into blocks outside the light ways
    (displays) become observer links."""
    m = Model(meta={k: v for k, v in graph.meta.items() if k not in ("paths", "branches")})
    for node in graph.nodes.values():
        m.add(replace(node, params=dict(node.params)))
    ports = {nid: {p.name: p for p in node_ports(n, registry)} for nid, n in graph.nodes.items()}

    def port_of(node_id: str, socket: str) -> str | None:
        return next(
            (p.name for p in ports[node_id].values() if socket in (p.in_socket, p.out_socket)), None
        )

    def is_light(node_id: str, socket: str) -> bool:
        return any(
            s.name == socket and s.kind is SocketKind.LIGHT
            for s in node_sockets(graph.nodes[node_id], registry)
        )

    light_links: list[Link] = []
    for lk in graph.links:
        if (
            is_light(*lk.source)
            and is_light(*lk.target)
            and _is_topology(registry, graph.nodes[lk.target[0]])
        ):
            light_links.append(lk)
        else:
            m.links.append(lk)
    for lk in light_links:
        a = (lk.source[0], port_of(*lk.source) or lk.source[1])
        b = (lk.target[0], port_of(*lk.target) or lk.target[1])
        existing = next((c for c in m.connections if {c.a, c.b} == {a, b}), None)
        if existing is None:
            if m.connection_at(a) is not None or m.connection_at(b) is not None:
                raise ValueError(
                    f"link {lk.source} -> {lk.target}: a port has one connection; split the "
                    "light in a block"
                )
            m.connect(a, b, lk.transport)
        elif existing.a == b and lk.transport != existing.transport:
            m.connections[m.connections.index(existing)] = replace(
                existing, transport_back=lk.transport
            )

    # branches: every way along links and through from a source to a detector
    def category(node_id: str) -> Category | None:
        t = registry.get(graph.nodes[node_id].type_id)
        return t.category if t is not None else None

    words: list[tuple[Step, ...]] = []

    def walk(prefix: tuple[Step, ...], node_id: str, in_socket: str) -> None:
        port_in = port_of(node_id, in_socket)
        if category(node_id) is Category.DETECTOR:
            words.append((*prefix, Step(node_id, port_in, None)))
            return
        for o, ins in node_through(graph.nodes[node_id], registry).items():
            if in_socket in ins and is_light(node_id, o):
                step = Step(node_id, port_in, port_of(node_id, o))
                for lk in light_links:
                    if lk.source == (node_id, o):
                        walk((*prefix, step), *lk.target)

    for nid in graph.nodes:
        if category(nid) is Category.SOURCE:
            for s in node_sockets(graph.nodes[nid], registry):
                if s.kind is SocketKind.LIGHT and s.direction is Direction.OUT:
                    first = Step(nid, None, port_of(nid, s.name))
                    for lk in light_links:
                        if lk.source == (nid, s.name):
                            walk((first,), *lk.target)
    # imported graphs: one branch per Raytatouille path; paths taking the same way (diffraction
    # orders – version 1 knows no channels) differ by the channel "path:<name>" at the detector
    named: dict[tuple[str, ...], list[str]] = {}
    for name, blocks in graph.meta.get("paths", {}).items():
        if isinstance(blocks, list):
            named.setdefault(tuple(blocks), []).append(str(name))
    if "paths" in graph.meta and graph.meta.get("imported_from") == "raytatouille":
        m.meta["raytatouille_paths"] = list(graph.meta["paths"])
    taken: set[str] = set()
    for k, w in enumerate(words):
        blocks = tuple(s.node for s in w if category(s.node) is not Category.SOURCE)
        names = named.get(blocks, [f"branch {k + 1}"])
        for name in names:
            steps = w
            if len(names) > 1:
                steps = (*w[:-1], replace(w[-1], channel=f"path:{name}"))
            unique = name if name not in taken else f"{name} ({k + 1})"
            taken.add(unique)
            m.branches.append(Branch(unique, steps))
    return m


# ---------------------------------------------------------------------------- format


def step_to_list(s: Step) -> list[Any]:
    out: list[Any] = [s.node, s.port_in, s.port_out]
    if s.channel:
        out.append(s.channel)
    return out


def branch_to_dict(b: Branch) -> dict[str, Any]:
    return {"name": b.name, "steps": [step_to_list(s) for s in b.steps]}


def branch_from_dict(d: dict[str, Any]) -> Branch:
    steps = []
    for s in d["steps"]:
        steps.append(
            Step(
                str(s[0]),
                None if s[1] is None else str(s[1]),
                None if s[2] is None else str(s[2]),
                str(s[3]) if len(s) > 3 else "",
            )
        )
    return Branch(str(d["name"]), tuple(steps))


__all__ = [
    "SUM_TYPE",
    "Branch",
    "Connection",
    "Group",
    "Model",
    "ModelError",
    "Step",
    "Variant",
    "connection_diagnostics",
    "derive_ports",
    "flatten",
    "model_from_graph",
    "node_ports",
    "node_scatter",
    "port_frames",
    "set_variant",
    "unroll",
    "validate_model",
]
