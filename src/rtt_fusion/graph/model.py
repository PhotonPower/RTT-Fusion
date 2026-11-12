"""Working graph of the node workflow (docs/konzept-node-workflow.md, sections 2 and 3).

A graph consists of blocks (`Node`) of a registered type (`NodeType`), sockets on them
(`SocketSpec`) and links between an output and an input socket (`Link`). Light sockets accept a set
of representations of the closed set of the Fusion ADR (R, H, G, S, T, P); data sockets carry one
named data type (N, Q, field_map, …; section 3.2 and L6). A link carries a transport: "auto" (the
planner chooses) or a fixed representation.

The graph is plain data: parameters are values, never code (Raytatouille decision O2/F3).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum
from typing import Any

#: Light representations (Fusion ADR point 3; P with M9).
LIGHT_REPRESENTATIONS = frozenset({"R", "H", "G", "S", "T", "P"})
#: Data types (Fusion ADR addendum R8–R10 and concept L6); new ones come with their first module.
DATA_TYPES = frozenset(
    {
        "N",
        "Q",
        "field_map",
        "occupancy",
        "voxel_materials",
        "shape",
        "height_map",
        "response_map",
        "scattering_stats",
        "result",
        "drive",
    }
)
#: Transports of a light link: "auto" or one fixed representation.
TRANSPORTS = frozenset({"auto"}) | LIGHT_REPRESENTATIONS


class SocketKind(str, Enum):
    LIGHT = "light"
    DATA = "data"


class Direction(str, Enum):
    IN = "in"
    OUT = "out"


@dataclass(frozen=True)
class SocketSpec:
    """A socket of a node type.

    Attributes:
        name: unique within the node type (and the node, for dynamic sockets).
        direction: input or output.
        kind: light or data.
        accepts: light: the representations the socket can take (subset of LIGHT_REPRESENTATIONS);
            data: exactly one data type.
        side: "front", "back" or "" (no side), for two-sided domains (concept section 4).
        label: display name.
    """

    name: str
    direction: Direction
    kind: SocketKind
    accepts: frozenset[str]
    side: str = ""
    label: str = ""

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("socket name must not be empty")
        allowed = LIGHT_REPRESENTATIONS if self.kind is SocketKind.LIGHT else DATA_TYPES
        if not self.accepts or not self.accepts <= allowed:
            raise ValueError(
                f"socket {self.name}: accepts {sorted(self.accepts)} not in {sorted(allowed)}"
            )
        if self.kind is SocketKind.DATA and len(self.accepts) != 1:
            raise ValueError(f"data socket {self.name} carries exactly one data type")
        if self.side not in ("", "front", "back"):
            raise ValueError(f"socket {self.name}: side must be 'front', 'back' or ''")

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {
            "name": self.name,
            "direction": self.direction.value,
            "kind": self.kind.value,
            "accepts": sorted(self.accepts),
        }
        if self.side:
            d["side"] = self.side
        if self.label:
            d["label"] = self.label
        return d

    @staticmethod
    def from_dict(d: dict[str, Any]) -> SocketSpec:
        return SocketSpec(
            name=str(d["name"]),
            direction=Direction(d["direction"]),
            kind=SocketKind(d["kind"]),
            accepts=frozenset(str(a) for a in d["accepts"]),
            side=str(d.get("side", "")),
            label=str(d.get("label", "")),
        )


def light_in(name: str, accepts: str = "RHGSTP", side: str = "", label: str = "") -> SocketSpec:
    return SocketSpec(name, Direction.IN, SocketKind.LIGHT, frozenset(accepts), side, label)


def light_out(name: str, accepts: str = "RHGSTP", side: str = "", label: str = "") -> SocketSpec:
    return SocketSpec(name, Direction.OUT, SocketKind.LIGHT, frozenset(accepts), side, label)


def data_in(name: str, data_type: str, label: str = "") -> SocketSpec:
    return SocketSpec(name, Direction.IN, SocketKind.DATA, frozenset({data_type}), "", label)


def data_out(name: str, data_type: str, label: str = "") -> SocketSpec:
    return SocketSpec(name, Direction.OUT, SocketKind.DATA, frozenset({data_type}), "", label)


@dataclass(frozen=True)
class PortSpec:
    """A light port (ADR ports, poses, branches, section 2): a reference surface of the node with
    two channels, *in* (light entering the node) and *out* (light leaving it). In the derived
    socket form each channel is a light socket; a port with one channel only (a source, a
    detector, a converter's sides) has the other name empty.

    Attributes:
        name: unique among the node's ports.
        in_socket: name of the light input socket of the in channel, or "".
        out_socket: name of the light output socket of the out channel, or "".
        side: "front", "back" or "".
        label: display name.
    """

    name: str
    in_socket: str = ""
    out_socket: str = ""
    side: str = ""
    label: str = ""

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("port name must not be empty")
        if not self.in_socket and not self.out_socket:
            raise ValueError(f"port {self.name} has no channel")

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {"name": self.name}
        for key, value in (
            ("in", self.in_socket),
            ("out", self.out_socket),
            ("side", self.side),
            ("label", self.label),
        ):
            if value:
                d[key] = value
        return d

    @staticmethod
    def from_dict(d: dict[str, Any]) -> PortSpec:
        return PortSpec(
            str(d["name"]),
            str(d.get("in", "")),
            str(d.get("out", "")),
            str(d.get("side", "")),
            str(d.get("label", "")),
        )


#: Entries of a scattering table (ADR section 3): what a node can deliver from one port to another.
SCATTER_KINDS = ("none", "oneway", "full", "ghost")


class Category(str, Enum):
    """Kinds of blocks (concept section 3.1)."""

    SOURCE = "source"
    DOMAIN = "domain"
    SPLITTER = "splitter"
    DETECTOR = "detector"
    CONVERSION = "conversion"
    DATA = "data"
    ANALYSIS = "analysis"
    GROUP = "group"
    META = "meta"


@dataclass(frozen=True)
class ParamSpec:
    """A parameter of a node type (A2): the same fields as nfcad.Param of Raytatouille Fusion."""

    name: str
    type: str  # number, integer, boolean, choice, text, file, files, numbers, point, entity
    default: Any = None
    unit: str = ""
    minimum: float | None = None
    maximum: float | None = None
    options: tuple[str, ...] = ()
    group: str = ""
    label: str = ""
    description: str = ""

    def to_dict(self) -> dict[str, Any]:
        d: dict[str, Any] = {"name": self.name, "type": self.type, "default": self.default}
        for key, value in (
            ("unit", self.unit),
            ("group", self.group),
            ("label", self.label),
            ("description", self.description),
        ):
            if value:
                d[key] = value
        if self.minimum is not None:
            d["min"] = self.minimum
        if self.maximum is not None:
            d["max"] = self.maximum
        if self.options:
            d["options"] = list(self.options)
        return d

    @staticmethod
    def from_dict(d: dict[str, Any]) -> ParamSpec:
        return ParamSpec(
            name=str(d["name"]),
            type=str(d["type"]),
            default=d.get("default"),
            unit=str(d.get("unit", "")),
            minimum=d.get("min"),
            maximum=d.get("max"),
            options=tuple(d.get("options", ())),
            group=str(d.get("group", "")),
            label=str(d.get("label", "")),
            description=str(d.get("description", "")),
        )


@dataclass(frozen=True)
class NodeType:
    """A registered block type.

    Attributes:
        type_id: e.g. "rtt.section", "fouriergrating.grating".
        category: kind of block.
        sockets: the fixed sockets; nodes may add dynamic ones if `dynamic_sockets`.
        params: parameter descriptions.
        through: internal dependencies of the light sockets, output name -> input names. A
            two-sided domain passed in both directions has front_in -> back_out and
            back_in -> front_out; the graph is acyclic on this level (concept section 4).
            Empty for a type means "every output depends on every input".
        dynamic_sockets: nodes of this type carry their own sockets (imported sections).
        extra_sockets: sockets added to the own sockets of such nodes (see the field).
        description: one line for the node library.
    """

    type_id: str
    category: Category
    sockets: tuple[SocketSpec, ...] = ()
    params: tuple[ParamSpec, ...] = ()
    through: tuple[tuple[str, tuple[str, ...]], ...] = ()
    dynamic_sockets: bool = False
    description: str = ""
    #: sockets every node of the type has in addition to its own (e.g. the "drive" input of an
    #: imported mirror, whose light sockets come from the import)
    extra_sockets: tuple[SocketSpec, ...] = ()
    #: light ports (ADR ports, section 2); empty: derived from the light sockets by their names
    #: (`topology.derive_ports`). Needed where the names do not tell, e.g. a mirror's in/out.
    ports: tuple[PortSpec, ...] = ()
    #: scattering table (port in, port out, kind, channels); empty: every pair of `through` is
    #: "oneway", every other pair "none" (`topology.node_scatter`)
    scatter: tuple[tuple[str, str, str, tuple[str, ...]], ...] = ()

    def socket(self, name: str) -> SocketSpec | None:
        return next((s for s in self.sockets if s.name == name), None)


@dataclass
class Node:
    """A block in the graph.

    `sockets`, `ports` and `through` override the type's (types with dynamic sockets; `through`
    also on the nodes of the socket form unrolled from branches). `position` is the place in the
    editor (no meaning for the computation).
    """

    id: str
    type_id: str
    name: str = ""
    params: dict[str, Any] = field(default_factory=dict)
    position: tuple[float, float] = (0.0, 0.0)
    muted: bool = False
    sockets: tuple[SocketSpec, ...] = ()
    through: tuple[tuple[str, tuple[str, ...]], ...] = ()
    ports: tuple[PortSpec, ...] = ()


@dataclass(frozen=True)
class Link:
    """A link from an output socket to an input socket, with a transport for light links."""

    source: tuple[str, str]  # (node id, socket name)
    target: tuple[str, str]
    transport: str = "auto"

    def __post_init__(self) -> None:
        if self.transport not in TRANSPORTS:
            raise ValueError(f"transport {self.transport!r} not in {sorted(TRANSPORTS)}")


class GraphError(ValueError):
    """An editing operation that cannot be applied (unknown node, socket, …)."""


@dataclass
class Graph:
    """The working graph: nodes by id, links, free metadata (e.g. the system it came from)."""

    nodes: dict[str, Node] = field(default_factory=dict)
    links: list[Link] = field(default_factory=list)
    meta: dict[str, Any] = field(default_factory=dict)

    # ------------------------------------------------------------------ queries
    def node(self, node_id: str) -> Node:
        try:
            return self.nodes[node_id]
        except KeyError:
            raise GraphError(f"unknown node {node_id!r}") from None

    def links_from(self, node_id: str, socket: str | None = None) -> list[Link]:
        return [
            lk
            for lk in self.links
            if lk.source[0] == node_id and (socket is None or lk.source[1] == socket)
        ]

    def links_to(self, node_id: str, socket: str | None = None) -> list[Link]:
        return [
            lk
            for lk in self.links
            if lk.target[0] == node_id and (socket is None or lk.target[1] == socket)
        ]

    # ------------------------------------------------------------------ editing
    def add(self, node: Node) -> Node:
        if node.id in self.nodes:
            raise GraphError(f"node id {node.id!r} exists")
        self.nodes[node.id] = node
        return node

    def unique_id(self, stem: str) -> str:
        if stem not in self.nodes:
            return stem
        k = 2
        while f"{stem}_{k}" in self.nodes:
            k += 1
        return f"{stem}_{k}"

    def remove(self, node_id: str) -> None:
        """Removes a node and its links."""
        self.node(node_id)
        del self.nodes[node_id]
        self.links = [lk for lk in self.links if node_id not in (lk.source[0], lk.target[0])]

    def connect(
        self, source: tuple[str, str], target: tuple[str, str], transport: str = "auto"
    ) -> Link:
        """Adds a link (checked by `validate`, not here, so a graph can be edited step by step)."""
        self.node(source[0])
        self.node(target[0])
        link = Link(source, target, transport)
        if link in self.links:
            raise GraphError(f"link {source} -> {target} exists")
        self.links.append(link)
        return link

    def disconnect(self, link: Link) -> None:
        try:
            self.links.remove(link)
        except ValueError:
            raise GraphError("unknown link") from None

    def set_transport(self, link: Link, transport: str) -> Link:
        """Fixes the transport of a link ("auto", "R", "H", …)."""
        index = self.links.index(link)
        self.links[index] = replace(link, transport=transport)
        return self.links[index]

    def insert_on_link(
        self, link: Link, node: Node, in_socket: str, out_socket: str
    ) -> tuple[Link, Link]:
        """Inserts `node` into `link`: source -> node.in_socket, node.out_socket -> target (the
        editor's drop on a link). Both new links keep the transport of the old one."""
        if node.id not in self.nodes:
            self.add(node)
        self.disconnect(link)
        a = self.connect(link.source, (node.id, in_socket), link.transport)
        b = self.connect((node.id, out_socket), link.target, link.transport)
        return a, b

    def replace_node(
        self, node_id: str, new: Node, socket_map: dict[str, str] | None = None
    ) -> Node:
        """Exchanges a block and keeps its links: each link of the old node's socket s moves to
        socket socket_map.get(s, s) of the new node (e.g. mirror -> FourierGrating, concept 1)."""
        old = self.node(node_id)
        mapping = socket_map or {}
        new.position = old.position if new.position == (0.0, 0.0) else new.position
        moved: list[Link] = []
        for lk in self.links:
            source, target = lk.source, lk.target
            if source[0] == node_id:
                source = (new.id, mapping.get(source[1], source[1]))
            if target[0] == node_id:
                target = (new.id, mapping.get(target[1], target[1]))
            moved.append(replace(lk, source=source, target=target))
        del self.nodes[node_id]
        if new.id in self.nodes:
            raise GraphError(f"node id {new.id!r} exists")
        self.nodes[new.id] = new
        self.links = moved
        return new

    def set_muted(self, node_id: str, muted: bool = True) -> None:
        """Mute / bypass (concept section 6): light passes unchanged or the block drops out."""
        self.node(node_id).muted = muted
