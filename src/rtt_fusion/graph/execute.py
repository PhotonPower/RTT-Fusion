"""Execution of a working graph without planner (stage R15, concept section 6).

Each output socket is evaluated on demand from the inputs it depends on (`through`), by an
evaluator registered for the block type. Results are cached under a content hash of the block
type, its parameters, the socket, the transports and the hashes of everything upstream, so
changing a parameter recomputes only what lies downstream (Blender-like). A muted block passes
its light through unchanged (light outputs) or yields nothing (data outputs).

Before running, the planner (stage R17, `rtt_fusion.graph.planner`) fixes the representation of
every light link from the representations the executor can run (`available`, default R: ray
bundles, `rtt_fusion.graph.light`; R18 adds H: fields, `rtt_fusion.graph.field`, for backends
whose evaluators produce ray bundles). "auto" picks R where both are possible, with the results of
R15; fixed transports outside `available` raise `exec.transport_unsupported`. The plan is kept in
`Executor.plan`.

Evaluators work in their native representation (attribute `representations` of the evaluator,
default R). The executor converts at the ports: a link planned H turns the rays of its source into
fields on a port plane (R -> H); if the evaluator at its end works with rays only, the fields
become rays again there (H -> R, after nothing but the conversion; propagation of fields between
blocks belongs to the evaluators that take fields, such as detectors).
"""

from __future__ import annotations

import hashlib
import json
from collections import OrderedDict
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any

from rtt_fusion.graph.field import FieldBundle, FieldError, field_to_rays, rays_to_field
from rtt_fusion.graph.light import RayBundle
from rtt_fusion.graph.model import Direction, Graph, Node, SocketKind
from rtt_fusion.graph.planner import Plan, Requirements, link_key, plan
from rtt_fusion.graph.pose import Isometry, bind_exported, placement, posed
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.graph.validate import errors, node_sockets, node_through, validate

#: evaluator(context, node, output socket, inputs) -> value of that output.
#: Light inputs and outputs are lists of RayBundle; data values are free objects.
Evaluator = Callable[["Context", Node, str, dict[str, Any]], Any]

#: Representations run by default (backends whose evaluators produce rays may add H, R18).
EXECUTABLE = frozenset({"R"})
WITH_FIELDS = frozenset({"R", "H"})
SUPPORTED_TRANSPORTS = frozenset({"auto"}) | EXECUTABLE


class ExecutionError(RuntimeError):
    def __init__(self, code: str, message: str, location: str = "") -> None:
        super().__init__(f"{code}: {message}" + (f" ({location})" if location else ""))
        self.code = code
        self.location = location


@dataclass
class Context:
    """What evaluators may need besides their inputs: the graph, shared objects (e.g. a compiled
    Raytatouille system), and a log."""

    graph: Graph
    shared: dict[str, Any] = field(default_factory=dict)
    log: list[str] = field(default_factory=list)
    plan: Plan | None = None  # the plan of the run (representation per link)
    evaluators: dict[str, Any] = field(default_factory=dict)  # for backends that call other blocks
    registry: Any = None
    #: value of another block's output (meta blocks, R20); set by the executor
    evaluate: Callable[[str, str], Any] | None = None
    #: progress of the run, part in [0, 1] of the whole and a text (GUI: status bar); may raise
    #: to cancel (the plugin's cancel check)
    progress: Callable[[float, str], None] | None = None
    #: called with the id of every block whose evaluator starts (GUI: the running block)
    on_node: Callable[[str], None] | None = None
    #: operator store (R26): the expensive, pose-independent parts of solver blocks in their
    #: local frame (a BEM factorisation, a built table), kept across runs; keys without the pose
    operators: Any = None
    _span: tuple[float, float] = (0.0, 1.0)

    def report(self, part: float, text: str = "") -> None:
        """Progress `part` (0 … 1) of the current step, mapped into the whole run."""
        if self.progress is None:
            return
        a, b = self._span
        self.progress(a + (b - a) * min(max(float(part), 0.0), 1.0), text)

    @contextmanager
    def span(self, part0: float, part1: float) -> Iterator[None]:
        """Nested progress: within the block, `report(0 … 1)` covers part0 … part1 of the
        current step (e.g. one value of a sweep)."""
        a, b = self._span
        self._span = (a + (b - a) * part0, a + (b - a) * part1)
        try:
            yield
        finally:
            self._span = (a, b)


class Executor:
    def __init__(
        self,
        graph: Graph,
        registry: NodeTypeRegistry,
        evaluators: dict[str, Evaluator],
        shared: dict[str, Any] | None = None,
        available: frozenset[str] = EXECUTABLE,
        cache: ResultCache | None = None,
        environment: object = None,
        operators: ResultCache | None = None,
    ) -> None:
        self.graph = graph
        self.registry = registry
        self.evaluators = evaluators
        self.context = Context(graph, dict(shared or {}), evaluators=evaluators, registry=registry)
        self._running: set[str] = set()
        self.context.evaluate = self._nested
        self.operators: ResultCache = operators if operators is not None else ResultCache()
        self.context.operators = self.operators
        self.cache: ResultCache = cache if cache is not None else ResultCache()
        self.environment = _hash(environment) if environment is not None else _environment(shared)
        self.calls: dict[tuple[str, str], int] = {}  # evaluator calls per output (tests, profiling)
        self.plan: Plan | None = None
        self.available = frozenset(available)

    # ------------------------------------------------------------------ public
    def check(self) -> None:
        """Raises for errors of the graph check, for transports that cannot run yet and for
        errors of the plan; keeps the plan in `self.plan`."""
        found = errors(validate(self.graph, self.registry))
        if found:
            d = found[0]
            raise ExecutionError(d.code, d.message, d.location)
        for lk in self.graph.links:
            if lk.transport != "auto" and lk.transport not in self.available:
                raise ExecutionError(
                    "exec.transport_unsupported",
                    f"transport {lk.transport} comes with a later stage (this backend runs "
                    f"{''.join(sorted(self.available))})",
                    f"{lk.source[0]}.{lk.source[1]} -> {lk.target[0]}.{lk.target[1]}",
                )
        self.plan = plan(self.graph, self.registry, Requirements(available=self.available))
        self.context.plan = self.plan
        for d in self.plan.diagnostics:
            if d.severity == "error":
                raise ExecutionError(d.code, d.message, d.location)

    def evaluate(self, node_id: str, socket: str) -> Any:
        """The value of an output socket (computing what it needs)."""
        self.check()
        return self._value(node_id, socket, ())[0]

    def results(self) -> dict[tuple[str, str], Any]:
        """All data outputs of detectors and analyses (what a run produces)."""
        self.check()
        out: dict[tuple[str, str], Any] = {}
        for node in self.graph.nodes.values():
            for s in node_sockets(node, self.registry):
                if s.direction is Direction.OUT and s.kind is SocketKind.DATA:
                    out[(node.id, s.name)] = self._value(node.id, s.name, ())[0]
        return out

    # ------------------------------------------------------------------ internals
    def _key(self, node_id: str, socket: str, stack: tuple[str, ...] = ()) -> str:
        """Content hash of an output without computing it: the environment, the block (type,
        parameters, socket, muted), the keys of everything it depends on (through its inputs and,
        for meta blocks, the output they observe) and the planned representation per link."""
        here = f"{node_id}.{socket}"
        if here in stack:
            raise ExecutionError("exec.cycle", " -> ".join((*stack, here)), here)
        node = self.graph.node(node_id)
        deps = node_through(node, self.registry).get(socket, ())
        input_keys: list[str] = []
        for in_name in deps:
            links = self.graph.links_to(node_id, in_name)
            if not links:
                input_keys.append(f"{in_name}=∅")
                continue
            lk = links[0]
            k = self._key(lk.source[0], lk.source[1], (*stack, here))
            rep = self.plan.representations.get(link_key(lk), "") if self.plan else ""
            input_keys.append(f"{in_name}={k}|{rep or lk.transport}")
        extra_keys: list[str] = []
        key_extra = getattr(self.evaluators.get(node.type_id), "key_extra", None)
        if key_extra is not None:
            for item in key_extra(self.context, node):
                if isinstance(item, tuple):  # another output the block depends on
                    extra_keys.append(self._key(item[0], item[1], (*stack, here)))
                else:  # a plain part of the key (e.g. the parameter a meta block varies)
                    extra_keys.append(str(item))
        return _hash(
            {
                "env": self.environment,
                "type": node.type_id,
                "params": node.params,
                "groups": [g.params for g in self._ancestors(node)],
                "socket": socket,
                "muted": node.muted,
                "inputs": input_keys,
                "extra": extra_keys,
            }
        )

    def _value(self, node_id: str, socket: str, stack: tuple[str, ...]) -> tuple[Any, str]:
        here = f"{node_id}.{socket}"
        if here in stack:
            raise ExecutionError("exec.cycle", " -> ".join((*stack, here)), here)
        node = self.graph.node(node_id)
        spec = next((s for s in node_sockets(node, self.registry) if s.name == socket), None)
        if spec is None or spec.direction is not Direction.OUT:
            raise ExecutionError("exec.socket_unknown", f"no output {socket!r}", node_id)
        key = self._key(node_id, socket, stack)
        evaluator = self.evaluators.get(node.type_id)
        uncached = bool(getattr(evaluator, "uncached", False))
        if not uncached and key in self.cache:
            return self.cache[key], key  # nothing in front of it is computed either
        deps = node_through(node, self.registry).get(socket, ())
        inputs: dict[str, Any] = {}
        for in_name in deps:
            links = self.graph.links_to(node_id, in_name)
            if not links:
                inputs[in_name] = []
                continue
            lk = links[0]
            value, _ = self._value(lk.source[0], lk.source[1], (*stack, here))
            rep = self.plan.representations.get(link_key(lk), "R") if self.plan else "R"
            inputs[in_name] = self._convert(value, rep, node, lk)
        if node.muted or node.type_id == "fusion.sum":  # a junction adds what arrives
            if spec.kind is SocketKind.LIGHT:
                value = [b for v in inputs.values() for b in (v or [])]
            else:
                value = None
        else:
            if evaluator is None:
                raise ExecutionError(
                    "exec.no_evaluator", f"no evaluator for {node.type_id}", node_id
                )
            self.calls[(node_id, socket)] = self.calls.get((node_id, socket), 0) + 1
            if self.context.on_node is not None:
                self.context.on_node(node_id)
            value = evaluator(self.context, self._posed(node), socket, inputs)
        if not uncached:
            self.cache[key] = value
        return value, key

    def _ancestors(self, node: Node) -> list[Node]:
        """The groups a block lies in, innermost first (R28: `_parent` of the unrolled blocks)."""
        out: list[Node] = []
        parent = node.params.get("_parent")
        while parent and parent in self.graph.nodes:
            out.append(self.graph.nodes[parent])
            parent = self.graph.nodes[parent].params.get("_parent")
        return out

    def _posed(self, node: Node) -> Node:
        """The node with its geometric parameters in world coordinates (pose ∘ offset, R26), the
        exported parameters of its groups applied and their placements composed (R28)."""
        t = self.registry.get(node.type_id)
        defaults = {p.name: p.default for p in t.params} if t is not None else {}
        ancestors = self._ancestors(node)
        if not ancestors:
            return posed(node, defaults)
        bound = bind_exported(node, ancestors)
        outer = Isometry.identity()
        for g in reversed(ancestors):  # outermost first
            outer = outer @ placement(g)
        return posed(bound, defaults, outer)

    def needed(self, wanted: list[tuple[str, str]]) -> dict[str, set[str]]:
        """The outputs per block that the outputs `wanted` depend on (through the inputs and the
        outputs meta blocks observe), `wanted` included – what a run of them computes."""
        out: dict[str, set[str]] = {}
        stack = list(wanted)
        while stack:
            node_id, socket = stack.pop()
            if socket in out.get(node_id, set()) or node_id not in self.graph.nodes:
                continue
            out.setdefault(node_id, set()).add(socket)
            node = self.graph.nodes[node_id]
            for in_name in node_through(node, self.registry).get(socket, ()):
                stack += [lk.source for lk in self.graph.links_to(node_id, in_name)]
            key_extra = getattr(self.evaluators.get(node.type_id), "key_extra", None)
            if key_extra is not None:
                stack += [i for i in key_extra(self.context, node) if isinstance(i, tuple)]
        return out

    def status(self, wanted: list[tuple[str, str]] | None = None) -> dict[str, bool]:
        """Per block: True if its outputs are in the cache (up to date), without computing
        anything (GUI: a green tick / a red cross). With `wanted` (e.g. the outputs of the
        display blocks) only the blocks those need are listed, each up to date if the outputs
        they need are – an output nobody needs (a splitter's way to a camera nobody shows) does
        not count."""
        self.check()
        if wanted is not None:
            result: dict[str, bool] = {}
            for node_id, sockets in self.needed(wanted).items():
                try:
                    result[node_id] = all(self._key(node_id, s) in self.cache for s in sockets)
                except ExecutionError:
                    result[node_id] = False
            return result
        out: dict[str, bool] = {}
        for node in self.graph.nodes.values():
            # outputs with a value (a "drive" carries none), and of those the linked ones (an
            # unused channel of a splitter is never computed); a block without linked outputs
            # (a display, a detector) counts all its outputs
            valued = [
                s.name
                for s in node_sockets(node, self.registry)
                if s.direction is Direction.OUT and "drive" not in s.accepts
            ]
            linked = [name for name in valued if self.graph.links_from(node.id, name)]
            outputs = linked or valued
            try:
                out[node.id] = bool(outputs) and all(
                    self._key(node.id, name) in self.cache for name in outputs
                )
            except ExecutionError:
                out[node.id] = False
        return out

    def _nested(self, node_id: str, socket: str) -> Any:
        """Evaluation for meta blocks (parameters may have changed; the cache keys contain them, the
        plan does not depend on them). A meta block may not collect its own output."""
        here = f"{node_id}.{socket}"
        if here in self._running:
            raise ExecutionError("exec.cycle", f"a meta block collects {here} inside itself", here)
        self._running.add(here)
        try:
            return self._value(node_id, socket, ())[0]
        finally:
            self._running.discard(here)

    def _convert(self, value: Any, rep: str, consumer: Node, lk: Any) -> Any:
        """Light on a link in its planned representation, as the consumer can take it."""
        if not isinstance(value, list):
            return value
        native = getattr(self.evaluators.get(consumer.type_id), "representations", frozenset("R"))
        where = f"{lk.source[0]}.{lk.source[1]} -> {lk.target[0]}.{lk.target[1]}"
        out: list[Any] = []
        try:
            for item in value:
                if rep == "H" and isinstance(item, RayBundle):
                    item = rays_to_field(item)
                # rays where the consumer cannot take fields, and where the user asked for rays
                # on the link (a block that got a field returns one)
                if isinstance(item, FieldBundle) and (
                    "H" not in native or consumer.muted or lk.transport == "R"
                ):
                    item = field_to_rays(item)
                out.append(item)
        except FieldError as e:
            raise ExecutionError("exec.field", str(e), where) from e
        return out


def _environment(shared: dict[str, Any] | None) -> str:
    """Fingerprint of the settings in `shared` that results depend on: plain values (numbers,
    text, lists) and the system's JSON; objects without a plain form (a compiled system, a
    material library) are left out – pass `environment` explicitly to cover them."""
    plain: dict[str, Any] = {}
    for k, v in sorted((shared or {}).items()):
        if isinstance(v, int | float | str | bool) or v is None:
            plain[k] = v
        elif isinstance(v, list | tuple) and all(isinstance(x, int | float | str) for x in v):
            plain[k] = list(v)
        elif k == "system" and hasattr(v, "to_json"):
            plain[k] = hashlib.sha256(v.to_json().encode("utf-8")).hexdigest()
    return _hash(plain)


def _nbytes(obj: Any, seen: set[int] | None = None, depth: int = 0) -> int:
    """Estimated memory of a result: numpy arrays inside dataclasses, lists, tuples and dicts."""
    seen = set() if seen is None else seen
    if id(obj) in seen or depth > 8:
        return 0
    seen.add(id(obj))
    nbytes = getattr(obj, "nbytes", None)
    if isinstance(nbytes, int) and hasattr(obj, "dtype"):
        return nbytes
    if isinstance(obj, dict):
        return sum(_nbytes(v, seen, depth + 1) for v in obj.values())
    if isinstance(obj, list | tuple | set | frozenset):
        return sum(_nbytes(v, seen, depth + 1) for v in obj)
    fields = getattr(obj, "__dict__", None)
    if fields is None and hasattr(obj, "__dataclass_fields__"):
        fields = {f: getattr(obj, f, None) for f in obj.__dataclass_fields__}
    if isinstance(fields, dict):
        return 64 + sum(_nbytes(v, seen, depth + 1) for v in fields.values())
    return 16


class ResultCache:
    """Results by content key, kept across runs (Blender-like: only what changed is computed
    again), at most `limit_bytes` (estimated); the least recently used results go first, with a
    line in `warnings`."""

    def __init__(self, limit_bytes: int = 1 << 30) -> None:
        self.limit_bytes = int(limit_bytes)
        self._items: OrderedDict[str, tuple[Any, int]] = OrderedDict()
        self.size_bytes = 0
        self.warnings: list[str] = []

    def __contains__(self, key: object) -> bool:
        return key in self._items

    def __getitem__(self, key: str) -> Any:
        value, _ = self._items[key]
        self._items.move_to_end(key)
        return value

    def __setitem__(self, key: str, value: Any) -> None:
        if key in self._items:
            self.size_bytes -= self._items.pop(key)[1]
        size = _nbytes(value)
        self._items[key] = (value, size)
        self.size_bytes += size
        evicted = 0
        while self.size_bytes > self.limit_bytes and len(self._items) > 1:
            _, (_, s) = self._items.popitem(last=False)
            self.size_bytes -= s
            evicted += 1
        if evicted:
            self.warnings.append(
                f"result cache over {self.limit_bytes / 2**30:.2g} GB: {evicted} older result(s) "
                "dropped (they are computed again when needed)"
            )

    def __len__(self) -> int:
        return len(self._items)

    def clear(self) -> None:
        self._items.clear()
        self.size_bytes = 0


def _hash(obj: object) -> str:
    text = json.dumps(obj, sort_keys=True, default=str, ensure_ascii=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:24]
