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
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any

from rtt_fusion.graph.field import FieldBundle, FieldError, field_to_rays, rays_to_field
from rtt_fusion.graph.light import RayBundle
from rtt_fusion.graph.model import Direction, Graph, Node, SocketKind
from rtt_fusion.graph.planner import Plan, Requirements, link_key, plan
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
    ) -> None:
        self.graph = graph
        self.registry = registry
        self.evaluators = evaluators
        self.context = Context(graph, dict(shared or {}), evaluators=evaluators, registry=registry)
        self._running: set[str] = set()
        self.context.evaluate = self._nested
        self.cache: dict[str, Any] = {}
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
    def _value(self, node_id: str, socket: str, stack: tuple[str, ...]) -> tuple[Any, str]:
        here = f"{node_id}.{socket}"
        if here in stack:
            raise ExecutionError("exec.cycle", " -> ".join((*stack, here)), here)
        node = self.graph.node(node_id)
        spec = next((s for s in node_sockets(node, self.registry) if s.name == socket), None)
        if spec is None or spec.direction is not Direction.OUT:
            raise ExecutionError("exec.socket_unknown", f"no output {socket!r}", node_id)
        deps = node_through(node, self.registry).get(socket, ())
        inputs: dict[str, Any] = {}
        input_keys: list[str] = []
        for in_name in deps:
            links = self.graph.links_to(node_id, in_name)
            if not links:
                inputs[in_name] = []
                input_keys.append(f"{in_name}=∅")
                continue
            lk = links[0]
            value, key = self._value(lk.source[0], lk.source[1], (*stack, here))
            rep = self.plan.representations.get(link_key(lk), "R") if self.plan else "R"
            value = self._convert(value, rep, node, lk)
            inputs[in_name] = value
            rep = self.plan.representations.get(link_key(lk), "") if self.plan else ""
            input_keys.append(f"{in_name}={key}|{rep or lk.transport}")
        key = _hash(
            {
                "type": node.type_id,
                "params": node.params,
                "socket": socket,
                "muted": node.muted,
                "inputs": input_keys,
            }
        )
        evaluator = self.evaluators.get(node.type_id)
        uncached = bool(getattr(evaluator, "uncached", False))
        if key in self.cache and not uncached:
            return self.cache[key], key
        if node.muted:
            if spec.kind is SocketKind.LIGHT:
                value = [b for v in inputs.values() for b in (v or [])]
            else:
                value = None
        else:
            evaluator = self.evaluators.get(node.type_id)
            if evaluator is None:
                raise ExecutionError(
                    "exec.no_evaluator", f"no evaluator for {node.type_id}", node_id
                )
            self.calls[(node_id, socket)] = self.calls.get((node_id, socket), 0) + 1
            if self.context.on_node is not None:
                self.context.on_node(node_id)
            value = evaluator(self.context, node, socket, inputs)
        if not uncached:
            self.cache[key] = value
        return value, key

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
                if isinstance(item, FieldBundle) and ("H" not in native or consumer.muted):
                    item = field_to_rays(item)
                out.append(item)
        except FieldError as e:
            raise ExecutionError("exec.field", str(e), where) from e
        return out


def _hash(obj: object) -> str:
    text = json.dumps(obj, sort_keys=True, default=str, ensure_ascii=False)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:24]
