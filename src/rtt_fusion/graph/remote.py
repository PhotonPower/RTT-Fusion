"""Blocks evaluated in another process (CAD ADR 0010, D2/D3, stage W2b).

Two sides of one small protocol; the transport (CAD's plugin channel) carries plain JSON:

- **Provider** (a plugin with blocks): `BlockProvider(types, evaluators)` describes its block types
  (`catalogue`, for the plugin's hello) and answers an evaluation request (`evaluate`): the inputs
  come as port data (`rtt_fusion.portdata`) in a folder, the outputs go as port data into another.
- **Engine**: `remote_types` turns the providers' catalogues into block types of its registry and
  `RemoteEvaluator` evaluates such a block by a request (`call(provider, request) -> response`).

A remote block gets its inputs, its parameters and the run's plain settings (numbers, strings,
lists – `shared` without objects); it cannot ask for other blocks' values (`ctx.evaluate`, meta
blocks stay in the engine) and has no port frames (code, not data).

Prepare (ADR 0010, D6, stage W3): blocks that need more than their inputs – Raytatouille blocks
make one system from the data of all `rtt.*` blocks and follow the branches of the whole graph –
get the run's unrolled graph once, before the first of their evaluations in that run
(`{"kind": "prepare", "run", "graph", "shared"}`); their evaluations then name the run and work in
its context (`Context.graph` the whole graph, `shared` the plain settings plus what the provider's
`prepare` hook adds, e.g. a material library).

Changes within a run (issue #70): meta blocks (`meta.sweep`, `meta.replicates`, `meta.optimise`)
change parameters of blocks in the engine's graph after the run was prepared – the drive of a
mirror, a pose, a parameter. Every evaluation therefore names the blocks whose parameters differ
from the prepared graph (`"changes": {id: params}`), and the provider sets every block of its copy
to the prepared parameters plus these changes before it evaluates (a value that goes back to the
prepared one is reset too).
"""

from __future__ import annotations

import itertools
import json
import shutil
import tempfile
import uuid
from collections import OrderedDict
from collections.abc import Callable
from pathlib import Path
from typing import Any

from rtt_fusion import portdata
from rtt_fusion.graph import io
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.model import Graph, Node, NodeType
from rtt_fusion.graph.registry import NodeTypeRegistry

#: request a provider answers: {"node", "socket", "inputs": {socket: handle}, "folders": {socket:
#: folder of that input}, "out": folder for the outputs, "shared": plain settings of the run}
#: response: {"ok": True, "outputs": handle, "log": [...]} or {"ok": False, "error": "..."}
Call = Callable[[str, dict[str, Any]], dict[str, Any]]


def _plain(shared: dict[str, Any]) -> dict[str, Any]:
    """The JSON-able part of the run's shared settings (numbers, strings, lists of them)."""
    out: dict[str, Any] = {}
    for k, v in shared.items():
        try:
            json.dumps(v, allow_nan=False)
        except (TypeError, ValueError):
            continue
        out[k] = v
    return out


#: what a provider's prepare hook adds to a run's shared settings (objects such as a library)
Prepare = Callable[[Graph, dict[str, Any]], dict[str, Any]]
#: prepared runs a provider keeps (the newest)
KEPT_RUNS = 4


class BlockProvider:
    """The blocks of a plugin: their types and the evaluators that run them in its process.
    `prepare`: what a run's context needs beyond the plain settings (D6), from its graph."""

    def __init__(
        self,
        types: list[NodeType],
        evaluators: dict[str, Evaluator],
        prepare: Prepare | None = None,
    ) -> None:
        missing = [t.type_id for t in types if t.type_id not in evaluators]
        if missing:
            raise ValueError(f"no evaluator for {', '.join(missing)}")
        self.types = {t.type_id: t for t in types}
        self.evaluators = evaluators
        self.prepare_hook = prepare
        #: kept across requests (operators, prepared solvers – ADR 0010 D5)
        self.operators: dict[str, Any] = {}
        self._runs: OrderedDict[str, tuple[Graph, dict[str, Any]]] = OrderedDict()
        #: the parameters of every block of a prepared run, as prepared (issue #70)
        self._prepared: dict[str, dict[str, str]] = {}

    def prepare(self, request: dict[str, Any]) -> dict[str, Any]:
        """Keep a run's graph and context (D6); the newest `KEPT_RUNS` runs stay."""
        try:
            graph = io.from_dict(request["graph"])
            shared = dict(request.get("shared", {}))
            if self.prepare_hook is not None:
                shared.update(self.prepare_hook(graph, shared))
            self._runs[str(request["run"])] = (graph, shared)
            self._prepared[str(request["run"])] = _snapshot(graph)
            while len(self._runs) > KEPT_RUNS:
                old, _ = self._runs.popitem(last=False)
                self._prepared.pop(old, None)
            return {"ok": True}
        except Exception as e:
            return {"ok": False, "error": f"{type(e).__name__}: {e}"}

    def handle(self, request: dict[str, Any]) -> dict[str, Any]:
        """One request of the engine: prepare a run or evaluate a block."""
        if request.get("kind") == "prepare":
            return self.prepare(request)
        return self.evaluate(request)

    def catalogue(self) -> list[dict[str, Any]]:
        """The block types for the hello: the static description, the representations each
        evaluator works in natively and whether a run is prepared for them (D6)."""
        out = []
        for type_id, t in self.types.items():
            reps = getattr(self.evaluators[type_id], "representations", frozenset("R"))
            entry = {**t.to_dict(), "representations": sorted(reps)}
            if self.prepare_hook is not None:
                entry["prepare"] = True
            out.append(entry)
        return out

    def evaluate(self, request: dict[str, Any]) -> dict[str, Any]:
        """Evaluate one block for a request of the engine; errors become {"ok": False}."""
        log: list[str] = []
        try:
            n = request["node"]
            node = Node(
                str(n["id"]),
                str(n["type"]),
                str(n.get("name", "")),
                params=dict(n.get("params", {})),
                muted=bool(n.get("muted", False)),
            )
            if node.type_id not in self.types:
                raise ValueError(f"this plugin provides no block type {node.type_id!r}")
            folders = request.get("folders", {})
            inputs = {k: portdata.read(h, folders[k]) for k, h in request.get("inputs", {}).items()}
            run = request.get("run")
            if run is not None and str(run) in self._runs:  # the prepared run's context (D6)
                graph, shared = self._runs[str(run)]
                _apply_changes(graph, self._prepared[str(run)], request.get("changes", {}))
                node = graph.nodes.get(node.id, node)
            else:
                graph, shared = Graph(), dict(request.get("shared", {}))
                graph.add(node)
            ctx = Context(graph, shared=dict(shared), operators=self.operators)
            value = self.evaluators[node.type_id](ctx, node, str(request["socket"]), inputs)
            log = list(ctx.log)
            handle = portdata.write(value, request["out"])
            return {"ok": True, "outputs": handle, "log": log}
        except Exception as e:
            return {"ok": False, "error": f"{type(e).__name__}: {e}", "log": log}


def _snapshot(graph: Graph) -> dict[str, str]:
    """The parameters of every block as JSON text (a comparable, independent copy)."""
    return {nid: json.dumps(n.params, sort_keys=True) for nid, n in graph.nodes.items()}


def _apply_changes(graph: Graph, prepared: dict[str, str], changes: dict[str, Any]) -> None:
    """Every block of a prepared copy to its prepared parameters plus the request's changes."""
    for nid, node in graph.nodes.items():
        if nid in changes:
            params = dict(changes[nid])
        elif nid in prepared:
            params = json.loads(prepared[nid])
        else:
            continue
        if params != node.params:
            node.params.clear()
            node.params.update(params)


def remote_types(catalogues: list[dict[str, Any]]) -> list[tuple[str, NodeType, frozenset[str]]]:
    """(provider, type, representations) of the providers' catalogue entries
    ({"provider", "type": {...}})."""
    out = []
    for entry in catalogues:
        t = entry["type"]
        out.append(
            (
                str(entry["provider"]),
                NodeType.from_dict(t),
                frozenset(t.get("representations", ["R"])),
            )
        )
    return out


class RemoteEvaluator:
    """An evaluator for a block of another plugin: inputs as port data into a folder below `root`,
    one request, the outputs read back. The provider's error – or a broken channel, a provider
    that stopped – is the block's ExecutionError "exec.remote"; the run goes on elsewhere.
    `prepare`: before the provider's first block of a run, its graph is sent (D6)."""

    _count = itertools.count()

    def __init__(
        self,
        provider: str,
        call: Call,
        representations: frozenset[str],
        root: str | Path,
        prepare: bool = False,
        runs: dict[tuple[str, int], tuple[str, dict[str, str]]] | None = None,
    ) -> None:
        self.provider = provider
        self.call = call
        self.representations = representations
        self.root = Path(root)
        self.prepare = prepare
        #: prepared runs by (provider, graph) with the parameters sent then: shared by the
        #: evaluators of one engine
        self.runs = runs if runs is not None else {}

    def _changes(self, ctx: Context, key: tuple[str, int]) -> dict[str, Any]:
        """The blocks whose parameters differ from those the run was prepared with (#70)."""
        prepared = self.runs[key][1]
        out: dict[str, Any] = {}
        for nid, n in ctx.graph.nodes.items():
            if nid in prepared and json.dumps(n.params, sort_keys=True) != prepared[nid]:
                out[nid] = n.params
        return out

    def _run(self, ctx: Context, node: Node) -> str | None:
        if not self.prepare:
            return None
        key = (self.provider, id(ctx.graph))
        if key not in self.runs:
            run = uuid.uuid4().hex
            request = {
                "kind": "prepare",
                "run": run,
                "graph": io.to_dict(ctx.graph),
                "shared": _plain(ctx.shared),
            }
            try:
                answer = self.call(self.provider, request)
            except Exception as e:
                raise ExecutionError("exec.remote", f"{self.provider}: {e}", node.id) from e
            if not answer.get("ok"):
                why = answer.get("error", "prepare failed")
                raise ExecutionError("exec.remote", f"{self.provider}: {why}", node.id)
            self.runs[key] = (run, _snapshot(ctx.graph))
        return self.runs[key][0]

    def __call__(self, ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
        run = self._run(ctx, node)
        self.root.mkdir(parents=True, exist_ok=True)
        folder = Path(tempfile.mkdtemp(prefix=f"{next(self._count)}-", dir=self.root))
        try:
            handles = {k: portdata.write(v, folder / "in" / k) for k, v in inputs.items()}
            request = {
                "node": {
                    "id": node.id,
                    "type": node.type_id,
                    "name": node.name,
                    "params": node.params,
                    "muted": node.muted,
                },
                "socket": socket,
                "inputs": handles,
                "folders": {k: str(folder / "in" / k) for k in handles},
                "out": str(folder / "out"),
                "shared": _plain(ctx.shared),
            }
            if run is not None:
                request["run"] = run
                request["changes"] = self._changes(ctx, (self.provider, id(ctx.graph)))
            try:
                response = self.call(self.provider, request)
            except Exception as e:
                raise ExecutionError("exec.remote", f"{self.provider}: {e}", node.id) from e
            ctx.log.extend(str(x) for x in response.get("log", []))
            if not response.get("ok"):
                raise ExecutionError(
                    "exec.remote",
                    f"{self.provider}: {response.get('error', 'evaluation failed')}",
                    node.id,
                )
            return portdata.read(response["outputs"], folder / "out")
        finally:
            shutil.rmtree(folder, ignore_errors=True)


def register_remote(
    registry: NodeTypeRegistry,
    catalogues: list[dict[str, Any]],
    call: Call,
    root: str | Path,
) -> dict[str, Evaluator]:
    """Adds the providers' block types to `registry`; their evaluators (by type id). A type the
    registry has already (e.g. `rtt.*`, built in) keeps its description and is evaluated by the
    provider. A catalogue entry with "prepare": true gets the run's graph first (D6)."""
    evaluators: dict[str, Evaluator] = {}
    runs: dict[tuple[str, int], tuple[str, dict[str, str]]] = {}
    for entry, (provider, t, reps) in zip(catalogues, remote_types(catalogues), strict=True):
        if registry.get(t.type_id) is None:
            registry.add(t)
        prepare = bool(entry["type"].get("prepare", False))
        evaluators[t.type_id] = RemoteEvaluator(provider, call, reps, root, prepare, runs)
    return evaluators


__all__ = ["BlockProvider", "Call", "RemoteEvaluator", "register_remote", "remote_types"]
