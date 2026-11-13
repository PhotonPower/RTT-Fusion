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
"""

from __future__ import annotations

import itertools
import json
import shutil
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import Any

from rtt_fusion import portdata
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


class BlockProvider:
    """The blocks of a plugin: their types and the evaluators that run them in its process."""

    def __init__(self, types: list[NodeType], evaluators: dict[str, Evaluator]) -> None:
        missing = [t.type_id for t in types if t.type_id not in evaluators]
        if missing:
            raise ValueError(f"no evaluator for {', '.join(missing)}")
        self.types = {t.type_id: t for t in types}
        self.evaluators = evaluators
        #: kept across requests (operators, prepared solvers – ADR 0010 D5)
        self.operators: dict[str, Any] = {}

    def catalogue(self) -> list[dict[str, Any]]:
        """The block types for the hello: the static description and the representations each
        evaluator works in natively."""
        out = []
        for type_id, t in self.types.items():
            reps = getattr(self.evaluators[type_id], "representations", frozenset("R"))
            out.append({**t.to_dict(), "representations": sorted(reps)})
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
            graph = Graph()
            graph.add(node)
            ctx = Context(graph, shared=dict(request.get("shared", {})), operators=self.operators)
            value = self.evaluators[node.type_id](ctx, node, str(request["socket"]), inputs)
            log = list(ctx.log)
            handle = portdata.write(value, request["out"])
            return {"ok": True, "outputs": handle, "log": log}
        except Exception as e:
            return {"ok": False, "error": f"{type(e).__name__}: {e}", "log": log}


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
    that stopped – is the block's ExecutionError "exec.remote"; the run goes on elsewhere."""

    _count = itertools.count()

    def __init__(
        self, provider: str, call: Call, representations: frozenset[str], root: str | Path
    ) -> None:
        self.provider = provider
        self.call = call
        self.representations = representations
        self.root = Path(root)

    def __call__(self, ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
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
            try:
                response = self.call(self.provider, request)
            except Exception as e:  # noqa: BLE001 - the channel: the provider stopped, …
                raise ExecutionError("exec.remote", f"{self.provider}: {e}", node.id) from e
            ctx.log.extend(str(x) for x in response.get("log", []))
            if not response.get("ok"):
                raise ExecutionError(
                    "exec.remote", f"{self.provider}: {response.get('error', 'evaluation failed')}", node.id
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
    """Adds the providers' block types to `registry`; their evaluators (by type id)."""
    evaluators: dict[str, Evaluator] = {}
    for provider, t, reps in remote_types(catalogues):
        if registry.get(t.type_id) is None:
            registry.add(t)
        evaluators[t.type_id] = RemoteEvaluator(provider, call, reps, root)
    return evaluators


__all__ = ["BlockProvider", "Call", "RemoteEvaluator", "register_remote", "remote_types"]
