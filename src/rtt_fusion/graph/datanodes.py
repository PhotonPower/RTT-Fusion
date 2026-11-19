"""Evaluators of the data blocks (D13, Fusion plan 6.2): blocks that deliver data, not light.

- `biooptomd.layer`: a BioOptoMD run directory (`adapters.biooptomd.read_layer`, manifest
  `biooptomd-layer.json`) → `profile` (N, `adapters.biooptomd.Profile`) and `dipoles` (Q,
  `adapters.biooptomd.DipoleSeries`: the tracker series per receptor state in BioOptoMD's sensor
  frame, unchanged; a node that computes light converts it with `DipoleSeries.in_rtt_frame`).

The files enter the cache key with their size and modification time, so a rewritten run is read
again although the block's parameters did not change.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from rtt_fusion.adapters import biooptomd as bmd
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.model import Node
from rtt_fusion.graph.registry import NodeTypeRegistry


def _directory(node: Node) -> str:
    path = str(node.params.get("directory", "") or "")
    if not path:
        raise ExecutionError("exec.param", "no BioOptoMD output directory", node.id)
    return path


def biooptomd_layer(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
    """The profile (N) or the dipole series (Q) of the run directory."""
    try:
        layer = bmd.read_layer(_directory(node))
    except ValueError as e:
        raise ExecutionError("exec.data", str(e), node.id) from e
    if socket == "profile":
        if layer.profile is None:
            raise ExecutionError("exec.data", "the run has no index profile", node.id)
        return layer.profile
    if socket == "dipoles":
        if layer.dipoles is None:
            raise ExecutionError("exec.data", "the run has no dipole series", node.id)
        return layer.dipoles
    raise ExecutionError("exec.socket_unknown", socket, node.id)


def _files_key(ctx: Context, node: Node) -> list[Any]:
    out: list[Any] = []
    path = str(node.params.get("directory", "") or "")
    if not path:
        return out
    for f in bmd.layer_files(path):
        try:
            st = Path(f).stat()
            out.append(f"{f}:{st.st_size}:{st.st_mtime_ns}")
        except OSError:
            out.append(f"{f}:missing")
    return out


biooptomd_layer.key_extra = _files_key  # type: ignore[attr-defined]
biooptomd_layer.representations = frozenset()  # type: ignore[attr-defined]  # no light

EVALUATORS: dict[str, Evaluator] = {"biooptomd.layer": biooptomd_layer}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    """The types are built in (`registry.builtin_registry`); this adds their evaluators."""
    return dict(EVALUATORS)
