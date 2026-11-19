"""Evaluators of the data blocks (D13, Fusion plan 6.2): blocks that deliver data, not light.

- `biooptomd.layer`: a BioOptoMD run directory (`adapters.biooptomd.read_layer`, manifest
  `biooptomd-layer.json`) → `profile` (N, `adapters.biooptomd.Profile`) and `dipoles` (Q,
  `adapters.biooptomd.DipoleSeries`: the tracker series per receptor state in BioOptoMD's sensor
  frame, unchanged; a node that computes light converts it with `DipoleSeries.in_rtt_frame`).

- `data.field_map`, `data.response_map`, `data.occupancy`: a file in the `raytatouille-result`
  envelope (types `FieldMap`, `EmitterResponseTable`, `Occupancy`; hp-FEM writes F and M,
  BindOptiX B) → its data port (D13), read by BindOptiX (`adapters.bindoptix.read`).
- `bindoptix.signal`: F, M, Q and B → the mean detected photon rate per receptor, spot, time and
  bound detector region (`adapters.bindoptix.signal`, a `SignalResult` on `rates`). Q carries all
  receptor states, so one Q input suffices (no socket per state).

The files enter the cache key with their size and modification time, so a rewritten file is read
again although the block's parameters did not change.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from rtt_fusion.adapters import bindoptix as bxa
from rtt_fusion.adapters import biooptomd as bmd
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.model import Category, Node, NodeType, ParamSpec, data_in, data_out
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


# ------------------------------------------------------------------ data files F, M, B (D13)

_PORT_OF = {
    "data.field_map": "field_map",
    "data.response_map": "response_map",
    "data.occupancy": "occupancy",
}


def _bindoptix_missing(node: Node, e: ImportError) -> ExecutionError:
    return ExecutionError(
        "exec.solver",
        f"BindOptiX is not installed ({e}); pip install "
        '"bindoptix @ git+https://github.com/PhotonPower/BindOptiX" '
        "-C cmake.define.BINDOPTIX_BUILD_CORE=OFF",
        node.id,
    )


def data_file(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
    """The object of a `raytatouille-result` file of a D13 data port."""
    path = str(node.params.get("file", "") or "")
    if not path:
        raise ExecutionError("exec.param", "no file", node.id)
    try:
        return bxa.read(_PORT_OF[node.type_id], path)
    except ImportError as e:
        raise _bindoptix_missing(node, e) from e
    except ValueError as e:
        raise ExecutionError("exec.data", str(e), node.id) from e


def _file_key(ctx: Context, node: Node) -> list[Any]:
    path = str(node.params.get("file", "") or "")
    try:
        st = Path(path).stat()
        return [f"{path}:{st.st_size}:{st.st_mtime_ns}"]
    except OSError:
        return [f"{path}:missing"]


data_file.key_extra = _file_key  # type: ignore[attr-defined]
data_file.representations = frozenset()  # type: ignore[attr-defined]


# ------------------------------------------------------------------ bindoptix.signal


def _number(node: Node, name: str) -> float:
    try:
        return float(node.params[name])
    except (KeyError, TypeError, ValueError) as e:
        raise ExecutionError("exec.param", f"{name}: a number is required", node.id) from e


def bindoptix_signal(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
    """Detected photon rate per receptor, spot, time and bound region (`SignalResult`)."""
    missing = [n for n in ("excitation", "response", "dipoles", "occupancy") if n not in inputs]
    if missing:
        raise ExecutionError("exec.data", f"inputs not connected: {missing}", node.id)
    q = inputs["dipoles"]
    if not isinstance(q, bmd.DipoleSeries):
        raise ExecutionError("exec.data", "dipoles: not a BioOptoMD dipole series (Q)", node.id)
    p = node.params
    dye = {
        "epsilon_per_molar_cm": _number(node, "epsilon_per_molar_cm"),
        "quantum_yield": _number(node, "quantum_yield"),
        "excitation_um": _number(node, "excitation_um"),
        "emission_um": _number(node, "emission_um"),
    }
    anchors = [float(x) for row in p.get("anchors_um", [[0.0, 0.0, 0.0]]) for x in row]
    if not anchors or len(anchors) % 3:
        raise ExecutionError("exec.param", "anchors_um: triples (x, y, z) in µm", node.id)
    irradiance = float(p.get("port_irradiance_w_m2_per_mw", 0.0) or 0.0)
    port = (
        (
            irradiance,
            str(p.get("port_irradiance_status", "P")),
            str(p.get("port_irradiance_reference", "")),
        )
        if irradiance > 0
        else None
    )
    try:
        return bxa.signal(
            inputs["excitation"],
            inputs["response"],
            q,
            inputs["occupancy"],
            dye,
            _number(node, "power_mw"),
            np.asarray(anchors).reshape(-1, 3),
            str(p.get("regions", "") or ""),
            port,
        )
    except ImportError as e:
        raise _bindoptix_missing(node, e) from e
    except ValueError as e:
        raise ExecutionError("exec.data", str(e), node.id) from e


bindoptix_signal.representations = frozenset()  # type: ignore[attr-defined]

_FILE = (ParamSpec("file", "file", "", description="raytatouille-result file"),)

TYPES = (
    NodeType(
        "data.field_map",
        Category.DATA,
        (data_out("field_map", "field_map"),),
        _FILE,
        description="Excitation field map F from a file (FieldMap; hp-FEM, BindOptiX)",
    ),
    NodeType(
        "data.response_map",
        Category.DATA,
        (data_out("response_map", "response_map"),),
        _FILE,
        description="Emitter response table M from a file (EmitterResponseTable; hp-FEM)",
    ),
    NodeType(
        "data.occupancy",
        Category.DATA,
        (data_out("occupancy", "occupancy"),),
        _FILE,
        description="Occupancy B of the receptor states over time from a file (BindOptiX)",
    ),
    NodeType(
        "bindoptix.signal",
        Category.ANALYSIS,
        (
            data_in("excitation", "field_map"),
            data_in("response", "response_map"),
            data_in("dipoles", "Q"),
            data_in("occupancy", "occupancy"),
            data_out("rates", "result"),
        ),
        (
            ParamSpec(
                "epsilon_per_molar_cm",
                "number",
                80000.0,
                unit="1/(M·cm)",
                minimum=0.0,
                group="dye",
                description="molar extinction coefficient at excitation",
            ),
            ParamSpec(
                "quantum_yield",
                "number",
                0.9,
                minimum=0.0,
                maximum=1.0,
                group="dye",
                description="intrinsic quantum yield q₀ in the host medium",
            ),
            ParamSpec("excitation_um", "number", 0.633, unit="µm", minimum=0.0, group="dye"),
            ParamSpec("emission_um", "number", 0.67, unit="µm", minimum=0.0, group="dye"),
            ParamSpec(
                "power_mw",
                "number",
                1.0,
                unit="mW",
                minimum=0.0,
                description="incident power at the input port",
            ),
            ParamSpec(
                "anchors_um",
                "numbers",
                [[0.0, 0.0, 0.0]],
                unit="µm",
                description="anchor (x, y, z) per spot of the occupancy, sensor frame; "
                "one triple for all spots",
            ),
            ParamSpec(
                "regions",
                "text",
                "",
                description="label=region, … : detector regions of the response table "
                "bound to detectors (empty: the only region)",
            ),
            ParamSpec(
                "port_irradiance_w_m2_per_mw",
                "number",
                0.0,
                unit="W/m² per mW",
                minimum=0.0,
                group="map normalised to I₀",
                description="irradiance at the spot per mW at the input port (0: the map "
                "is per mW already)",
            ),
            ParamSpec(
                "port_irradiance_status",
                "choice",
                "P",
                options=("L", "K", "P"),
                group="map normalised to I₀",
                description="L literature, K calibration, P placeholder",
            ),
            ParamSpec(
                "port_irradiance_reference",
                "text",
                "",
                group="map normalised to I₀",
                description="source of L or K",
            ),
        ),
        description="Detected photon rate per receptor, spot, time and detector region "
        "(BindOptiX signal model from F, M, Q, B)",
    ),
)

EVALUATORS: dict[str, Evaluator] = {
    "biooptomd.layer": biooptomd_layer,
    "data.field_map": data_file,
    "data.response_map": data_file,
    "data.occupancy": data_file,
    "bindoptix.signal": bindoptix_signal,
}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    """Adds the D13 data types (`biooptomd.layer` is built in) and returns the evaluators."""
    for t in TYPES:
        registry.add(t)
    return dict(EVALUATORS)
