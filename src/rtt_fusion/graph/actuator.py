"""Actuators: one value that sets a parameter of the Raytatouille system and moves blocks along.

An actuator (`fusion.actuator`) stands for a stage of the setup – a translation table under the
exposure head, a rotation mount, a focus drive. Its value (`value`, plus `drive` when a meta block
drives it through its `drive` input) does two things:

- `parameter`: the free parameter of that name of the Raytatouille system (the `rtt_system` of
  the `rtt.source`) is set to the value, within the range (`min`, `max`) its row gives. Every
  Raytatouille element whose pose refers to the parameter moves, and expressions over it
  (`Z_G2 = Z_G1 + DP`) follow.
- `moves`: blocks the stage carries, rigidly: entries "block.component" or
  "block.component*factor" (component dx, dy, dz in mm or rx, ry, rz), separated by commas. Each
  gets the offset factor·(value − `zero`) in that component as its own pose offset
  `pose_by_<actuator>` (drivers add, R26); a rotation takes the value in rad, or in degrees when
  the actuator's `unit` is "deg". `zero` is the value at which the carried blocks stand where the
  model puts them.

A block an actuator carries keeps its own pose (`keep_pose`): a block that replaced a Raytatouille
element is otherwise placed where that element is, and the move would count twice. This makes
the result the same whether the Raytatouille blocks run in the engine's process or behind a
provider (Raytatouille Fusion), where the placement does not reach the other blocks.

`apply_actuators(graph)` writes all of it into the graph (the unrolled working graph, never the
model); the executor calls it before every run (`Executor.check`) and a meta block after each
value it sets on an actuator. A muted actuator does nothing.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any

from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.model import Category, Graph, Node, NodeType, ParamSpec, data_in
from rtt_fusion.graph.pose import COMPONENTS
from rtt_fusion.graph.registry import NodeTypeRegistry

TYPE_ID = "fusion.actuator"
UNITS = ("mm", "deg", "rad", "µm", "")

_MOVE = re.compile(
    r"^(?P<block>[^\s.*]+)\.(?P<component>dx|dy|dz|rx|ry|rz)"
    r"(\s*\*\s*(?P<factor>[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?))?$"
)


@dataclass(frozen=True)
class Move:
    """One block the actuator carries: the offset component and the factor on the value."""

    block: str
    component: str
    factor: float = 1.0


def parse_moves(text: str, where: str = "") -> list[Move]:
    """The entries of `moves`: "block.component" or "block.component*factor", comma-separated."""
    out: list[Move] = []
    for entry in (e.strip() for e in re.split(r"[,;\n]", str(text or ""))):
        if not entry:
            continue
        m = _MOVE.match(entry)
        if m is None:
            raise ExecutionError(
                "exec.param",
                f"moves: {entry!r} – write block.component or block.component*factor "
                f"(component {', '.join(COMPONENTS)})",
                where,
            )
        out.append(Move(m["block"], m["component"], float(m["factor"] or 1.0)))
    return out


def _value(node: Node) -> float:
    return float(node.params.get("value", 0.0) or 0.0) + float(node.params.get("drive", 0.0) or 0.0)


def _blocks(graph: Graph, block: str) -> list[Node]:
    """The blocks of the working graph that stand for `block` of the model: its copies ("x@2")
    and, for a group, the blocks inside it ("g/x")."""
    return [
        n
        for nid, n in graph.nodes.items()
        if nid.split("@")[0] == block or nid.split("@")[0].startswith(block + "/")
    ]


def _set_parameter(graph: Graph, name: str, value: float, where: str) -> None:
    sources = [
        n
        for n in graph.nodes.values()
        if n.type_id == "rtt.source" and isinstance(n.params.get("rtt_system"), dict)
    ]
    found = False
    for src in sources:
        system = src.params["rtt_system"]
        rows = list(system.get("parameters") or [])
        for k, row in enumerate(rows):
            if not isinstance(row, dict) or row.get("name") != name:
                continue
            found = True
            if "expression" in row:
                raise ExecutionError(
                    "exec.param",
                    f"parameter {name} of the system is an expression ({row['expression']}); "
                    "an actuator sets a free parameter",
                    where,
                )
            low, high = row.get("min"), row.get("max")
            below = low is not None and value < float(low)
            if below or (high is not None and value > float(high)):
                raise ExecutionError(
                    "exec.param",
                    f"{name} = {value:g} lies outside its range {low} … {high} in the system",
                    where,
                )
            if row.get("value") == value:
                continue
            rows[k] = {**row, "value": value}
            src.params["rtt_system"] = {**system, "parameters": rows}  # the model keeps its own
    if not found:
        names = sorted(
            str(r.get("name"))
            for s in sources
            for r in s.params["rtt_system"].get("parameters") or []
            if isinstance(r, dict) and "expression" not in r
        )
        raise ExecutionError(
            "exec.param",
            f"the system has no parameter {name!r}"
            + (f" (free parameters: {', '.join(names)})" if names else ""),
            where,
        )


def apply_actuators(graph: Graph) -> None:
    """Sets the system parameters and the pose offsets of the carried blocks from the actuators
    of `graph` (see the module); raises `exec.param` for what does not exist."""
    for act in [n for n in graph.nodes.values() if n.type_id == TYPE_ID]:
        key = f"pose_by_{act.id}"
        moves = parse_moves(str(act.params.get("moves", "")), act.id)
        if act.muted:
            for n in graph.nodes.values():
                n.params.pop(key, None)
            continue
        value = _value(act)
        name = str(act.params.get("parameter", "") or "").strip()
        if name:
            _set_parameter(graph, name, value, act.id)
        delta = value - float(act.params.get("zero", 0.0) or 0.0)
        per_block: dict[str, list[float]] = {}
        for mv in moves:
            targets = _blocks(graph, mv.block)
            if not targets:
                raise ExecutionError("exec.param", f"moves: no block {mv.block!r}", act.id)
            k = COMPONENTS.index(mv.component)
            amount = mv.factor * delta
            if k >= 3 and str(act.params.get("unit", "")) == "deg":
                amount = math.radians(amount)
            for n in targets:
                per_block.setdefault(n.id, [0.0] * 6)[k] += amount
        for n in graph.nodes.values():
            if n.id in per_block:
                if n.params.get(key) != per_block[n.id]:
                    n.params[key] = per_block[n.id]
                if (
                    "rtt" in n.params
                    and not n.type_id.startswith("rtt.")
                    and not n.params.get("keep_pose")
                ):
                    n.params["keep_pose"] = True  # carried: placed by the stage, not the system
            elif key in n.params:  # no longer carried
                del n.params[key]


def actuator(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> Any:
    """An actuator has no output: the executor applies it to the graph before the run."""
    raise ExecutionError("exec.socket_unknown", "an actuator has no output", node.id)


TYPES = (
    NodeType(
        TYPE_ID,
        Category.META,
        (data_in("drive", "drive", "drive"),),
        (
            ParamSpec("value", "number", 0.0, description="position of the stage"),
            ParamSpec(
                "unit",
                "choice",
                "mm",
                options=UNITS,
                description="unit of the value; deg: rotations of the carried blocks in degrees",
            ),
            ParamSpec(
                "parameter",
                "text",
                "",
                description="free parameter of the Raytatouille system set to the value, e.g. DX",
            ),
            ParamSpec(
                "moves",
                "text",
                "",
                description="blocks carried along: block.component or block.component*factor, "
                "comma-separated (dx dy dz in mm, rx ry rz), e.g. G1.dx, G2.dx",
            ),
            ParamSpec(
                "zero",
                "number",
                0.0,
                description="value at which the carried blocks stand where the model puts them",
            ),
            ParamSpec(
                "drive", "number", 0.0, description="added to the value (set by a meta block)"
            ),
        ),
        description="Actuator: a stage that sets a system parameter and carries blocks along",
    ),
)

actuator.representations = frozenset()  # type: ignore[attr-defined]  # no light

EVALUATORS: dict[str, Evaluator] = {TYPE_ID: actuator}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    for t in TYPES:
        registry.add(t)
    return dict(EVALUATORS)


__all__ = ["TYPE_ID", "Move", "apply_actuators", "parse_moves", "register"]
