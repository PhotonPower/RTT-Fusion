"""Poses of nodes (ADR ports, poses, branches, section 4; stage R26).

Convention of Raytatouille (`docs/architecture.md`, "Transformationen"): translation, then rotation
about the pivot; rotations intrinsic X -> Y -> Z, R = Rx·Ry·Rz;
p_parent = position + pivot + R·(p_child − pivot). Internally isometries (R, t); the six numbers
are only the form at interfaces.

A node without geometry of its own carries its pose as the parameter `pose`
({"position": [mm], "rotation": [rad], "pivot": [mm]}, default the identity); a node with geometry
(a Raytatouille block) takes its pose from the system. Both have the offset input `pose`: six
values (dx, dy, dz in mm; rx, ry, rz in rad) in the node's local system at its zero position,
rotation about its pivot (the reference point of the node: `point` or `center`). The offset is the
parameter `pose_offset` plus, per driver (sweep, replicates, optimisation linked to the input), a
parameter `pose_by_<driver>` – drivers add. Location = zero position ∘ offset.

`posed(node)` gives the node with its geometric parameters in world coordinates; the executor
hands that to the evaluators, so every backend computes in the world.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from rtt_fusion._types import FloatArray
from rtt_fusion.graph.model import Node

#: geometric parameters: points move with the full isometry, directions with its rotation
POINT_PARAMS = ("point", "center")
DIRECTION_PARAMS = ("normal", "direction", "axis", "up", "grating_vector", "polarisation")
#: names of the offset components (the `component` of a driver)
COMPONENTS = ("dx", "dy", "dz", "rx", "ry", "rz")


def rotation_xyz(angles: Any) -> FloatArray:
    """R = Rx(a)·Ry(b)·Rz(c) (intrinsic X -> Y -> Z), angles in rad."""
    a, b, c = (float(v) for v in angles)
    ca, sa, cb, sb, cc, sc = np.cos(a), np.sin(a), np.cos(b), np.sin(b), np.cos(c), np.sin(c)
    rx = np.array([[1, 0, 0], [0, ca, -sa], [0, sa, ca]])
    ry = np.array([[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]])
    rz = np.array([[cc, -sc, 0], [sc, cc, 0], [0, 0, 1]])
    return np.asarray(rx @ ry @ rz, dtype=float)


def angles_xyz(r: Any) -> FloatArray:
    """The angles (a, b, c) with R = Rx(a)·Ry(b)·Rz(c); b in [−π/2, π/2]."""
    m = np.asarray(r, dtype=float)
    b = float(np.arcsin(np.clip(m[0, 2], -1.0, 1.0)))
    if abs(m[0, 2]) < 1.0 - 1e-12:
        a = float(np.arctan2(-m[1, 2], m[2, 2]))
        c = float(np.arctan2(-m[0, 1], m[0, 0]))
    else:  # gimbal lock: put everything into a
        a = float(np.arctan2(m[2, 1], m[1, 1]))
        c = 0.0
    return np.array([a, b, c])


@dataclass(frozen=True)
class Isometry:
    """p -> R·p + t."""

    rotation: FloatArray
    translation: FloatArray

    @staticmethod
    def identity() -> Isometry:
        return Isometry(np.eye(3), np.zeros(3))

    @staticmethod
    def from_pose(
        position: Any = (0, 0, 0), rotation: Any = (0, 0, 0), pivot: Any = (0, 0, 0)
    ) -> Isometry:
        """Raytatouille's pose: p -> position + pivot + R·(p − pivot)."""
        r = rotation_xyz(rotation)
        c = np.asarray(pivot, dtype=float)
        return Isometry(r, np.asarray(position, dtype=float) + c - r @ c)

    def __matmul__(self, other: Isometry) -> Isometry:
        """Composition: (self ∘ other)(p) = self(other(p))."""
        return Isometry(
            self.rotation @ other.rotation, self.rotation @ other.translation + self.translation
        )

    def inverse(self) -> Isometry:
        rt = self.rotation.T
        return Isometry(rt, -rt @ self.translation)

    def point(self, p: Any) -> FloatArray:
        return np.asarray(self.rotation @ np.asarray(p, dtype=float) + self.translation)

    def vector(self, v: Any) -> Any:
        return self.rotation @ np.asarray(v)

    def is_identity(self, tol: float = 0.0) -> bool:
        return bool(
            np.abs(self.rotation - np.eye(3)).max() <= tol and np.abs(self.translation).max() <= tol
        )


def offset_of(node: Node) -> FloatArray:
    """The offset of a node: `pose_offset` plus every `pose_by_<driver>` (six values)."""
    total = np.zeros(6)
    for name, value in node.params.items():
        if name == "pose_offset" or name.startswith("pose_by_"):
            v = np.asarray(([float(x) for x in (value or [])] + [0.0] * 6)[:6])
            total = total + v
    return total


def pivot_of(node: Node, defaults: dict[str, Any] | None = None) -> FloatArray:
    """The reference point of the node in its local system (`point`, `center`, else 0)."""
    for name in POINT_PARAMS:
        if name in node.params:
            return np.asarray(node.params[name], dtype=float)
        if defaults and name in defaults:
            return np.asarray(defaults[name], dtype=float)
    return np.zeros(3)


def placement(node: Node, defaults: dict[str, Any] | None = None) -> Isometry:
    """World <- local of the node: its pose ∘ its offset (about its pivot)."""
    p = node.params.get("pose") or {}
    pose = Isometry.from_pose(
        p.get("position", (0, 0, 0)), p.get("rotation", (0, 0, 0)), p.get("pivot", (0, 0, 0))
    )
    o = offset_of(node)
    if not np.any(o):
        return pose
    return pose @ Isometry.from_pose(o[:3], o[3:], pivot_of(node, defaults))


def posed(
    node: Node, defaults: dict[str, Any] | None = None, outer: Isometry | None = None
) -> Node:
    """The node with its geometric parameters in world coordinates (`defaults`: the type's
    defaults of parameters the node does not set, e.g. a detector's `up`; `outer`: the placement
    of the groups it lies in, R28)."""
    iso = placement(node, defaults)
    if outer is not None:
        iso = outer @ iso
    if iso.is_identity():
        return node
    params = dict(node.params)
    for name in (*POINT_PARAMS, *DIRECTION_PARAMS):
        value = params.get(name, (defaults or {}).get(name))
        if value is None:
            continue
        if name in POINT_PARAMS:
            params[name] = [float(x) for x in iso.point(value)]
        elif name == "polarisation":
            v = np.asarray(value, dtype=complex)
            w = iso.vector(v.real) + 1j * iso.vector(v.imag)
            params[name] = [complex(x) if abs(complex(x).imag) > 0 else float(x.real) for x in w]
        else:
            params[name] = [float(x) for x in iso.vector(np.asarray(value, dtype=float))]
    return replace(node, params=params)


def bind_exported(node: Node, groups: list[Node]) -> Node:
    """A block inside groups with the exported parameters set (R28): `_bind` maps an own
    parameter ("focal_length" or "point[2]") to the name of the exported parameter of its group
    (`_parent`); a value set on the group wins over the block's own."""
    binding = node.params.get("_bind") or {}
    if not binding or not groups:
        return node
    params = dict(node.params)
    group = groups[0]
    for target, exported in binding.items():
        if exported not in group.params:
            continue
        value = group.params[exported]
        name, _, index = str(target).partition("[")
        if index:
            current = list(params.get(name) or [])
            k = int(index.rstrip("]"))
            while len(current) <= k:
                current.append(0.0)
            current[k] = value
            params[name] = current
        else:
            params[name] = value
    return replace(node, params=params)


__all__ = [
    "COMPONENTS",
    "DIRECTION_PARAMS",
    "POINT_PARAMS",
    "Isometry",
    "angles_xyz",
    "bind_exported",
    "offset_of",
    "pivot_of",
    "placement",
    "posed",
    "rotation_xyz",
]
