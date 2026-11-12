"""R29: condensing a group into plane-wave tables (ADR ports, poses, branches, section 7;
acceptance PA5 condensed).

Expected values are independent of the code under test: the same group resolved into its blocks
(the plate of W50 from two interfaces with three round trips as inner branches) – equal up to
rounding where the incident direction lies on the table's grid, converging with the grid between
its points (the tolerance of the table, bilinear: error ~ step²); the table is built once and kept
when the group moves (operator store, as PA9); exactly the inner branches the model names are
applied; a group that is not laterally invariant (a lens) is refused.
"""

from __future__ import annotations

import json
import warnings
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import (
    Group,
    Model,
    Node,
    Variant,
    builtin_registry,
    condense,
    model_from_dict,
    model_to_dict,
    optics,
    unroll,
)
from rtt_fusion.graph.execute import ExecutionError, Executor, ResultCache
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.registry import NodeTypeRegistry

N_PLATE, D_MM, LAM_UM = 1.5, 0.0123456, 0.6


def _registry() -> tuple[NodeTypeRegistry, dict[str, Any]]:
    reg = builtin_registry()
    ev = dict(optics.register(reg))
    ev.update(condense.register(reg))
    return reg, ev


def plate_model(
    condensed: dict[str, Any] | None,
    direction: tuple[float, float, float] = (0.0, 0.0, 1.0),
    polarisation: tuple[float, float, float] = (0.0, 1.0, 0.0),
    trips: tuple[int, ...] = (0, 1, 2),
    pose: dict[str, Any] | None = None,
) -> Model:
    m = Model()
    m.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "center": [0, 0, -5],
                "direction": list(direction),
                "radius": 1.0,
                "rings": 6,
                "wavelengths_um": [LAM_UM],
                "polarisation": list(polarisation),
            },
        )
    )
    inner = Model()
    inner.add(
        Node(
            "I1",
            "optics.interface",
            params={"point": [0, 0, 0], "normal": [0, 0, -1], "n_front": 1.0, "n_back": N_PLATE},
        )
    )
    inner.add(
        Node(
            "I2",
            "optics.interface",
            params={"point": [0, 0, D_MM], "normal": [0, 0, -1], "n_front": N_PLATE, "n_back": 1.0},
        )
    )
    inner.connect(("I1", "back"), ("I2", "front"))
    for k in range(3):
        inner.branch(
            f"T{k}",
            ("I1", "front", "back"),
            *[("I2", "front", "front"), ("I1", "back", "back")] * k,
            ("I2", "front", "back"),
        )
    Group(
        {"plate": Variant(inner, {"in": ("I1", "front"), "out": ("I2", "back")})},
        "plate",
        condensed,
    ).add_to(m, Node("plate", "group", params={"pose": pose} if pose else {}))
    m.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, 10],
                "normal": [0, 0, -1],
                "half_size": [3, 3],
                "pixels": [24, 24],
            },
        )
    )
    m.connect(("src", "out"), ("plate", "in"))
    m.connect(("plate", "out"), ("cam", "in"))
    for k in trips:
        m.branch(
            f"k{k}", ("src", None, "out"), ("plate", "in", "out", f"T{k}"), ("cam", "in", None)
        )
    return m


def _image(m: Model, operators: ResultCache | None = None) -> tuple[Image, Executor]:
    reg, ev = _registry()
    ex = Executor(unroll(m, reg), reg, ev, operators=operators)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        im: Image = ex.evaluate("cam", "image")
    return im, ex


TABLE = {"sigma_max": 0.2, "count": 5}
ON_GRID = [  # incident directions on the grid of TABLE; s or p (the chain of interfaces is
    # exact for pure s or p light, the table for any – mixed light would compare two methods)
    ((0.0, 0.0, 1.0), "s"),  # normal incidence
    ((0.1, 0.0, float(np.sqrt(0.99))), "s"),
    ((0.1, 0.0, float(np.sqrt(0.99))), "p"),
    ((0.2, -0.1, float(np.sqrt(0.95))), "s"),
    ((0.2, -0.1, float(np.sqrt(0.95))), "p"),
]


@pytest.mark.parametrize(("direction", "kind"), ON_GRID)
def test_w55_condensed_equals_open_on_the_grid(
    direction: tuple[float, float, float], kind: str
) -> None:
    k = np.asarray(direction)
    s = np.cross(k, [0.0, 0.0, 1.0])
    s = s / np.linalg.norm(s) if np.linalg.norm(s) > 0 else np.array([0.0, 1.0, 0.0])
    vec = s if kind == "s" else np.cross(s, k)
    pol = (float(vec[0]), float(vec[1]), float(vec[2]))
    open_, _ = _image(plate_model(None, direction, pol))
    closed, ex = _image(plate_model(TABLE, direction, pol))
    assert np.abs(closed.values - open_.values).max() < 1e-10 * open_.values.max()
    assert sum("built" in line for line in ex.context.log) == 3  # one table per inner branch


def test_w55_between_the_grid_points_the_table_converges() -> None:
    d = (0.137, 0.05, float(np.sqrt(1 - 0.137**2 - 0.05**2)))
    open_, _ = _image(plate_model(None, d))
    lit = open_.values > 0.5 * open_.values.max()  # inside the beam, away from its edge
    centre = open_.values[lit].mean()
    errors = []
    for count in (5, 9, 17, 33):
        closed, _ = _image(plate_model({"sigma_max": 0.2, "count": count}, d))
        errors.append(abs(closed.values[lit].mean() - centre) / centre)
    assert errors[-1] < 1e-4
    assert errors[3] < errors[2] < errors[1]  # bilinear: the error falls with the step


def test_w55_table_kept_when_the_group_moves() -> None:
    store = ResultCache()
    _image(plate_model(TABLE), store)
    pose = {"position": [0.3, -0.2, 0.5], "rotation": [0.0, 0.0, 0.4]}  # a rotation about z too
    moved, ex = _image(plate_model(TABLE, pose=pose), store)
    assert not any("built" in line for line in ex.context.log)  # the operator store
    open_moved, _ = _image(plate_model(None, pose=pose))
    assert np.abs(moved.values - open_moved.values).max() < 1e-10 * open_moved.values.max()


def test_w55_exactly_the_named_inner_branches() -> None:
    open_, _ = _image(plate_model(None, trips=(0, 2)))
    closed, ex = _image(plate_model(TABLE, trips=(0, 2)))
    assert np.abs(closed.values - open_.values).max() < 1e-10 * open_.values.max()
    assert sum("built" in line for line in ex.context.log) == 2  # T1 is never used


def test_w55_a_lens_is_not_condensed() -> None:
    m = plate_model(None)
    inner = Model()
    inner.add(Node("L", "optics.thin_lens", params={"point": [0, 0, 0], "focal_length": 50.0}))
    inner.branch("T0", ("L", "in", "out"))
    m.groups["plate"] = Group(
        {"lens": Variant(inner, {"in": ("L", "in"), "out": ("L", "out")})}, "lens", TABLE
    )
    m.branches = [m.branches[0]]
    with pytest.raises(ExecutionError) as e:
        _image(m)
    assert e.value.code == "exec.condense"


def test_w55_file_round_trip() -> None:
    m = plate_model(TABLE)
    reg, _ = _registry()
    again = model_from_dict(json.loads(json.dumps(model_to_dict(m))), reg)
    assert again.groups["plate"].condense == TABLE
    np.testing.assert_array_equal(_image(again)[0].values, _image(m)[0].values)
