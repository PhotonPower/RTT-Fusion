"""R26: poses of nodes, the offset input `pose`, port frames and the checks of connections, the
operator store (ADR ports, poses, branches; acceptance PA7, PA8, PA9 and the diagnostics of
section 5).

Expected values are independent of the code under test: Raytatouille's compiled layout for the
pose convention, invariance under a rigid motion of the whole set-up (PA7), the fringe period
λ / sin 2α of two plane waves after a mirror tilt α (as W4) and the equality of a pose shift along
the normal with the drive of W35 (PA8), the count of operator preparations (PA9), and the
geometry of the connection checks constructed by hand.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from rtt_fusion.graph import (
    Graph,
    Model,
    Node,
    builtin_registry,
    graph_from_raytatouille,
    meta,
    optics,
)
from rtt_fusion.graph.execute import ExecutionError, Executor, ResultCache
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.pose import Isometry, angles_xyz, offset_of, posed, rotation_xyz
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.graph.topology import connection_diagnostics
from rtt_fusion.graph.validate import node_sockets

sys.path.insert(0, str(Path(__file__).parent))
from test_execute import fringe_period

DATA = Path(__file__).parent / "data"
WAVELENGTH_MM = 0.6328e-3


def _registry() -> tuple[NodeTypeRegistry, dict[str, Any]]:
    reg = builtin_registry()
    ev = dict(optics.register(reg))
    ev.update(meta.register(reg))
    return reg, ev


def michelson(theta: float = 0.0) -> Graph:
    """W4's Michelson with plane elements (mirror B tilted by θ about z)."""
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"center": [0, 0, -10], "direction": [0, 0, 1], "radius": 2.0, "rings": 30},
        )
    )
    g.add(Node("bs", "optics.beam_splitter", params={"point": [0, 0, 0], "normal": [-1, 0, 1]}))
    g.add(Node("mA", "optics.plane_mirror", params={"point": [0, 0, 20.0], "normal": [0, 0, -1]}))
    g.add(
        Node(
            "mB",
            "optics.plane_mirror",
            params={"point": [20.0, 0, 0], "normal": [-np.cos(theta), np.sin(theta), 0]},
        )
    )
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [-30, 0, 0],
                "normal": [1, 0, 0],
                "half_size": [2.5, 2.5],
                "pixels": [256, 256],
            },
        )
    )
    g.connect(("src", "out"), ("bs", "in"))
    g.connect(("bs", "t"), ("mA", "in"))
    g.connect(("mA", "out"), ("bs", "t_in"))
    g.connect(("bs", "r"), ("mB", "in"))
    g.connect(("mB", "out"), ("bs", "r_in"))
    g.connect(("bs", "merged"), ("cam", "in"))
    return g


def _image(g: Graph, node: str = "cam", socket: str = "image") -> Any:
    reg, ev = _registry()
    return Executor(g, reg, ev).evaluate(node, socket)


# ---------------------------------------------------------------------------- W45 convention


def test_w45_rotation_convention() -> None:
    rng = np.random.default_rng(45)
    for _ in range(20):
        a = rng.uniform(-np.pi, np.pi, 3)
        a[1] = rng.uniform(-1.4, 1.4)
        r = rotation_xyz(a)
        # intrinsic X -> Y -> Z (scipy: upper case = intrinsic)
        np.testing.assert_allclose(r, Rotation.from_euler("XYZ", a).as_matrix(), atol=1e-14)
        np.testing.assert_allclose(angles_xyz(r), a, atol=1e-12)
    p, q = (rng.uniform(-5, 5, 3) for _ in range(2))
    iso = Isometry.from_pose(rng.uniform(-5, 5, 3), rng.uniform(-1, 1, 3), rng.uniform(-1, 1, 3))
    np.testing.assert_allclose(iso.inverse().point(iso.point(p)), p, atol=1e-12)
    other = Isometry.from_pose(q, (0.1, 0.2, 0.3))
    np.testing.assert_allclose((iso @ other).point(p), iso.point(other.point(p)), atol=1e-12)


def test_w45_pose_as_raytatouille_compiles_it() -> None:
    rt = pytest.importorskip("raytatouille")
    d = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    angles, position, pivot = [17.0, -23.0, 41.0], [1.0, 2.0, 100.0], [0.3, -0.2, 0.5]
    m = next(c for c in d["root"]["children"] if c["name"] == "reference mirror")
    m["pose"] = {"position": position, "rotation_deg": angles, "pivot": pivot}
    compiled = rt.compile(rt.System.from_json(json.dumps(d)), rt.MaterialLibrary())
    layout = next(x for x in rt.layout.surfaces(compiled) if x.id == "M1")
    iso = Isometry.from_pose(position, np.radians(angles), pivot)
    np.testing.assert_allclose(np.asarray(layout.rotation), iso.rotation, atol=1e-14)
    np.testing.assert_allclose(np.asarray(layout.translation), iso.point(np.zeros(3)), atol=1e-12)


def test_w45_every_light_block_has_a_pose_input() -> None:
    reg, _ = _registry()
    for type_id in ("optics.plane_mirror", "optics.detector", "rtt.section", "rtt.mirror"):
        names = {s.name for s in node_sockets(Node("n", type_id), reg)}
        assert "pose" in names
    assert "pose" not in {s.name for s in node_sockets(Node("n", "meta.sweep"), reg)}


# ---------------------------------------------------------------------------- W46 (PA7)


@pytest.mark.parametrize("seed", range(3))
def test_w46_rigid_motion_of_the_whole_set_up(seed: int) -> None:
    rng = np.random.default_rng(seed)
    pose = {
        "position": [float(x) for x in rng.uniform(-50, 50, 3)],
        "rotation": [float(x) for x in rng.uniform(-np.pi, np.pi, 3)],
        "pivot": [float(x) for x in rng.uniform(-5, 5, 3)],
    }
    still = _image(michelson(theta=1e-3))
    g = michelson(theta=1e-3)
    for n in g.nodes.values():
        n.params["pose"] = pose
    moved: Image = _image(g)
    assert np.abs(moved.values - still.values).max() < 1e-9 * still.values.max()
    assert moved.power == pytest.approx(still.power, rel=1e-9)  # rounding of ~50 mm moves


# ---------------------------------------------------------------------------- W47 (PA8)


def test_w47_mirror_tilted_through_the_pose_input() -> None:
    """A sweep drives the pose input of mirror B (component rz): fringe period λ / sin 2α."""
    alpha = 1e-3
    g = michelson()
    g.add(Node("tilt", "meta.sweep", params={"values": [alpha], "component": "rz"}))
    g.connect(("tilt", "drive"), ("mB", "pose"))
    g.connect(("cam", "image"), ("tilt", "observe"))
    result = _image(g, "tilt", "results")
    im: Image = result.results[0]
    period = fringe_period(im.values, im.extent, 1.5)
    assert period == pytest.approx(WAVELENGTH_MM / np.sin(2 * alpha), rel=1e-6)
    assert "pose_by_tilt" not in g.nodes["mB"].params  # restored after the sweep
    # the same with the offset as a parameter, about the mirror's point (its pivot)
    g2 = michelson()
    g2.nodes["mB"].params["pose_offset"] = [0, 0, 0, 0, 0, alpha]
    np.testing.assert_allclose(_image(g2).values, im.values, atol=1e-12)


def test_w47_drivers_add() -> None:
    node = Node("m", "optics.plane_mirror", params={"pose_offset": [0, 0, 1e-4, 0, 0, 0]})
    node.params["pose_by_a"] = [0, 0, 2e-4, 0, 0, 0]
    np.testing.assert_allclose(offset_of(node), [0, 0, 3e-4, 0, 0, 0])
    one = posed(
        Node(
            "m",
            "optics.plane_mirror",
            params={"point": [0, 0, 0], "pose_offset": [0, 0, 3e-4, 0, 0, 0]},
        )
    )
    two = posed(Node("m", "optics.plane_mirror", params={"point": [0, 0, 0], **node.params}))
    np.testing.assert_allclose(one.params["point"], two.params["point"], atol=1e-18)


def _rtt_michelson(**params: Any) -> Image:
    rt = pytest.importorskip("raytatouille")
    from rtt_fusion.graph import backend_rtt

    d = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    g = graph_from_raytatouille(d)
    g.nodes["reference_mirror"].params.update(params)
    ex = Executor(
        g,
        builtin_registry(),
        backend_rtt.EVALUATORS,
        {"system": rt.System.from_json(json.dumps(d)), "rings": 30, "pixels": [320, 240]},
    )
    image: Image = ex.evaluate("camera", "image")
    return image


def test_w47_raytatouille_mirror_tilted_by_its_pose() -> None:
    alpha = 1e-3
    im = _rtt_michelson(pose_offset=[0, 0, 0, alpha, 0, 0])
    period = fringe_period(im.values, im.extent, 2.0)
    assert period == pytest.approx(WAVELENGTH_MM / np.sin(2 * alpha), rel=1e-2)


def test_w47_pose_dz_is_the_drive() -> None:
    # the reference mirror's local z is its normal: dz of the pose = W35's drive
    delta = WAVELENGTH_MM / 8
    by_pose = _rtt_michelson(pose_offset=[0, 0, delta, 0, 0, 0])
    by_drive = _rtt_michelson(drive=delta)
    np.testing.assert_allclose(by_pose.values, by_drive.values, atol=1e-12 * by_drive.values.max())
    still = _rtt_michelson()
    assert np.abs(by_pose.values - still.values).max() > 0.1 * still.values.max()


def test_w47_pose_of_a_partial_section_is_refused() -> None:
    rt = pytest.importorskip("raytatouille")
    from rtt_fusion.graph import backend_rtt

    d = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    g = graph_from_raytatouille(d)
    mirror = g.nodes["reference_mirror"]
    mirror.params["surfaces"] = [*mirror.params["surfaces"], "M2"]  # claims more than its node
    mirror.params["pose_offset"] = [0, 0, 1e-3, 0, 0, 0]
    ex = Executor(
        g,
        builtin_registry(),
        backend_rtt.EVALUATORS,
        {"system": rt.System.from_json(json.dumps(d))},
    )
    with pytest.raises(ExecutionError) as e:
        ex.evaluate("camera", "image")
    assert e.value.code == "exec.pose"


# ---------------------------------------------------------------------------- W48 (PA9)


def test_w48_moved_solver_keeps_its_operator() -> None:
    pytest.importorskip("cliffordbem")
    from rtt_fusion.graph import particle

    def run(dz: float, operators: ResultCache) -> tuple[Any, Executor]:
        reg = builtin_registry()
        ev = optics.register(reg)
        ev.update(particle.register(reg))
        g = Graph()
        g.add(
            Node(
                "src",
                "optics.collimated_source",
                params={"radius": 1.6, "rings": 8, "center": [0, 0, -1], "wavelengths_um": [0.6]},
            )
        )
        g.add(
            Node(
                "lens",
                "optics.aplanatic_lens",
                params={"focal_length": 2.0, "aperture_radius": 1.6},
            )
        )
        g.add(
            Node(
                "p",
                "fusion.particle",
                params={
                    "point": [0, 0, 2.0],
                    "radius_um": 0.2,
                    "refinement": 2,
                    "n_theta": 12,
                    "n_phi": 24,
                    "pose_offset": [0, 0, dz, 0, 0, 0],
                },
            )
        )
        g.connect(("src", "out"), ("lens", "in"))
        g.connect(("lens", "out"), ("p", "in"))
        ex = Executor(g, reg, ev, operators=operators)
        return ex.evaluate("p", "out"), ex

    store = ResultCache()
    run(0.0, store)
    prepared = len(list(store._items))
    moved, ex = run(2e-4, store)
    assert len(list(store._items)) == prepared  # no new preparation (operator store)
    assert not any("operator prepared" in line for line in ex.context.log)
    assert ex.calls[("p", "out")] == 1  # the result itself is computed again
    fresh, ex2 = run(2e-4, ResultCache())  # the same with a fresh operator
    assert any("operator prepared" in line for line in ex2.context.log)
    a, b = np.asarray(moved[0].amplitude), np.asarray(fresh[0].amplitude)
    # equal up to rounding (CliffordBem's parallel sums differ by ~1e-14 between preparations)
    assert np.abs(a - b).max() <= 1e-12 * np.abs(b).max()


# ---------------------------------------------------------------------------- W49 connections


def _chain(**detector: Any) -> Model:
    m = Model()
    m.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"center": [0, 0, -5], "radius": 1.0, "wavelengths_um": [0.6]},
        )
    )
    m.add(
        Node(
            "lens",
            "optics.aplanatic_lens",
            params={"point": [0, 0, 0], "focal_length": 10.0, "aperture_radius": 2.0},
        )
    )
    params = {"point": [0, 0, 10.0], "normal": [0, 0, -1], "half_size": [0.5, 0.5]}
    params.update(detector)
    m.add(Node("cam", "optics.detector", params=params))
    m.connect(("src", "out"), ("lens", "in"))
    m.connect(("lens", "out"), ("cam", "in"))
    m.branch("main", ("src", None, "out"), ("lens", "in", "out"), ("cam", "in", None))
    return m


def _diagnostics(m: Model) -> set[str]:
    reg, _ = _registry()
    return {d.code for d in connection_diagnostics(m, reg)}


def test_w49_connection_checks() -> None:
    assert _diagnostics(_chain()) == set()
    assert _diagnostics(_chain(normal=[0, 0, 1])) == {"graph.port_facing"}
    assert _diagnostics(_chain(point=[0, 0, -3.0])) == {"graph.port_behind"}
    assert _diagnostics(_chain(point=[5.0, 0, 10.0])) == {"graph.port_aperture"}
    assert _diagnostics(_chain(point=[0, 0, 0.002])) == {"graph.port_near_field"}  # 2 µm
    tilted = _chain(normal=[0, np.sin(1e-3), -np.cos(1e-3)])
    assert _diagnostics(tilted) == set()  # rays: no limit
    tilted.connections[1] = type(tilted.connections[1])(("lens", "out"), ("cam", "in"), "H")
    assert _diagnostics(tilted) == {"graph.port_tilt"}


def test_w49_frames_move_with_the_pose() -> None:
    m = _chain()
    pose = {"position": [3.0, -1.0, 7.0], "rotation": [0.3, -0.2, 1.1]}
    for n in m.nodes.values():
        n.params["pose"] = pose
    assert _diagnostics(m) == set()  # moved together: still facing
    m.nodes["cam"].params["pose"] = {"position": [0, 0, 0], "rotation": [np.pi, 0, 0]}
    assert "graph.port_facing" in _diagnostics(m) or "graph.port_behind" in _diagnostics(m)
