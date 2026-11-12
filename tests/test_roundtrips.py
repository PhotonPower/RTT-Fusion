"""R27: interfaces, the general splitter with four ports and round trips the modeller enters as
branches (ADR ports, poses, branches; acceptance PA3, PA4, ghosts in Raytatouille sections).

Expected values are independent of the code under test: the partial sums of the Airy series of a
plate (Fresnel coefficients at normal incidence, phase 2 k₀ n d per round trip), the properties of
a lossless reciprocal splitter (S unitary and symmetric), and Raytatouille's own trace of the same
ghost as an explicit path of the system (tests/data/ghost_plates.rtt.json, Raytatouille M4) with
its power ρ = R(B.S1)·R(A.S2).
"""

from __future__ import annotations

import copy
import json
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion._types import ComplexArray
from rtt_fusion.graph import (
    Model,
    Node,
    builtin_registry,
    graph_from_raytatouille,
    model_from_raytatouille,
    optics,
    unroll,
    validate_model,
)
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.registry import NodeTypeRegistry

DATA = Path(__file__).parent / "data"
N_PLATE, D_MM, LAM_UM = 1.5, 0.0123456, 0.6
R_F = (N_PLATE - 1) / (N_PLATE + 1)  # |r| of both faces at normal incidence
DELTA = 2 * np.pi / (LAM_UM * 1e-3) * N_PLATE * D_MM


def _registry() -> tuple[NodeTypeRegistry, dict[str, Any]]:
    reg = builtin_registry()
    return reg, dict(optics.register(reg))


def _source(center: list[float], direction: list[float]) -> Node:
    return Node(
        "src",
        "optics.collimated_source",
        params={
            "center": center,
            "direction": direction,
            "radius": 1.0,
            "rings": 6,
            "wavelengths_um": [LAM_UM],
            "polarisation": [0.0, 1.0, 0.0],  # transverse to every beam here
        },
    )


def _detector(node_id: str, point: list[float], normal: list[float], half: float = 0.5) -> Node:
    return Node(
        node_id,
        "optics.detector",
        params={"point": point, "normal": normal, "half_size": [half, half], "pixels": [8, 8]},
    )


def plate(terms: list[int], side: str) -> Model:
    """Plate n = 1.5, thickness d, from two interface blocks; a splitter in front sends the
    reflected light to a detector. Branch k: k round trips in the plate, entered by the modeller."""
    m = Model()
    m.add(_source([0, 0, -10], [0, 0, 1]))
    m.add(Node("bs", "optics.splitter", params={"point": [0, 0, -5], "normal": [1, 0, -1]}))
    m.add(
        Node(
            "I1",
            "optics.interface",
            params={"point": [0, 0, 0], "normal": [0, 0, -1], "n_front": 1.0, "n_back": N_PLATE},
        )
    )
    m.add(
        Node(
            "I2",
            "optics.interface",
            params={"point": [0, 0, D_MM], "normal": [0, 0, -1], "n_front": N_PLATE, "n_back": 1.0},
        )
    )
    m.add(_detector("T", [0, 0, 10], [0, 0, -1]))
    m.add(_detector("R", [-10, 0, -5], [1, 0, 0]))
    m.connect(("src", "out"), ("bs", "front_a"))
    m.connect(("bs", "back_a"), ("I1", "front"))
    m.connect(("I1", "back"), ("I2", "front"))
    m.connect(("I2", "back"), ("T", "in"))
    m.connect(("bs", "back_b"), ("R", "in"))
    into = [("src", None, "out"), ("bs", "front_a", "back_a")]
    trip = [("I2", "front", "front"), ("I1", "back", "back")]
    for k in terms:
        if side == "T":
            steps = [
                *into,
                ("I1", "front", "back"),
                *trip * k,
                ("I2", "front", "back"),
                ("T", "in", None),
            ]
        elif k == 0:
            steps = [*into, ("I1", "front", "front"), ("bs", "back_a", "back_b"), ("R", "in", None)]
        else:
            steps = [
                *into,
                ("I1", "front", "back"),
                *trip * (k - 1),
                ("I2", "front", "front"),
                ("I1", "back", "front"),
                ("bs", "back_a", "back_b"),
                ("R", "in", None),
            ]
        m.branch(f"{side}{k}", *steps)
    other = "R" if side == "T" else "T"  # the other detector sees the plain way (one branch)
    if other == "T":
        m.branch("T0", *into, ("I1", "front", "back"), ("I2", "front", "back"), ("T", "in", None))
    else:
        m.branch(
            "R0", *into, ("I1", "front", "front"), ("bs", "back_a", "back_b"), ("R", "in", None)
        )
    return m


def _irradiance(m: Model, detector: str) -> float:
    reg, ev = _registry()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        im: Image = Executor(unroll(m, reg), reg, ev).evaluate(detector, "image")
    return float(im.values[3:5, 3:5].mean())


def transmitted(terms: list[int]) -> complex:
    """Airy terms of the transmitted amplitude (power-normalised): t t' (r² e^{2iδ})^k e^{iδ}."""
    return complex(
        sum((1 - R_F**2) * (R_F**2 * np.exp(2j * DELTA)) ** k * np.exp(1j * DELTA) for k in terms)
    )


def reflected(terms: list[int]) -> complex:
    """−r, then t t' r e^{2iδ} (r² e^{2iδ})^(k−1)."""
    total = 0j
    for k in terms:
        total += (
            -R_F
            if k == 0
            else (1 - R_F**2) * R_F * np.exp(2j * DELTA) * (R_F**2 * np.exp(2j * DELTA)) ** (k - 1)
        )
    return complex(total)


# ---------------------------------------------------------------------------- W50 (PA3)


@pytest.mark.parametrize("n", [0, 1, 2, 5])
def test_w50_airy_partial_sums(n: int) -> None:
    terms = list(range(n + 1))
    # transmitted: the splitter passes |t|² = 1/2 of the irradiance 1
    assert _irradiance(plate(terms, "T"), "T") == pytest.approx(
        0.5 * abs(transmitted(terms)) ** 2, rel=1e-9
    )
    # reflected: the splitter passes t, then r' on the way back: 1/4
    assert _irradiance(plate(terms, "R"), "R") == pytest.approx(
        0.25 * abs(reflected(terms)) ** 2, rel=1e-9
    )


def test_w50_without_branch_n_exactly_term_n_less() -> None:
    terms = [0, 1, 2, 4]  # the modeller left out the second round trip
    assert _irradiance(plate(terms, "T"), "T") == pytest.approx(
        0.5 * abs(transmitted(terms)) ** 2, rel=1e-9
    )
    assert abs(transmitted(terms)) != pytest.approx(abs(transmitted([0, 1, 2, 3, 4])), rel=1e-6)


def test_w50_energy_of_many_round_trips() -> None:
    # 25 round trips entered as branches: R + T of the plate approaches 1 (the rest is the tail
    # (r²)^26 of the series, not summed by the framework)
    terms = list(range(26))
    t = _irradiance(plate(terms, "T"), "T") / 0.5
    r = _irradiance(plate(terms, "R"), "R") / 0.25
    assert t + r == pytest.approx(
        abs(transmitted(terms)) ** 2 + abs(reflected(terms)) ** 2, rel=1e-9
    )
    assert abs(t + r - 1.0) < 1e-11  # the tail (r²)^26 ~ 1e−37; the rest is rounding


# ---------------------------------------------------------------------------- W51 (PA4)

_PLACES: dict[
    str, tuple[list[float], list[float]]
] = {  # where light enters each port, and its direction (normal (1, 0, −1)/√2 at 0)
    "front_a": ([0, 0, -10], [0, 0, 1]),
    "front_b": ([10, 0, 0], [-1, 0, 0]),
    "back_a": ([0, 0, 10], [0, 0, -1]),
    "back_b": ([-10, 0, 0], [1, 0, 0]),
}
_OUT: dict[str, tuple[str, str]] = {  # light into a port leaves reflected and transmitted at
    "front_a": ("front_b", "back_a"),
    "front_b": ("front_a", "back_b"),
    "back_a": ("back_b", "front_a"),
    "back_b": ("back_a", "front_b"),
}


def _s_matrix(r: complex, t: complex) -> tuple[ComplexArray, dict[tuple[str, str], float]]:
    ports = list(_PLACES)
    s = np.zeros((4, 4), dtype=complex)
    powers: dict[tuple[str, str], float] = {}
    reg, ev = _registry()
    for j, port in enumerate(ports):
        m = Model()
        m.add(_source(*_PLACES[port]))
        m.add(
            Node(
                "bs",
                "optics.splitter",
                params={"point": [0, 0, 0], "normal": [1, 0, -1], "r": r, "t": t},
            )
        )
        m.connect(("src", "out"), ("bs", port))
        for out in _OUT[port]:
            where, direction = _PLACES[out]
            m.add(_detector(f"d_{out}", where, [-x for x in direction], half=1.2))  # whole beam
            m.connect(("bs", out), (f"d_{out}", "in"))
            m.branch(out, ("src", None, "out"), ("bs", port, out), (f"d_{out}", "in", None))
        g = unroll(m, reg)
        ex = Executor(g, reg, ev)
        for out in _OUT[port]:
            b = ex.evaluate("bs", f"{out}_out")[0]
            ray = int(np.argmin(np.linalg.norm(np.asarray(b.source_xy), axis=1)))
            k0 = 2 * np.pi / (LAM_UM * 1e-3)
            # amplitude relative to the light arriving at the splitter (same OPL for r and t)
            s[ports.index(out), j] = np.asarray(b.amplitude)[ray] * np.exp(
                1j * k0 * (np.asarray(b.opl)[ray] - 10.0)
            )
            # input power: irradiance 1 on the source grid, a hexapolar polygon of 6 rings
            # (36 corners on the edge): area 18 sin(2π/36) R²
            area = 18 * np.sin(2 * np.pi / 36) * 1.0**2
            powers[(port, out)] = ex.evaluate(f"d_{out}", "image").ray_power / area
    return s, powers


@pytest.mark.parametrize(
    ("r", "t"), [(np.sqrt(0.5), np.sqrt(0.5)), (0.6j, 0.8), (0.3, 0.4 + 0.866j)]
)
def test_w51_lossless_reciprocal_splitter(r: complex, t: complex) -> None:
    t = np.sqrt(1 - abs(r) ** 2) * t / abs(t)  # lossless: |r|² + |t|² = 1, phase of t kept
    s, powers = _s_matrix(r, t)
    np.testing.assert_allclose(s.conj().T @ s, np.eye(4), atol=1e-12)  # energy: S unitary
    np.testing.assert_allclose(s, s.T, atol=1e-12)  # reciprocity
    for port in _PLACES:  # all four ports lit in turn: the power leaves through two ports
        assert sum(p for (i, _), p in powers.items() if i == port) == pytest.approx(1.0, rel=1e-12)
    # phase relation of the lossless splitter: r t* + t r'* = 0
    r_front, t_ab = s[1, 0], s[2, 0]
    r_back = s[3, 2]
    assert abs(r_front * np.conj(t_ab) + t_ab * np.conj(r_back)) < 1e-12


# ---------------------------------------------------------------------------- W52 ghosts


def _rtt_image(
    d: dict[str, Any], model: Any = None, detector: str = "image", **shared: Any
) -> Image:
    rt = pytest.importorskip("raytatouille")
    from rtt_fusion.graph import backend_rtt

    reg = builtin_registry()
    g = unroll(model, reg) if model is not None else graph_from_raytatouille(d)
    ex = Executor(
        g,
        reg,
        backend_rtt.EVALUATORS,
        {"system": rt.System.from_json(json.dumps(d)), "rings": 12, "pixels": [32, 32], **shared},
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        im: Image = ex.evaluate(detector, "image")
    return im


def test_w52_ghost_branch_equals_raytatouille_ghost_path() -> None:
    d = json.loads((DATA / "ghost_plates.rtt.json").read_text(encoding="utf-8"))
    m = model_from_raytatouille(d)
    m.branches = []
    m.branch(
        "ghost",
        ("source", None, "stop"),
        ("stop", "source", "A"),
        ("A", "stop", "B"),
        ("B", "A", "A", "ghost:B.S1"),
        ("A", "B", "B", "ghost:A.S2"),
        ("B", "A", "image"),
        ("image", "B", None),
    )
    assert not validate_model(m, builtin_registry())
    node_wise = _rtt_image(d, m)
    explicit = copy.deepcopy(d)
    explicit["paths"] = [
        {
            "name": "ghost",
            "events": [
                {"surface": "STO"},
                {"surface": "A.S1"},
                {"surface": "A.S2"},
                {"surface": "B.S1", "kind": "reflect"},
                {"surface": "A.S2", "kind": "reflect"},
                {"surface": "B.S1"},
                {"surface": "B.S2"},
                {"surface": "IMG"},
            ],
        }
    ]
    whole = _rtt_image(explicit, rtt_execution="paths")
    assert np.abs(node_wise.values - whole.values).max() <= 1e-9 * whole.values.max()
    main = _rtt_image(d)
    rho = ((2.0 - 1) / 3.0) ** 2 * ((1.5 - 1) / 2.5) ** 2  # Raytatouille M4: Fresnel products
    assert node_wise.power / main.power == pytest.approx(rho, rel=1e-9)


def test_w52_ghost_checks() -> None:
    d = json.loads((DATA / "ghost_plates.rtt.json").read_text(encoding="utf-8"))
    m = model_from_raytatouille(d)
    m.branch(
        "odd but forward",
        ("source", None, "stop"),
        ("stop", "source", "A"),
        ("A", "stop", "B", "ghost:A.S2"),  # one reflection leaves back through "stop"
        ("B", "A", "image"),
        ("image", "B", None),
    )
    m.branch(
        "foreign surface",
        ("source", None, "stop"),
        ("stop", "source", "A"),
        ("A", "stop", "stop", "ghost:B.S1"),
        ("stop", "A", "source"),
        ("source", "stop", None),
    )
    codes = [x.code for x in validate_model(m, builtin_registry())]
    assert codes.count("graph.branch_backend") >= 2


def test_w52_ghost_through_curved_surfaces() -> None:
    """A ghost of a plano-convex singlet (Raytatouille M4 ghost_singlet): reflected at L.S2, then
    at L.S1; the surfaces passed on the way refract – equal to Raytatouille's own path with these
    events (a straight pass would not focus the ghost)."""
    d = json.loads((DATA / "raytatouille" / "ghost_singlet.rtt.json").read_text(encoding="utf-8"))
    m = model_from_raytatouille(d)
    lens = next(n for n in m.nodes.values() if "L.S1" in n.params.get("surfaces", []))
    stop = next(n for n in m.nodes.values() if "STO" in n.params.get("surfaces", []))
    image = next(n for n in m.nodes.values() if n.type_id == "rtt.detector")
    m.branches = []
    m.branch(
        "ghost",
        ("source", None, stop.id),
        (stop.id, "source", lens.id),
        (lens.id, stop.id, image.id, "ghost:L.S2,L.S1"),
        (image.id, lens.id, None),
    )
    assert not validate_model(m, builtin_registry())
    node_wise = _rtt_image(d, m, detector=image.id)
    explicit = copy.deepcopy(d)
    explicit["paths"] = [
        {
            "name": "ghost",
            "events": [
                {"surface": "STO"},
                {"surface": "L.S1"},
                {"surface": "L.S2", "kind": "reflect"},
                {"surface": "L.S1", "kind": "reflect"},
                {"surface": "L.S2"},
                {"surface": "IMG"},
            ],
        }
    ]
    whole = _rtt_image(explicit, detector=image.id, rtt_execution="paths")
    assert np.abs(node_wise.values - whole.values).max() <= 1e-9 * whole.values.max()
