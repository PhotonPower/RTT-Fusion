"""S2 of ADR Nachtrag 1: the Raytatouille backend computes with the system the workflow makes.

Expected values are independent of the generation: the images computed with the imported file
itself (the graph without the kept data falls back to shared["system"]); an edit in the workflow
must give the image of the same edit made in the system file; a branch of the modeller made of
known traversals becomes a path of the generated system and runs as a whole path.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.special import jv

from rtt_fusion.graph import (
    Node,
    builtin_registry,
    graph_from_raytatouille,
    grating,
    model_from_raytatouille,
    replace_block,
    unroll,
)
from rtt_fusion.graph.execute import WITH_FIELDS, Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.model import Graph

rt = pytest.importorskip("raytatouille")

DATA = Path(__file__).parent / "data"
KEPT = ("rtt", "runs", "rtt_system", "rtt_tree", "rtt_paths")


def _image(
    g: Graph, d: dict[str, Any], detector: str, transport: str = "auto", **shared: Any
) -> tuple[Image, Executor]:
    from rtt_fusion.graph import backend_rtt

    reg = builtin_registry()
    ev = dict(backend_rtt.EVALUATORS)
    ev.update(grating.register(reg))
    for link in g.links_to(detector):
        g.set_transport(link, transport)
    system = rt.System.from_json(json.dumps(d))
    ex = Executor(g, reg, ev, {"system": system, **shared}, available=WITH_FIELDS)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        im: Image = ex.evaluate(detector, "image")
    return im, ex


def _from_file(g: Graph) -> Graph:
    """The graph without the data the import keeps: the backend uses the file (before S2)."""
    for n in g.nodes.values():
        for key in KEPT:
            n.params.pop(key, None)
    return g


def _same(a: Image, b: Image, tol: float = 1e-12) -> None:
    assert np.abs(a.values - b.values).max() <= tol * max(a.values.max(), 1e-300)
    assert b.power == pytest.approx(a.power, rel=tol, abs=1e-300)


def _michelson(tilt: float = 2e-3) -> dict[str, Any]:
    d: dict[str, Any] = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    mirror = next(c for c in d["root"]["children"] if c["name"] == "test mirror")
    mirror["pose"]["rotation_deg"][0] = -90.0 + np.degrees(tilt)
    d["aperture"] = {"type": "epd", "value": 1.0}
    return d


W18: dict[str, Any] = {"rings": 20, "pixels": [48, 48], "gaussian_waist": 0.5}


# ---------------------------------------------------------------------------- W58


@pytest.mark.parametrize("transport", ["auto", "H"])
def test_w58_michelson_from_the_workflow(transport: str) -> None:
    if transport == "H":
        pytest.importorskip("fieldtracing")
    d = _michelson()
    images = []
    for g in (graph_from_raytatouille(d), _from_file(graph_from_raytatouille(d))):
        g.nodes["camera"].params["half_size"] = [0.3, 0.3]
        images.append(_image(g, d, "camera", transport, **W18)[0])
    _same(images[0], images[1])


@pytest.mark.parametrize(
    ("name", "detector", "transport", "shared"),
    [
        ("paraboloid_f50.rtt.json", "focus", "H", {"rings": 16, "pixels": [21, 21]}),
        ("linnik_mini.rtt.json", "camera", "auto", {"rings": 12, "pixels": [32, 32]}),
        ("ghost_plates.rtt.json", "image", "auto", {"rings": 12, "pixels": [24, 24]}),
    ],
)
def test_w58_systems_from_the_workflow(
    name: str, detector: str, transport: str, shared: dict[str, Any]
) -> None:
    if transport == "H":
        pytest.importorskip("fieldtracing")
    d = json.loads((DATA / name).read_text(encoding="utf-8"))
    images = []
    for g in (graph_from_raytatouille(d), _from_file(graph_from_raytatouille(d))):
        if detector == "focus":
            g.nodes[detector].params["half_size"] = [0.1, 0.1]
        images.append(_image(g, d, detector, transport, **shared)[0])
    _same(images[0], images[1])


def test_w58_grating_sample_from_the_workflow(tmp_path: Path) -> None:
    """R19b: the replaced sample keeps its element's data, the grating acts there."""
    sys.path.insert(0, str(Path(__file__).parent))
    from test_grating_rtt import PERIOD_UM, PHI0, load_system, write_table

    table = write_table(
        tmp_path / "grating.json", {m: float(jv(m, PHI0)) for m in (-2, -1, 0, 1, 2)}, PERIOD_UM
    )
    d = load_system()
    reg = builtin_registry()
    grating.register(reg)
    images = []
    for strip in (False, True):
        g = graph_from_raytatouille(d)
        if strip:
            _from_file(g)
        replace_block(g, reg, "sample", Node("sample", "fusion.grating", params={"table": table}))
        images.append(_image(g, d, "camera", rings=15, pixels=[128, 8])[0])
    _same(images[0], images[1])


def test_w58_an_edit_in_the_workflow_is_computed() -> None:
    """The test mirror tilted in the workflow gives the image of the same tilt in the file."""
    d = _michelson(tilt=0.0)
    m = model_from_raytatouille(d)
    mirror = m.nodes["test_mirror"]
    ((location, element),) = mirror.params["rtt"].items()
    element["pose"]["rotation_deg"][0] = -90.0 + np.degrees(2e-3)  # edited in the node's data
    g = unroll(m, builtin_registry())
    g.nodes["camera"].params["half_size"] = [0.3, 0.3]
    edited, _ = _image(g, d, "camera", **W18)  # the file still has the untilted mirror
    reference = _from_file(graph_from_raytatouille(_michelson(tilt=2e-3)))
    reference.nodes["camera"].params["half_size"] = [0.3, 0.3]
    expected, _ = _image(reference, _michelson(tilt=2e-3), "camera", **W18)
    _same(edited, expected)
    untilted = _from_file(graph_from_raytatouille(d))
    untilted.nodes["camera"].params["half_size"] = [0.3, 0.3]
    plain, _ = _image(untilted, d, "camera", **W18)
    assert np.abs(plain.values - edited.values).max() > 0.1 * plain.values.max()  # it matters


def test_w58_own_branch_is_a_path() -> None:
    """A branch of the modeller made of known traversals (a renamed arm) is a path of the
    generated system and runs as a whole path (not block by block)."""
    from rtt_fusion.graph import backend_rtt

    d = _michelson()
    m = model_from_raytatouille(d)
    arm = next(b for b in m.branches if b.name == "test arm")
    m.branches = [type(arm)("my arm", arm.steps)]
    g = unroll(m, builtin_registry())
    g.nodes["camera"].params["half_size"] = [0.3, 0.3]
    mine, ex = _image(g, d, "camera", **W18)
    assert not backend_rtt._nodewise(ex.context)
    assert [p["name"] for p in json.loads(ex.context.shared["generated_system"][1])["paths"]] == [
        "my arm"
    ]
    whole = model_from_raytatouille(d)
    whole.branches = [b for b in whole.branches if b.name == "test arm"]
    g2 = unroll(whole, builtin_registry())
    g2.nodes["camera"].params["half_size"] = [0.3, 0.3]
    _same(_image(g2, d, "camera", **W18)[0], mine)
