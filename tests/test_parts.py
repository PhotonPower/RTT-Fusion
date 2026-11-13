"""S4 of ADR Nachtrag 1: new Raytatouille elements from the library compute like imported ones.

Expected values are the images of the imported systems themselves (the graph without the kept
data computes with the file): an imported element swapped for a part of the library with the same
values (splitter, mirror, camera of the Michelson; the lens of a singlet) gives the same image; an
edit of the part (a tilt) gives the image of the same edit in the file. Plus the traversals of the
parts as their ports give them.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import Model, Node, builtin_registry, model_from_raytatouille, unroll
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.model import Graph
from rtt_fusion.graph.parts import part_events
from rtt_fusion.graph.topology import Branch, Connection, Step

rt = pytest.importorskip("raytatouille")

DATA = Path(__file__).parent / "data"
KEPT = ("rtt", "runs", "rtt_system", "rtt_tree", "rtt_paths")
SHARED: dict[str, Any] = {"rings": 20, "pixels": [48, 48], "gaussian_waist": 0.5}


def _image(g: Graph, d: dict[str, Any], detector: str) -> Image:
    from rtt_fusion.graph import backend_rtt

    g.nodes[detector].params["half_size"] = [0.3, 0.3]
    ex = Executor(
        g,
        builtin_registry(),
        backend_rtt.EVALUATORS,
        {"system": rt.System.from_json(json.dumps(d)), **SHARED},
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        im: Image = ex.evaluate(detector, "image")
    return im


def _from_file(d: dict[str, Any], detector: str) -> Image:
    m = model_from_raytatouille(d)
    g = unroll(m, builtin_registry())
    for n in g.nodes.values():
        for key in KEPT:
            n.params.pop(key, None)
    return _image(g, d, detector)


def swap(
    m: Model, old: str, kind: str, ports: dict[str, str], element: dict[str, Any] | None = None
) -> str:
    """Replaces an imported block by a part of the library with the element's values (as a user
    enters them): connections and branch steps move to the part's ports."""
    block = m.nodes.pop(old)
    ((_, original),) = block.params["rtt"].items()
    new = f"new_{old}"
    m.add(Node(new, f"rtt.part.{kind}", params={"rtt_new": element or original}))
    m.connections = [
        Connection(
            (new, ports[c.a[1]]) if c.a[0] == old else c.a,
            (new, ports[c.b[1]]) if c.b[0] == old else c.b,
            c.transport,
            c.transport_back,
        )
        for c in m.connections
    ]
    m.branches = [
        Branch(
            b.name,
            tuple(
                Step(
                    new,
                    ports.get(s.port_in or "", None) if s.port_in else None,
                    ports.get(s.port_out or "", None) if s.port_out else None,
                    s.channel,
                )
                if s.node == old
                else s
                for s in b.steps
            ),
        )
        for b in m.branches
    ]
    return new


def _michelson(tilt: float = 2e-3) -> dict[str, Any]:
    d: dict[str, Any] = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    mirror = next(c for c in d["root"]["children"] if c["name"] == "test mirror")
    mirror["pose"]["rotation_deg"][0] = -90.0 + np.degrees(tilt)
    d["aperture"] = {"type": "epd", "value": 1.0}
    return d


def _same(a: Image, b: Image, tol: float = 1e-12) -> None:
    assert np.abs(a.values - b.values).max() <= tol * a.values.max()
    assert b.power == pytest.approx(a.power, rel=tol)


# ---------------------------------------------------------------------------- W59


def test_w59_splitter_mirror_and_camera_from_the_library() -> None:
    d = _michelson()
    m = model_from_raytatouille(d)
    swap(
        m,
        "beam_splitter",
        "splitter",
        {
            "source": "front_a",
            "reference_mirror": "back_a",
            "test_mirror": "front_b",
            "camera": "back_b",
        },
    )
    swap(m, "test_mirror", "mirror", {"beam_splitter": "front"})
    camera = swap(m, "camera", "detector", {"beam_splitter": "in"})
    g = unroll(m, builtin_registry())
    assert {g.nodes[n].type_id for n in g.nodes} >= {
        "rtt.part.splitter",
        "rtt.part.mirror",
        "rtt.part.detector",
    }
    _same(_from_file(d, "camera"), _image(g, d, camera))


def test_w59_a_lens_from_the_library() -> None:
    base = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    d: dict[str, Any] = {
        **{k: v for k, v in base.items() if k not in ("root", "paths", "name", "aperture")},
        "name": "singlet",
        "aperture": {"type": "epd", "value": 4.0},
        "root": {
            "type": "assembly",
            "name": "singlet",
            "children": [
                {
                    "type": "stop",
                    "name": "stop",
                    "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 2.0}}],
                },
                {
                    "type": "lens",
                    "name": "L",
                    "pose": {"position": [0.0, 0.0, 5.0]},
                    "material": "CONST:1.5",
                    "surfaces": [
                        {
                            "id": "L.S1",
                            "shape": {"base": {"type": "conic", "radius": 60.0}},
                            "aperture": {"type": "circular", "radius": 6.0},
                        },
                        {
                            "id": "L.S2",
                            "pose": {"position": [0.0, 0.0, 3.0]},
                            "shape": {"base": {"type": "conic", "radius": -60.0}},
                            "aperture": {"type": "circular", "radius": 6.0},
                        },
                    ],
                },
                {
                    "type": "detector",
                    "name": "image",
                    "pose": {"position": [0.0, 0.0, 70.0]},
                    "surfaces": [
                        {
                            "id": "IMG",
                            "aperture": {
                                "type": "rectangular",
                                "half_width_x": 2.0,
                                "half_width_y": 2.0,
                            },
                        }
                    ],
                },
            ],
        },
        "paths": [
            {
                "name": "main",
                "events": [
                    {"surface": "STO"},
                    {"surface": "L.S1"},
                    {"surface": "L.S2"},
                    {"surface": "IMG"},
                ],
            }
        ],
    }
    m = model_from_raytatouille(d)
    swap(m, "L", "lens", {"stop": "front", "image": "back"})
    g = unroll(m, builtin_registry())
    _same(_from_file(d, "image"), _image(g, d, "image"))


def test_w59_an_edit_of_a_part_is_computed() -> None:
    d = _michelson(tilt=0.0)
    m = model_from_raytatouille(d)
    new = swap(m, "test_mirror", "mirror", {"beam_splitter": "front"})
    m.nodes[new].params["rtt_new"]["pose"]["rotation_deg"][0] = -90.0 + np.degrees(2e-3)
    edited = _image(unroll(m, builtin_registry()), d, "camera")
    _same(_from_file(_michelson(tilt=2e-3), "camera"), edited)


def test_w59_traversals_of_the_parts() -> None:
    assert part_events("mirror", ["m.M"], "front", "front") == [
        {"surface": "m.M", "kind": "reflect"}
    ]
    assert part_events("mirror", ["m.M"], "front", "") is None
    lens = ["l.S1", "l.S2"]
    assert [e["surface"] for e in part_events("lens", lens, "back", "front") or []] == [
        "l.S2",
        "l.S1",
    ]
    assert part_events("splitter", ["b.BS"], "front_a", "front_b") == [
        {"surface": "b.BS", "kind": "reflect"}
    ]
    assert part_events("splitter", ["b.BS"], "front_a", "back_a") == [
        {"surface": "b.BS", "kind": "transmit"}
    ]
    assert part_events("splitter", ["b.BS"], "front_a", "back_b") is None


def test_w59_a_new_part_is_in_the_exported_system() -> None:
    from rtt_fusion.graph.importer import system_from_model
    from rtt_fusion.graph.parts import PARTS

    m = model_from_raytatouille(_michelson())
    m.add(Node("extra", "rtt.part.lens", params={"rtt_new": PARTS["lens"][1]}))
    exported = system_from_model(m)
    assert exported["root"]["children"][-1]["name"] == "Lens"
    assert [s["id"] for s in exported["root"]["children"][-1]["surfaces"]] == [
        "extra.S1",
        "extra.S2",
    ]


# ---------------------------------------------------------------------------- W64


def test_w64_a_new_system_without_import_computes() -> None:
    from rtt_fusion.graph import backend_rtt, new_system, validate_model
    from rtt_fusion.graph.importer import system_from_model

    m = new_system("bench")
    reg = builtin_registry()
    assert not [x for x in validate_model(m, reg) if x.severity == "error"]
    d = system_from_model(m)
    assert d["name"] == "bench" and [c["type"] for c in d["root"]["children"]] == [
        "stop",
        "detector",
    ]
    assert not rt.validate(rt.System.from_json(json.dumps(d)))  # a valid Raytatouille file
    ex = Executor(unroll(m, reg), reg, backend_rtt.EVALUATORS, {"rings": 12, "pixels": [32, 32]})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        im: Image = ex.evaluate("detector", "image")
    # the collimated beam lies inside the window (the pixels at the rim of the disc miss ~4 %)
    assert im.power > 0.9 * im.ray_power > 0


def test_w64_a_lens_inserted_equals_the_saved_file() -> None:
    from rtt_fusion.graph import io, new_system
    from rtt_fusion.graph.importer import system_from_model
    from rtt_fusion.graph.parts import PARTS

    m = new_system()
    lens: dict[str, Any] = json.loads(json.dumps(PARTS["lens"][1]))
    lens["pose"] = {"position": [0.0, 0.0, 5.0], "rotation_deg": [0.0, 0.0, 0.0]}
    m.add(Node("lens", "rtt.part.lens", params={"rtt_new": lens}))
    # as in the editor: the connection stop–detector replaced by stop–lens–detector, the branch
    # recorded again through the lens
    m.connections = [c for c in m.connections if {c.a[0], c.b[0]} != {"stop", "detector"}]
    m.connect(("stop", "back"), ("lens", "front"))
    m.connect(("lens", "back"), ("detector", "in"))
    m.branches = [
        Branch(
            "main",
            (
                Step("source", None, "out"),
                Step("stop", "front", "back"),
                Step("lens", "front", "back"),
                Step("detector", "in", None),
            ),
        )
    ]
    m.nodes["source"].params["rtt_system"]["wavelengths"] = [{"um": 0.6328, "reference": True}]
    again = io.model_from_dict(json.loads(json.dumps(io.model_to_dict(m))), builtin_registry())
    d = system_from_model(again)  # the system as saved
    assert d["wavelengths"] == [{"um": 0.6328, "reference": True}]
    assert [e["surface"] for e in d["paths"][0]["events"]] == [
        "stop.STO",
        "lens.S1",
        "lens.S2",
        "detector.IMG",
    ]
    from_model = _image(unroll(again, builtin_registry()), d, "detector")
    imported = model_from_raytatouille(d)  # the saved file, imported again: ids from the names
    (camera,) = [n.id for n in imported.nodes.values() if n.type_id == "rtt.detector"]
    _same(_from_file(d, camera), from_model)
    assert from_model.power > 0
