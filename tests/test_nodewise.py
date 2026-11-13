"""R25: node-wise execution in the Raytatouille backend (ADR ports, poses, branches; stage R25).

Each Raytatouille block passes the light of a branch along its own events; the detector images
what arrives. Expected values: the whole-path execution of R15/R18b/R19b (W9, W17, W18, W22),
itself checked against analytic references – node-wise must give the same image (up to rounding)
and the same power. Plus the choice of the mode and the checks of branches the system does not
know.
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
    validate_model,
)
from rtt_fusion.graph.execute import WITH_FIELDS, Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.model import Graph
from rtt_fusion.graph.topology import Step

rt = pytest.importorskip("raytatouille")

DATA = Path(__file__).parent / "data"


def _image(
    g: Graph, system: Any, detector: str, mode: str, transport: str = "auto", **shared: Any
) -> tuple[Image, Executor]:
    from rtt_fusion.graph import backend_rtt

    reg = builtin_registry()
    ev = dict(backend_rtt.EVALUATORS)
    ev.update(grating.register(reg))
    for link in g.links_to(detector):
        g.set_transport(link, transport)
    ex = Executor(
        g, reg, ev, {"system": system, "rtt_execution": mode, **shared}, available=WITH_FIELDS
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        im: Image = ex.evaluate(detector, "image")
    return im, ex


def _same(a: Image, b: Image, tol: float = 1e-9) -> None:
    assert np.abs(a.values - b.values).max() <= tol * a.values.max()
    assert b.power == pytest.approx(a.power, rel=tol)


def _michelson() -> dict[str, Any]:
    """The Michelson of W18: test mirror tilted by 2 mrad, EPD 1 mm."""
    d: dict[str, Any] = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    mirror = next(c for c in d["root"]["children"] if c["name"] == "test mirror")
    mirror["pose"]["rotation_deg"][0] = -90.0 + np.degrees(2e-3)
    d["aperture"] = {"type": "epd", "value": 1.0}
    return d


SHARED_W18: dict[str, Any] = {"rings": 20, "pixels": [64, 64], "gaussian_waist": 0.5}


# ---------------------------------------------------------------------------- W43


@pytest.mark.parametrize("transport", ["auto", "H"])
def test_w43_michelson_node_wise_is_whole_path(transport: str) -> None:
    pytest.importorskip("fieldtracing")
    d = _michelson()
    system = rt.System.from_json(json.dumps(d))
    images = []
    for mode in ("paths", "nodes"):
        g = graph_from_raytatouille(d)
        g.nodes["camera"].params["half_size"] = [0.3, 0.3]
        images.append(_image(g, system, "camera", mode, transport, **SHARED_W18)[0])
    _same(images[0], images[1])


def test_w43_michelson_with_vignetted_rays() -> None:
    # W9's Michelson: part of the start grid lies behind the tilted splitter and is lost along
    # the whole path – node-wise as well
    d = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    system = rt.System.from_json(json.dumps(d))
    a, _ = _image(graph_from_raytatouille(d), system, "camera", "paths", rings=30)
    b, ex = _image(graph_from_raytatouille(d), system, "camera", "nodes", rings=30)
    _same(a, b, 1e-12)
    assert len(b.groups) == 2  # node-wise: two input polarisations, same image
    assert ex.calls[("reference_mirror", "to_beam_splitter")] == 1


def test_w43_paraboloid_focus_node_wise() -> None:
    pytest.importorskip("fieldtracing")
    path = DATA / "paraboloid_f50.rtt.json"  # a path with "auto" events
    system = rt.load(str(path))
    images = []
    for mode in ("paths", "nodes"):
        g = graph_from_raytatouille(str(path))
        g.nodes["focus"].params["half_size"] = [0.1, 0.1]
        images.append(_image(g, system, "focus", mode, "H", rings=20, pixels=[31, 31])[0])
    _same(images[0], images[1])


@pytest.mark.parametrize("transport", ["auto", "H"])
def test_w43_linnik_node_wise(transport: str) -> None:
    if transport == "H":
        pytest.importorskip("fieldtracing")
    d = json.loads((DATA / "linnik_mini.rtt.json").read_text(encoding="utf-8"))
    system = rt.System.from_json(json.dumps(d))
    images = [
        _image(
            graph_from_raytatouille(d), system, "camera", mode, transport, rings=12, pixels=[48, 48]
        )[0]
        for mode in ("paths", "nodes")
    ]
    _same(images[0], images[1])


def test_w43_rays_lost_at_the_last_event_of_a_run() -> None:
    """Rays vignetted at the objective – a run of one event – are lost node-wise as along the
    whole path (Raytatouille: a ray lost at event j has count = j + 2; for the last event of a
    run that is count = slots, so "all slots" alone does not mean "passed")."""
    d = json.loads((DATA / "linnik_mini.rtt.json").read_text(encoding="utf-8"))
    objective = next(c for c in d["root"]["children"] if c["name"] == "objective")
    objective["surfaces"][0]["aperture"]["radius"] = 0.05  # cuts most of the beam
    system = rt.System.from_json(json.dumps(d))
    images = [
        _image(
            graph_from_raytatouille(d), system, "camera", mode, "auto", rings=12, pixels=[48, 48]
        )[0]
        for mode in ("paths", "nodes")
    ]
    full = json.loads((DATA / "linnik_mini.rtt.json").read_text(encoding="utf-8"))
    open_ = _image(
        graph_from_raytatouille(full),
        rt.System.from_json(json.dumps(full)),
        "camera",
        "paths",
        "auto",
        rings=12,
        pixels=[48, 48],
    )[0]
    assert images[0].values.sum() < 0.9 * open_.values.sum()  # the aperture does cut
    _same(images[0], images[1])


def test_w43_grating_sample_node_wise(tmp_path: Path) -> None:
    """W22/W23's Linnik with a phase grating as the sample (R19b split path) node-wise."""
    sys.path.insert(0, str(Path(__file__).parent))
    from test_grating_rtt import PERIOD_UM, PHI0, load_system, write_table

    table = write_table(
        tmp_path / "grating.json", {m: float(jv(m, PHI0)) for m in (-2, -1, 0, 1, 2)}, PERIOD_UM
    )
    d = load_system()
    system = rt.System.from_json(json.dumps(d))
    images = []
    for mode in ("paths", "nodes"):
        g = graph_from_raytatouille(d)
        replace_block(
            g,
            builtin_registry_with_grating(),
            "sample",
            Node("sample", "fusion.grating", params={"table": table}),
        )
        images.append(_image(g, system, "camera", mode, rings=15, pixels=[256, 8])[0])
    _same(images[0], images[1])


def builtin_registry_with_grating() -> Any:
    reg = builtin_registry()
    grating.register(reg)
    return reg


# ---------------------------------------------------------------------------- W44


def test_w44_branch_of_the_model_runs_node_wise() -> None:
    """A branch that is no path of the system (a renamed arm, only one arm) is computed
    node-wise automatically: the arm alone, as the whole path of that arm gives it."""
    d = _michelson()
    system = rt.System.from_json(json.dumps(d))
    m = model_from_raytatouille(d)
    arm = next(b for b in m.branches if b.name == "test arm")
    m.branches = [type(arm)("my arm", arm.steps)]
    assert not validate_model(m, builtin_registry())
    g = unroll(m, builtin_registry())
    g.nodes["camera"].params["half_size"] = [0.3, 0.3]
    auto, ex = _image(g, system, "camera", "auto", **SHARED_W18)
    assert ex.calls.get(("test_mirror", "to_beam_splitter")) == 1  # light passed block by block
    whole = model_from_raytatouille(d)
    whole.branches = [b for b in whole.branches if b.name == "test arm"]
    g2 = unroll(whole, builtin_registry())
    g2.nodes["camera"].params["half_size"] = [0.3, 0.3]
    path_image, _ = _image(g2, system, "camera", "paths", **SHARED_W18)
    _same(path_image, auto)


def test_w44_unknown_traversal() -> None:
    """A branch through a block in a way no system path takes (here: the source light straight
    to the camera through the splitter) is reported before anything is computed."""
    d = _michelson()
    m = model_from_raytatouille(d)
    m.branch(
        "ghost",
        Step("source", None, "beam_splitter"),
        Step("beam_splitter", "source", "camera"),
        Step("camera", "beam_splitter", None),
    )
    codes = {x.code for x in validate_model(m, builtin_registry())}
    assert "graph.branch_backend" in codes
