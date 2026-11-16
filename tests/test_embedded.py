"""Blocks hold their data instead of naming files (W75): a workflow document computes on its own.

Grating and stack used to name a FourierGrating configuration or a Jones table by path; a
document taken to another machine then failed with "cannot read …". They now hold the content
(`config`, `table_data`, `rtt_fusion.graph.embedded`). References:
- W75a the grating with its Jones table embedded gives exactly the orders of the same block
  reading the table file (whose values W19/W20 check against the grating equation) – no solver;
- W75b the grating with its FourierGrating configuration embedded gives exactly the order powers
  of the block reading the file (checked against FourierGrating's direct solution, W26);
- W75c the stack with its configuration embedded gives exactly the Fresnel powers of the block
  reading the file (ST0);
- W75d embedding: a readable file goes into the block and its path is cleared, also inside a
  group's variant, and the document keeps it; an unreadable one is reported and kept; relative
  CAD files are made absolute and reported; without file and copy the run says how to fix it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import Graph, Node, builtin_registry, embedded, grating, io, optics
from rtt_fusion.graph.execute import ExecutionError, Executor
from rtt_fusion.graph.topology import Group, Model, Variant
from rtt_fusion.results import load_json

REFERENCE = Path(__file__).parent.parent / "reference" / "jones_table_lamellar.result.json"
LAMBDA_UM = 0.6


def _orders(params: dict[str, Any], alpha_deg: float = 2.0) -> dict[int, tuple[Any, ...]]:
    """Directions, amplitudes and OPL of every order behind a grating with `params`."""
    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(grating.register(reg))
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": 0.01,
                "rings": 6,
                "center": [0, 0, -0.1],
                "wavelengths_um": [LAMBDA_UM],
                "field_angles_deg": [[alpha_deg, 0.0]],
                "polarisation": [0, 1, 0],
            },
        )
    )
    g.add(Node("gr", "fg.grating", params=params))
    g.connect(("src", "out"), ("gr", "in"))
    out = Executor(g, reg, ev).evaluate("gr", "out")
    return {
        int(b.path[-1].split("[")[1].rstrip("]")): (
            np.asarray(b.direction),
            np.asarray(b.amplitude),
            np.asarray(b.opl),
            np.asarray(b.prt),
        )
        for b in out
    }


def _same(a: dict[int, tuple[Any, ...]], b: dict[int, tuple[Any, ...]]) -> None:
    assert a.keys() == b.keys() and a
    for m in a:
        for x, y in zip(a[m], b[m], strict=True):
            np.testing.assert_array_equal(x, y)


def _fg_config() -> dict[str, Any]:
    _, d = load_json(REFERENCE.read_text(encoding="utf-8"))
    cfg: dict[str, Any] = json.loads(d["parameters"]["fouriergrating_config"])
    return cfg


# ------------------------------------------------------------------ evaluation from held data
def test_w75a_embedded_table_equals_the_table_file() -> None:
    doc = json.loads(REFERENCE.read_text(encoding="utf-8"))
    from_file = _orders({"table": str(REFERENCE)})
    held = _orders({"table_data": doc})
    _same(from_file, held)


def test_w75b_embedded_configuration_equals_the_file(tmp_path: Path) -> None:
    pytest.importorskip("fouriergrating")
    cfg = _fg_config()
    path = tmp_path / "lamellar.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    common = {"side": "reflected", "orders": [-1, 0, 1]}
    _same(
        _orders({**common, "fouriergrating": str(path)}),
        _orders({**common, "config": cfg}),
    )


def test_w75c_embedded_stack_equals_the_file(tmp_path: Path) -> None:
    pytest.importorskip("fouriergrating")
    from test_stack import _rays, tmm

    from rtt_fusion.graph import stack
    from rtt_fusion.graph.execute import Context

    cfg = {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 0},
        "period_um": {"x": 1.0},
        "materials": {"m_in": {"n": 1.0}, "m_out": {"n": 1.5}, "l0": {"n": 2.0}},
        "superstrate": "m_in",
        "substrate": "m_out",
        "layers": [{"type": "homogeneous", "thickness_um": 0.1, "material": "l0"}],
    }
    path = tmp_path / "stack.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    theta = np.radians(30.0)
    powers = []
    for params in ({"fouriergrating": str(path)}, {"config": cfg}):
        node = Node("st", stack.TYPE_ID, params=params)
        rays = _rays(theta, "s", 0.6)
        out = stack.stack(Context(Graph()), node, "front_out", {"front_in": [rays]})
        b = out[0]
        f = np.asarray(b.amplitude)[:, None] * b.field_direction()
        powers.append(float(np.sum(np.abs(f[0]) ** 2)))
    r_ref, _ = tmm(1.0, [(2.0, 0.1)], 1.5, theta, 0.6, "s")
    assert powers[0] == powers[1]
    assert powers[1] == pytest.approx(r_ref, abs=1e-11)


# ------------------------------------------------------------------ embedding
def test_w75d_embedding_moves_the_files_into_the_blocks(tmp_path: Path) -> None:
    cfg = _fg_config()
    path = tmp_path / "lamellar.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    inner = Model()
    inner.add(Node("inner_gr", "fg.grating", params={"table": str(REFERENCE)}))
    m = Model()
    m.add(Node("gr", "fg.grating", params={"fouriergrating": str(path), "side": "reflected"}))
    m.groups["g"] = Group({"a": Variant(inner, {})}, "a")
    done = embedded.embed_files(m)
    assert done.diagnostics == []
    assert len(done.stored) == 2
    assert m.nodes["gr"].params["config"] == cfg
    assert m.nodes["gr"].params["fouriergrating"] == ""
    held = inner.nodes["inner_gr"].params
    assert held["table"] == "" and held["table_data"]["format"] == "raytatouille-result"
    # a second time: nothing to do
    assert embedded.embed_files(m).stored == []


def test_w75d_the_document_keeps_the_data_without_paths(tmp_path: Path) -> None:
    cfg = _fg_config()
    path = tmp_path / "lamellar.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    reg = builtin_registry()
    grating.register(reg)
    m = Model()
    m.add(Node("gr", "fg.grating", params={"fouriergrating": str(path)}))
    embedded.embed_files(m)
    text = json.dumps(io.model_to_dict(m))
    path.unlink()
    assert str(tmp_path) not in text
    back = io.model_from_dict(json.loads(text), reg)
    assert back.nodes["gr"].params["config"] == cfg


def test_w75d_workflow_embeds_also_in_a_socket_graph(tmp_path: Path) -> None:
    from rtt_fusion.graph import workflow

    g = Graph()
    g.add(Node("gr", "fg.grating", params={"table": str(REFERENCE)}))
    done = workflow.embed_files(g)
    assert len(done.stored) == 1 and g.nodes["gr"].params["table"] == ""


def test_w75d_unreadable_file_is_reported_and_kept(tmp_path: Path) -> None:
    missing = str(tmp_path / "nowhere" / "grating.fg.json")
    m = Model()
    m.add(Node("gr", "fg.grating", params={"fouriergrating": missing}))
    done = embedded.embed_files(m)
    assert done.stored == []
    assert [(d.code, d.severity, d.location) for d in done.diagnostics] == [
        ("graph.file_unreadable", "warning", "gr")
    ]
    assert m.nodes["gr"].params == {"fouriergrating": missing}


def test_w75d_relative_cad_file_made_absolute_and_reported(tmp_path: Path) -> None:
    cfg = {**_fg_config(), "layers": [{"type": "cad", "file": "body.step", "thickness_um": 0.2}]}
    (tmp_path / "cad.json").write_text(json.dumps(cfg), encoding="utf-8")
    node = Node("gr", "fg.grating", params={"fouriergrating": str(tmp_path / "cad.json")})
    done = embedded.embed_node(node)
    assert node.params["config"]["layers"][0]["file"] == str(tmp_path / "body.step")
    assert [d.code for d in done.diagnostics] == ["graph.file_external"]


def test_w75d_without_file_and_copy_the_run_says_how_to_fix_it(tmp_path: Path) -> None:
    missing = str(tmp_path / "grating.fg.json")
    with pytest.raises(ExecutionError) as e:
        _orders({"fouriergrating": missing})
    assert e.value.code == "exec.param" and "Check workflow" in str(e.value)
