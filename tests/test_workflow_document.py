"""The workflow document of the node editor (CAD ADR 0010, W1): model with sockets and ports,
library, diagnostics and plan; the editable fields from the providing plugin, kept by the engine."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from rtt_fusion.graph import model_from_raytatouille, section_to_group
from rtt_fusion.graph.workflow import (
    block_states,
    catalogue,
    kept_fields,
    model_document,
    model_id,
    registry,
)

DATA = Path(__file__).parent / "data"


def _linnik() -> Any:
    return model_from_raytatouille(json.loads((DATA / "linnik_mini.rtt.json").read_text("utf-8")))


def _describe(node: dict[str, Any]) -> list[dict[str, Any]]:
    return (
        [{"name": "/location", "type": "text", "label": "where"}]
        if "location" in node.get("params", {})
        else []
    )


def test_w66_document_of_a_model() -> None:
    reg, _ = registry()
    doc = model_document(_linnik(), reg, _describe)
    assert set(doc) == {"model", "types", "diagnostics", "plan"}
    assert not [d for d in doc["diagnostics"] if d["severity"] == "error"]
    assert doc["plan"]["links"] and doc["plan"]["connections"]  # planned per link and connection
    nodes = {n["id"]: n for n in doc["model"]["nodes"]}
    assert nodes["objective"]["sockets"] and nodes["objective"]["ports"]
    assert nodes["objective"]["fields"] == [{"name": "/location", "type": "text", "label": "where"}]
    assert "fields" not in nodes["source"]
    assert {t["type"] for t in doc["types"]} >= {"rtt.section", "rtt.source"}
    assert doc["types"] == catalogue(reg)
    json.dumps(doc)  # plain JSON for the channel to CAD


def test_w66_fields_are_kept_also_inside_groups() -> None:
    reg, _ = registry()
    d = json.loads((DATA / "raytatouille" / "feature_tour.rtt.json").read_text("utf-8"))
    m = model_from_raytatouille(d)
    section_to_group(m, "polarization_optics")
    first = model_document(m, reg, _describe)
    kept = kept_fields(first)
    assert {"polarizer", "quarter_wave_plate"} <= set(kept)  # blocks inside the group
    again = model_document(m, reg, None, kept)  # the engine: no provider at hand, fields kept
    assert again["model"] == first["model"]


def test_w66_block_states_per_model_block() -> None:
    raw = {"a@1": True, "a@2": False, "g/x": True, "g/y": True, "b": True, "other": True}
    assert block_states(raw, {"a", "g", "b"}) == {"a": False, "g": True, "b": True}
    assert block_states(raw, None) == raw
    assert model_id("g/x@3") == "g"


def test_w66_glass_catalogues_are_a_parameter_of_the_source() -> None:
    from rtt_fusion.graph import validate

    reg, _ = registry()
    m = _linnik()
    m.nodes["source"].params["rtt_catalogs"] = ["/somewhere/patent.agf"]
    from rtt_fusion.graph import unroll

    assert not [d for d in validate(unroll(m, reg), reg) if d.code == "graph.param_unknown"]
