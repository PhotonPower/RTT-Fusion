"""The workflow document of the node editor (CAD ADR 0010, W1): model with sockets and ports,
library, diagnostics and plan; the editable fields from the providing plugin, kept by the engine."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from rtt_fusion.graph import model_from_raytatouille, section_to_group
from rtt_fusion.graph.model import Category, NodeType, light_in, light_out
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
    assert set(doc) == {"model", "types", "diagnostics", "plan", "providers"}
    assert doc["providers"] == {}  # every block type is built in
    assert not [d for d in doc["diagnostics"] if d["severity"] == "error"]
    assert doc["plan"]["links"] and doc["plan"]["connections"]  # planned per link and connection
    nodes = {n["id"]: n for n in doc["model"]["nodes"]}
    assert nodes["objective"]["sockets"] and nodes["objective"]["ports"]
    assert nodes["objective"]["fields"] == [{"name": "/location", "type": "text", "label": "where"}]
    assert "fields" not in nodes["source"]
    assert {t["type"] for t in doc["types"]} >= {"rtt.section", "rtt.source", "optics.thin_lens"}
    assert "fusion.condensed" in registry()[1]  # the evaluator of condensed groups (R29)
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


def test_w70_provider_records_and_a_missing_plugin() -> None:
    """CAD ADR 0010, D7: the document records which plugin provides a block type its model uses;
    without that plugin the block stays with its data and the check says which plugin is missing."""
    reg, _ = registry()
    m = _linnik()
    grating = NodeType(
        "fg.grating",
        Category.DOMAIN,
        (light_in("in", "R"), light_out("out", "R")),
        description="a grating",
    )
    reg.add(grating)
    m.nodes["sample"].type_id = "fg.grating"
    now = {"fg.grating": {"plugin": "fouriergrating", "version": "0.2"}}
    doc = model_document(m, reg, None, None, now)
    assert doc["providers"] == now
    assert (
        next(t for t in doc["types"] if t["type"] == "fg.grating")["provider"] == "fouriergrating"
    )
    assert "provider" not in next(t for t in doc["types"] if t["type"] == "rtt.section")
    # the plugin is gone: its type is unknown, the record of the document says whose it was
    bare, _ = registry()
    again = model_document(m, bare, None, None, {}, doc["providers"])
    missing = [d for d in again["diagnostics"] if d["code"] == "graph.provider_missing"]
    assert (
        len(missing) == 1
        and missing[0]["location"] == "sample"
        and missing[0]["severity"] == "error"
    )
    assert "fouriergrating" in missing[0]["message"] and "0.2" in missing[0]["message"]
    assert not [d for d in again["diagnostics"] if d["code"] == "graph.node_type_unknown"]
    assert again["providers"] == now  # kept for the next time
    assert next(n for n in again["model"]["nodes"] if n["id"] == "sample")["type"] == "fg.grating"
    # never recorded: an unknown type as before
    plain = model_document(m, bare)
    assert [d["code"] for d in plain["diagnostics"]] == ["graph.node_type_unknown"]
