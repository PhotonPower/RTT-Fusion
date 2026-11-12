"""S1 of ADR Nachtrag 1: after the import the model carries the whole system.

Expected values are the imported files themselves: `system_from_model` of a model fresh from
the import is the imported system (canonical JSON), for every test system of RTT-Fusion and three
reference systems of Raytatouille (tests/data/raytatouille: `feature_tour` with every model
element of the schema, `calcite_walkoff` with ordinary/extraordinary channels, `fresnel_bk7` with
an "auto" path); an edit in a block's data appears in the system; a branch the modeller removed is
no path any more.
"""

from __future__ import annotations

import glob
import json
from pathlib import Path
from typing import Any

import pytest

from rtt_fusion.graph import model_from_raytatouille
from rtt_fusion.graph.importer import event_channel, system_from_model

DATA = Path(__file__).parent / "data"
SYSTEMS = sorted(glob.glob(str(DATA / "*.rtt.json"))) + sorted(
    glob.glob(str(DATA / "raytatouille" / "*.rtt.json"))
)


def _canonical(d: object) -> str:
    return json.dumps(d, sort_keys=True)


@pytest.mark.parametrize("path", SYSTEMS, ids=lambda p: Path(p).name)
def test_w57_import_keeps_everything(path: str) -> None:
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    m = model_from_raytatouille(d)
    assert _canonical(system_from_model(m)) == _canonical(d)


def test_w57_every_tree_node_has_one_owner() -> None:
    d = json.loads((DATA / "raytatouille" / "feature_tour.rtt.json").read_text(encoding="utf-8"))
    m = model_from_raytatouille(d)
    places = [loc for n in m.nodes.values() for loc in n.params.get("rtt", {})]
    places += list(m.nodes["source"].params["rtt_tree"])
    assert len(places) == len(set(places))  # nothing twice


def test_w57_an_edit_in_a_block_is_in_the_system() -> None:
    d = json.loads((DATA / "linnik_mini.rtt.json").read_text(encoding="utf-8"))
    m = model_from_raytatouille(d)
    sample = next(n for n in m.nodes.values() if "sample" in n.id.lower())
    location, element = next(iter(sample.params["rtt"].items()))
    element["pose"] = {"position": [0.0, 0.0, 0.123]}
    exported = system_from_model(m)
    node: Any = exported["root"]
    for part in location.split("/")[2:]:  # /root/children/i/…
        node = node[int(part)] if part.isdigit() else node[part]
    assert node["pose"] == {"position": [0.0, 0.0, 0.123]}
    m.branches = []  # the modeller removed every branch: no path is left
    assert system_from_model(m)["paths"] == []


def test_w57_channels() -> None:
    assert event_channel({"surface": "G", "order": -1}) == "G:m-1"
    assert event_channel({"surface": "P.S1", "kind": "extraordinary"}) == "P.S1:extraordinary"
    assert event_channel({"surface": "BS", "kind": "reflect"}) == ""  # the ports tell it
