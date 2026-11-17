"""A system that does not compile is an execution error with what it lacks, not a traceback: an
unknown glass names its catalogue and where to add it (a document whose catalogue file is
missing on this machine)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from rtt_fusion.graph.execute import ExecutionError

rt = pytest.importorskip("raytatouille")

from rtt_fusion.graph.backend_rtt import _compile  # noqa: E402

DATA = Path(__file__).parent / "data"


def _with_material(material: str) -> object:
    text = (DATA / "ghost_plates.rtt.json").read_text(encoding="utf-8")
    text = text.replace('"CONST:1.5"', json.dumps(material), 1)
    return rt.System.from_json(text)


def test_an_unknown_catalogue_names_it_and_the_fix() -> None:
    with pytest.raises(ExecutionError) as e:
        _compile(_with_material("DE10316415:N-BK7"), rt.MaterialLibrary())
    assert e.value.code == "exec.compile"
    assert "DE10316415" in str(e.value) and "Glass catalogues" in str(e.value)


def test_a_compiling_system_is_unchanged() -> None:
    compiled = _compile(_with_material("CONST:1.5"), rt.MaterialLibrary())
    assert compiled is not None
