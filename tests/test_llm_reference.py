"""The reference for agents (src/rtt_fusion/REFERENCE.md, shipped in the wheel).

- The block catalogue in it is the one generated from the registry now
  (tools/gen_llm_reference.py), and every registered block type says what it needs.
- Every example marked ```python run``` runs; one with `needs=a,b` is skipped without those
  modules (RTT_FUSION_REQUIRE_MODULES=1: fails instead).
- The file is package data (importlib.resources), `rtt_fusion.reference()` and
  `python -m rtt_fusion.reference` give it or a section of it.
"""

from __future__ import annotations

import importlib.util
import os
import re
import subprocess
import sys
from importlib import resources
from pathlib import Path

import pytest

import rtt_fusion
from rtt_fusion import _reference

FENCE = re.compile(r"^```python run(?P<options>[^\n]*)\n(?P<code>.*?)^```", re.M | re.S)
KNOWN_MODULES = {"fieldtracing", "fouriergrating", "cliffordbem", "raytatouille", "specklebem"}


def _examples() -> list[tuple[int, tuple[str, ...], str]]:
    doc = _reference.text()
    out = []
    for m in FENCE.finditer(doc):
        line = doc.count("\n", 0, m.start()) + 1
        needs = tuple(
            name
            for option in m.group("options").split()
            if option.startswith("needs=")
            for name in option.removeprefix("needs=").split(",")
            if name
        )
        out.append((line, needs, m.group("code")))
    return out


EXAMPLES = _examples()


def test_block_catalogue_is_up_to_date() -> None:
    assert (
        _reference.generated() == _reference.catalogue_markdown()
    ), "REFERENCE.md: block catalogue out of date – run python tools/gen_llm_reference.py"


def test_every_block_type_says_what_it_needs() -> None:
    from rtt_fusion.graph import workflow

    reg, _ = workflow.registry()
    assert {t.type_id for t in reg.types()} == set(_reference.NEEDS)


def test_examples_are_found() -> None:
    assert len(EXAMPLES) >= 10
    assert {n for _, needs, _ in EXAMPLES for n in needs} <= KNOWN_MODULES
    assert any(not needs for _, needs, _ in EXAMPLES)  # some run with numpy/scipy alone


@pytest.mark.parametrize(
    ("line", "needs", "code"), EXAMPLES, ids=[f"line {line}" for line, _, _ in EXAMPLES]
)
def test_example_runs(
    line: int, needs: tuple[str, ...], code: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    missing = [name for name in needs if importlib.util.find_spec(name) is None]
    if missing:
        if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
            raise RuntimeError(f"REFERENCE.md line {line} needs {', '.join(missing)}")
        pytest.skip(f"needs {', '.join(missing)}")
    monkeypatch.chdir(tmp_path)  # examples may write files
    exec(compile(code, f"REFERENCE.md:{line}", "exec"), {"__name__": "__reference__"})


def test_reference_is_package_data() -> None:
    text = resources.files("rtt_fusion").joinpath("REFERENCE.md").read_text(encoding="utf-8")
    assert text.startswith("# RTT-Fusion reference for agents")
    assert text == rtt_fusion.reference()
    assert "rtt_fusion.reference()" in (rtt_fusion.__doc__ or "").splitlines()[0]


def test_sections() -> None:
    assert rtt_fusion.reference("1.4").startswith("### 1.4 Planner")
    assert rtt_fusion.reference("planner").startswith("### 1.4 Planner")
    assert rtt_fusion.reference("2").startswith("## 2 Recipes")
    assert "### 2.11" in rtt_fusion.reference("2")  # a section includes its subsections
    entry = rtt_fusion.reference("fg.grating")
    assert entry.startswith("### fg.grating") and "\n### " not in entry
    with pytest.raises(KeyError, match="sections"):
        rtt_fusion.reference("no such section")
    # comments in code blocks are no headings
    titles = [title for _, _, title in _reference.headings()]
    assert not any(t.startswith(("Raytatouille Fusion document", "FourierGrating")) for t in titles)


def test_command_line() -> None:
    def run(*args: str) -> str:
        done = subprocess.run(
            [sys.executable, "-m", "rtt_fusion.reference", *args],
            capture_output=True,
            text=True,
            encoding="utf-8",
            check=True,
        )
        return done.stdout

    assert run() == rtt_fusion.reference()
    assert "1.4 Planner" in run("--list")
    assert run("--section", "1.4") == rtt_fusion.reference("1.4")
    blocks = run("--blocks", "optics.detector", "fusion.convert")
    assert blocks.startswith("### optics.detector") and "### fusion.convert" in blocks
    assert run("--blocks") == _reference.generated()
    assert run("--blocks", "--live") == _reference.catalogue_markdown()


def test_explicitly_imported_module_is_callable() -> None:
    import rtt_fusion.reference as module

    assert module("1.4") == _reference.get("1.4")  # type: ignore[operator]
    assert rtt_fusion.reference("1.4") == _reference.get("1.4")
