"""Write the generated block catalogue of src/rtt_fusion/REFERENCE.md from the registry.

python tools/gen_llm_reference.py          # rewrite the section between the markers
python tools/gen_llm_reference.py --check  # exit 1 if it is out of date

The catalogue comes from `rtt_fusion._reference.catalogue_markdown` (the workflow's registry, its
evaluators and those of the Raytatouille backend); it needs none of the optional solver modules.
"""

from __future__ import annotations

import sys
from pathlib import Path

from rtt_fusion._reference import catalogue_markdown, generated, replace_generated

PATH = Path(__file__).resolve().parents[1] / "src" / "rtt_fusion" / "REFERENCE.md"


def main(argv: list[str]) -> int:
    doc = PATH.read_text(encoding="utf-8")
    body = catalogue_markdown()
    if "--check" in argv:
        if generated(doc) != body:
            print(f"{PATH.name}: block catalogue out of date – run python {__file__}")
            return 1
        return 0
    new = replace_generated(doc, body)
    if new != doc:
        PATH.write_text(new, encoding="utf-8")
    print(f"{PATH}: {len(new.splitlines())} lines, {len(new.encode()) / 1024:.1f} kB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
