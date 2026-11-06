"""Write the reference files for the C++ port: reference/<case>.result.json.

python tools/export_reference.py              # all cases
python tools/export_reference.py hemisphere_quadrature
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from reference_cases import CASES

from rtt_fusion.results import FIELD_REFERENCE, to_json

OUT = Path(__file__).resolve().parents[1] / "reference"


def main(names: list[str]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name in names or list(CASES):
        path = OUT / f"{name}.result.json"
        path.write_text(to_json(FIELD_REFERENCE, CASES[name](), indent=1) + "\n", encoding="utf-8")
        print(f"{path.relative_to(OUT.parent)}: {path.stat().st_size / 1024:.1f} kB")


if __name__ == "__main__":
    main(sys.argv[1:])
