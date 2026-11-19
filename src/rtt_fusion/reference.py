"""The reference for agents on the command line: `python -m rtt_fusion.reference`.

    python -m rtt_fusion.reference                  # the whole REFERENCE.md
    python -m rtt_fusion.reference --list           # section headings
    python -m rtt_fusion.reference --section 3      # one section (number, title or word)
    python -m rtt_fusion.reference --blocks         # the block catalogue
    python -m rtt_fusion.reference --blocks fg.grating optics.detector
    python -m rtt_fusion.reference --blocks --live  # catalogue from the installed registry

In Python: `rtt_fusion.reference()` or `rtt_fusion.reference("planner")`. After an explicit
`import rtt_fusion.reference` the name is this module, which can be called the same way.
"""

from __future__ import annotations

import sys
import types

from rtt_fusion._reference import get, main


class _CallableModule(types.ModuleType):
    def __call__(self, section: str | None = None) -> str:
        return get(section)


sys.modules[__name__].__class__ = _CallableModule

__all__ = ["get", "main"]

if __name__ == "__main__":
    raise SystemExit(main())
