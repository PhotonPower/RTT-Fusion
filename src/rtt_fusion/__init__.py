"""RTT-Fusion. Agents: read rtt_fusion.reference() or run python -m rtt_fusion.reference.

Prototype of the Raytatouille Fusion: wave-optical solvers coupled through ports in a node
workflow (`rtt_fusion.graph`). The compact reference for language models and agents ships in
the package as REFERENCE.md: `print(rtt_fusion.reference())`, one section with
`rtt_fusion.reference("planner")`, on the command line `python -m rtt_fusion.reference --list`.

Conventions: mm, µm (wavelengths), rad; fields ∝ exp(i(k·r − ωt)); Fresnel Convention A.
"""

__version__ = "0.3.0"


def reference(section: str | None = None) -> str:
    """The reference for agents (REFERENCE.md in the package), or one section of it: a number
    ("3"), a title, or a block type id ("fg.grating"). Raises KeyError listing the sections."""
    from rtt_fusion._reference import get

    return get(section)
