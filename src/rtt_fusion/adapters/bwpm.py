"""Adapter for the BWPM (stage R2): `bwpm.adapter` in the common port types.

The BWPM uses the conventions of Raytatouille already (BWPM docs/konventionen.md), so this module
only translates the types: `rtt_fusion.ports.GriddedField` (representation G) in and out. Physics
and conformance of the BWPM itself: BWPM tests K1–K16.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from rtt_fusion.ports import GriddedField, Port


def _bw() -> Any:
    from bwpm import adapter

    return adapter


def describe() -> dict[str, Any]:
    d = dict(_bw().describe())
    d["representation"] = "G"
    d["adapter"] = "rtt_fusion.adapters.bwpm"
    d["conventions_of_solver"] = ["none (BWPM uses the conventions of Raytatouille)"]
    return d


def structure(*args: Any, **kwargs: Any) -> Any:
    """`bwpm.adapter.Structure`: nx, ny, dx, dy, dz, indices, labels, index_before, z0, mask."""
    return _bw().Structure(*args, **kwargs)


def prepare(structure: Any, wavelengths_um: list[float] | tuple[float, ...]) -> Any:
    return _bw().prepare(structure, wavelengths_um)


def front_port(handle: Any, wavelength_um: float, direction: int = 1) -> Port:
    s = handle.structure
    return Port(s.z0, direction, complex(s.index_before), wavelength_um)


def apply(
    handle: Any, incident: GriddedField, iterations: int = 3, backend: str = "auto"
) -> dict[str, GriddedField]:
    """Gridded field on the front port (+z) → reflected (front, −z) and transmitted (back, +z)."""
    bw = _bw()
    s = handle.structure
    p = incident.port
    if (incident.dx, incident.dy) != (s.dx, s.dy) or incident.shape != (s.ny, s.nx):
        raise ValueError("incident grid differs from the structure grid")
    bw_port = bw.Port(p.z, p.direction, complex(p.index), p.wavelength_um, s.nx, s.ny, s.dx, s.dy)
    out = bw.apply(handle, bw.PortField(bw_port, np.asarray(incident.field)), iterations, backend)
    result = {}
    for side in ("front", "back"):
        q = out[side].port
        result[side] = GriddedField(
            Port(q.z, q.direction, complex(q.index), q.wavelength_um), q.dx, q.dy, out[side].field
        )
    return result


def apply_many(
    handle: Any, incidents: list[GriddedField], iterations: int = 3, backend: str = "auto"
) -> list[dict[str, GriddedField]]:
    """`apply` for many incident fields (L4); one march each – BWPM has no assembly to share (a
    batch of fields in BWPM itself: plan L4, stage O3)."""
    return [apply(handle, inc, iterations, backend) for inc in incidents]


def estimate(handle: Any, iterations: int = 3, backend: str = "auto") -> dict[str, float]:
    return dict(_bw().estimate(handle, iterations, backend))
