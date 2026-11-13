"""R2: BWPM adapter in the common port types.

Checks: identical to `bwpm.adapter` (the wrapper adds nothing), power through the common
`rtt_fusion.ports.power` conserved for a lossless stack, and the cross-check against the
FourierGrating adapter through the common types (G → S on the BWPM side).
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pytest

try:
    import bwpm  # noqa: F401
except ImportError:
    if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
        raise
    pytest.skip("bwpm not installed", allow_module_level=True)

from rtt_fusion.adapters import bwpm as bwa
from rtt_fusion.ports import GriddedField, power, to_plane_waves

LAM = 0.6
DZ = LAM * 1e-3 / 40
N = 80  # |k⊥| = k0·sqrt(mx² + my²)·40/N


def _handle(media: list[complex]) -> Any:
    unique = sorted(set(media), key=lambda c: (c.real, c.imag))
    labels = np.array([np.full((N, N), unique.index(m)) for m in media])
    return bwa.prepare(bwa.structure(N, N, DZ, DZ, DZ, tuple(unique), labels, 1.0 + 0j), [LAM])


def _plane_wave(handle: Any, mx: int, my: int, e: tuple[complex, complex]) -> GriddedField:
    x = (np.arange(N) - N // 2) * DZ
    kx, ky = 2 * np.pi * mx / (N * DZ), 2 * np.pi * my / (N * DZ)
    phase = np.exp(1j * (kx * x[None, :] + ky * x[:, None]))
    field = np.stack([e[0] * phase, e[1] * phase], axis=-1)
    return GriddedField(bwa.front_port(handle, LAM), DZ, DZ, field)


MEDIA = [2.1 + 0j] * 6 + [1.45 + 0j] * 10 + [1.5 + 0j]


def test_describe() -> None:
    d = bwa.describe()
    assert d["representation"] == "G" and d["name"] == "BWPM"


def test_wrapper_equals_bwpm_adapter() -> None:
    from bwpm import adapter as direct

    handle = _handle(MEDIA)
    inc = _plane_wave(handle, 1, 1, (0.5 + 0.1j, -0.2 + 0.4j))
    out = bwa.apply(handle, inc, iterations=10, backend="python")
    port = direct.Port(0.0, 1, 1.0, LAM, N, N, DZ, DZ)
    ref = direct.apply(handle, direct.PortField(port, inc.field), 10, "python")
    for side in ("front", "back"):
        assert np.array_equal(out[side].field, ref[side].field)
        assert out[side].port.direction == ref[side].port.direction


def test_energy_through_common_power() -> None:
    handle = _handle(MEDIA)
    inc = _plane_wave(handle, 1, 0, (0.3, 0.7j))
    out = bwa.apply(handle, inc, iterations=60)
    balance = (power(out["front"]) + power(out["back"])) / power(inc)
    assert abs(balance - 1.0) < 1e-12


def test_cross_check_with_fouriergrating_through_common_types() -> None:
    fga = pytest.importorskip("rtt_fusion.adapters.fouriergrating")
    pytest.importorskip("fouriergrating")
    handle = _handle(MEDIA)
    e = (0.5 + 0.1j, -0.2 + 0.4j)
    bw_out = bwa.apply(handle, _plane_wave(handle, 1, 1, e), iterations=60)
    bw_r = to_plane_waves(bw_out["front"], threshold=1e-9)
    assert len(bw_r.kx) == 1
    cfg = {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 2, "Nz": 0},
        "period_um": {"x": 0.7, "z": 0.0},
        "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}, "hi": {"n": 2.1}, "mid": {"n": 1.45}},
        "superstrate": "air",
        "substrate": "glass",
        "layers": [
            {"type": "homogeneous", "thickness_um": 10 * DZ * 1e3, "material": "mid"},
            {"type": "homogeneous", "thickness_um": 6 * DZ * 1e3, "material": "hi"},
        ],
    }
    fg_handle = fga.prepare(fga.Structure(cfg), [LAM])
    fg_out = fga.apply(
        fg_handle,
        bw_r.__class__(
            fga.front_port(fg_handle, LAM), bw_r.kx, bw_r.ky, np.array([e], dtype=np.complex128)
        ),
    )
    assert fg_out["front"].orders is not None
    j = next(i for i, o in enumerate(fg_out["front"].orders) if tuple(o) == (0, 0))
    assert abs(fg_out["front"].kx[j] - bw_r.kx[0]) < 1e-6
    assert np.max(np.abs(fg_out["front"].field[j] - bw_r.field[0])) < 1e-9


def test_l4_apply_many_equals_apply() -> None:
    handle = _handle(MEDIA)
    incidents = [_plane_wave(handle, mx, my, (0.5, 0.2j)) for mx, my in ((0, 0), (1, 1), (2, -1))]
    many = bwa.apply_many(handle, incidents, iterations=10, backend="python")
    for inc, b in zip(incidents, many, strict=True):
        a = bwa.apply(handle, inc, iterations=10, backend="python")
        for side in ("front", "back"):
            assert np.array_equal(a[side].field, b[side].field)
