"""R9: conformance of the hp-FEM grating adapter, and the cross-check against FourierGrating.

hp-FEM and FourierGrating solve the same grating with entirely different methods (finite
elements with PML and Bloch faces vs. Fourier modal method) and different conventions (s/p
signs, power-normalised Jones matrices, references at the cell edge and the stack interfaces).
References: vector Fresnel (flat interface), the FourierGrating adapter (lamellar grating of
pilot K1, normal, oblique and conical incidence) and energy conservation.
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pytest

try:
    import hpfem
    from hpfem import meshing, units
except ImportError:
    if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
        raise
    pytest.skip("hpfem not installed", allow_module_level=True)

from test_adapter_fouriergrating import _fresnel_vector
from test_pilot_k1 import CONFIG

from rtt_fusion.adapters import hpfem as hpa
from rtt_fusion.ports import PlaneWaves

LAM = 0.6
K0 = 2 * np.pi / (LAM * 1e-3)
NM = units.nm
AIR, HI, SUB = 1, 2, 4


def _incident(port: Any, theta: float, phi: float, e: tuple[complex, complex]) -> PlaneWaves:
    k = float(np.real(port.index)) * K0
    return PlaneWaves(
        port,
        np.array([k * np.sin(theta) * np.cos(phi)]),
        np.array([k * np.sin(theta) * np.sin(phi)]),
        np.array([e], dtype=np.complex128),
    )


def _by_order(waves: PlaneWaves) -> dict[int, np.ndarray]:  # type: ignore[type-arg]
    assert waves.orders is not None
    return {int(o[0]): np.asarray(waves.field[i]) for i, o in enumerate(waves.orders)}


@pytest.fixture(scope="module")
def grating() -> Any:
    omega = units.angular_frequency(wavelength=LAM * 1e-6)
    air, hi = hpfem.Material.dielectric(1.0), hpfem.Material.dielectric(2.1)
    mats = {AIR: air, HI: hi, SUB: air}
    cell = meshing.UnitCell(
        period=700 * NM,
        y_bottom=-900 * NM,
        y_top=1150 * NM,
        slabs=[meshing.Slab(SUB, -900 * NM, 0.0)],
        shapes=[
            meshing.Shape("rectangle", HI, dict(x=0.0, y=0.0, width=300 * NM, height=250 * NM))
        ],
        background_tag=AIR,
    )
    mesh, _ = meshing.unit_cell_mesh(cell, meshing.element_sizes(mats, omega, p=2))
    stack = hpfem.LayerStack2D(air, [], air, 0.0)
    # ports: front at the ridge top (y = 250 nm), back at its foot (y = 0), as FourierGrating
    return hpa.prepare(hpa.Structure(mesh, mats, stack, 1.0, 1.0, 250 * NM, 0.0, 5), [LAM])


def test_flat_interface_equals_vector_fresnel() -> None:
    omega = units.angular_frequency(wavelength=LAM * 1e-6)
    air, glass = hpfem.Material.dielectric(1.0), hpfem.Material.dielectric(1.5)
    mats = {AIR: air, SUB: glass}
    cell = meshing.UnitCell(
        period=500 * NM,
        y_bottom=-900 * NM,
        y_top=900 * NM,
        slabs=[meshing.Slab(SUB, -900 * NM, 0.0)],
        shapes=[],
        background_tag=AIR,
    )
    mesh, _ = meshing.unit_cell_mesh(cell, meshing.element_sizes(mats, omega, p=2))
    h = hpa.prepare(
        hpa.Structure(mesh, mats, hpfem.LayerStack2D(air, [], glass, 0.0), 1.0, 1.5, 0.0, 0.0, 5),
        [LAM],
    )
    e = (0.6 + 0.2j, -0.3 + 0.5j)
    for theta, phi in ((0.0, 0.0), (0.5, 0.0), (0.5, 0.8)):
        inc = _incident(hpa.front_port(h, LAM), theta, phi, e)
        out = hpa.apply(h, inc)
        t_ref, r_ref = _fresnel_vector(np.array(e), float(inc.kx[0]), float(inc.ky[0]), 1.0, 1.5)
        assert np.max(np.abs(_by_order(out["front"])[0] - r_ref)) < 1e-12
        assert np.max(np.abs(_by_order(out["back"])[0] - t_ref)) < 1e-12


@pytest.mark.parametrize(("theta", "phi"), [(0.0, 0.0), (0.08, 0.0), (0.08, 0.7)])
def test_cross_check_with_fouriergrating(grating: Any, theta: float, phi: float) -> None:
    fga = pytest.importorskip("rtt_fusion.adapters.fouriergrating")
    pytest.importorskip("fouriergrating")
    fg = fga.prepare(fga.Structure(CONFIG), [LAM])
    e = (0.6 + 0.2j, -0.3 + 0.5j)
    a = hpa.apply(grating, _incident(hpa.front_port(grating, LAM), theta, phi, e))
    b = fga.apply(fg, _incident(fga.front_port(fg, LAM), theta, phi, e))
    for side in ("front", "back"):
        mine, theirs = _by_order(a[side]), _by_order(b[side])
        assert set(mine) == set(theirs) == {-1, 0, 1}
        for m, value in theirs.items():
            # strong orders (T, R0) on this mesh within 0.5 %, the weak R±1 (3 %) within 3 %
            tol = 0.03 if (side == "front" and m != 0) else 0.006
            assert np.max(np.abs(mine[m] - value)) < tol * np.max(np.abs(value))


def test_energy_over_orders(grating: Any) -> None:
    inc = _incident(hpa.front_port(grating, LAM), 0.08, 0.7, (0.8, 0.3j))
    out = hpa.apply(grating, inc)
    p = sum(float(np.sum(out[s].power_density())) for s in ("front", "back"))
    assert abs(p / float(inc.power_density()[0]) - 1) < 5e-3  # lossless; FEM residual ~1e-3
