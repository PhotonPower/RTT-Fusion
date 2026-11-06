"""Pilot K3 (R4): particle in the focus of an objective, backscatter collected by the objective.

Chain: focused beam as plane waves (`rtt_fusion.couplings.focus.richards_wolf`, NA 0.8,
x-polarised) → CliffordBem adapter (sphere R = 0.25 µm, n = 1.5, in the focus) → backscattered
spectrum at the front port → power within the NA. References (Fusion plan, K3):

- the focal field equals CliffordBem's own Richards–Wolf beam (`BeamField.focused`, f₀ → ∞) up to
  its power normalisation;
- the collected power fraction equals the one computed directly in CliffordBem with
  `BeamField.focused` and its far field (no adapter, no port spectra).
"""

from __future__ import annotations

import os
from typing import Any

import numpy as np
import pytest

try:
    import cliffordbem as cb
except ImportError:
    if os.environ.get("RTT_FUSION_REQUIRE_MODULES") == "1":
        raise
    pytest.skip("cliffordbem not installed", allow_module_level=True)

from rtt_fusion.adapters import cliffordbem as cba
from rtt_fusion.convert import hemisphere
from rtt_fusion.couplings.focus import richards_wolf
from rtt_fusion.ports import spectrum_power

LAM, R, N_P, NA = 0.6, 0.25, 1.5, 0.8
K_UM = 2 * np.pi / LAM
CENTER = (0.0, 0.0, 0.002)


@pytest.fixture(scope="module")
def setup() -> Any:
    m0 = cb.make_icosphere(6)
    mesh = cb.TriangleMesh(np.asarray(m0.points) * R, m0.triangles)
    handle = cba.prepare(cba.Structure(mesh, N_P**2, 1.0, CENTER, 0.0, 0.004), [LAM])
    beam = richards_wolf(cba.front_port(handle, LAM), NA, CENTER, (1.0, 0.0), 40, 80)
    reference = cb.BeamField.focused(
        cb.Medium(eps=1.0), K_UM, np.zeros(3), NA, 1e6, np.array([1, 0, 0], complex), 40, 80
    )
    return mesh, handle, beam, reference


def _beam_field(beam: Any, points_mm: np.ndarray) -> np.ndarray:  # type: ignore[type-arg]
    kz = np.real(beam.port.kz(beam.kx, beam.ky))
    ez = beam.ez()
    out = []
    for p in points_mm:
        ph = np.exp(1j * (beam.kx * p[0] + beam.ky * p[1] + kz * (p[2] - beam.port.z)))
        out.append([np.sum(beam.field[:, 0] * ph), np.sum(beam.field[:, 1] * ph), np.sum(ez * ph)])
    return np.array(out)


def test_focal_field_equals_cliffordbem_beam(setup: Any) -> None:
    _, _, beam, reference = setup
    pts_um = np.array([[0, 0, 0], [0.3, 0, 0], [0, 0.3, 0], [0, 0, 0.5], [0.2, -0.1, -0.3]])
    ours = _beam_field(beam, np.asarray(CENTER) + pts_um * 1e-3)
    theirs = np.asarray(reference.eval(pts_um)[0])
    scale = theirs[0, 0] / ours[0, 0]
    assert abs(scale.imag) < 1e-12 * abs(scale)  # same phase convention
    assert np.max(np.abs(ours * scale - theirs)) < 1e-10 * abs(theirs[0, 0])


def test_collected_backscatter_fusion_equals_direct(setup: Any) -> None:
    mesh, handle, beam, reference = setup
    out = cba.apply(handle, beam, 40, 80)
    fusion = spectrum_power(out["front"], NA) / spectrum_power(beam)
    # direct: CliffordBem only (power of BeamField = 1/2 (2π/k)² ∫|a|² dΩ, ε = μ = 1)
    medium = cb.Medium(eps=1.0)
    b = np.asarray(reference.project(mesh, medium))
    r = cba._problem(handle, LAM).solve_rhs(b, cb.SolveOptions(tol=1e-8))
    h_scat = np.asarray(r.h) - b
    dirs, w = hemisphere(40, 80, False)
    e = np.asarray(cb.far_field_E(mesh, h_scat, K_UM, 1.0, dirs))
    e = e - np.sum(e * dirs, axis=1, keepdims=True) * dirs
    inside = np.sqrt(1 - dirs[:, 2] ** 2) <= NA
    direct = np.sum(w[inside] * np.sum(np.abs(e[inside]) ** 2, axis=1)) / (2 * reference.power)
    assert abs(fusion / direct - 1) < 1e-6
    assert 0.01 < fusion < 0.1  # 3.7 % of the beam power returns into the objective
