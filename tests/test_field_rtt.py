"""Fields with the Raytatouille backend (R18b): H on the link into a detector of an imported system.

References, independent of the code under test:
- W17 paraboloid f = 100 mm, D = 2 mm (Raytatouille reference m6/paraboloid_f50): on the axis a
  paraboloid is free of aberrations, so its focus is the Airy pattern [2 J₁(v)/v]²,
  v = k·NA·r, NA = sin(2·arctan(D/4f)), with the peak P·π NA²/λ² for the power P through the
  pupil (unit irradiance over the unit pupil: P = π). The ray image (R) cannot resolve it.
- W18 Michelson (Raytatouille reference m0, test mirror tilted by 2 mrad) with a Gaussian beam:
  over the H segment (splitter -> camera, 50 mm ≪ z_R ≈ 310 mm) the field image equals the ray
  image (W9) up to diffraction, with the same power.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pytest
from scipy.special import j1

from rtt_fusion.graph import builtin_registry, graph_from_raytatouille
from rtt_fusion.graph.execute import WITH_FIELDS, Executor
from rtt_fusion.graph.imaging import Image

rt = pytest.importorskip("raytatouille")
pytest.importorskip("fieldtracing")

DATA = Path(__file__).parent / "data"


def run(graph: object, system: object, node: str, **shared: object) -> tuple[Image, Executor]:
    from rtt_fusion.graph import backend_rtt

    ex = Executor(
        graph,  # type: ignore[arg-type]
        builtin_registry(),
        backend_rtt.EVALUATORS,
        {"system": system, **shared},
        available=WITH_FIELDS,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        image: Image = ex.evaluate(node, "image")
    return image, ex


def test_w17_paraboloid_focus_is_airy() -> None:
    path = DATA / "paraboloid_f50.rtt.json"
    system = rt.load(str(path))
    g = graph_from_raytatouille(str(path))
    half, n, wavelength_mm = 0.1, 61, 0.5876e-3
    g.nodes["focus"].params["half_size"] = [half, half]
    rays, _ = run(g, system, "focus", rings=30, pixels=[n, n])
    assert rays.values.max() == 0.0  # the rays meet in less than a pixel
    g.set_transport(g.links[-1], "H")
    image, ex = run(g, system, "focus", rings=30, pixels=[n, n])
    assert len(image.groups) == 2  # unpolarised source: two polarisations, incoherent
    u = (np.arange(n) + 0.5) / n * 2 * half - half
    uu, vv = np.meshgrid(u, u)
    na = np.sin(2 * np.arctan(2.0 / (4 * 100.0)))
    v = 2 * np.pi / wavelength_mm * na * np.hypot(uu, vv)
    airy = np.where(v > 0, (2 * j1(v) / np.where(v > 0, v, 1.0)) ** 2, 1.0)
    assert np.abs(image.values / image.values.max() - airy).max() < 2e-3
    peak = np.pi * (np.pi * na**2 / wavelength_mm**2)
    assert image.values.max() == pytest.approx(peak, rel=1e-2)
    assert any("H propagation" in line for line in ex.context.log)


def test_w18_michelson_fields_like_rays() -> None:
    d = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    mirror = next(c for c in d["root"]["children"] if c["name"] == "test mirror")
    mirror["pose"]["rotation_deg"][0] = -90.0 + np.degrees(2e-3)
    d["aperture"] = {"type": "epd", "value": 1.0}
    g = graph_from_raytatouille(d)
    system = rt.System.from_json(json.dumps(d))
    g.nodes["camera"].params["half_size"] = [0.3, 0.3]
    shared = {"rings": 20, "pixels": [64, 64], "gaussian_waist": 0.5}
    rays, _ = run(g, system, "camera", **shared)
    for link in g.links_to("camera"):
        g.set_transport(link, "H")
    field, _ = run(g, system, "camera", **shared)
    lit = rays.values > 0.1 * rays.values.max()
    assert np.abs(field.values - rays.values)[lit].max() < 0.12 * rays.values.max()
    assert field.power == pytest.approx(rays.power, rel=2e-2)


def test_auto_link_into_a_focus_computes_fields() -> None:
    """With `auto_fields` a detector on an "auto" link notices that the rays meet near it (W17's
    focus) and computes fields from the last surface on: the image of the link set to H."""
    path = DATA / "paraboloid_f50.rtt.json"
    system = rt.load(str(path))
    g = graph_from_raytatouille(str(path))
    half, n = 0.1, 31
    g.nodes["focus"].params["half_size"] = [half, half]
    auto, ex = run(g, system, "focus", rings=30, pixels=[n, n], auto_fields=True)
    assert any("focuses near the detector" in line for line in ex.context.log)
    g.set_transport(g.links[-1], "H")
    fields, _ = run(g, system, "focus", rings=30, pixels=[n, n])
    assert auto.values.max() > 0.0
    assert np.abs(auto.values - fields.values).max() <= 1e-12 * fields.values.max()


def test_auto_link_away_from_a_focus_keeps_rays() -> None:
    """A collimated beam on the camera (W18's Michelson): no caustic, the auto link stays R."""
    d = json.loads((DATA / "michelson.rtt.json").read_text(encoding="utf-8"))
    d["aperture"] = {"type": "epd", "value": 1.0}
    g = graph_from_raytatouille(d)
    system = rt.System.from_json(json.dumps(d))
    g.nodes["camera"].params["half_size"] = [0.3, 0.3]
    shared = {"rings": 12, "pixels": [32, 32], "gaussian_waist": 0.5}
    rays, _ = run(g, system, "camera", **shared)
    auto, ex = run(g, system, "camera", auto_fields=True, **shared)
    assert not any("focuses near the detector" in line for line in ex.context.log)
    assert np.array_equal(auto.values, rays.values)


def test_a_window_wider_than_rayleigh_sommerfeld_takes_debye() -> None:
    """Rayleigh–Sommerfeld computes 64² points λ/(8 NA) apart (0.47 mm for W17's focus): for a
    2 mm window the plan is Debye, with its estimated error in the log, and the Airy pattern's
    rings reach far beyond that grid instead of being cut off."""
    path = DATA / "paraboloid_f50.rtt.json"
    system = rt.load(str(path))
    g = graph_from_raytatouille(str(path))
    g.nodes["focus"].params["half_size"] = [1.0, 1.0]
    g.set_transport(g.links[-1], "H")
    image, ex = run(g, system, "focus", rings=30, pixels=[41, 41])
    assert any("Debye instead" in line for line in ex.context.log)
    assert any("by debye" in line for line in ex.context.log)
    assert image.values[20, 20] == image.values.max()
    # 0.78 mm off the axis, far beyond the ±0.24 mm of the Rayleigh–Sommerfeld grid
    assert image.values[20, 4] > 0.0
