"""Representation H on links (R18): fields of Fieldtracing between blocks.

References, independent of the code under test:
- W13 focus of an aplanatic lens: Richards and Wolf (1959), `fieldtracing.focus.richards_wolf`
  (Bessel-integral form, not the FFT/Debye path used here).
- W14 Gaussian beam: I(r, z) = I₀ (w₀/w)² exp(−2r²/w²), w = w₀ sqrt(1 + (z/z_R)²),
  z_R = π w₀²/λ; power π w₀² I₀/2.
- W15 Michelson with H: two-beam interference as in W5 (half the intensity at δ = λ/8) and the
  on-axis decay of the Gaussian beams over the H segment, 1/(1 + (L/z_R)²).
- W16 R -> H -> R: the ray image must not change (the conversion alone propagates nothing).
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

pytest.importorskip("fieldtracing")

from fieldtracing.focus import richards_wolf  # noqa: E402

from rtt_fusion.graph import Graph, Node, builtin_registry, optics
from rtt_fusion.graph.execute import WITH_FIELDS, ExecutionError, Executor
from rtt_fusion.graph.field import field_to_rays, rays_to_field
from rtt_fusion.graph.imaging import Image, coherent_image
from rtt_fusion.graph.light import Plane

WAVELENGTH_UM = 0.6328
K = 2 * np.pi / (WAVELENGTH_UM * 1e-3)


def pixel_grid(half: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    u = (np.arange(n) + 0.5) / n * 2 * half - half
    return np.meshgrid(u, u)


# ------------------------------------------------------------------------------ W13 focus


def test_w13_focus_against_richards_wolf() -> None:
    f, a, half, n = 2.0, 1.0, 2.5e-3, 41  # NA 0.5, ±3 Airy radii
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": a, "rings": 30, "center": [0, 0, -1]},
        )
    )
    g.add(Node("lens", "optics.aplanatic_lens", params={"focal_length": f, "aperture_radius": a}))
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, f],
                "normal": [0, 0, -1],
                "half_size": [half, half],
                "pixels": [n, n],
            },
        )
    )
    g.connect(("src", "out"), ("lens", "in"))
    g.connect(("lens", "out"), ("cam", "in"), "H")
    ex = Executor(g, reg, ev, available=WITH_FIELDS)
    im: Image = ex.evaluate("cam", "image")
    assert ex.plan is not None and ex.plan.representations["lens.out -> cam.in"] == "H"
    assert any("debye" in line for line in ex.context.log)  # Fieldtracing chose Debye
    plane = Plane.make([0, 0, f], [0, 0, -1])
    uu, vv = pixel_grid(half, n)
    pts = uu.reshape(-1, 1) * plane.u[:2] + vv.reshape(-1, 1) * plane.v[:2]
    ref = np.sum(np.abs(richards_wolf(pts, f, K, a / f)) ** 2, axis=1).reshape(n, n)
    assert im.values.max() == pytest.approx(ref.max(), rel=2e-3)
    assert np.abs(im.values - ref).max() / ref.max() < 5e-3


# ------------------------------------------------------------------------------ W14 Gaussian beam


@pytest.mark.parametrize("z_over_zr", [0.5, 1.0, 2.0])
def test_w14_gaussian_beam(z_over_zr: float) -> None:
    w0 = 0.05
    zr = np.pi * w0**2 / (WAVELENGTH_UM * 1e-3)
    z = z_over_zr * zr
    w = w0 * np.sqrt(1 + z_over_zr**2)
    half, n = 2.0 * w, 41
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": 2.5 * w0, "rings": 25, "gaussian_waist": w0},
        )
    )
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, z],
                "normal": [0, 0, -1],
                "half_size": [half, half],
                "pixels": [n, n],
            },
        )
    )
    g.connect(("src", "out"), ("cam", "in"), "H")
    im: Image = Executor(g, reg, ev, available=WITH_FIELDS).evaluate("cam", "image")
    uu, vv = pixel_grid(half, n)
    ref = (w0 / w) ** 2 * np.exp(-2 * (uu**2 + vv**2) / w**2)
    assert np.abs(im.values - ref).max() / ref.max() < 2e-2
    assert im.ray_power == pytest.approx(np.pi * w0**2 / 2, rel=1e-2)  # power at the source


# ------------------------------------------------------------------------------ W15 Michelson


def michelson_h(transports: dict[int, str], delta_mm: float = 0.0, theta: float = 0.0) -> Image:
    """The Michelson of W4 made small: Gaussian beam w₀ = 0.1 mm, arms 3 mm, camera 6 mm behind
    the splitter."""
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"center": [0, 0, -10], "radius": 0.25, "rings": 20, "gaussian_waist": 0.1},
        )
    )
    g.add(Node("bs", "optics.beam_splitter", params={"point": [0, 0, 0], "normal": [-1, 0, 1]}))
    g.add(Node("mA", "optics.plane_mirror", params={"point": [0, 0, 3.0], "normal": [0, 0, -1]}))
    g.add(
        Node(
            "mB",
            "optics.plane_mirror",
            params={"point": [3.0 + delta_mm, 0, 0], "normal": [-np.cos(theta), np.sin(theta), 0]},
        )
    )
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [-6, 0, 0],
                "normal": [1, 0, 0],
                "half_size": [0.2, 0.2],
                "pixels": [48, 48],
            },
        )
    )
    g.connect(("src", "out"), ("bs", "in"))
    g.connect(("bs", "t"), ("mA", "in"))
    g.connect(("mA", "out"), ("bs", "t_in"))
    g.connect(("bs", "r"), ("mB", "in"))
    g.connect(("mB", "out"), ("bs", "r_in"))
    g.connect(("bs", "merged"), ("cam", "in"))
    for i, t in transports.items():
        g.set_transport(g.links[i], t)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result: Image = Executor(g, reg, ev, available=WITH_FIELDS).evaluate("cam", "image")
    return result


def test_w15_michelson_with_fields() -> None:
    zr = np.pi * 0.1**2 / (WAVELENGTH_UM * 1e-3)
    decay = 1.0 / (1.0 + (6.0 / zr) ** 2)  # H segment: splitter (port plane) -> camera, 6 mm
    c = 24
    rays = michelson_h({})
    field = michelson_h({5: "H"})
    assert field.values[c, c] / rays.values[c, c] == pytest.approx(decay, abs=3e-3)
    # two-beam interference: δ = λ/8 halves the intensity (W5)
    eighth = michelson_h({5: "H"}, delta_mm=WAVELENGTH_UM * 1e-3 / 8)
    assert eighth.values[c, c] / field.values[c, c] == pytest.approx(0.5, abs=2e-3)
    # H on every link: the arms go R -> H -> R at the mirrors, the image does not change
    every = michelson_h(dict.fromkeys(range(6), "H"))
    assert np.abs(every.values - field.values).max() < 1e-6 * field.values.max()


def test_w15_tilted_arm() -> None:
    # arm B tilted by 2 mrad: its field plane is tilted against the camera and is resampled
    rays = michelson_h({}, theta=2e-3)
    field = michelson_h({5: "H"}, theta=2e-3)
    assert np.abs(field.values - rays.values).max() / rays.values.max() < 3e-2


# ------------------------------------------------------------------------------ W16 round trip


def test_w16_rays_field_rays() -> None:
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": 0.25,
                "rings": 20,
                "gaussian_waist": 0.1,
                "field_angles_deg": [[0.05, 0.0]],
            },
        )
    )
    bundle = Executor(g, reg, ev).evaluate("src", "out")[0]
    plane = Plane.make([0, 0, 5.0], [0, 0, -1])
    direct = coherent_image([bundle], plane, 0.2, 0.2, 48, 48)
    back = field_to_rays(rays_to_field(bundle))
    via = coherent_image([back], plane, 0.2, 0.2, 48, 48)
    assert np.abs(via.values - direct.values).max() < 2e-3 * direct.values.max()


def test_field_errors() -> None:
    # H needs a backend that runs fields; the default executor runs R only
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(Node("src", "optics.collimated_source"))
    g.add(Node("cam", "optics.detector", params={"point": [0, 0, 5], "normal": [0, 0, -1]}))
    g.connect(("src", "out"), ("cam", "in"), "H")
    with pytest.raises(ExecutionError, match="exec.transport_unsupported"):
        Executor(g, reg, ev).evaluate("cam", "image")
