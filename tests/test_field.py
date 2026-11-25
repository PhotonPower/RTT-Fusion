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
from dataclasses import replace
from types import SimpleNamespace
from typing import Any, cast

import numpy as np
import pytest

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.graph import Graph, Node, builtin_registry, optics
from rtt_fusion.graph.execute import WITH_FIELDS, ExecutionError, Executor
from rtt_fusion.graph.field import field_to_rays, rays_to_field
from rtt_fusion.graph.imaging import Image, coherent_image
from rtt_fusion.graph.light import Plane, RayBundle

# fieldtracing is optional (requirements-modules.txt); without it these tests are skipped
richards_wolf = pytest.importorskip("fieldtracing.focus").richards_wolf

WAVELENGTH_UM = 0.6328
K = 2 * np.pi / (WAVELENGTH_UM * 1e-3)


def pixel_grid(half: float, n: int) -> tuple[FloatArray, FloatArray]:
    u = (np.arange(n) + 0.5) / n * 2 * half - half
    uu, vv = np.meshgrid(u, u)
    return np.asarray(uu, dtype=np.float64), np.asarray(vv, dtype=np.float64)


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


@pytest.mark.parametrize("given", [None, 0.05])
def test_field_tolerance_reaches_the_planner(given: float | None, monkeypatch: Any) -> None:
    """ε of the run (`shared["field_tolerance"]`, the Workflow panel) is the planner's tolerance;
    without it Fieldtracing's default holds."""
    import fieldtracing.chain as chain

    seen: list[float] = []
    real = chain.plan_propagation

    def spy(*args: Any, **kwargs: Any) -> Any:
        seen.append(kwargs.get("tolerance", 0.02))
        return real(*args, **kwargs)

    monkeypatch.setattr(chain, "plan_propagation", spy)
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(Node("src", "optics.collimated_source", params={"radius": 1.0, "rings": 10}))
    g.add(
        Node("lens", "optics.aplanatic_lens", params={"focal_length": 2.0, "aperture_radius": 1.0})
    )
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, 2.0],
                "normal": [0, 0, -1],
                "half_size": [2e-3, 2e-3],
                "pixels": [9, 9],
            },
        )
    )
    g.connect(("src", "out"), ("lens", "in"))
    g.connect(("lens", "out"), ("cam", "in"), "H")
    shared = {} if given is None else {"field_tolerance": given}
    Executor(g, reg, ev, shared, available=WITH_FIELDS).evaluate("cam", "image")
    assert seen and set(seen) == {0.02 if given is None else given}


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


# ---------------------------------------------------------------------------- W65


def _spot_image(half: float, pixels: int = 32, behind: float = 0.05) -> Image:
    """A ray image (R) of a focus: aplanatic lens f = 2 mm, NA 0.5, detector `behind` the focus
    (the spot is then about 2·0.5·behind wide)."""
    f, a = 2.0, 1.0
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": a, "rings": 12, "center": [0, 0, -1]},
        )
    )
    g.add(Node("lens", "optics.aplanatic_lens", params={"focal_length": f, "aperture_radius": a}))
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, f + behind],
                "normal": [0, 0, -1],
                "half_size": [half, half],
                "pixels": [pixels, pixels],
            },
        )
    )
    g.connect(("src", "out"), ("lens", "in"))
    g.connect(("lens", "out"), ("cam", "in"), "R")
    im: Image = Executor(g, reg, ev).evaluate("cam", "image")
    return im


def test_w65_a_spot_smaller_than_a_pixel_is_noted() -> None:
    big = _spot_image(5.0)  # 0.31 mm pixels, a spot of about 0.05 mm
    assert big.ray_power > 0 and big.power < 0.5 * big.ray_power
    assert "pixel centre" in big.note and "half size to about" in big.note
    suggested = float(big.note.split("half size to about ")[1].split(" mm")[0])
    fitted = _spot_image(suggested)  # the window the note suggests
    assert fitted.note == ""
    assert fitted.power == pytest.approx(fitted.ray_power, rel=0.1)


def test_w65_a_fitting_window_and_interference_are_no_reason() -> None:
    assert _spot_image(0.05).note == ""
    # two coherent copies in antiphase: the image is dark by interference, well sampled
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": 1.0, "rings": 12, "center": [0, 0, -1]},
        )
    )
    rays = Executor(g, reg, ev).evaluate("src", "out")[0]
    opposite = replace(rays, amplitude=-np.asarray(rays.amplitude))
    dark = coherent_image([rays, opposite], Plane.make([0, 0, 0], [0, 0, -1]), 1.5, 1.5, 32, 32)
    assert dark.ray_power > 0 and dark.power < 1e-9 * dark.ray_power
    assert dark.note == ""


def test_sampling_a_grid_equals_the_cubic_interpolation_of_the_whole_grid() -> None:
    """`_sample` fits the spline only where the pixels lie: as close to the sampled function as
    scipy's cubic interpolation over the whole grid (which took 26 s per field on a 2016² grid),
    0 outside the grid."""
    from types import SimpleNamespace

    from scipy.interpolate import RegularGridInterpolator

    from rtt_fusion.graph.field import _sample

    def field(u: FloatArray, v: FloatArray) -> ComplexArray:
        phase = np.exp(1j * (5.0 * u - 2.0 * v)) * np.exp(-(u**2 + 2 * v**2) / 2)
        return np.stack([phase, 0.5 * phase.conj(), 0.1j * phase], axis=-1)

    x = np.linspace(-3.0, 3.0, 301)
    xx, yy = np.meshgrid(x, x)
    values = field(xx, yy)
    f = SimpleNamespace(x=x, physical=lambda: values)
    rng = np.random.default_rng(3)
    q = np.concatenate([rng.uniform(-1.0, 0.5, (500, 2)), [[3.5, 0.0], [0.0, -3.2]]])
    got = _sample(cast(Any, f), q)
    want = np.zeros_like(got)
    for c in range(3):
        for part, unit in ((np.real, 1.0), (np.imag, 1j)):
            interp = RegularGridInterpolator(
                (x, x), part(values[:, :, c]), method="cubic", bounds_error=False, fill_value=0.0
            )
            want[:, c] += unit * interp(np.stack([q[:, 1], q[:, 0]], axis=-1))
    exact = field(q[:-2, 0], q[:-2, 1])
    error = np.abs(got[:-2] - exact).max()
    assert error < 1e-6
    assert error <= 1.5 * np.abs(want[:-2] - exact).max()  # the same order (other end conditions)
    assert np.abs(got - want).max() < 1e-6
    assert np.all(got[-2:] == 0)


def test_a_method_without_operators_falls_back(monkeypatch: Any) -> None:
    """RTT-Fusion#119: the planner may choose a method Fieldtracing cannot compute yet ("mixed",
    F_min ≤ 1 and F_max ≥ 1/ε, e.g. many rings of a focus). The propagation takes the method the
    next selection rule gives (here Debye, the focus of W13) and says so, instead of a ValueError
    out of Fieldtracing."""
    import fieldtracing.chain as chain

    real = chain.plan_propagation

    def mixed(*args: Any, **kwargs: Any) -> Any:
        plan = real(*args, **kwargs)
        if kwargs.get("method") is None:
            return replace(plan, method="mixed", reason="mixed (test)", operators=())
        return plan

    monkeypatch.setattr(chain, "plan_propagation", mixed)
    f, a, half, n = 2.0, 1.0, 2.5e-3, 21
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(Node("src", "optics.collimated_source", params={"radius": a, "rings": 20}))
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
    assert any("cannot compute yet – debye instead" in line for line in ex.context.log)
    plane = Plane.make([0, 0, f], [0, 0, -1])
    uu, vv = pixel_grid(half, n)
    pts = uu.reshape(-1, 1) * plane.u[:2] + vv.reshape(-1, 1) * plane.v[:2]
    ref = np.sum(np.abs(richards_wolf(pts, f, K, a / f)) ** 2, axis=1).reshape(n, n)
    assert np.abs(im.values - ref).max() / ref.max() < 1e-2


# ------------------------------------------------- #88: off-axis foci (Debye on the field's axis)


def _converging(x_c: float, x_f: float, length: float, radius: float, rings: int = 16) -> Any:
    """An ideally converging bundle: rays from a disc of `radius` around (x_c, 0, 0) to the focus
    (x_f, 0, length)."""
    from rtt_fusion.graph.light import CoherenceKey

    px, py = [0.0], [0.0]
    for r in range(1, rings + 1):
        for k in range(6 * r):
            a = 2 * np.pi * k / (6 * r)
            px.append(r / rings * np.cos(a))
            py.append(r / rings * np.sin(a))
    start = np.stack(
        [x_c + radius * np.array(px), radius * np.array(py), np.zeros(len(px))], axis=-1
    )
    d = np.array([x_f, 0.0, length]) - start
    dist = np.linalg.norm(d, axis=1)
    return RayBundle(
        key=CoherenceKey("source"),
        wavelength_um=0.54607,
        position=start,
        direction=d / dist[:, None],
        opl=dist.max() - dist,
        amplitude=np.ones(len(px), dtype=complex),
        source_xy=start[:, :2].copy(),
        alive=np.ones(len(px), dtype=bool),
        source_frame=np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    )


@pytest.mark.parametrize(
    ("x_c", "x_f", "origin"),
    [(0.0, 1.5, 0.0), (0.0, 1.8, 0.0), (-1.2, 1.2, -1.2), (1.2, -1.2, 1.2)],
)
def test_off_axis_focus_lands_at_the_geometric_focus(x_c: float, x_f: float, origin: float) -> None:
    """#88 finding 1: a bundle converging beside the axis of its field coordinates (the port
    plane's origin lies in its centroid) focused onto a periodic copy of the Debye integral. The
    peak must lie at the geometric focus – also beyond half a Debye period (2.5 mm here)."""
    from rtt_fusion.graph.field import field_image

    length, half, n = 58.0, 2.0, 201
    fb = rays_to_field(
        _converging(x_c, x_f, length, 1.2), Plane.make([origin, 0.0, 0.0], [0.0, 0.0, 1.0])
    )
    detector = Plane.make([0.0, 0.0, length], [0.0, 0.0, 1.0])
    log: list[str] = []
    v = field_image([fb], detector, half, half, n, n, log).values
    row, col = np.unravel_index(np.argmax(v), v.shape)
    pixel = 2 * half / n
    assert abs((col + 0.5) * pixel - half - x_f) <= pixel, log
    assert abs((row + 0.5) * pixel - half) <= pixel


def test_two_coherent_beams_interfere_in_their_common_focus() -> None:
    """#88 case A: two parallel coherent beams (±1 orders of a grating pair, 0.98 mm off the axis)
    meet in the focus of a lens: a spot in the middle with fringes of period λ·f / (2·offset)
    (34.0 µm), contrast ≈ 1. Each bundle's field has its origin in its own centroid, as a port
    plane puts it; before the fix the image was empty in the middle."""
    from rtt_fusion.graph.field import field_image

    length, offset, lam_mm = 164.5, 0.98, 0.405e-3
    period = lam_mm * length / (2 * offset)
    fields = []
    for x_c in (offset, -offset):
        b = _converging(x_c, 0.0, length, 0.15, rings=8)
        b = replace(b, wavelength_um=0.405)
        fields.append(rays_to_field(b, Plane.make([x_c, 0.0, 0.0], [0.0, 0.0, 1.0])))
    detector = Plane.make([0.0, 0.0, length], [0.0, 0.0, 1.0])
    half, n = 0.06, 241  # 0.5 µm pixels
    log: list[str] = []
    v = field_image(fields, detector, half, half, n, n, log).values
    x = (np.arange(n) + 0.5) / n * 2 * half - half
    row = v[n // 2]
    centre = float(row[np.argmin(np.abs(x))])
    assert centre >= 0.97 * float(v.max()), log
    dark = float(row[np.argmin(np.abs(x - period / 2))])
    bright = float(row[np.argmin(np.abs(x - period))])
    assert dark < 0.03 * centre
    assert bright > 0.5 * centre


def test_a_steep_sub_beam_in_its_focus_is_not_dark() -> None:
    """#88 case B: a sub-beam converging at 37° (sin α = 0.6) onto a detector at its focus: the
    field image (Debye) must show the spot where the rays meet and hold the beam's power – before
    the fix it was empty (power 0, no note)."""
    from rtt_fusion.graph.field import field_image

    length = 1.2 / np.tan(np.arcsin(0.6))
    b = _converging(1.2, 0.0, length, 0.15, rings=12)
    b = replace(b, wavelength_um=0.405)
    fb = rays_to_field(b, Plane.make([1.2, 0.0, 0.0], [0.0, 0.0, 1.0]))
    detector = Plane.make([0.0, 0.0, length], [0.0, 0.0, 1.0])
    half, n = 0.006, 97
    log: list[str] = []
    h = field_image([fb], detector, half, half, n, n, log)
    assert h.values.max() > 0, log
    row, col = np.unravel_index(np.argmax(h.values), h.values.shape)
    assert abs(int(row) - n // 2) <= 1 and abs(int(col) - n // 2) <= 1
    # the spot (FWHM ≈ 2.7 µm) lies in the ±6 µm window: nearly all of the power arrives
    pixel_area = (2 * half / n) ** 2
    assert 0.9 * h.power <= float(h.values.sum()) * pixel_area <= 1.05 * h.power


# ------------------------------------------- #135: scattered rays, Debye period, auto at a window


def _two_beams_with_scattered_rays(offset: float = 1.0834, spread: float = 8.0) -> list[Any]:
    """The ±1 orders of the FrObSy head with fine fringes, as two coherent beams converging into
    a common focus 164.5 mm on: every 7th ray strongly aberrated (spread across `spread` times the
    pupil at the port plane, aimed 8 mm beside the focus, a tenth of the amplitude)."""
    length = 164.5
    rng = np.random.default_rng(1)
    fields = []
    for x_c in (offset, -offset):
        b = replace(_converging(x_c, 0.0, length, 0.3, rings=10), wavelength_um=0.405)
        idx = rng.choice(len(b.direction), size=len(b.direction) // 7, replace=False)
        p, d, a = np.array(b.position), np.array(b.direction), np.array(b.amplitude)
        centre = p.mean(axis=0)
        p[idx] = centre + (p[idx] - centre) * spread
        aim = np.array([8.0 * np.sign(x_c), 0.0, length]) - p[idx]
        d[idx] = aim / np.linalg.norm(aim, axis=1)[:, None]
        a[idx] *= 0.1
        b = replace(b, position=p, direction=d, amplitude=a)
        fields.append(rays_to_field(b, Plane.make([x_c, 0.0, 0.0], [0.0, 0.0, 1.0])))
    return fields


def test_scattered_rays_do_not_take_the_focus_off_the_window() -> None:
    """#135 case 2: a few strongly aberrated rays set the field's extent (pupil 5 mm instead of
    0.3 mm) and the grid of the propagation, which put the focus beside the camera's window –
    an empty image. Samples landing more than two window widths off are left out (logged): the
    two beams meet in the middle with fringes of λ·L/(2·offset)."""
    from rtt_fusion.graph.field import field_image

    length, offset, half, n = 164.5, 1.0834, 0.1523, 201
    period = 0.405e-3 * length / (2 * offset)
    log: list[str] = []
    detector = Plane.make([0.0, 0.0, length], [0.0, 0.0, 1.0])
    v = field_image(_two_beams_with_scattered_rays(offset), detector, half, half, n, n, log).values
    x = (np.arange(n) + 0.5) / n * 2 * half - half
    row = v[n // 2]
    centre = float(row[np.argmin(np.abs(x))])
    assert v.max() > 0 and centre >= 0.9 * float(v.max()), log
    assert float(row[np.argmin(np.abs(x - period / 2))]) < 0.05 * centre
    assert any("samples left out" in line for line in log), log


def test_debye_about_its_axis_needs_a_period_holding_the_window(monkeypatch: Any) -> None:
    """#135 case 2: Debye about the field's own axis gives one period n·λ/(2 sin θ) around it
    (zero beyond); a spot outside the inner half of it shows shifted by a period. Then the
    propagation goes about the focus, as for a spot far off the axis."""
    from rtt_fusion.graph import field

    fb = _two_beams_with_scattered_rays(spread=1.0)[0]  # converging 1.08 mm beside its axis
    detector = Plane.make([0.0, 0.0, 164.5], [0.0, 0.0, 1.0])
    used: list[Any] = []
    calls: list[float] = []

    def plan(fb_: Any, distance: float, window: Any, tolerance: Any, notes: list[str]) -> Any:
        calls.append(float(window))  # first about the focus, then about the field's own axis
        own = len(calls) == 2
        return SimpleNamespace(
            method="debye",
            reason="test",
            # about its own axis a period shorter than the window it needs (2.47 mm)
            parameters={"length_u_mm": 5.1, "samples": 238.0, "period_mm": 1.55 if own else 1.0},
            centred=not own,
        )

    monkeypatch.setattr(field, "plan_onto", plan)

    def run(fb_: Any, d: float, p: Any, notes: list[str]) -> Any:
        used.append(p)
        return fb_, notes

    monkeypatch.setattr(field, "_run", run)
    _, notes = field.propagate_to(fb, detector, window=0.3046)
    assert calls[1] > 2.0  # the window about its own axis holds the spot 1.08 mm off
    assert used[-1].centred, notes
    assert any("Debye period 1.55 mm" in line for line in notes), notes


def test_auto_judges_only_light_near_the_window() -> None:
    """#135 case 3: a bundle that focuses millimetres beside the detector's window does not
    decide about it; one focusing in the window asks for fields."""
    from rtt_fusion.graph.field import near_focus

    length = 164.5
    detector = Plane.make([0.0, 0.0, length], [0.0, 0.0, 1.0])
    beside = replace(_converging(0.0, 5.0, length, 0.3, rings=8), wavelength_um=0.405)
    inside = replace(_converging(1.0, 0.0, length, 0.3, rings=8), wavelength_um=0.405)
    assert near_focus(rays_to_field(beside), detector, 0.3046, 0.2) is None
    assert near_focus(rays_to_field(inside), detector, 0.3046, 0.2) is not None


# ------------------------------------------------- #135: a folded ray map, one field per sheet


def _aberrated(c: float = 0.5) -> Any:
    """A bundle converging with strong spherical aberration: the ray from pupil radius ρ meets
    the axis at f·(1 − c·ρ²/a²) (f = 50 mm, a = 2 mm); beyond the marginal rays' caustic the map
    from the pupil onto a plane folds."""
    f0, a = 50.0, 2.0
    b = replace(_converging(0.0, 0.0, f0, a, rings=40), wavelength_um=0.5)
    p = np.asarray(b.position)
    rho2 = np.sum(p[:, :2] ** 2, axis=1)
    focus = f0 * (1.0 - c * rho2 / a**2)
    d = np.stack([-p[:, 0], -p[:, 1], focus], axis=-1)
    d /= np.linalg.norm(d, axis=1)[:, None]
    length = np.sqrt(rho2 + focus**2)
    return replace(b, direction=d, opl=length.max() - length)


def test_an_unfolded_ray_map_gives_one_field() -> None:
    from rtt_fusion.graph.field import rays_to_fields

    b = _converging(0.3, 0.0, 58.0, 1.2)
    plane = Plane.make([0.0, 0.0, 10.0], [0.0, 0.0, 1.0])
    one = rays_to_field(b, plane)
    many = rays_to_fields(b, plane)
    assert len(many) == 1
    np.testing.assert_array_equal(many[0].field.points, one.field.points)
    np.testing.assert_array_equal(many[0].field.amplitude, one.field.amplitude)


def test_a_folded_ray_map_gives_a_field_per_sheet() -> None:
    """#135 case 2: behind an aberrated lens the rays of the rim cross the others; one field of
    all samples counted the overlap once (half the power of a beam behind the FrObSy objective).
    Here: the image 45 mm on, made from a plane where the map folds (32 mm, three sheets) and,
    as the reference route, from the start plane where it does not. One field per sheet comes
    close to the reference (samples at the folds are left out, the geometric amplitude does not
    hold there); one field of all samples has a fifth of the power."""
    from rtt_fusion.graph.field import field_image, rays_to_fields

    b = _aberrated()
    detector = Plane.make([0.0, 0.0, 45.0], [0.0, 0.0, 1.0])
    half, n = 0.4, 81
    folded = Plane.make([0.0, 0.0, 32.0], [0.0, 0.0, 1.0])
    sheets = rays_to_fields(b, folded)
    assert len(sheets) >= 2

    def image(fields: list[Any]) -> FloatArray:
        return np.asarray(field_image(fields, detector, half, half, n, n, []).values)

    reference = image([rays_to_field(b, Plane.make([0.0, 0.0, 1e-3], [0.0, 0.0, 1.0]))])
    by_sheet = image(sheets)
    one = image([rays_to_field(b, folded)])

    def agreement(v: FloatArray) -> tuple[float, float]:
        return float(np.corrcoef(reference.ravel(), v.ravel())[0, 1]), float(
            v.sum() / reference.sum()
        )

    corr_sheets, power_sheets = agreement(by_sheet)
    corr_one, power_one = agreement(one)
    assert corr_sheets > 0.75 and power_sheets > 0.45, (corr_sheets, power_sheets)
    assert corr_one < 0.65 and power_one < 0.3, (corr_one, power_one)
