"""Issue #69: fields (H) of off-axis field points focus where the rays do.

Raytatouille starts the rays of a field point at infinity with OPL 0 on its start plane, which is
the wavefront only on the axis; `backend_rtt._start_opl` adds the plane wave's own phase n·d̂·(p −
p_ref). An ideal lens (f = 10 mm, entrance pupil 2 mm at the lens, NA 0.1; no stop, so that the
backend starts the rays itself, as for the folded Linnik of #69) images field angles of 0°, 0.5°
and 1° onto a camera in its focus through transport H: the spot sits at f·tan θ (independent
reference; the ray image cannot be one: in a perfect focus all rays meet in one point) with the
on-axis peak, and ψ = k₀·OPL agrees with the rays' directions already behind the source."""

from __future__ import annotations

import json
import warnings
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import backend_rtt, builtin_registry, graph_from_raytatouille
from rtt_fusion.graph.execute import WITH_FIELDS, Executor

rt = pytest.importorskip("raytatouille")
pytest.importorskip("fieldtracing")

F_MM, LAM = 10.0, 0.55


def _system(field_deg: float) -> dict[str, Any]:
    return {
        "schema_version": "0.9.0",
        "name": "ideal lens, field angle",
        "units": {"length": "mm", "wavelength": "um"},
        "wavelengths": [{"um": LAM, "reference": True}],
        # no stop, an entrance pupil instead: Raytatouille cannot aim the rays, the backend starts
        # them over the EPD on a plane ⟂ the axis – the start rays of #69 (a folded Linnik)
        "aperture": {"type": "epd", "value": 2.0},
        "fields": {"type": "angle_deg", "points": [{"y": field_deg}] if field_deg else [{}]},
        "root": {
            "type": "assembly",
            "name": "imaging",
            "children": [
                {
                    "type": "thin_element",
                    "name": "lens",
                    "pose": {"position": [0.0, 0.0, 0.0]},
                    "surfaces": [
                        {
                            "id": "L",
                            "aperture": {"type": "circular", "radius": 2.0},
                            "interaction": {"type": "ideal_lens", "focal_length": F_MM},
                        }
                    ],
                },
                {
                    "type": "detector",
                    "name": "camera",
                    "pose": {"position": [0.0, 0.0, F_MM]},
                    "surfaces": [
                        {
                            "id": "CAM",
                            "aperture": {
                                "type": "rectangular",
                                "half_width_x": 1.0,
                                "half_width_y": 1.0,
                            },
                        }
                    ],
                },
            ],
        },
        "paths": [
            {
                "name": "image",
                "events": [
                    {"surface": "L", "kind": "transmit"},
                    {"surface": "CAM", "kind": "transmit"},
                ],
            }
        ],
    }


def _image(field_deg: float, transport: str = "H") -> Any:
    d = _system(field_deg)
    g = graph_from_raytatouille(d)
    camera = next(n for n in g.nodes.values() if n.type_id == "rtt.detector")
    camera.params["half_size"] = [0.25, 0.25]
    for lk in g.links_to(camera.id):
        g.set_transport(lk, transport)
    ex = Executor(
        g,
        builtin_registry(),
        dict(backend_rtt.EVALUATORS),
        {"system": rt.System.from_json(json.dumps(d)), "rings": 10, "pixels": [201, 201]},
        available=WITH_FIELDS,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ex.evaluate(camera.id, "image"), ex


def _peak(image: Any) -> tuple[float, float, float]:
    v = np.asarray(image.values)
    x0, x1, y0, y1 = image.extent
    j, i = np.unravel_index(v.argmax(), v.shape)
    ny, nx = v.shape
    return x0 + (i + 0.5) * (x1 - x0) / nx, y0 + (j + 0.5) * (y1 - y0) / ny, float(v.max())


def _encircled(image: Any, centre: tuple[float, float], radius: float) -> float:
    v = np.asarray(image.values)
    x0, x1, y0, y1 = image.extent
    ny, nx = v.shape
    xs = x0 + (np.arange(nx) + 0.5) * (x1 - x0) / nx
    ys = y0 + (np.arange(ny) + 0.5) * (y1 - y0) / ny
    xx, yy = np.meshgrid(xs, ys)
    inside = np.hypot(xx - centre[0], yy - centre[1]) <= radius
    return float(v[inside].sum() / v.sum())


def test_issue_69_off_axis_field_points_focus_at_f_tan_theta() -> None:
    """The light of 0.5° and 1° gathers at f·tan θ as on the axis: the share inside two Airy radii
    (6.7 µm) about the image point is the on-axis one (before the fix the spot lay 22 µm away and
    the share was a few per cent); the power is that of the axis (nothing vignettes).

    The peak is checked at field angles whose image point is a pixel centre (the Airy disc, radius
    3.4 µm, is sampled by 2.5 µm pixels; half a pixel off lowers the peak a lot): within 3 % of the
    axis. The start term needs the index of object space – an index estimated from the first
    segment (here the ideal lens, which adds its phase there) gave 0.987 and left a tilt that moved
    the field's focus by about 1 µm."""
    airy = 0.61 * LAM * 1e-3 / 0.1
    image0, _ = _image(0.0)
    share0 = _encircled(image0, (0.0, 0.0), 2 * airy)
    assert share0 > 0.8  # content guard: an Airy disc
    for field_deg in (0.5, 1.0):
        image, _ = _image(field_deg)
        expected = F_MM * np.tan(np.radians(field_deg))
        x, y, _ = _peak(image)
        centre = (0.0, float(np.sign(y)) * expected)
        assert _encircled(image, centre, 2 * airy) == pytest.approx(share0, abs=0.03), field_deg
        assert image.power == pytest.approx(image0.power, rel=1e-3)
    _, _, peak0 = _peak(image0)
    pixel = 0.5 / 201
    for near in (0.5, 1.0):
        on_pixel = round(F_MM * np.tan(np.radians(near)) / pixel) * pixel  # 201 pixels: one at 0
        image, _ = _image(float(np.degrees(np.arctan(on_pixel / F_MM))))
        _, y, peak = _peak(image)
        assert abs(abs(y) - on_pixel) < 1e-9 and peak == pytest.approx(peak0, rel=0.03), near


def test_issue_69_start_opl_agrees_with_the_directions() -> None:
    """Behind the source of a tilted field point ψ = k₀·OPL changes between neighbouring rays as
    their common direction says (it was off by about 90 rad across the pupil in #69)."""
    from scipy.spatial import Delaunay

    from rtt_fusion.graph.field import rays_to_field

    d = _system(1.0)
    g = graph_from_raytatouille(d)
    ex = Executor(
        g,
        builtin_registry(),
        dict(backend_rtt.EVALUATORS),
        {"system": rt.System.from_json(json.dumps(d)), "rings": 10, "rtt_execution": "nodes"},
        available=WITH_FIELDS,
    )
    source = next(n for n in g.nodes.values() if n.type_id == "rtt.source")
    link = next(lk for lk in g.links if lk.source[0] == source.id)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bundle = ex.evaluate(source.id, link.source[1])[0]
    from fieldtracing.field import GridlessField

    f = rays_to_field(bundle).field
    assert isinstance(f, GridlessField)
    p = np.asarray(f.points)
    tri = Delaunay(p).simplices
    edges = np.unique(
        np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [0, 2]]]), axis=1), axis=0
    )
    a, b = edges[:, 0], edges[:, 1]
    grad, phase = np.asarray(f.gradient), np.asarray(f.phase)
    mismatch = np.abs(
        (phase[b] - phase[a]) - 0.5 * np.sum((grad[a] + grad[b]) * (p[b] - p[a]), axis=1)
    )
    assert float(mismatch.max()) < 1e-6, float(mismatch.max())
