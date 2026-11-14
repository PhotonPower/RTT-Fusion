"""Issue #69, side finding: a field reaches only as far as its outermost samples. Where an aperture
clips the sampled pupil, the start rays get points at the clipping edge (bisection between
neighbours that pass and are lost, `backend_rtt._edge_refined`), so that the field's power no
longer depends on how the hexapolar grid falls on the edge.

System: an entrance pupil of 2 mm (the backend starts the rays over it), an aperture of radius
0.7 mm behind it that clips the pupil (off-centre for a field angle), an ideal lens, the camera in
its focus; transport H. Reference: the clipped area times the source's irradiance – the field's
power is the same for 6 and 12 rings and equals the area of the aperture's opening as seen in
the pupil (an off-centre circle inside the pupil: π·0.7² of the pupil's area π·1²).
"""

from __future__ import annotations

import json
import warnings
from typing import Any

import pytest

from rtt_fusion.graph import backend_rtt, builtin_registry, graph_from_raytatouille
from rtt_fusion.graph.execute import WITH_FIELDS, Executor
from rtt_fusion.graph.field import field_image, field_power

rt = pytest.importorskip("raytatouille")
pytest.importorskip("fieldtracing")


def _system(field_deg: float, clip: float) -> dict[str, Any]:
    return {
        "schema_version": "0.9.0",
        "name": "clipped pupil",
        "units": {"length": "mm", "wavelength": "um"},
        "wavelengths": [{"um": 0.55, "reference": True}],
        "aperture": {"type": "epd", "value": 2.0},
        "fields": {"type": "angle_deg", "points": [{"y": field_deg}] if field_deg else [{}]},
        "root": {
            "type": "assembly",
            "name": "imaging",
            "children": [
                {
                    "type": "thin_element",
                    "name": "aperture",
                    "pose": {"position": [0.0, 0.0, 2.0]},
                    "surfaces": [{"id": "AP", "aperture": {"type": "circular", "radius": clip}}],
                },
                {
                    "type": "thin_element",
                    "name": "lens",
                    "pose": {"position": [0.0, 0.0, 3.0]},
                    "surfaces": [
                        {
                            "id": "L",
                            "aperture": {"type": "circular", "radius": 3.0},
                            "interaction": {"type": "ideal_lens", "focal_length": 10.0},
                        }
                    ],
                },
                {
                    "type": "detector",
                    "name": "camera",
                    "pose": {"position": [0.0, 0.0, 13.0]},
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
                    {"surface": "AP", "kind": "transmit"},
                    {"surface": "L", "kind": "transmit"},
                    {"surface": "CAM", "kind": "transmit"},
                ],
            }
        ],
    }


def _field_power(field_deg: float, clip: float, rings: int, monkeypatch: Any) -> float:
    seen: list[float] = []

    def spy(fields: Any, *args: Any, **kwargs: Any) -> Any:
        seen.append(sum(field_power(f) for f in fields))
        return field_image(fields, *args, **kwargs)

    monkeypatch.setattr("rtt_fusion.graph.backend_rtt.field_image", spy)
    d = _system(field_deg, clip)
    g = graph_from_raytatouille(d)
    camera = next(n for n in g.nodes.values() if n.type_id == "rtt.detector")
    camera.params["half_size"] = [0.05, 0.05]
    for lk in g.links_to(camera.id):
        g.set_transport(lk, "H")
    ex = Executor(
        g,
        builtin_registry(),
        dict(backend_rtt.EVALUATORS),
        {"system": rt.System.from_json(json.dumps(d)), "rings": rings, "pixels": [11, 11]},
        available=WITH_FIELDS,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ex.evaluate(camera.id, "image")
    return seen[-1]


@pytest.mark.parametrize("field_deg", [0.0, 1.0])
def test_issue_69_field_power_does_not_depend_on_the_ring_count(
    field_deg: float, monkeypatch: Any
) -> None:
    unclipped = _field_power(field_deg, 5.0, 8, monkeypatch)  # the aperture does not clip
    coarse = _field_power(field_deg, 0.7, 6, monkeypatch)
    fine = _field_power(field_deg, 0.7, 12, monkeypatch)
    # 1 %: the support is the polygon through the edge points, which lies inside the circle (about
    # 0.5 % of the area with the ~36 edge points of 6 rings); before the refinement both lay ~10 %
    # short of the true area
    assert coarse == pytest.approx(fine, rel=1e-2)
    assert fine / unclipped == pytest.approx(0.7**2, rel=1e-2)  # the clipped share of the pupil
