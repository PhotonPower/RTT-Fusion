"""Start rays of a path that make_rays cannot aim (#125): the fallback puts them on the first
surface and sends them along the field direction in that surface's frame – a stop turned out of
the axis (a laser exit along +x before a fold) is hit; a path no start ray passes says so in the
log instead of an image of zero power without a word."""

from __future__ import annotations

import json
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import Graph
from rtt_fusion.graph.execute import Context

rt = pytest.importorskip("raytatouille")

from rtt_fusion.graph import backend_rtt  # noqa: E402


def _system(rotation_y: float) -> Any:
    d = {
        "schema_version": "0.9.0",
        "name": "turned stop",
        "units": {"length": "mm", "wavelength": "um"},
        "wavelengths": [{"um": 0.55, "reference": True}],
        "aperture": {"type": "epd", "value": 1.0},
        "fields": {"points": [{}]},
        "root": {
            "type": "assembly",
            "name": "root",
            "children": [
                {
                    "type": "stop",
                    "name": "exit",
                    "pose": {"rotation_deg": [0.0, rotation_y, 0.0]},
                    "surfaces": [{"id": "S", "aperture": {"type": "circular", "radius": 0.5}}],
                },
                {
                    "type": "detector",
                    "name": "cam",
                    "pose": {"position": [10.0, 0.0, 0.0], "rotation_deg": [0.0, 90.0, 0.0]},
                    "surfaces": [
                        {
                            "id": "D",
                            "aperture": {
                                "type": "rectangular",
                                "half_width_x": 2.0,
                                "half_width_y": 2.0,
                            },
                        }
                    ],
                },
            ],
        },
        "paths": [{"name": "main", "events": [{"surface": "S"}, {"surface": "D"}]}],
    }
    return rt.System.from_json(json.dumps(d))


@pytest.fixture
def no_aiming(monkeypatch: pytest.MonkeyPatch) -> None:
    def refuse(*args: Any, **kwargs: Any) -> Any:
        raise rt.RaytatouilleError("not aimable")

    monkeypatch.setattr(rt.trace, "make_rays", refuse)


def _trace(rotation_y: float) -> tuple[Any, Any]:
    system = _system(rotation_y)
    ctx = Context(Graph(), shared={"system": system})
    compiled = rt.compile(system, rt.MaterialLibrary())
    rays = backend_rtt._start_rays(ctx, compiled, "main", 3, 0, 0)
    _, rec = rt.trace.trace(
        compiled, rays, path="main", record_path=True, max_recorded_rays=rays.size
    )
    return rays, rec


def test_start_rays_run_along_the_normal_of_a_turned_stop(no_aiming: None) -> None:
    rays, rec = _trace(90.0)
    d = np.stack([np.asarray(rays.dir_x), np.asarray(rays.dir_y), np.asarray(rays.dir_z)], -1)
    assert np.allclose(d, [1.0, 0.0, 0.0], atol=1e-12)
    assert np.all(backend_rtt._passed(rec))


def test_an_axial_stop_is_unchanged(no_aiming: None) -> None:
    rays, _ = _trace(0.0)
    d = np.stack([np.asarray(rays.dir_x), np.asarray(rays.dir_y), np.asarray(rays.dir_z)], -1)
    assert np.allclose(d, [0.0, 0.0, 1.0], atol=1e-12)


def test_a_dark_path_is_logged_once() -> None:
    from rtt_fusion.graph import Node

    ctx = Context(Graph())
    node = Node("cam", "rtt.detector")
    for _ in range(2):
        backend_rtt._dark(ctx, node, "main", 0, 0)
    assert len([x for x in ctx.log if "no start ray passes path 'main'" in x]) == 1
