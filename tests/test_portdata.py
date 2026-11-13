"""Port data (CAD ADR 0010, D4, W2a): the values on the links through files, exactly – floats bit
for bit (NaN and infinities too), arrays as .npy; only allowed classes are created."""

from __future__ import annotations

import dataclasses
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion import portdata
from rtt_fusion.graph import Graph, Node, builtin_registry, optics
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.light import CoherenceKey, RayBundle
from rtt_fusion.ports import GriddedField, PlaneWaves, Port


def same(a: Any, b: Any) -> bool:
    """Equal bit for bit (NaN equals NaN, arrays by dtype, shape and content)."""
    if type(a) is not type(b):
        return False
    if isinstance(a, np.ndarray):
        return bool(a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes())
    if isinstance(a, float):
        return a.hex() == b.hex() if not math.isnan(a) else math.isnan(b)
    if isinstance(a, complex):
        return same(a.real, b.real) and same(a.imag, b.imag)
    if dataclasses.is_dataclass(a):
        return all(same(getattr(a, f.name), getattr(b, f.name)) for f in dataclasses.fields(a))
    if isinstance(a, list | tuple):
        return len(a) == len(b) and all(same(x, y) for x, y in zip(a, b, strict=True))
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)
    return bool(a == b)


def round_trip(value: Any, tmp: Path) -> Any:
    handle = portdata.write(value, tmp)
    text = json.dumps(handle, allow_nan=False)  # plain JSON: CAD's parser knows no NaN
    return portdata.read(json.loads(text), tmp)


def _rays() -> RayBundle:
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"radius": 1.0, "rings": 4, "center": [0, 0, -1]},
        )
    )
    bundle: RayBundle = Executor(g, reg, ev).evaluate("src", "out")[0]
    return bundle


def test_w68_ray_bundles_exactly(tmp_path: Path) -> None:
    rays = _rays()
    assert same(round_trip([rays], tmp_path / "a"), [rays])  # an evaluator's output: a list
    odd = dataclasses.replace(
        rays,
        opl=np.where(np.arange(len(rays.opl)) % 2 == 0, np.nan, np.inf),
        irradiance=float("-inf"),
        index=float("nan"),
        branch=("main", 2),
        jones=np.array([1.0 + 2.0j, -0.0, 3e-310j]),
    )
    back = round_trip(odd, tmp_path / "b")
    assert same(back, odd)
    assert isinstance(back.key, CoherenceKey)


def test_w68_ports_images_and_containers(tmp_path: Path) -> None:
    port = Port(0.002, -1, 1.5 + 1e-3j, 0.6328)
    waves = PlaneWaves(
        port,
        np.array([0.0, 1.5]),
        np.array([2.0, -1.0]),
        np.array([[1 + 1j, 0], [0.5, -2j]]),
        np.array([[0, 0], [1, -1]], dtype=np.int64),
    )
    grid = GriddedField(port, 1e-4, 2e-4, np.ones((3, 4, 2), dtype=complex) * (1 - 1j))
    image = Image(
        np.arange(12.0).reshape(3, 4),
        (-1.0, 1.0, -0.5, 0.5),
        2.5,
        [CoherenceKey("s", 0, 1, 0)],
        "a note",
    )
    value = {"waves": waves, "grid": grid, "image": image, "n": 7, "flag": True, "none": None}
    assert same(round_trip(value, tmp_path), value)


def test_w68_refusals(tmp_path: Path) -> None:
    handle = portdata.write(_rays(), tmp_path)
    bad = json.loads(json.dumps(handle))
    bad["root"]["cls"] = "subprocess.Popen"
    with pytest.raises(portdata.PortDataError, match="outside"):
        portdata.read(bad, tmp_path)
    bad = json.loads(json.dumps(handle))
    bad["root"]["fields"]["position"] = {"t": "array", "file": "../elsewhere.npy"}
    with pytest.raises(portdata.PortDataError, match="outside the folder"):
        portdata.read(bad, tmp_path)
    with pytest.raises(portdata.PortDataError, match="not rtt-fusion-portdata"):
        portdata.read({"format": "other", "version": 1, "root": {"t": "none"}}, tmp_path)
    with pytest.raises(portdata.PortDataError, match="cannot write"):
        portdata.write(object(), tmp_path)


def test_w68_arrays_mapped_read_only(tmp_path: Path) -> None:
    rays = _rays()
    handle = portdata.write(rays, tmp_path)
    mapped = portdata.read(handle, tmp_path, mmap=True)
    assert isinstance(mapped.position, np.memmap) and not mapped.position.flags.writeable
    assert same(np.asarray(mapped.position).copy(), rays.position)


def test_w68_field_bundles_with_fieldtracing_types(tmp_path: Path) -> None:
    pytest.importorskip("fieldtracing")
    from rtt_fusion.graph.field import field_power, rays_to_field

    fb = rays_to_field(_rays())
    back = round_trip(fb, tmp_path)
    assert type(back.field) is type(fb.field)  # Fieldtracing's own class
    assert same(back, fb)
    assert field_power(back) == field_power(fb)
