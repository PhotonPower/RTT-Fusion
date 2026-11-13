"""Layer stack `fg.stack` with fields H (grating ADR, stage F4b; ST3, GF6 for the stack).

References independent of FourierGrating: the transfer matrix of `test_stack` (Macleod) per plane
wave of the field's own spectrum, the ray path of the same block (exact per ray, ST0–ST2) and the
in-process evaluation for the field across a process boundary. Powers are measured spectrally
(Σ Re k_z·|E|² per plane wave), which holds for gridded and gridless fields alike.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.graph import Graph, Node, builtin_registry, optics, spectral, stack
from rtt_fusion.graph.execute import WITH_FIELDS, Context, Executor
from rtt_fusion.graph.field import FieldBundle, rays_to_field
from rtt_fusion.graph.light import Plane, RayBundle

pytest.importorskip("fouriergrating")
pytest.importorskip("fieldtracing")

from test_stack import _config, tmm  # noqa: E402

LAM = 0.55
BRAGG = [(2.4, LAM / (4 * 2.4)), (1.46, LAM / (4 * 1.46))] * 2


def _collimated(theta_deg: float, rings: int = 12) -> RayBundle:
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": 0.5,
                "rings": rings,
                "center": [0, 0, 1.0],
                "direction": [0, 0, -1],
                "wavelengths_um": [LAM],
                "field_angles_deg": [[theta_deg, 0.0]],
            },
        )
    )
    bundle: RayBundle = Executor(g, reg, ev).evaluate("src", "out")[0]
    return bundle


def _focused(na: float) -> RayBundle:
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    f = 2.0
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": na * f,
                "rings": 16,
                "center": [0, 0, 3.0],
                "direction": [0, 0, -1],
                "wavelengths_um": [LAM],
            },
        )
    )
    g.add(
        Node(
            "lens",
            "optics.aplanatic_lens",
            params={
                "point": [0, 0, f],
                "axis": [0, 0, -1],
                "focal_length": f,
                "aperture_radius": 1,
            },
        )
    )
    g.connect(("src", "out"), ("lens", "in"))
    bundle: RayBundle = Executor(g, reg, ev).evaluate("lens", "out")[0]
    return bundle


def _power(field: Any, n: float) -> float:
    """Σ Re k_z·|E|² over the plane waves of a field (same arbitrary factor for all fields on the
    same wavelength)."""
    s = spectral.sample(field, 4096)
    sp = spectral.spectrum(s.residual)
    kxx, kyy = s.grid.frequencies()
    k = 2 * np.pi / (LAM * 1e-3) * n
    kx, ky = kxx + s.grid.carrier[0], kyy + s.grid.carrier[1]
    kz = np.sqrt((k**2 - kx**2 - ky**2).astype(complex)).real
    return float(np.sum(kz[..., None] * np.abs(sp) ** 2))


def _both(path: str, fb: FieldBundle) -> tuple[FieldBundle, FieldBundle, FieldBundle, Context]:
    node = Node("st", stack.TYPE_ID, params={"fouriergrating": path})
    ctx = Context(Graph())
    plane = stack._plane(node)
    reflected = stack._apply_field(ctx, node, fb, plane, False, True)
    transmitted = stack._apply_field(ctx, node, fb, plane, False, False)
    incident = stack._to_plane(Context(Graph()), node, fb, plane, False)
    return incident, reflected, transmitted, ctx


@pytest.mark.parametrize(("layers", "n_sub", "theta_deg"), [([], 1.5, 20.0), (BRAGG, 1.52, 25.0)])
def test_st3_collimated_field_against_the_transfer_matrix(
    tmp_path: Path, layers: list[tuple[float, float]], n_sub: float, theta_deg: float
) -> None:
    """A collimated beam (p polarised in the plane of incidence) as a field: reflected and
    transmitted power as the transfer matrix at its angle (its angular spread, λ/D ≈ 1e−3 rad,
    changes R by less than the tolerance); energy kept."""
    path = _config(tmp_path, 1.0, layers, n_sub, "stack")
    fb = rays_to_field(_collimated(theta_deg))
    incident, r_field, t_field, ctx = _both(path, fb)
    assert not any("left out" in line for line in ctx.log)
    p0 = _power(incident.field, 1.0)
    r, t = _power(r_field.field, 1.0) / p0, _power(t_field.field, n_sub) / p0
    r_ref, t_ref = tmm(1.0, layers, n_sub, np.radians(theta_deg), LAM, "p")
    assert r == pytest.approx(r_ref, abs=2e-5) and t == pytest.approx(t_ref, abs=2e-5)
    assert r + t == pytest.approx(1.0, abs=2e-5)
    assert r_field.plane.normal @ np.array([0, 0, 1.0]) > 0  # reflected: back up
    assert t_field.field.refractive_index == pytest.approx(n_sub)


def test_st3_focus_on_the_stack_equals_the_ray_path(tmp_path: Path) -> None:
    """A focus (NA 0.1) on a Bragg stack: the field reaches it by wave propagation (gridded); its
    reflectance equals the ray path's (exact per ray) and the transfer matrix applied to every
    plane wave of the field's spectrum."""
    path = _config(tmp_path, 1.0, BRAGG, 1.52, "bragg")
    rays = _focused(0.1)
    fb = rays_to_field(rays, Plane.make([0, 0, 0.5], [0, 0, -1]))
    incident, r_field, t_field, ctx = _both(path, fb)
    assert any("wave propagation" in line for line in ctx.log)
    p0 = _power(incident.field, 1.0)
    r, t = _power(r_field.field, 1.0) / p0, _power(t_field.field, 1.52) / p0
    node = Node("st", stack.TYPE_ID, params={"fouriergrating": path})
    plane = stack._plane(node)

    def ray_power(b: RayBundle | None) -> float:
        assert b is not None
        f = np.asarray(b.amplitude)[:, None] * b.field_direction()
        return float(np.sum(np.where(np.asarray(b.alive), np.sum(np.abs(f) ** 2, axis=1), 0.0)))

    c = Context(Graph())
    p_rays = ray_power(rays)
    r_rays = ray_power(stack._apply(c, node, rays, plane, False, True)) / p_rays
    assert r == pytest.approx(r_rays, abs=5e-5) and r + t == pytest.approx(1.0, abs=1e-4)
    assert 0.5 < r < 0.9  # content guard: the stack reflects a good part


def test_st3_field_link_in_a_graph(tmp_path: Path) -> None:
    """With transport H on the link the block gets and returns fields (no conversion to rays);
    the camera behind the AR-coated window sees the light (its link carries rays: a field image
    would propagate the field onto the pixels, which is not what this test is about)."""
    n_c = float(np.sqrt(1.52))
    path = _config(tmp_path, 1.0, [(n_c, LAM / (4 * n_c))], 1.52, "ar")
    reg = builtin_registry()
    ev = optics.register(reg)
    ev.update(stack.register(reg))
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": 0.5,
                "rings": 10,
                "center": [0, 0, 1.0],
                "direction": [0, 0, -1],
                "wavelengths_um": [LAM],
            },
        )
    )
    g.add(Node("st", stack.TYPE_ID, params={"fouriergrating": path}))
    g.add(
        Node(
            "cam",
            "optics.detector",
            params={
                "point": [0, 0, -0.5],
                "normal": [0, 0, 1],
                "half_size": [0.8, 0.8],
                "pixels": [16, 16],
            },
        )
    )
    link = g.connect(("src", "out"), ("st", "front_in"))
    g.set_transport(link, "H")
    g.set_transport(g.connect(("st", "back_out"), ("cam", "in")), "R")
    ex = Executor(g, reg, ev, available=WITH_FIELDS)
    out = ex.evaluate("st", "back_out")
    assert out and all(isinstance(b, FieldBundle) for b in out)
    image = ex.evaluate("cam", "image")
    assert float(np.asarray(image.values).sum()) > 0


def test_st3_a_field_running_away_is_refused(tmp_path: Path) -> None:
    from rtt_fusion.graph.execute import ExecutionError

    path = _config(tmp_path, 1.0, [], 1.5, "bare")
    fb = rays_to_field(_collimated(0.0, rings=6))
    node = Node("st", stack.TYPE_ID, params={"fouriergrating": path})
    with pytest.raises(ExecutionError, match="does not run into the stack"):
        stack.stack(Context(Graph()), node, "front_out", {"back_in": [fb]})


def test_gf6_stack_field_in_another_process(tmp_path: Path) -> None:
    """The field path of the stack behind a JSON channel (a block provider, as the FourierGrating
    plugin of Raytatouille Fusion runs it): the field crosses the process border as port data and
    comes back equal to the in-process one."""
    from rtt_fusion.graph.remote import BlockProvider, register_remote

    path = _config(tmp_path, 1.0, BRAGG, 1.52, "bragg")
    fb = rays_to_field(_collimated(15.0, rings=8))
    node = Node("st", stack.TYPE_ID, params={"fouriergrating": path})
    local = stack.stack(Context(Graph()), node, "front_out", {"front_in": [fb]})[0]
    provider = BlockProvider([stack.TYPE], dict(stack.EVALUATORS))
    calls: list[str] = []

    def call(name: str, request: dict[str, Any]) -> dict[str, Any]:
        wire = json.loads(json.dumps(request, allow_nan=False))
        calls.append(name)
        answer: dict[str, Any] = json.loads(json.dumps(provider.handle(wire), allow_nan=False))
        return answer

    reg = builtin_registry()
    catalogue = [{"provider": "fg", "type": e} for e in provider.catalogue()]
    remote = register_remote(reg, catalogue, call, tmp_path / "portdata")[stack.TYPE_ID]
    again = remote(Context(Graph()), node, "front_out", {"front_in": [fb]})[0]
    assert calls == ["fg"]
    from fieldtracing.field import GridlessField

    assert isinstance(again, FieldBundle) and isinstance(again.field, GridlessField)
    assert isinstance(local.field, GridlessField)  # a collimated beam: the geometric route
    np.testing.assert_allclose(
        np.asarray(again.field.amplitude), np.asarray(local.field.amplitude), rtol=0, atol=1e-12
    )
    np.testing.assert_allclose(np.asarray(again.field.points), np.asarray(local.field.points))
