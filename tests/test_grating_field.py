"""Grating as an operator on fields (grating ADR, stage F1: GF1, GF2, interpolated transfer).

References: `fouriergrating.solve` for the plane wave at the beam's direction (efficiency per
order, the GUI's and the plugin's numbers), the grating's own exact path for the interpolated one.
Powers are measured spectrally (Σ Re k_z·|E|² per plane wave of each order's field).
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from rtt_fusion.couplings import grating_field as gf
from rtt_fusion.graph import Graph, Node, builtin_registry, optics, spectral
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.field import rays_to_field
from rtt_fusion.graph.light import Plane, RayBundle

fg = pytest.importorskip("fouriergrating")
pytest.importorskip("fieldtracing")

LAM = 0.55
LAMELLAR = {
    "dimension": "1D",
    "solver": "rcwa",
    "truncation": {"Nx": 8},
    "period_um": {"x": 0.7},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}},
    "superstrate": "air",
    "substrate": "glass",
    "layers": [
        {
            "type": "lamellar1d",
            "thickness_um": 0.25,
            "breakpoints_um": [0, 0.35, 0.7],
            "materials": ["glass", "air"],
        }
    ],
}


def _field(theta_deg: float, radius: float, rings: int) -> Any:
    """A collimated beam at θ in the x-z plane (polarised in that plane), as a field in the
    grating's plane z = 0 and frame (x = grating vector)."""
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": radius,
                "rings": rings,
                "center": [0, 0, 1.0],
                "direction": [0, 0, -1],
                "wavelengths_um": [LAM],
                "field_angles_deg": [[theta_deg, 0.0]],
            },
        )
    )
    rays: RayBundle = Executor(g, reg, ev).evaluate("src", "out")[0]
    log: list[str] = []
    target = gf.onto_plane(rays_to_field(rays), Plane.make([0, 0, 0], [0, 0, 1]), False, log)
    return gf.in_frame(target.field, gf.frame(target, [1.0, 0.0, 0.0]).q)


def _power(field: Any, n: float) -> float:
    s = spectral.sample(field, 4096)
    sp = spectral.spectrum(s.residual)
    kxx, kyy = s.grid.frequencies()
    k = 2 * np.pi / (LAM * 1e-3) * n
    kx, ky = kxx + s.grid.carrier[0], kyy + s.grid.carrier[1]
    kz = np.sqrt((k**2 - kx**2 - ky**2).astype(complex)).real
    return float(np.sum(kz[..., None] * np.abs(sp) ** 2))


def _plane_wave_efficiencies(field: Any) -> dict[tuple[str, tuple[int, int]], float]:
    """FourierGrating's efficiencies for the plane wave at the field's mean direction (TM)."""
    sin_x = float(np.asarray(field.gradient).mean(axis=0)[0] / (2 * np.pi / (LAM * 1e-3)))
    cfg = copy.deepcopy(LAMELLAR)
    cfg["incidence"] = {
        "lambda_um": LAM,
        "phi_deg": float(np.degrees(np.arcsin(abs(sin_x)))),
        "psi_deg": 0.0 if sin_x >= 0 else 180.0,
        "delta_deg": 90.0,
    }
    res = fg.solve(cfg)
    out = {("front", (int(o["n"]), 0)): float(o["eff"]) for o in res.reflected}
    out |= {("back", (int(o["n"]), 0)): float(o["eff"]) for o in res.transmitted}
    return out


def _efficiencies(field: Any, applied: gf.Applied) -> dict[tuple[str, tuple[int, int]], float]:
    p0 = _power(field, 1.0)
    return {
        key: _power(f, 1.0 if key[0] == "front" else 1.5) / p0 for key, f in applied.fields.items()
    }


def _handle() -> Any:
    from rtt_fusion.adapters import fouriergrating as fga

    return fga.prepare(fga.Structure(LAMELLAR), (LAM,))


def test_gf1_orders_of_a_beam_are_the_plane_wave_efficiencies() -> None:
    """GF1: a collimated beam (radius 0.15 mm, λ/D ≈ 2e−3) on a lamellar grating, exact path: the
    power of every order is FourierGrating's efficiency of the plane wave; energy kept; every
    order a field with the carrier k_c + nK."""
    field = _field(10.0, 0.15, 5)
    applied = gf.apply(_handle(), field)
    eff = _efficiencies(field, applied)
    ref = _plane_wave_efficiencies(field)
    assert set(eff) == set(ref)  # the same propagating orders
    worst = max(abs(eff[k] - ref[k]) for k in ref)
    assert worst < 2e-5, (worst, eff, ref)
    assert sum(eff.values()) == pytest.approx(1.0, abs=3e-5)
    k_step = 2e3 * np.pi / 0.7
    carrier = np.asarray(field.gradient).mean(axis=0)
    for (side, (n, _)), f in applied.fields.items():
        g = np.asarray(f.gradient).mean(axis=0)
        expected = carrier[0] + n * k_step
        assert g[0] == pytest.approx(expected if side == "back" else expected, rel=1e-9)


def test_gf2_a_wider_beam_comes_closer_to_the_plane_wave() -> None:
    """GF2: the deviation from the plane wave's efficiencies shrinks as the beam gets wider
    (0.15 → 0.3 mm radius). Interpolated transfer: the wide beam's 128² spectrum is verified on a
    33² sub-grid; the narrow one's 32² spectrum is too coarse for 1e−8 and takes the exact path."""
    handle = _handle()
    errors = []
    for radius, rings in ((0.15, 5), (0.3, 12)):
        field = _field(10.0, radius, rings)
        log: list[str] = []
        applied = gf.apply(handle, field, transfer="interpolated", log=log)
        assert any(("verified" if radius > 0.2 else "exact") in line for line in log), log
        eff = _efficiencies(field, applied)
        ref = _plane_wave_efficiencies(field)
        errors.append(max(abs(eff[k] - ref[k]) for k in ref))
    assert errors[1] < errors[0] < 2e-5, errors


def test_interpolated_transfer_equals_the_exact_one() -> None:
    """The interpolated transfer, verified to a tolerance of |T|, gives the exact path's fields
    within it (tolerance 1e−6: a 17² sub-grid of the 32² spectrum)."""
    handle = _handle()
    field = _field(10.0, 0.15, 5)
    exact = gf.apply(handle, field)
    log: list[str] = []
    interp = gf.apply(handle, field, transfer="interpolated", tolerance=1e-6, log=log)
    assert any("verified" in line for line in log), log
    assert set(exact.fields) == set(interp.fields)
    scale = max(float(np.max(np.abs(f.amplitude))) for f in exact.fields.values())
    for key, f in exact.fields.items():
        diff = np.abs(np.asarray(interp.fields[key].amplitude) - np.asarray(f.amplitude))
        assert float(diff.max()) < 2e-6 * scale, key


def test_interpolation_falls_back_near_a_rayleigh_anomaly() -> None:
    """At 12.36° the reflected order +1 grazes (sin θ + λ/Λ = 1): it appears inside the beam's
    spectrum, the interpolation is refused and the exact path taken."""
    handle = _handle()
    field = _field(float(np.degrees(np.arcsin(1 - LAM / 0.7))), 0.15, 5)
    log: list[str] = []
    applied = gf.apply(handle, field, transfer="interpolated", log=log)
    assert any("exact" in line for line in log), log
    assert ("back", (0, 0)) in applied.fields


def test_orders_can_be_chosen() -> None:
    field = _field(10.0, 0.15, 5)
    applied = gf.apply(_handle(), field, orders=lambda n, m: n == 0, transfer="interpolated")
    assert set(applied.fields) == {("front", (0, 0)), ("back", (0, 0))}


# ------------------------------------------------------------------ in the graph (stage F2)
def _graph_beam(
    theta_deg: float, radius: float = 0.15, rings: int = 5, waist: float = 0.0, lam: float = LAM
) -> RayBundle:
    reg = builtin_registry()
    ev = optics.register(reg)
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={
                "radius": radius,
                "rings": rings,
                "center": [0, 0, 1.0],
                "direction": [0, 0, -1],
                "wavelengths_um": [lam],
                "field_angles_deg": [[theta_deg, 0.0]],
                "gaussian_waist": waist,
            },
        )
    )
    beam: RayBundle = Executor(g, reg, ev).evaluate("src", "out")[0]
    return beam


def _config_file(tmp_path: Any) -> str:
    import json

    path = tmp_path / "lamellar.json"
    path.write_text(json.dumps(LAMELLAR), encoding="utf-8")
    return str(path)


def _ray_power(b: RayBundle) -> float:
    f = np.asarray(b.amplitude)[:, None] * b.field_direction()
    return float(np.sum(np.where(np.asarray(b.alive), np.sum(np.abs(f) ** 2, axis=1), 0.0)))


@pytest.mark.parametrize("side", ["reflected", "transmitted"])
def test_gf4_field_path_equals_the_ray_path(tmp_path: Any, side: str) -> None:
    """GF4: the same beam through `fg.grating` with a FourierGrating configuration, as rays (Jones
    table built around them) and as a field: the same power per order."""
    from rtt_fusion.graph import grating
    from rtt_fusion.graph.execute import Context

    node = Node(
        "gr",
        "fg.grating",
        params={
            "fouriergrating": _config_file(tmp_path),
            "side": side,
            "orders": [-1, 0, 1, 2],  # the ray path builds (−1, 0, 1) unless told otherwise
            "field_transfer": "exact",
        },
    )
    rays = _graph_beam(10.0)
    by_rays = grating.grating(Context(Graph()), node, "out", {"in": [rays]})
    p_rays = _ray_power(rays)
    ray_eff = {b.path[-1]: _ray_power(b) / p_rays for b in by_rays}
    ctx = Context(Graph())
    fb = rays_to_field(rays)
    by_field = grating.grating(ctx, node, "out", {"in": [fb]})
    n = 1.0 if side == "reflected" else 1.5
    p0 = _power(gf.onto_plane(fb, Plane.make([0, 0, 0], [0, 0, 1]), False, []).field, 1.0)
    field_eff = {b.path[-1]: _power(b.field, n) / p0 for b in by_field}
    assert set(field_eff) == set(ray_eff), (field_eff, ray_eff)
    for k, v in ray_eff.items():
        assert field_eff[k] == pytest.approx(v, abs=1e-4), (k, field_eff[k], v)
    assert all(b.field.refractive_index == pytest.approx(n) for b in by_field)


def test_a_field_at_a_table_grating_becomes_rays(tmp_path: Any) -> None:
    """With only a Jones table (no configuration) a field is converted to rays at the grating
    (the table describes rays, W21) and the orders leave as rays."""
    from rtt_fusion.graph import grating
    from rtt_fusion.graph.execute import Context

    table = str(Path(__file__).parent.parent / "reference" / "jones_table_lamellar.result.json")
    node = Node("gr", "fg.grating", params={"table": table})
    ctx = Context(Graph())
    out = grating.grating(
        ctx, node, "out", {"in": [rays_to_field(_graph_beam(0.0, 0.01, 12, lam=0.6))]}
    )
    assert out and all(isinstance(b, RayBundle) for b in out)
    assert any("converted to rays" in line for line in ctx.log)


def test_gf6_grating_field_in_another_process(tmp_path: Any) -> None:
    """GF6: the grating's field path behind a JSON channel (as in the FourierGrating plugin of
    Raytatouille Fusion): every order comes back equal to the in-process one."""
    import json

    from rtt_fusion.graph import grating
    from rtt_fusion.graph.execute import Context
    from rtt_fusion.graph.field import FieldBundle
    from rtt_fusion.graph.remote import BlockProvider, register_remote

    node = Node(
        "gr",
        "fg.grating",
        params={
            "fouriergrating": _config_file(tmp_path),
            "side": "transmitted",
            "field_transfer": "interpolated",
        },
    )
    fb = rays_to_field(_graph_beam(10.0, 0.15, 5))
    local = grating.grating(Context(Graph()), node, "out", {"in": [fb]})
    provider = BlockProvider([grating.TYPE], dict(grating.EVALUATORS))

    def call(name: str, request: dict[str, Any]) -> dict[str, Any]:
        wire = json.loads(json.dumps(request, allow_nan=False))
        answer: dict[str, Any] = json.loads(json.dumps(provider.handle(wire), allow_nan=False))
        return answer

    reg = builtin_registry()
    catalogue = [{"provider": "fg", "type": e} for e in provider.catalogue()]
    remote = register_remote(reg, catalogue, call, tmp_path / "portdata")["fg.grating"]
    again = remote(Context(Graph()), node, "out", {"in": [fb]})
    assert [b.path for b in again] == [b.path for b in local] and len(local) >= 3
    from fieldtracing.field import GridlessField

    for a, b in zip(again, local, strict=True):
        assert isinstance(a, FieldBundle) and isinstance(a.field, GridlessField)
        assert isinstance(b.field, GridlessField)
        np.testing.assert_allclose(
            np.asarray(a.field.amplitude), np.asarray(b.field.amplitude), rtol=0, atol=1e-12
        )


def test_gf3_orders_superpose_with_the_right_phases(tmp_path: Any) -> None:
    """GF3: the coherent sum of the orders −1, 0, 1 just behind a dielectric lamellar grating
    (Λ = 10 µm, normal incidence, Gaussian beam) is, at the beam's centre and divided by its
    envelope, the superposition Σ E_n·e^(inKx) of the orders of one plane wave from the adapter –
    the phase reference x = 0 and the bookkeeping ψ + nK·x of every order. (The Talbot revival
    itself needs a beam much wider than the orders' walk-off z_T·K/k ≈ 20 µm – millimetres, too
    large for a test – and is not checked here.)"""
    import json

    from rtt_fusion.adapters import fouriergrating as fga
    from rtt_fusion.graph import grating
    from rtt_fusion.graph.execute import Context
    from rtt_fusion.ports import PlaneWaves

    cfg = {
        "dimension": "1D",
        "solver": "rcwa",
        "truncation": {"Nx": 10},
        "period_um": {"x": 10.0},
        "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}},
        "superstrate": "air",
        "substrate": "air",
        "layers": [
            {
                "type": "lamellar1d",
                "thickness_um": 0.4,
                "breakpoints_um": [0, 5.0, 10.0],
                "materials": ["glass", "air"],
            }
        ],
    }
    path = tmp_path / "phase.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    node = Node(
        "gr",
        "fg.grating",
        params={
            "fouriergrating": str(path),
            "side": "transmitted",
            "orders": [-1, 0, 1],
            "field_transfer": "interpolated",
        },
    )
    waist = 0.05
    fb = rays_to_field(_graph_beam(0.0, radius=0.15, rings=30, waist=waist))
    orders = [b.field for b in grating.grating(Context(Graph()), node, "out", {"in": [fb]})]
    assert len(orders) == 3
    # every order is a field on the incident samples: their sum there needs no interpolation
    target = gf.onto_plane(fb, Plane.make([0, 0, 0], [0, 0, 1]), False, [])
    incident = gf.in_frame(target.field, gf.frame(target, [1.0, 0.0, 0.0]).q)  # grating's frame

    def total(f: Any) -> Any:
        return np.asarray(f.amplitude) * np.exp(1j * np.asarray(f.phase))[:, None]

    pts = np.asarray(incident.points)
    centre = np.hypot(pts[:, 0], pts[:, 1]) < 0.03  # well inside the beam
    assert all(np.allclose(np.asarray(f.points), pts) for f in orders)
    e_in = total(incident)[centre]
    near = sum(total(f)[centre] for f in orders) / e_in[:, :1]  # per unit incident E_x
    # reference: one plane wave (unit x polarisation) through the adapter, orders superposed
    handle = fga.prepare(fga.Structure(cfg), (LAM,))
    port = fga.front_port(handle, LAM)
    unit = PlaneWaves(port, np.array([0.0]), np.array([0.0]), np.array([[1.0, 0.0]], dtype=complex))
    waves = fga.apply_many(handle, [unit])[0]["back"]
    ez = waves.ez()
    assert waves.orders is not None
    x = pts[centre, 0]
    ref = np.zeros((len(x), 3), dtype=complex)
    for i, o in enumerate(waves.orders):
        if int(o[0]) in (-1, 0, 1):
            e = np.array([waves.field[i][0], waves.field[i][1], ez[i]])
            ref += np.exp(1j * waves.kx[i] * x)[:, None] * e[None, :]
    scale = float(np.max(np.abs(ref)))
    # 3e−3: the beam is no plane wave – its angular spread λ/(πw) ≈ 3.5e−3 rad meets transfers of
    # the orders ±1 that change linearly with the angle; a wrong sign of nK or a shifted phase
    # reference gives deviations of order 1
    assert len(x) > 20 and float(np.max(np.abs(near - ref))) < 3e-3 * scale
    intensity = np.sum(np.abs(ref) ** 2, axis=1)
    assert float(intensity.max() - intensity.min()) > 0.05 * float(intensity.max())  # a pattern
