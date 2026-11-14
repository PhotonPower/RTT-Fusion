"""Periodic structures as operators on fields (grating ADR, stages F1 and F4b).

A layer stack (FourierGrating with N = 0) and a grating act on every plane wave of an incident
field: the wave k∥ leaves in order (n, m) as k∥ + n·K_x + m·K_z with the amplitude FourierGrating
gives for it. This module applies that to a field (H) in the block's plane, without the graph:

1. `onto_plane`: the field into the plane – geometrically (field → rays → field, keeping ψ) unless
   the rays have a caustic there (a focus on the sample), then by Fieldtracing's wave propagation;
2. `in_frame`: into the block's frame (x = grating vector, z along the incident light);
3. `apply`: the field on a grid with its carrier k_c removed (`graph.spectral`), every
   significant, propagating plane wave through FourierGrating on its own (`apply_many`, exact per
   direction; on its own, so that the order of each output is known), and per side and order one
   outgoing field on the incident field's samples.

**Orders as fields of their own.** Order (n, m) of the plane wave k_c + κ has the transverse wave
vector (k_c + nK) + κ: with the carrier k_c + nK its residual spectrum sits on the same bins κ as
the incident one. Each order therefore becomes a field with the phase ψ + nK·x (gradient + nK,
same Hessian) on the same samples – no window commensurate with the period is needed (an
addendum to the ADR, section 2), and the orders are channels as in the ray path of `fg.grating`.
They share the coherence key of the incident field and interfere downstream.

**Phase reference.** FourierGrating refers its amplitudes to x = 0 of the grating; the field's
coordinates are taken from the block's point, so x = 0 is the grating's origin.

**Transfer per plane wave.** `transfer="exact"` (the ADR's default): every significant,
propagating plane wave through FourierGrating on its own. `transfer="interpolated"`: the transfer
T(κ) of every order (3 × 2) from FourierGrating on a coarse sub-grid of the spectrum, cubic in
between, **verified** at bins between the nodes against FourierGrating; above `tolerance` (relative
to the largest |T|) the sub-grid is refined, down to the exact case. A narrow spectrum (a collimated
beam spans a degree or so) needs only a few nodes; near a Rayleigh anomaly (an order appearing or
vanishing inside the spectrum) the verification fails and the exact path is taken.

**Exactness at the samples.** The transfer t0 of the carrier wave (per order, 3 × 2) is applied
pointwise on the samples and only (T(κ) − t0)·Ẽ(κ) goes through the grid, which keeps hard beam
edges free of the reconstruction's ripple.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np

from rtt_fusion._types import ComplexArray
from rtt_fusion.graph import spectral
from rtt_fusion.graph.field import (
    FieldBundle,
    FieldError,
    field_to_rays,
    propagate_to,
    rays_to_field,
)
from rtt_fusion.graph.light import Plane


def onto_plane(fb: FieldBundle, plane: Plane, from_back: bool, log: list[str]) -> FieldBundle:
    """The field in `plane` (its normal points to the front side), light running in from the front
    (or the back): geometrically unless the rays have a caustic there, then by wave propagation.
    FieldError if the field does not run into the plane."""
    from fieldtracing.field import GridlessField

    light = np.asarray(plane.normal) * (1.0 if from_back else -1.0)
    ahead = float((np.asarray(plane.origin) - np.asarray(fb.plane.origin)) @ light)
    if float(np.asarray(fb.plane.normal) @ light) <= 0 or ahead < -1e-12:
        raise FieldError("the field does not run into the block (wrong side or direction)")
    on = Plane.make(plane.origin, light, plane.v)
    if isinstance(fb.field, GridlessField):
        try:
            return rays_to_field(field_to_rays(fb), on)
        except FieldError as e:
            log.append(f"{e} – wave propagation onto the block instead")
    target, notes = propagate_to(fb, on)
    log.extend(notes)
    return target


def in_frame(field: Any, q: Any) -> Any:
    """A field with its in-plane coordinates turned by the 2×2 rotation q (new = q·old); FieldError
    for a gridded field that would need turning (its grid stays axis-parallel)."""
    from fieldtracing.field import GridlessField

    if np.allclose(q, np.eye(2), atol=1e-12):
        return field
    if not isinstance(field, GridlessField):
        raise FieldError(
            "a gridded field whose axes are not those of the block (a grating or anisotropic "
            "layers need the block's frame): turn the block's axes onto the field's"
        )
    u = np.asarray(field.amplitude)
    u = np.concatenate([u[:, :2] @ q.T, u[:, 2:]], axis=1)
    return GridlessField(
        np.asarray(field.points) @ q.T,
        np.asarray(field.phase),
        np.asarray(field.gradient) @ q.T,
        np.einsum("ij,njk,lk->nil", q, np.asarray(field.hessian), q),
        np.asarray(u, dtype=np.complex128),
        float(field.wavelength_um),
        float(field.refractive_index),
    )


@dataclass(frozen=True)
class Frame:
    """The block's frame in the plane: unit x (grating vector), y, z (along the incident light),
    and the rotation q from the field's plane axes to (x, y)."""

    x: Any
    y: Any
    z: Any
    q: Any


def frame(target: FieldBundle, plane_x: Any, keep_field_axes: bool = False) -> Frame:
    """The block's frame for a field in its plane: x = `plane_x` projected; `keep_field_axes`: use
    the field's own axes (isotropic structures, where the orientation does not matter)."""
    z = np.asarray(target.plane.normal, dtype=float)
    if keep_field_axes:
        return Frame(np.asarray(target.plane.u), np.asarray(target.plane.v), z, np.eye(2))
    x = np.asarray(plane_x, dtype=float) - z * float(np.asarray(plane_x) @ z)
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    q = np.array(
        [[x @ target.plane.u, x @ target.plane.v], [y @ target.plane.u, y @ target.plane.v]]
    )
    return Frame(x, y, z, q)


@dataclass
class Applied:
    """The result of `apply` for one field: outgoing fields by (side, (n, m)), with the share of
    the incident spectrum that was left out."""

    fields: dict[tuple[str, tuple[int, int]], Any]
    dropped: float


def apply(
    handle: Any,
    field: Any,
    *,
    orders: Callable[[int, int], bool] | None = None,
    threshold: float = 1e-10,
    grid_max: int = 1024,
    transfer: str = "exact",
    tolerance: float = 1e-8,
    log: list[str] | None = None,
    cache: dict[tuple[float, float], Any] | None = None,
) -> Applied:
    """Every significant, propagating plane wave of `field` (in the block's frame, light along +z)
    through FourierGrating; per side ("front": reflected, "back": transmitted) and order (n, m) an
    outgoing field object (GridlessField or GriddedField) in the frame of its side. `orders`
    selects orders (default: all); `cache` keeps the transfer per direction across calls with the
    same handle and order selection. SpectralError, ValueError from FourierGrating pass through."""
    from rtt_fusion.adapters import fouriergrating as fga

    lam = float(field.wavelength_um)
    sampled = spectral.sample(field, grid_max)
    grid = sampled.grid
    spec = spectral.spectrum(sampled.residual)
    kxx, kyy = grid.frequencies()
    kx, ky = kxx + grid.carrier[0], kyy + grid.carrier[1]
    port = fga.front_port(handle, lam)
    n_in = handle.front_index.real
    keep = np.hypot(kx, ky) < n_in * port.k0 * (1.0 - 1e-12)
    keep &= spectral.significant(spec[..., :2], threshold)
    total = float(np.sum(np.abs(spec[..., :2]) ** 2)) or 1.0
    dropped = float(np.sum(np.abs(spec[..., :2][~keep]) ** 2)) / total
    rows, cols = np.nonzero(keep)
    c = grid.carrier
    if transfer == "interpolated":
        done = _interpolated(
            handle, port, spec, kx, ky, keep, c, orders, tolerance, log if log is not None else []
        )
        if done is not None:
            return Applied(_back(handle, sampled, done[0], done[1], n_in), dropped)
    # exact: the transfer T (3 × 2 per side and order) of every bin from FourierGrating – one
    # S-matrix per direction serves both polarisations – applied to the spectrum; T depends on the
    # direction only, so `cache` keeps it for further fields of the same geometry (the other
    # polarisation, the next run)
    cache = cache if cache is not None else {}
    keys = [
        (round(float(kx[r, k]), 9), round(float(ky[r, k]), 9))
        for r, k in zip(rows, cols, strict=True)
    ]
    carrier_key = (round(float(c[0]), 9), round(float(c[1]), 9))
    missing = sorted({*keys, carrier_key} - set(cache))
    if missing:
        found = _transfer_at(
            handle,
            port,
            np.array([m[0] for m in missing]),
            np.array([m[1] for m in missing]),
            orders,
        )
        cache.update(zip(missing, found, strict=True))
    elif log is not None:
        log.append("transfer of every plane wave from the cache")
    specs: dict[tuple[str, tuple[int, int]], ComplexArray] = {}
    for key, r, k in zip(keys, rows, cols, strict=True):
        for so, t in cache[key].items():
            s_ = specs.setdefault(so, np.zeros_like(spec))
            s_[r, k] = t @ spec[r, k, :2]
    t0 = cache[carrier_key]
    return Applied(_back(handle, sampled, specs, t0, n_in), dropped)


def _back(
    handle: Any,
    sampled: Any,
    specs: dict[tuple[str, tuple[int, int]], ComplexArray],
    t0: dict[tuple[str, tuple[int, int]], ComplexArray],
    n_in: float,
) -> dict[tuple[str, tuple[int, int]], Any]:
    """The outgoing fields from the outgoing spectra, per side and order."""
    period = handle.structure.config.get("period_um", {})
    px, pz = float(period.get("x", 0.0) or 0.0), float(period.get("z", 0.0) or 0.0)
    kx_step = 2e3 * np.pi / px if px else 0.0
    kz_step = 2e3 * np.pi / pz if pz else 0.0
    fields: dict[tuple[str, tuple[int, int]], Any] = {}
    n_out = {"front": n_in, "back": handle.back_index.real}
    for (side, order), s in specs.items():
        shift = np.array([order[0] * kx_step, order[1] * kz_step])
        fields[(side, order)] = sampled.back(
            spectral.inverse(s), n_out[side], side == "front", t0.get((side, order)), shift
        )
    return fields


def _transfer_at(
    handle: Any, port: Any, kx: Any, ky: Any, orders: Callable[[int, int], bool] | None
) -> list[dict[tuple[str, tuple[int, int]], ComplexArray]]:
    """T (3 × 2) of every side and order at the points (kx, ky), from FourierGrating (two unit
    incidences per point; one S-matrix serves both)."""
    from rtt_fusion.adapters import fouriergrating as fga
    from rtt_fusion.ports import PlaneWaves

    units = [
        PlaneWaves(port, np.array([a]), np.array([b]), np.array([u], dtype=complex))
        for a, b in zip(kx, ky, strict=True)
        for u in ((1.0, 0.0), (0.0, 1.0))
    ]
    results = fga.apply_many(handle, units)
    out: list[dict[tuple[str, tuple[int, int]], ComplexArray]] = []
    for p in range(len(kx)):
        t: dict[tuple[str, tuple[int, int]], ComplexArray] = {}
        for j in (0, 1):
            for side, waves in results[2 * p + j].items():
                if waves.orders is None:
                    continue
                ez = waves.ez()
                for i, o in enumerate(waves.orders):
                    order = (int(o[0]), int(o[1]))
                    if orders is not None and not orders(*order):
                        continue
                    m = t.setdefault((side, order), np.zeros((3, 2), dtype=complex))
                    m[:2, j] = np.asarray(waves.field)[i]
                    m[2, j] = ez[i]
        out.append(t)
    return out


def _interpolated(
    handle: Any,
    port: Any,
    spec: ComplexArray,
    kx: Any,
    ky: Any,
    keep: Any,
    carrier: Any,
    orders: Callable[[int, int], bool] | None,
    tolerance: float,
    log: list[str],
) -> tuple[dict[tuple[str, tuple[int, int]], ComplexArray], dict[Any, ComplexArray]] | None:
    """Outgoing spectra and carrier transfers with T interpolated from a verified sub-grid (module
    docstring); None if the exact path has to be taken."""
    ny, nx = keep.shape
    # monotonic axes of the spectrum (fftshift order) and the bins in that order
    ax = np.fft.fftshift(kx[0])
    ay = np.fft.fftshift(ky[:, 0])
    rng = np.random.default_rng(0)
    t0 = _transfer_at(handle, port, np.array([carrier[0]]), np.array([carrier[1]]), orders)[0]
    previous = None
    for stride in (16, 8, 4, 2):
        ix = np.unique(np.r_[np.arange(0, nx, stride), nx - 1])
        iy = np.unique(np.r_[np.arange(0, ny, stride), ny - 1])
        if len(ix) < 4 or len(iy) < 4:
            continue
        gx, gy = np.meshgrid(ax[ix], ay[iy])
        if np.any(np.hypot(gx, gy) >= handle.front_index.real * port.k0 * (1.0 - 1e-12)):
            log.append("interpolated transfer: the spectrum reaches evanescent directions – exact")
            return None
        nodes = _transfer_at(handle, port, gx.ravel(), gy.ravel(), orders)
        keys = set(nodes[0])
        if any(set(n) != keys for n in nodes):
            log.append("interpolated transfer: an order appears inside the spectrum – exact")
            return None
        tables = {
            key: np.stack([n[key] for n in nodes]).reshape(len(iy), len(ix), 3, 2) for key in keys
        }

        def at(key: Any, qx: Any, qy: Any, tables: Any = tables, ix: Any = ix, iy: Any = iy) -> Any:
            return spectral.bicubic(ax[ix], ay[iy], tables[key], qx, qy)

        # verification at bins between the nodes
        cand = np.argwhere(keep)
        pick = cand[rng.choice(len(cand), size=min(12, len(cand)), replace=False)]
        qx, qy = kx[pick[:, 0], pick[:, 1]], ky[pick[:, 0], pick[:, 1]]
        exact = _transfer_at(handle, port, qx, qy, orders)
        scale = max(float(np.max(np.abs(tables[k]))) for k in keys) or 1.0
        error = (
            max(
                float(np.max(np.abs(at(k, qx, qy) - np.stack([e[k] for e in exact])))) for k in keys
            )
            / scale
        )
        if error <= tolerance:
            log.append(
                f"interpolated transfer: {len(ix)}×{len(iy)} nodes, verified to {error:.1e} "
                f"(tolerance {tolerance:g})"
            )
            rows, cols = np.nonzero(keep)
            specs: dict[tuple[str, tuple[int, int]], ComplexArray] = {}
            for key in keys:
                t = at(key, kx[rows, cols], ky[rows, cols])
                s = np.zeros_like(spec)
                s[rows, cols] = np.einsum("nij,nj->ni", t, spec[rows, cols, :2])
                specs[key] = s
            return specs, t0
        if previous is not None and error > previous / 2:
            # a smooth T loses ~16× per halving (cubic, h⁴), a kink ~1–2× (random check points:
            # give up below 2×): a kink – an order appearing or a Wood
            # anomaly inside the spectrum – does not converge; the exact path at once
            log.append(
                f"interpolated transfer: no convergence ({previous:.1e} → {error:.1e}), T has a "
                "kink inside the spectrum (Rayleigh/Wood anomaly) – exact"
            )
            return None
        previous = error
    log.append("interpolated transfer: not within the tolerance on a stride-2 sub-grid – exact")
    return None
