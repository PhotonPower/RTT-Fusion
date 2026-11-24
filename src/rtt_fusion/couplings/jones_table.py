"""Tabulated diffraction interaction (stage R6): Jones matrix per order over (λ, incidence).

Prototype of the interaction proposed in docs/adr-entwurf-beugungstabelle.md.

**Basis.** The table works in the transverse Cartesian components of the local surface frame
(x = grating vector, y = grooves, z = surface normal, light towards +z): an incident wave with
tangential direction cosines σ = (σ_x, σ_y) = k⊥/k₀ and transverse field (E_x, E_y) leaves in
order m with σ_m = σ + m·λ/Λ·x̂ and transverse field J_m·(E_x, E_y); E_z follows from
transversality. This basis has no degeneracy at normal incidence – unlike an s/p basis, whose
s is undefined there (Raytatouille resolves it with an axis rule, rtt/polar/prt.hpp, which
would make the result depend on that choice for an anisotropic grating).

**Normalisation.** J_m is power normalised: the field Jones matrix times
sqrt(Re k_z,m / Re k_z,in), so that for a unit field (|E| = 1, with E_z) the power fraction of the
order is |E_m|² with E_m the full field built from J_m·E⊥ and transversality – the semantics of
Raytatouille's P (ADR 0021).

**PRT matrix.** For unit wave vectors k_in, k_m (global = local here):
P_m = L_m J_m Π (I − k_in k_inᵀ) + k_m k_inᵀ, with Π the projection on (x, y) and L_m the lift
(E_x, E_y) → (E_x, E_y, −(σ_x,m E_x + σ_y,m E_y)/(n cos θ_m)) in direction-cosine form. Then
P_m k_in = k_m and P_m E = E_m for transverse E.

**Phase reference.** J_m refers to the local origin x = 0 of the grating (FourierGrating's x
origin = Raytatouille's φ(0, 0) = 0, ADR 0025 point 1); the phase m·φ(x, y) at other hit points is
the phase layer's own contribution.

**Solvers.** The table is built from any prepared grating whose adapter gives port waves:
FourierGrating (`adapters.fouriergrating`) and, since stage R35, hp-FEM (`adapters.hpfem`, which
refers its Jones matrices to x = 0 as well).
"""

from __future__ import annotations

import importlib
from dataclasses import dataclass
from typing import Any

import numpy as np

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.ports import PlaneWaves

#: a prepared grating: FourierGrating's handle (`adapters.fouriergrating`) or hp-FEM's
#: (`adapters.hpfem`) – anything whose adapter module has `front_port` and `apply_many` and that
#: has `front_index`, `back_index` and a period
GratingHandle = Any


def _adapter(handle: GratingHandle) -> Any:
    """The adapter module a handle comes from (its `front_port`, `apply_many`)."""
    return importlib.import_module(type(handle).__module__)


def _period_um(handle: GratingHandle) -> float:
    period = getattr(handle, "period_um", None)
    if period is not None:
        return float(period)
    return float(handle.structure.config["period_um"]["x"])


@dataclass(frozen=True)
class JonesTable:
    """Power-normalised Jones matrices J[order, λ, σ_x, σ_y, 2, 2] (NaN where not propagating).

    Attributes:
        orders: order numbers m (1D gratings, grating vector along local x).
        side: "transmitted" or "reflected".
        wavelengths_um, sigma_x, sigma_y: grid axes (σ = k⊥/k₀ of the incident wave).
        period_um: grating period Λ.
        n_in, n_out: real indices of the incident medium and of the medium of the order.
        jones: complex array (orders, λ, σ_x, σ_y, 2, 2).
    """

    orders: tuple[int, ...]
    side: str
    wavelengths_um: FloatArray
    sigma_x: FloatArray
    sigma_y: FloatArray
    period_um: float
    n_in: float
    n_out: float
    jones: ComplexArray


def direct_jones(
    handle: GratingHandle,
    wavelength_um: float,
    sigma_x: float,
    sigma_y: float,
    side: str = "transmitted",
) -> dict[int, ComplexArray]:
    """J_m for every propagating order at one point, directly from the solver (FourierGrating:
    2 solves; hp-FEM: one `grating.jones`, i.e. 2 solves)."""
    pair = _adapter(handle).apply_many(
        handle, _unit_incidents(handle, wavelength_um, sigma_x, sigma_y)
    )
    return _jones_of(handle, wavelength_um, sigma_x, sigma_y, side, pair)


def _unit_incidents(
    handle: GratingHandle, wavelength_um: float, sigma_x: float, sigma_y: float
) -> list[PlaneWaves]:
    """The two unit incidences (x, y polarised) at one point of the table."""
    port = _adapter(handle).front_port(handle, wavelength_um)
    kx, ky = sigma_x * port.k0, sigma_y * port.k0
    return [
        PlaneWaves(port, np.array([kx]), np.array([ky]), np.array([unit], dtype=np.complex128))
        for unit in ((1.0, 0.0), (0.0, 1.0))
    ]


def _jones_of(
    handle: GratingHandle,
    wavelength_um: float,
    sigma_x: float,
    sigma_y: float,
    side: str,
    pair: list[dict[str, PlaneWaves]],
) -> dict[int, ComplexArray]:
    """J_m from the results of the two unit incidences."""
    port = _adapter(handle).front_port(handle, wavelength_um)
    k0 = port.k0
    kx, ky = sigma_x * k0, sigma_y * k0
    columns: dict[int, list[ComplexArray]] = {}
    out_side = "back" if side == "transmitted" else "front"
    n_out = (handle.back_index if side == "transmitted" else handle.front_index).real
    kz_in = float(np.real(port.kz(np.array(kx), np.array(ky))))
    for result in pair:
        waves = result[out_side]
        assert waves.orders is not None
        for j, order in enumerate(waves.orders):
            m = int(order[0])
            kz_m = np.sqrt(complex((n_out * k0) ** 2 - waves.kx[j] ** 2 - waves.ky[j] ** 2))
            scale = np.sqrt(np.real(kz_m) / kz_in)
            columns.setdefault(m, []).append(np.asarray(waves.field[j]) * scale)
    return {m: np.column_stack(c) for m, c in columns.items() if len(c) == 2}


def jones_many(
    handle: GratingHandle,
    wavelength_um: float,
    sigmas: list[tuple[float, float]],
    side: str = "transmitted",
) -> list[dict[int, ComplexArray]]:
    """J_m of every propagating order at many points (σ_x, σ_y) of one wavelength, directly from
    the solver in one `apply_many` (two unit incidences per point)."""
    incidents = [inc for sx, sy in sigmas for inc in _unit_incidents(handle, wavelength_um, sx, sy)]
    results = _adapter(handle).apply_many(handle, incidents)
    return [
        _jones_of(handle, wavelength_um, sx, sy, side, results[2 * j : 2 * j + 2])
        for j, (sx, sy) in enumerate(sigmas)
    ]


def build_table(
    handle: GratingHandle,
    wavelengths_um: FloatArray,
    sigma_x: FloatArray,
    sigma_y: FloatArray,
    orders: tuple[int, ...],
    side: str = "transmitted",
) -> JonesTable:
    """Tabulate J_m on a regular grid (NaN where the order does not propagate)."""
    shape = (len(orders), len(wavelengths_um), len(sigma_x), len(sigma_y), 2, 2)
    jones = np.full(shape, np.nan + 0j, dtype=np.complex128)
    for a, lam in enumerate(wavelengths_um):
        # every point of a wavelength in one apply_many (L4): with an operator from
        # FourierGrating (plan L4, O3) one eigenproblem per direction instead of two solves
        points = [
            (b, c, float(sx), float(sy))
            for b, sx in enumerate(sigma_x)
            for c, sy in enumerate(sigma_y)
        ]
        incidents = [
            inc for _, _, sx, sy in points for inc in _unit_incidents(handle, float(lam), sx, sy)
        ]
        results = _adapter(handle).apply_many(handle, incidents)
        for j, (b, c, sx, sy) in enumerate(points):
            d = _jones_of(handle, float(lam), sx, sy, side, results[2 * j : 2 * j + 2])
            for o, m in enumerate(orders):
                if m in d:
                    jones[o, a, b, c] = d[m]
    n_out = (handle.back_index if side == "transmitted" else handle.front_index).real
    return JonesTable(
        tuple(orders),
        side,
        np.asarray(wavelengths_um, dtype=np.float64),
        np.asarray(sigma_x, dtype=np.float64),
        np.asarray(sigma_y, dtype=np.float64),
        _period_um(handle),
        handle.front_index.real,
        n_out,
        jones,
    )


def point_table(
    handle: GratingHandle,
    wavelength_um: float,
    sigma_x: float,
    sigma_y: float,
    orders: tuple[int, ...],
    side: str = "transmitted",
    half_width: float = 1e-7,
) -> JonesTable:
    """J_m exactly at one incidence (λ, σ) – a plane wave, e.g. one field point of a Köhler
    illumination – as a table: every node of a tiny 4×4×4 grid (± `half_width`) holds the direct
    solution, so the interpolation returns it unchanged (4 solves instead of a table's 4·n² nodes,
    and no stencil that reaches over a nearby Rayleigh anomaly). `handle` must be prepared for
    `wavelength_um`."""
    d = direct_jones(handle, wavelength_um, sigma_x, sigma_y, side)
    offsets = np.array([-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0]) * half_width
    jones = np.full((len(orders), 4, 4, 4, 2, 2), np.nan + 0j, dtype=np.complex128)
    for o, m in enumerate(orders):
        if m in d:
            jones[o] = d[m]
    n_out = (handle.back_index if side == "transmitted" else handle.front_index).real
    return JonesTable(
        tuple(orders),
        side,
        wavelength_um + offsets,
        sigma_x + offsets,
        sigma_y + offsets,
        _period_um(handle),
        handle.front_index.real,
        n_out,
        jones,
    )


def _cubic_weights(axis: FloatArray, x: float) -> tuple[np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    """Indices and weights of 4-point Lagrange (cubic) interpolation on a uniform axis."""
    h = float(axis[1] - axis[0])
    i = int(np.clip(np.floor((x - axis[0]) / h), 1, len(axis) - 3))
    idx = np.arange(i - 1, i + 3)
    nodes = axis[idx]
    w = np.array(
        [
            np.prod([(x - nodes[k]) / (nodes[j] - nodes[k]) for k in range(4) if k != j])
            for j in range(4)
        ]
    )
    return idx, w


def _linear_weights(axis: FloatArray, x: float) -> tuple[np.ndarray, np.ndarray]:  # type: ignore[type-arg]
    h = float(axis[1] - axis[0])
    i = int(np.clip(np.floor((x - axis[0]) / h), 0, len(axis) - 2))
    t = (x - axis[i]) / h
    return np.array([i, i + 1]), np.array([1 - t, t])


def interpolate(
    table: JonesTable,
    order: int,
    wavelength_um: float,
    sigma_x: float,
    sigma_y: float,
    method: str = "cubic",
) -> ComplexArray:
    """J_m at (λ, σ) by tensor-product interpolation of the complex entries (uniform axes).

    Returns NaN if a contributing grid node is not propagating (order appears or vanishes in the
    stencil: Rayleigh anomaly, see `near_rayleigh_anomaly`).
    """
    weights = _cubic_weights if method == "cubic" else _linear_weights
    o = table.orders.index(order)
    ia, wa = weights(table.wavelengths_um, wavelength_um)
    ib, wb = weights(table.sigma_x, sigma_x)
    ic, wc = weights(table.sigma_y, sigma_y)
    block = table.jones[o][np.ix_(ia, ib, ic)]
    return np.asarray(np.einsum("a,b,c,abcij->ij", wa, wb, wc, block))


def near_rayleigh_anomaly(
    table: JonesTable,
    wavelength_um: float,
    sigma_x: float,
    sigma_y: float,
    method: str = "cubic",
    max_order: int = 8,
) -> bool:
    """True if a Rayleigh anomaly crosses the interpolation stencil of the query point: some order
    m′ (|m′| ≤ max_order) is propagating at one node of the stencil box and evanescent at another,
    in the incident or the exit medium (|σ + m′λ/Λ x̂| − n changes sign). There every order's J
    has a square-root branch point and the interpolation is not reliable."""
    weights = _cubic_weights if method == "cubic" else _linear_weights
    ia, _ = weights(table.wavelengths_um, wavelength_um)
    ib, _ = weights(table.sigma_x, sigma_x)
    ic, _ = weights(table.sigma_y, sigma_y)
    lam = table.wavelengths_um[[ia[0], ia[-1]]]
    sx = table.sigma_x[[ib[0], ib[-1]]]
    sy = table.sigma_y[[ic[0], ic[-1]]]
    ll, xx, yy = np.meshgrid(lam, sx, sy, indexing="ij")
    for m in range(-max_order, max_order + 1):
        rho = np.hypot(xx + m * ll / table.period_um, yy)
        # also the interior: rho is monotone along each axis except where σ + mλ/Λ crosses 0
        for n in {table.n_in, table.n_out}:
            d = rho - n
            if d.min() < 0 < d.max():
                return True
    return False


#: use the C++ extension `rtt_fusion._diffraction` where it is installed (tests switch it off to
#: compare with the Python loop)
USE_CPP = True


def prt_many(
    table: JonesTable, order: int, wavelength_um: float, sigma: FloatArray, alive: Any
) -> tuple[ComplexArray, Any] | None:
    """The local PRT matrices (n, 3, 3) of `order` for n incidences `sigma` (n, 2) and where they
    are valid (n,): alive, σ on the table's axes, no Rayleigh anomaly over the stencil (orders up
    to ±8), every stencil node propagating – `near_rayleigh_anomaly`, `interpolate` and
    `prt_matrix` per ray in one call of the C++ table (stage R12). None without the extension
    (or with `USE_CPP` off)."""
    if not USE_CPP:
        return None
    try:
        _diffraction = importlib.import_module("rtt_fusion._diffraction")
    except ImportError:
        return None
    prt, ok = _diffraction.prt_many(
        list(table.orders),
        table.side == "reflected",
        np.ascontiguousarray(table.wavelengths_um, dtype=np.float64),
        np.ascontiguousarray(table.sigma_x, dtype=np.float64),
        np.ascontiguousarray(table.sigma_y, dtype=np.float64),
        float(table.period_um),
        float(table.n_in),
        float(table.n_out),
        np.ascontiguousarray(table.jones, dtype=np.complex128),
        int(order),
        float(wavelength_um),
        np.ascontiguousarray(sigma, dtype=np.float64),
        np.ascontiguousarray(alive, dtype=bool),
    )
    return np.asarray(prt), np.asarray(ok)


def prt_matrix(
    jones: ComplexArray,
    sigma_in: tuple[float, float],
    sigma_out: tuple[float, float],
    n_in: float,
    n_out: float,
    reflected: bool = False,
) -> ComplexArray:
    """3×3 PRT matrix of the order (module docstring): P k_in = k_m, P E = E_m."""
    k_in = (
        np.array([sigma_in[0], sigma_in[1], np.sqrt(n_in**2 - sigma_in[0] ** 2 - sigma_in[1] ** 2)])
        / n_in
    )
    cz = np.sqrt(n_out**2 - sigma_out[0] ** 2 - sigma_out[1] ** 2)
    sign = -1.0 if reflected else 1.0  # reflected orders run in −z (n_out = n_in)
    k_m = np.array([sigma_out[0], sigma_out[1], sign * cz]) / n_out
    lift = np.array(
        [[1.0, 0.0], [0.0, 1.0], [-sign * sigma_out[0] / cz, -sign * sigma_out[1] / cz]]
    )
    proj = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    transverse = np.eye(3) - np.outer(k_in, k_in)
    return np.asarray(lift @ jones @ proj @ transverse + np.outer(k_m, k_in))
