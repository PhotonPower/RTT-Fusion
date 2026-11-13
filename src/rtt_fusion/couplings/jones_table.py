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
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.adapters import fouriergrating as fga
from rtt_fusion.ports import PlaneWaves


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
    handle: fga.Handle,
    wavelength_um: float,
    sigma_x: float,
    sigma_y: float,
    side: str = "transmitted",
) -> dict[int, ComplexArray]:
    """J_m for every propagating order at one point, directly from FourierGrating (4 solves)."""
    pair = fga.apply_many(handle, _unit_incidents(handle, wavelength_um, sigma_x, sigma_y))
    return _jones_of(handle, wavelength_um, sigma_x, sigma_y, side, pair)


def _unit_incidents(
    handle: fga.Handle, wavelength_um: float, sigma_x: float, sigma_y: float
) -> list[PlaneWaves]:
    """The two unit incidences (x, y polarised) at one point of the table."""
    port = fga.front_port(handle, wavelength_um)
    kx, ky = sigma_x * port.k0, sigma_y * port.k0
    return [
        PlaneWaves(port, np.array([kx]), np.array([ky]), np.array([unit], dtype=np.complex128))
        for unit in ((1.0, 0.0), (0.0, 1.0))
    ]


def _jones_of(
    handle: fga.Handle,
    wavelength_um: float,
    sigma_x: float,
    sigma_y: float,
    side: str,
    pair: list[dict[str, PlaneWaves]],
) -> dict[int, ComplexArray]:
    """J_m from the results of the two unit incidences."""
    port = fga.front_port(handle, wavelength_um)
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
    handle: fga.Handle,
    wavelength_um: float,
    sigmas: list[tuple[float, float]],
    side: str = "transmitted",
) -> list[dict[int, ComplexArray]]:
    """J_m of every propagating order at many points (σ_x, σ_y) of one wavelength, directly from
    FourierGrating in one `apply_many` (two unit incidences per point)."""
    incidents = [inc for sx, sy in sigmas for inc in _unit_incidents(handle, wavelength_um, sx, sy)]
    results = fga.apply_many(handle, incidents)
    return [
        _jones_of(handle, wavelength_um, sx, sy, side, results[2 * j : 2 * j + 2])
        for j, (sx, sy) in enumerate(sigmas)
    ]


def build_table(
    handle: fga.Handle,
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
        results = fga.apply_many(handle, incidents)
        for j, (b, c, sx, sy) in enumerate(points):
            d = _jones_of(handle, float(lam), sx, sy, side, results[2 * j : 2 * j + 2])
            for o, m in enumerate(orders):
                if m in d:
                    jones[o, a, b, c] = d[m]
    cfg = handle.structure.config
    n_out = (handle.back_index if side == "transmitted" else handle.front_index).real
    return JonesTable(
        tuple(orders),
        side,
        np.asarray(wavelengths_um, dtype=np.float64),
        np.asarray(sigma_x, dtype=np.float64),
        np.asarray(sigma_y, dtype=np.float64),
        float(cfg["period_um"]["x"]),
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
