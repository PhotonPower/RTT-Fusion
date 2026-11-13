"""Adapter for FourierGrating (Fourier modal method / RCWA), stage R1.

Interface: docs/schnittstelle.md; conventions and their conversion: docs/konventionen.md.

FourierGrating (FG) layers in y with light from above (−y). This adapter maps the frames by a
rotation about x: x′ = x, y′ = z, z′ = −y (RTT frame primed). The front port of the structure is
the top interface of the stack (z′ = 0), the back port the bottom interface (z′ = total
thickness), matching FG's local phase references (reflected orders at the top, transmitted
orders at the bottom interface). Fields cross the boundary as full Cartesian vectors, so FG's
s/p bases (fg-conv: s = (γ, 0, −α)/|κ|, s = (0, 0, 1) at κ = 0; incident
s_in = (sin ψ, 0, −cos ψ), p_in = (cos φ cos ψ, sin φ, cos φ sin ψ)) never leak out. An
arbitrary incident polarisation is the superposition of the runs δ = 0° and δ = 90°.

Representation: S (plane waves, docs/konventionen.md). Only propagating orders are returned
(FG reports only those); superstrate and substrate must be non-absorbing isotropic media.
"""

from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from rtt_fusion import __version__
from rtt_fusion._types import ComplexArray, FloatArray
from rtt_fusion.ports import PlaneWaves, Port

ADAPTER_VERSION = "0.1.0"


@dataclass(frozen=True)
class Structure:
    """FG configuration without incidence: dimension, solver, truncation, period_um, materials,
    superstrate, substrate, layers and solver settings (FG `docs/config_format.md`), lengths in
    µm as FG expects. The x origin of the FG geometry is x′ = 0 of the ports."""

    config: dict[str, Any]


@dataclass
class Handle:
    structure: Structure
    wavelengths_um: tuple[float, ...]
    content_hash: str
    front_index: complex
    back_index: complex
    thickness_mm: float
    _sims: dict[int, Any] = field(default_factory=dict, repr=False)


def _fg() -> Any:
    import fouriergrating

    return fouriergrating


def _index(config: dict[str, Any], name: str) -> complex:
    material = config["materials"][name]
    if "n" not in material:
        raise ValueError(f"medium {name!r} at a port must be isotropic (key 'n')")
    n = material["n"]
    value = complex(n[0], n[1]) if isinstance(n, list) else complex(n)
    if value.imag != 0.0:
        raise ValueError(f"medium {name!r} at a port must be non-absorbing")
    return value


def _thickness_um(layers: list[dict[str, Any]]) -> float:
    total = 0.0
    for layer in layers:
        total += float(layer.get("thickness_um", layer.get("depth_um", 0.0)))
    return total


def describe() -> dict[str, Any]:
    fg = _fg()
    return {
        "name": "FourierGrating",
        "version": str(fg.__version__),
        "adapter_version": ADAPTER_VERSION,
        "method": "Fourier modal method: RCWA or differential method (implicit Runge-Kutta)",
        "representation": "S",
        "ports": [
            {"name": "front", "plane": "top interface, z' = 0", "in": "+z", "out": "-z"},
            {"name": "back", "plane": "bottom interface, z' = thickness", "out": "+z"},
        ],
        "linear": True,
        "reciprocal": True,
        "operator": False,
        "deterministic": {"cpu": "bitwise", "cuda": "in tolerance"},
        "conventions_of_solver": [
            "layering in y, light from +y: frame x' = x, y' = z, z' = -y",
            "order amplitudes (E.s, E.p), s = (gamma, 0, -alpha)/|kappa|, at kappa = 0 s = z",
            "incident A = cos(delta) s_in + sin(delta) p_in, s_in = (sin psi, 0, -cos psi)",
            "phase reference: top interface (reflected), bottom interface (transmitted)",
            "units um and degrees",
        ],
        "validity": [
            "only propagating orders are returned",
            "superstrate and substrate non-absorbing and isotropic",
            "convergence in the truncation (Nx, Nz) is the user's responsibility",
        ],
    }


def prepare(structure: Structure, wavelengths_um: list[float] | tuple[float, ...]) -> Handle:
    """Check the configuration, compute the content hash; no computation."""
    cfg = structure.config
    for key in ("materials", "superstrate", "substrate", "layers", "period_um", "dimension"):
        if key not in cfg:
            raise ValueError(f"configuration lacks {key!r}")
    if not wavelengths_um:
        raise ValueError("at least one wavelength is required")
    front, back = _index(cfg, cfg["superstrate"]), _index(cfg, cfg["substrate"])
    h = hashlib.sha256()
    h.update(f"rtt_fusion {__version__} fg-adapter {ADAPTER_VERSION}".encode())
    h.update(f"fouriergrating {_fg().__version__}".encode())
    clean = {k: v for k, v in cfg.items() if k not in ("incidence", "output", "sweep")}
    h.update(json.dumps(clean, sort_keys=True).encode())
    h.update(np.asarray(wavelengths_um, dtype=np.float64).tobytes())
    return Handle(
        structure,
        tuple(float(w) for w in wavelengths_um),
        h.hexdigest(),
        front,
        back,
        _thickness_um(cfg["layers"]) * 1e-3,
    )


def front_port(handle: Handle, wavelength_um: float, direction: int = 1) -> Port:
    return Port(0.0, direction, handle.front_index, wavelength_um)


def back_port(handle: Handle, wavelength_um: float) -> Port:
    return Port(handle.thickness_mm, 1, handle.back_index, wavelength_um)


def _simulation(handle: Handle, wavelength_um: float) -> Any:
    key = int(np.argmin([abs(w - wavelength_um) for w in handle.wavelengths_um]))
    if abs(handle.wavelengths_um[key] - wavelength_um) > 1e-12:
        raise ValueError("wavelength was not prepared")
    if key not in handle._sims:
        cfg = copy.deepcopy(handle.structure.config)
        cfg["incidence"] = {
            "lambda_um": wavelength_um,
            "phi_deg": 0.0,
            "psi_deg": 0.0,
            "delta_deg": 0.0,
        }
        cfg.pop("sweep", None)
        cfg.setdefault("output", {"dir": "out/rtt_fusion"})
        handle._sims[key] = _fg().Simulation(cfg)
    return handle._sims[key]


def _order_basis(
    alpha: float, beta: float, gamma: float, up: bool
) -> tuple[FloatArray, FloatArray]:
    """FG's s/p basis of an order (fg-conv, efficiency.cpp): s = (γ, 0, −α)/κ, (0, 0, 1) at κ = 0;
    p = k̂ × s with k = (α, ±β, γ) (+ for upward, reflected orders)."""
    kappa = float(np.hypot(alpha, gamma))
    s = np.array([0.0, 0.0, 1.0]) if kappa < 1e-12 else np.array([gamma, 0.0, -alpha]) / kappa
    k = np.array([alpha, beta if up else -beta, gamma])
    k = k / np.linalg.norm(k)
    return s, np.cross(k, s)


def _to_rtt(e_fg: ComplexArray) -> ComplexArray:
    """Cartesian FG field (x, y, z) → RTT frame (x′, y′, z′) = (x, z, −y)."""
    return np.array([e_fg[0], e_fg[2], -e_fg[1]])


def apply(handle: Handle, incident: PlaneWaves) -> dict[str, PlaneWaves]:
    """Plane waves incident on the front port → reflected (front, −z) and transmitted (back, +z)
    orders as plane waves in the RTT frame; amplitudes of all incident waves superpose linearly.
    """
    p = incident.port
    if p.direction != 1 or p.z != 0.0 or complex(p.index) != handle.front_index:
        raise ValueError("incident waves must be on the front port with direction +1")
    sim = _simulation(handle, p.wavelength_um)
    cfg = handle.structure.config
    lam = p.wavelength_um
    nu = handle.front_index.real
    n_sub = handle.back_index.real
    period_x = float(cfg["period_um"]["x"])
    period_z = float(cfg["period_um"].get("z", 0.0) or 0.0)
    k0 = p.k0
    outputs: dict[str, dict[tuple[float, float], ComplexArray]] = {"front": {}, "back": {}}
    order_of: dict[str, dict[tuple[float, float], tuple[int, int]]] = {"front": {}, "back": {}}
    ez = incident.ez()
    for j in range(len(incident.kx)):
        kx, ky = float(incident.kx[j]), float(incident.ky[j])
        sin_phi = np.hypot(kx, ky) / (nu * k0)
        if sin_phi >= 1.0:
            raise ValueError("incident wave must propagate in the superstrate")
        phi = float(np.arcsin(sin_phi))
        psi = float(np.arctan2(ky, kx)) if sin_phi > 0 else 0.0
        e_rtt = np.array([incident.field[j, 0], incident.field[j, 1], ez[j]])
        e_fg = np.array([e_rtt[0], -e_rtt[2], e_rtt[1]])  # inverse of _to_rtt
        s_in = np.array([np.sin(psi), 0.0, -np.cos(psi)])
        p_in = np.array([np.cos(phi) * np.cos(psi), np.sin(phi), np.cos(phi) * np.sin(psi)])
        weights = (complex(e_fg @ s_in), complex(e_fg @ p_in))
        alpha0, gamma0 = nu * np.sin(phi) * np.cos(psi), nu * np.sin(phi) * np.sin(psi)
        for weight, delta in zip(weights, (0.0, 90.0), strict=True):
            if weight == 0:
                continue
            res = sim.solve(
                lambda_um=lam, phi_deg=np.degrees(phi), psi_deg=np.degrees(psi), delta_deg=delta
            )
            for side, table, medium, up in (
                ("front", res.reflected, nu, True),
                ("back", res.transmitted, n_sub, False),
            ):
                for row in table:
                    n, m = int(row["n"]), int(row["m"])
                    alpha = alpha0 + n * lam / period_x
                    gamma = gamma0 + (m * lam / period_z if period_z > 0 else 0.0)
                    beta = float(np.sqrt(max(medium * medium - alpha * alpha - gamma * gamma, 0.0)))
                    s_vec, p_vec = _order_basis(alpha, beta, gamma, up)
                    e = _to_rtt(complex(row["amp_s"]) * s_vec + complex(row["amp_p"]) * p_vec)
                    key = (round(alpha * k0, 9), round(gamma * k0, 9))
                    outputs[side][key] = outputs[side].get(key, np.zeros(3, complex)) + weight * e
                    order_of[side][key] = (n, m)
    result: dict[str, PlaneWaves] = {}
    for side, port in (("front", front_port(handle, lam, -1)), ("back", back_port(handle, lam))):
        keys = sorted(outputs[side])
        kxs = np.array([k[0] for k in keys])
        kys = np.array([k[1] for k in keys])
        fields = np.array([outputs[side][k][:2] for k in keys]).reshape(len(keys), 2)
        orders = np.array([order_of[side][k] for k in keys], dtype=np.int64).reshape(len(keys), 2)
        result[side] = PlaneWaves(port, kxs, kys, fields.astype(np.complex128), orders)
    return result


def apply_many(
    handle: Handle, incidents: list[PlaneWaves], bloch: bool = False
) -> list[dict[str, PlaneWaves]]:
    """`apply` for many incidences (L4, stage O3). With FourierGrating's S-matrix in the s/p basis
    of the orders (`Simulation.smatrix`, FG Plan eq:Ssp) one eigen-solve per wavelength and
    direction serves both polarisations (`apply`: one solve each); results equal `apply` to
    rounding.

    `bloch`: one eigen-solve per Bloch family (α₀ + n λ/Λ) – a wave of the family is the column of
    its order. FourierGrating truncates the orders around the incident one, so an off-centre column
    sees a shifted truncation: the result differs from `apply` by the truncation error (about 7e−3
    at N = 6 in the tests), not by rounding. Without `smatrix` (an older FourierGrating) or for an
    order outside the truncation: `apply`."""
    families: dict[tuple[int, float, float], Any] = {}
    return [_apply_by_smatrix(handle, inc, families, bloch) for inc in incidents]


def _apply_by_smatrix(
    handle: Handle,
    incident: PlaneWaves,
    families: dict[tuple[int, float, float], Any],
    bloch: bool,
) -> dict[str, PlaneWaves]:
    p = incident.port
    if p.direction != 1 or p.z != 0.0 or complex(p.index) != handle.front_index:
        raise ValueError("incident waves must be on the front port with direction +1")
    sim = _simulation(handle, p.wavelength_um)
    if not hasattr(sim, "smatrix"):
        return apply(handle, incident)
    cfg = handle.structure.config
    lam = p.wavelength_um
    nu = handle.front_index.real
    n_sub = handle.back_index
    period_x = float(cfg["period_um"]["x"])
    period_z = float(cfg["period_um"].get("z", 0.0) or 0.0)
    step_x, step_z = lam / period_x, (lam / period_z if period_z > 0 else 0.0)
    k0 = p.k0
    wl_key = int(np.argmin([abs(w - lam) for w in handle.wavelengths_um]))
    outputs: dict[str, dict[tuple[float, float], ComplexArray]] = {"front": {}, "back": {}}
    order_of: dict[str, dict[tuple[float, float], tuple[int, int]]] = {"front": {}, "back": {}}
    ez = incident.ez()
    for j in range(len(incident.kx)):
        kx, ky = float(incident.kx[j]), float(incident.ky[j])
        sin_phi = np.hypot(kx, ky) / (nu * k0)
        if sin_phi >= 1.0:
            raise ValueError("incident wave must propagate in the superstrate")
        phi = float(np.arcsin(sin_phi))
        psi = float(np.arctan2(ky, kx)) if sin_phi > 0 else 0.0
        alpha0, gamma0 = nu * np.sin(phi) * np.cos(psi), nu * np.sin(phi) * np.sin(psi)
        if bloch:  # the Bloch family: α₀ modulo λ/Λx (and γ₀ modulo λ/Λz for 2D gratings)
            fx = (alpha0 / step_x) % 1.0
            fz = (gamma0 / step_z) % 1.0 if step_z > 0 else gamma0
            key = (wl_key, round(fx % (1.0 - 1e-10), 9), round(fz, 9))
        else:  # the direction itself: the incident order is the centre of the truncation
            key = (wl_key, round(alpha0, 12), round(gamma0, 12))
        if key not in families:
            families[key] = (
                sim.smatrix(lambda_um=lam, phi_deg=np.degrees(phi), psi_deg=np.degrees(psi)),
                alpha0,
                gamma0,
            )
        smat, alpha_r, gamma_r = families[key]
        n = round((alpha0 - alpha_r) / step_x)
        m = round((gamma0 - gamma_r) / step_z) if step_z > 0 else 0
        try:
            r = smat.index(n, m)
        except KeyError:
            r = -1
        if (
            r < 0
            or abs(complex(smat.alpha[r]) - alpha0) > 1e-9
            or abs(complex(smat.gamma[r]) - gamma0) > 1e-9
        ):
            single = PlaneWaves(
                p, incident.kx[j : j + 1], incident.ky[j : j + 1], incident.field[j : j + 1]
            )
            for side, waves in apply(handle, single).items():
                assert waves.orders is not None
                for i in range(len(waves.kx)):
                    e3 = np.array([waves.field[i, 0], waves.field[i, 1], 0.0], dtype=complex)
                    kxy = (float(waves.kx[i]), float(waves.ky[i]))
                    outputs[side][kxy] = outputs[side].get(kxy, np.zeros(3, complex)) + e3
                    order_of[side][kxy] = (int(waves.orders[i][0]), int(waves.orders[i][1]))
            continue
        nf = len(smat.orders)
        e_rtt = np.array([incident.field[j, 0], incident.field[j, 1], ez[j]])
        e_fg = np.array([e_rtt[0], -e_rtt[2], e_rtt[1]])  # inverse of _to_rtt
        beta0 = float(np.sqrt(max(nu * nu - alpha0 * alpha0 - gamma0 * gamma0, 0.0)))
        s_in, p_in = _order_basis(alpha0, beta0, gamma0, False)
        a = np.zeros(2 * nf, dtype=complex)
        a[r], a[nf + r] = complex(e_fg @ s_in), complex(e_fg @ p_in)
        lossless_sub = abs(complex(n_sub).imag) < 1e-14
        for side, amps, betas, up in (
            ("front", smat.R_top @ a, smat.beta_sup, True),
            ("back", smat.T_down @ a, smat.beta_sub, False),
        ):
            if side == "back" and not lossless_sub:
                continue  # FourierGrating lists no transmitted orders then (as apply)
            for q in range(nf):
                beta = complex(betas[q])
                if abs(beta.imag) >= 1e-12 or beta.real <= 0:
                    continue  # evanescent
                alpha, gamma = (
                    float(complex(smat.alpha[q]).real),
                    float(complex(smat.gamma[q]).real),
                )
                s_vec, p_vec = _order_basis(alpha, beta.real, gamma, up)
                e = _to_rtt(amps[q] * s_vec + amps[nf + q] * p_vec)
                kxy = (round(alpha * k0, 9), round(gamma * k0, 9))
                outputs[side][kxy] = outputs[side].get(kxy, np.zeros(3, complex)) + e
                # orders counted from the wave's own incidence, as apply does
                order_of[side][kxy] = (int(smat.orders[q][0]) - n, int(smat.orders[q][1]) - m)
    result: dict[str, PlaneWaves] = {}
    for side, port in (("front", front_port(handle, lam, -1)), ("back", back_port(handle, lam))):
        keys = sorted(outputs[side])
        kxs = np.array([k[0] for k in keys])
        kys = np.array([k[1] for k in keys])
        fields = np.array([outputs[side][k][:2] for k in keys]).reshape(len(keys), 2)
        orders = np.array([order_of[side][k] for k in keys], dtype=np.int64).reshape(len(keys), 2)
        result[side] = PlaneWaves(port, kxs, kys, fields.astype(np.complex128), orders)
    return result


def estimate(handle: Handle) -> dict[str, float]:
    """Rough cost: FG reports no model yet; matrix size (2 n_f)³ per slice as a proxy."""
    cfg = handle.structure.config
    nx = int(cfg.get("truncation", {}).get("Nx", 0))
    nz = int(cfg.get("truncation", {}).get("Nz", 0))
    n_f = (2 * nx + 1) * (2 * nz + 1)
    slices = max(len(cfg["layers"]), 1) * int(cfg.get("rcwa", {}).get("slices", 1) or 1)
    return {"seconds": 1e-9 * (4 * n_f) ** 3 * slices, "memory_bytes": 16.0 * (4 * n_f) ** 2 * 8}
