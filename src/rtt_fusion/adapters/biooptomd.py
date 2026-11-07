"""Adapter for BioOptoMD data (stage R10): refractive-index profiles n(z) and transition dipoles.

BioOptoMD is no light solver: it runs molecular dynamics of receptors on sensor surfaces (GPU)
and its optics observers produce the inputs of light solvers (Fusion plan: data ports N and Q).
It needs CUDA; this adapter therefore reads its **output formats** and needs no import of
BioOptoMD:

- n(z) (`RefractiveIndexProfiler`, BioOptoMD docs/PHYSICS_MODEL.md §6.3, Lorentz–Lorenz with the
  solvent as background): CSV of `biooptomd.analysis.export_refractive_index` – columns
  `z_nm, n_replica_0, …`, z = bin centres above the surface, bin width Δz; real n at the
  wavelength of the polarisabilities of the run;
- dipoles (`DipoleAndDistanceTracker`, §6.1): time series (n_replicas, n_frames, 7) =
  (R_F [nm], μ̂, d [nm]).

BioOptoMD conventions (docs/UNITS_AND_CONVENTIONS.md, commit `4a4c492`): lengths in nm, sensor
surface at z = 0, solution at z > 0, surface normal +ẑ. Here: lengths in mm, light along +z of
the RTT frame; `light_from` says from which side light meets the layer ("substrate": the layer
lies after the surface along +z, as in a Kretschmann prism; "solution": before it, mirrored).

The merging rule of `profile_layers` is the one of `biooptomd.analysis.refractive_index_layers`
(neighbouring bins with |Δn| ≤ tol merged, thickness-weighted; trailing solvent bins removed),
reimplemented so that this module needs no CUDA build.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from rtt_fusion._types import FloatArray

MM_PER_NM = 1e-6


@dataclass(frozen=True)
class Profile:
    """n(z) of one run: bin centres z_nm (n_bins,), n (n_replicas, n_bins), bin width dz_nm,
    solvent index, wavelength of the polarisabilities and a free-text source."""

    z_nm: FloatArray
    n: FloatArray
    dz_nm: float
    n_solvent: float
    wavelength_um: float
    source: str = ""

    def __post_init__(self) -> None:
        n = np.atleast_2d(np.asarray(self.n, dtype=np.float64))
        if n.shape[-1] != np.asarray(self.z_nm).shape[0]:
            raise ValueError("n must have shape (n_replicas, n_bins) matching z_nm")
        if not np.isfinite(n).all():
            raise ValueError("n contains NaN (L >= 1 in BioOptoMD: unphysical parameters)")
        if np.any(np.abs(np.diff(self.z_nm) - self.dz_nm) > 1e-9 * self.dz_nm):
            raise ValueError("bins must be uniform with width dz_nm")


def describe() -> dict[str, Any]:
    return {
        "name": "BioOptoMD",
        "version": "0.0.1.dev0 (output formats of commit 4a4c492)",
        "adapter_version": "0.1.0",
        "method": "molecular dynamics with optics observers (no light solver)",
        "representation": "N (index profile), Q (dipoles): data ports",
        "ports": [],
        "linear": None,
        "reciprocal": None,
        "operator": False,
        "deterministic": {"data": "files are fixed results; replicas give the spread"},
        "conventions_of_solver": [
            "nm; surface at z = 0, solution at z > 0, normal +z",
            "n real at the wavelength of the run (Lorentz-Lorenz, solvent background)",
        ],
        "validity": [
            "n(z) only at the wavelength of the run (dispersion via the run's Cauchy term)",
            "layer model: laterally averaged over the box (Lx, Ly); no lateral structure",
        ],
    }


def read_profile_csv(
    path: str | Path, dz_nm: float, n_solvent: float, wavelength_um: float
) -> Profile:
    """Reads the CSV of `biooptomd.analysis.export_refractive_index`."""
    data = np.loadtxt(Path(path), delimiter=",", comments="#", ndmin=2)
    return Profile(data[:, 0], data[:, 1:].T, dz_nm, n_solvent, wavelength_um, str(path))


def profile_layers(
    profile: Profile, replica: int | str = "mean", tol: float = 1e-9
) -> tuple[FloatArray, FloatArray]:
    """n(z) → layers (thickness [nm], n) from the surface outward (BioOptoMD's rule).

    replica: index, or "mean" for the mean profile over the replicas.
    """
    n_all = np.atleast_2d(np.asarray(profile.n, dtype=np.float64))
    nn = n_all.mean(axis=0) if replica == "mean" else n_all[int(replica)]
    thick, vals = [profile.dz_nm], [float(nn[0])]
    for v in nn[1:]:
        if abs(v - vals[-1]) <= tol:
            vals[-1] = (vals[-1] * thick[-1] + v * profile.dz_nm) / (thick[-1] + profile.dz_nm)
            thick[-1] += profile.dz_nm
        else:
            thick.append(profile.dz_nm)
            vals.append(float(v))
    t, n = np.array(thick), np.array(vals)
    while len(n) and abs(n[-1] - profile.n_solvent) <= max(tol, 1e-7):
        t, n = t[:-1], n[:-1]
    return t, n


def stack_along_z(
    thickness_nm: FloatArray, n: FloatArray, light_from: str = "substrate"
) -> list[tuple[float, float]]:
    """Layers in the order light meets them along +z of the RTT frame, thickness in mm."""
    if light_from not in ("substrate", "solution"):
        raise ValueError("light_from must be 'substrate' or 'solution'")
    layers = [(float(t) * MM_PER_NM, float(v)) for t, v in zip(thickness_nm, n, strict=True)]
    return layers if light_from == "substrate" else layers[::-1]


def excess_optical_thickness_nm(profile: Profile, replica: int | str = "mean") -> float:
    """∫ (n(z) − n_solvent) dz [nm]: the quantity a thin adlayer shifts a resonance by to first
    order (de Feijter); independent of the layering."""
    n_all = np.atleast_2d(np.asarray(profile.n, dtype=np.float64))
    nn = n_all.mean(axis=0) if replica == "mean" else n_all[int(replica)]
    return float(np.sum(nn - profile.n_solvent) * profile.dz_nm)


@dataclass(frozen=True)
class Dipoles:
    """Emitters in the RTT frame: positions (N, 3) in mm, unit orientations (N, 3)."""

    positions_mm: FloatArray
    orientations: FloatArray


def dipoles_from_tracker(
    series: FloatArray, surface_z_mm: float, light_from: str = "substrate"
) -> Dipoles:
    """(n_replicas, n_frames, 7) = (R_F [nm], μ̂, d [nm]) → positions and orientations in the RTT
    frame; the surface lies at z = surface_z_mm. "substrate": BioOptoMD's +z is our +z;
    "solution": mirrored (z → −z; μ̂_z → −μ̂_z), x and y kept (a rotation by π about x would flip
    y as well; a mirror keeps the handedness of the lateral coordinates of the box)."""
    a = np.asarray(series, dtype=np.float64).reshape(-1, 7)
    sign = 1.0 if light_from == "substrate" else -1.0
    pos = np.column_stack(
        [a[:, 0] * MM_PER_NM, a[:, 1] * MM_PER_NM, surface_z_mm + sign * a[:, 6] * MM_PER_NM]
    )
    mu = a[:, 3:6] * np.array([1.0, 1.0, sign])
    return Dipoles(pos, mu / np.linalg.norm(mu, axis=1, keepdims=True))
