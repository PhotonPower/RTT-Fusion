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
  (R_F [nm], μ̂, d [nm]);
- a run directory (D13): the manifest `biooptomd-layer.json` naming the profile CSV with its bin
  width, solvent index and wavelength, and one .npy of the tracker series per receptor state
  (`read_layer`, `write_layer`). BioOptoMD has no export of the series yet; `write_layer` is the
  reference for one. Q (`DipoleSeries`) keeps the series in BioOptoMD's frame; only a node that
  computes light converts it (`DipoleSeries.in_rtt_frame`, `dipoles_from_tracker`).

BioOptoMD conventions (docs/UNITS_AND_CONVENTIONS.md, commit `4a4c492`): lengths in nm, sensor
surface at z = 0, solution at z > 0, surface normal +ẑ. Here: lengths in mm, light along +z of
the RTT frame; `light_from` says from which side light meets the layer ("substrate": the layer
lies after the surface along +z, as in a Kretschmann prism; "solution": before it, mirrored).

The merging rule of `profile_layers` is the one of `biooptomd.analysis.refractive_index_layers`
(neighbouring bins with |Δn| ≤ tol merged, thickness-weighted; trailing solvent bins removed),
reimplemented so that this module needs no CUDA build.
"""

from __future__ import annotations

import json
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


# ------------------------------------------------------------------ directory of a run (D13)

#: the manifest of a BioOptoMD run directory (format "biooptomd-layer", version 1)
MANIFEST = "biooptomd-layer.json"
LAYER_FORMAT, LAYER_VERSION = "biooptomd-layer", 1


@dataclass(frozen=True)
class DipoleSeries:
    """Data port Q (D13, Fusion plan 6.2): BioOptoMD's tracker series per receptor state,
    unchanged – series[state] (n_replicas, n_frames, 7) = (R_F [nm], μ̂, d [nm]) in BioOptoMD's
    sensor frame (surface z = 0, solution z > 0), frames `dt_ps` apart. Only a node that computes
    light converts it into the RTT frame (`in_rtt_frame`, the mirror rule DP8 in one place)."""

    series: dict[str, FloatArray]
    dt_ps: float | None = None
    source: str = ""

    def __post_init__(self) -> None:
        if not self.series:
            raise ValueError("no receptor state")
        for state, s in self.series.items():
            a = np.asarray(s)
            if a.ndim != 3 or a.shape[-1] != 7:
                raise ValueError(
                    f"state {state!r}: shape {a.shape}, want (n_replicas, n_frames, 7)"
                )
            if not np.isfinite(a).all():
                raise ValueError(f"state {state!r}: NaN or infinity")
            if np.any(np.linalg.norm(a[..., 3:6], axis=-1) < 1e-12):
                raise ValueError(f"state {state!r}: zero dipole orientation")

    @property
    def states(self) -> list[str]:
        return list(self.series)

    def in_rtt_frame(
        self, state: str, surface_z_mm: float, light_from: str = "substrate"
    ) -> Dipoles:
        """The emitters of one state (all replicas and frames) in the RTT frame."""
        return dipoles_from_tracker(self.series[state], surface_z_mm, light_from)


@dataclass(frozen=True)
class Layer:
    """A BioOptoMD run directory: the index profile (data port N) and/or the dipoles (Q)."""

    profile: Profile | None
    dipoles: DipoleSeries | None
    source: str = ""


def _manifest_path(path: str | Path) -> Path:
    p = Path(path)
    return p if p.suffix == ".json" else p / MANIFEST


def layer_files(path: str | Path) -> list[Path]:
    """The manifest and the files it names (for cache keys); missing files are left out."""
    manifest = _manifest_path(path)
    try:
        doc = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return [manifest]
    base = manifest.parent
    names = [doc.get("profile", {}).get("file")]
    names += list(doc.get("dipoles", {}).get("states", {}).values())
    return [manifest] + [base / n for n in names if isinstance(n, str)]


def read_layer(path: str | Path) -> Layer:
    """Reads a run directory (or its manifest `biooptomd-layer.json`):

    {"format": "biooptomd-layer", "version": 1, "source": "...",
     "profile": {"file": "n_profile.csv", "dz_nm": 0.05, "n_solvent": 1.333,
                 "wavelength_um": 0.633},
     "dipoles": {"dt_ps": 10.0, "states": {"free": "dipoles_free.npy", ...}}}

    The profile is the CSV of `export_refractive_index`, each state an .npy of the tracker series
    (n_replicas, n_frames, 7) in nm. Either part may be missing, not both."""
    manifest = _manifest_path(path)
    try:
        doc = json.loads(manifest.read_text(encoding="utf-8"))
    except OSError as e:
        raise ValueError(f"no BioOptoMD manifest {manifest}: {e.strerror}") from e
    except ValueError as e:
        raise ValueError(f"{manifest}: not JSON ({e})") from e
    if doc.get("format") != LAYER_FORMAT:
        raise ValueError(f"{manifest}: format {doc.get('format')!r}, want {LAYER_FORMAT!r}")
    if doc.get("version") != LAYER_VERSION:
        raise ValueError(f"{manifest}: version {doc.get('version')!r}, want {LAYER_VERSION}")
    base, source = manifest.parent, str(doc.get("source", ""))
    profile = dipoles = None
    if "profile" in doc:
        p = doc["profile"]
        try:
            profile = read_profile_csv(
                base / p["file"],
                float(p["dz_nm"]),
                float(p["n_solvent"]),
                float(p["wavelength_um"]),
            )
        except KeyError as e:
            raise ValueError(f"{manifest}: profile without {e.args[0]!r}") from e
        except OSError as e:
            raise ValueError(f"{manifest}: profile {p['file']}: {e.strerror}") from e
    if "dipoles" in doc:
        q = doc["dipoles"]
        series = {}
        for state, name in dict(q.get("states", {})).items():
            try:
                series[str(state)] = np.load(base / name, allow_pickle=False).astype(np.float64)
            except OSError as e:
                raise ValueError(f"{manifest}: dipoles of {state!r} ({name}): {e}") from e
        dt = q.get("dt_ps")
        dipoles = DipoleSeries(series, None if dt is None else float(dt), source)
    if profile is None and dipoles is None:
        raise ValueError(f"{manifest}: neither 'profile' nor 'dipoles'")
    return Layer(profile, dipoles, source)


def write_layer(
    directory: str | Path,
    profile: Profile | None = None,
    dipoles: DipoleSeries | None = None,
    source: str = "",
) -> Path:
    """Writes a run directory that `read_layer` reads (the reference for an export in BioOptoMD);
    returns the manifest."""
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    doc: dict[str, Any] = {"format": LAYER_FORMAT, "version": LAYER_VERSION, "source": source}
    if profile is not None:
        n = np.atleast_2d(np.asarray(profile.n))
        np.savetxt(
            d / "n_profile.csv",
            np.column_stack([profile.z_nm, n.T]),
            delimiter=",",
            header="z_nm," + ",".join(f"n_replica_{i}" for i in range(n.shape[0])),
            comments="# ",
        )
        doc["profile"] = {
            "file": "n_profile.csv",
            "dz_nm": profile.dz_nm,
            "n_solvent": profile.n_solvent,
            "wavelength_um": profile.wavelength_um,
        }
    if dipoles is not None:
        states = {}
        for i, (state, s) in enumerate(dipoles.series.items()):
            name = f"dipoles_{i}.npy"
            np.save(d / name, np.asarray(s, dtype=np.float64), allow_pickle=False)
            states[state] = name
        doc["dipoles"] = {"dt_ps": dipoles.dt_ps, "states": states}
    manifest = d / MANIFEST
    manifest.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    return manifest
