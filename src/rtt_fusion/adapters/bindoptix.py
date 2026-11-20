"""Adapter for BindOptiX (D13 step 3, RTT-Fusion#106): the signal of a fluorescence biosensor from
data ports, no light solver.

BindOptiX (`PhotonPower/BindOptiX`, MIT) is the time, transport and biology layer of the
biosensor pilots (K8). Here it serves as a data node: the excitation field map F (`field_map`),
the emitter response table M (`response_map`), the dipoles Q per receptor state
(`adapters.biooptomd.DipoleSeries`) and the occupancy B (`occupancy`) give the mean detected
photon rate per receptor, spot, time and detector region (`bindoptix.signal.spot_rates_s`).

Formats: F, M and B are files in the `raytatouille-result` envelope with the types `FieldMap`,
`EmitterResponseTable` and `Occupancy` (BindOptiX `docs/adr-entwurf-datenports.md`; hp-FEM writes
F and M). BindOptiX reads them; it is imported only when needed (like `adapters.hpfem`). Without
its C++ core it needs NumPy and SciPy only:
`pip install "bindoptix @ git+https://github.com/PhotonPower/BindOptiX" -C
cmake.define.BINDOPTIX_BUILD_CORE=OFF`.

Frames: everything here stays in BindOptiX's and BioOptoMD's sensor frame (µm resp. nm, surface
z = 0, solution z > 0); Q goes in unconverted (`EmitterEnsemble.from_tracker`, no mirror) – only
nodes that compute light convert into the RTT frame (D13, DP8).

Detector regions (D13): the regions of a response table (e.g. `substrate_na_lt_1.45`) are bound
to detectors of the graph by a parameter "label=region, …"; a region the table lacks, or a label
or region given twice, is an error. A table with exactly one region needs no binding.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from rtt_fusion._types import FloatArray
from rtt_fusion.adapters.biooptomd import DipoleSeries

#: the envelope types of the three data ports
TYPES = {"field_map": "FieldMap", "response_map": "EmitterResponseTable", "occupancy": "Occupancy"}


def _bx() -> Any:
    """The BindOptiX modules used here (signal, ports, binding); ImportError if not installed."""
    import bindoptix.binding
    import bindoptix.ports
    import bindoptix.signal

    return bindoptix


def describe() -> dict[str, Any]:
    return {
        "name": "BindOptiX",
        "version": _version(),
        "adapter_version": "0.1.0",
        "method": "signal model of a fluorescence biosensor from data ports (no light solver)",
        "representation": "field_map (F), response_map (M), Q, occupancy (B) → result",
        "ports": [],
        "linear": None,
        "reciprocal": None,
        "operator": False,
        "deterministic": {"data": "rates are expectations; no sampling here"},
        "conventions_of_solver": [
            "µm, s, mW; sensor frame: surface z = 0, solution z > 0; F with e^(−iωt), per mW",
        ],
        "validity": [
            "linear excitation (no saturation, no triplet blinking); bleaching not applied here",
        ],
    }


def _version() -> str:
    try:
        return str(_bx().__version__)
    except ImportError:
        return "not installed"


def _envelope_type(path: Path) -> str:
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except OSError as e:
        raise ValueError(f"{path}: {e.strerror}") from e
    except ValueError as e:
        raise ValueError(f"{path}: not JSON ({e})") from e
    kind = doc.get("type") if isinstance(doc, dict) else None
    return str(kind) if kind is not None else ""


def read(port: str, path: str | Path) -> Any:
    """A data file of one of the ports `field_map`, `response_map`, `occupancy` as BindOptiX's
    object (`FieldMap`, `EmitterResponseTable`, `Occupancy`). ValueError for a missing file,
    another envelope type or invalid content; ImportError without BindOptiX."""
    want = TYPES[port]
    p = Path(path)
    kind = _envelope_type(p)
    if kind and kind != want:
        raise ValueError(f"{p}: a {kind}, not a {want} (port {port})")
    bx = _bx()
    cls = {
        "field_map": bx.ports.FieldMap,
        "response_map": bx.ports.EmitterResponseTable,
        "occupancy": bx.binding.Occupancy,
    }[port]
    try:
        return cls.from_json(p.read_text(encoding="utf-8"))
    except (KeyError, TypeError) as e:
        raise ValueError(f"{p}: invalid {want} ({e})") from e


def bind_regions(text: str, regions: Sequence[str]) -> list[tuple[str, int]]:
    """(label, index into `regions`) from "label=region, …" (D13). Empty text: the only region
    under its own name; with several regions the binding is required."""
    regions = list(regions)
    if not text.strip():
        if len(regions) == 1:
            return [(regions[0], 0)]
        raise ValueError(
            f"the response table has {len(regions)} detector regions {regions}: bind them, "
            "e.g. 'camera=" + (regions[0] if regions else "region") + "'"
        )
    out: list[tuple[str, int]] = []
    for item in (s.strip() for s in text.split(",")):
        if not item:
            continue
        label, sep, region = (s.strip() for s in item.partition("="))
        if not sep or not label or not region:
            raise ValueError(f"'{item}': want label=region")
        if region not in regions:
            raise ValueError(f"region {region!r} is not in the response table {regions}")
        if any(label == lb for lb, _ in out):
            raise ValueError(f"label {label!r} given twice")
        if any(regions.index(region) == i for _, i in out):
            raise ValueError(f"region {region!r} bound twice")
        out.append((label, regions.index(region)))
    if not out:
        raise ValueError("no region bound")
    return out


@dataclass(frozen=True)
class SignalResult:
    """Mean detected photon rate per receptor: rates[t, spot, label] [1/s] at `times_s`, for the
    spots of the occupancy and the bound detector labels."""

    times_s: FloatArray
    spots: tuple[int, ...]
    labels: tuple[str, ...]
    rates: FloatArray
    source: dict[str, str] = field(default_factory=dict)

    def to_curves(self) -> tuple[list[tuple[str, list[float], list[float]]], str, str]:
        """Curves for `display.plot`: one per spot and label over time."""
        t = [float(x) for x in self.times_s]
        curves = [
            (f"spot {s} – {label}", t, [float(x) for x in self.rates[:, i, j]])
            for i, s in enumerate(self.spots)
            for j, label in enumerate(self.labels)
        ]
        return curves, "time (s)", "detected rate per receptor (1/s)"


def ensembles(q: DipoleSeries, brightness: dict[str, float] | None = None) -> dict[str, Any]:
    """BindOptiX ensembles per state from Q, unconverted (sensor frame; `from_tracker`)."""
    bx = _bx()
    b = brightness or {}
    return {
        s: bx.signal.EmitterEnsemble.from_tracker(q.series[s], brightness=float(b.get(s, 1.0)))
        for s in q.states
    }


def signal(
    excitation: Any,
    table: Any,
    q: DipoleSeries,
    occupancy: Any,
    dye: dict[str, float],
    power_mw: float,
    anchors_um: FloatArray,
    regions: str = "",
    port_irradiance: tuple[float, str, str] | None = None,
) -> SignalResult:
    """`bindoptix.signal.spot_rates_s` for the bound regions. `dye`: epsilon_per_molar_cm,
    quantum_yield, excitation_um, emission_um (bleach_yield optional). `port_irradiance`
    (W/m² per mW, status L/K/P, reference) converts a map normalised to I₀ (as solvers write it)
    to one per mW at the input port."""
    bx = _bx()
    if not excitation.is_port_normalised:
        if port_irradiance is None:
            raise ValueError(
                "the field map is normalised to the incident irradiance I₀: give the irradiance "
                "per mW at the spot (port_irradiance_w_m2_per_mw) with its status"
            )
        excitation = excitation.per_port_power(*port_irradiance)
    bound = bind_regions(regions, table.detectors)
    missing = [s for s in occupancy.states if s not in q.series]
    if missing:
        raise ValueError(f"Q has no ensemble for the occupancy states {missing} (Q: {q.states})")
    anchors = np.atleast_2d(np.asarray(anchors_um, dtype=np.float64))
    if anchors.size == 3 and len(occupancy.spots) > 1:
        anchors = np.tile(anchors.reshape(1, 3), (len(occupancy.spots), 1))
    rates = bx.signal.spot_rates_s(
        occupancy,
        ensembles(q),
        anchors,
        excitation,
        table,
        bx.signal.Fluorophore(**dye),
        float(power_mw),
    )
    idx = [i for _, i in bound]
    src = {"bindoptix": _version(), "q": q.source}
    src.update({f"table_{k}": str(v) for k, v in getattr(table, "source", {}).items()})
    return SignalResult(
        np.asarray(occupancy.times_s, dtype=np.float64),
        tuple(int(s) for s in occupancy.spots),
        tuple(label for label, _ in bound),
        np.asarray(rates[..., idx], dtype=np.float64),
        src,
    )


def parse_regions(text: str) -> list[tuple[str, str, float, float]]:
    """Detector regions for `response_table` from "name=side:na_lo:na_hi, …" (side "substrate"
    or "solution"; the band na_lo ≤ n sin θ < na_hi of that half-space). Bands of one side must
    not overlap (A_rad − Σ A_det ≥ 0, DP6)."""
    out: list[tuple[str, str, float, float]] = []
    for item in (s.strip() for s in text.split(",")):
        if not item:
            continue
        name, sep, rest = (s.strip() for s in item.partition("="))
        parts = [s.strip() for s in rest.split(":")]
        if not sep or not name or len(parts) != 3 or parts[0] not in ("substrate", "solution"):
            raise ValueError(f"'{item}': want name=substrate|solution:na_lo:na_hi")
        if any(name == r[0] for r in out):
            raise ValueError(f"region {name!r} given twice")
        try:
            lo, hi = float(parts[1]), float(parts[2])
        except ValueError as e:
            raise ValueError(f"'{item}': NA bounds must be numbers") from e
        if not 0 <= lo < hi:
            raise ValueError(f"'{item}': want 0 <= na_lo < na_hi")
        for other, side, lo2, hi2 in out:
            if side == parts[0] and lo < hi2 and lo2 < hi:
                raise ValueError(
                    f"regions {other!r} and {name!r} overlap in the {side}: each power counts once"
                )
        out.append((name, parts[0], lo, hi))
    if not out:
        raise ValueError("no detector region")
    return out


def response_tensors(
    stack: Any, wavelength_um: float, heights_um: Any, regions: list[tuple[str, str, float, float]]
) -> tuple[FloatArray, FloatArray, FloatArray, float]:
    """A_tot, A_rad (nz, 3, 3) and A_det (nz, D, 3, 3) of a planar stack (`couplings.
    dipole_layers`, sensor frame, heights in µm) and the index of the emitters' medium. The stack
    is invariant about z, so every tensor is diag(a∥, a∥, a⊥) from an x and a z dipole."""
    from rtt_fusion.couplings import dipole_layers as dl

    z = np.asarray(heights_um, dtype=float).ravel() * 1e-3  # mm
    media = {stack.medium_at(float(v)) for v in z}
    if len(media) != 1:
        raise ValueError("all heights must lie in one medium of the stack")
    n_host = float(np.real(stack.indices[media.pop()]))
    tot = np.zeros((len(z), 3, 3))
    rad = np.zeros((len(z), 3, 3))
    det = np.zeros((len(z), len(regions), 3, 3))
    for i, z0 in enumerate(z):
        for axes, p in (((0, 1), [1.0, 0.0, 0.0]), ((2,), [0.0, 0.0, 1.0])):
            r = dl.rates(stack, wavelength_um, float(z0), p)
            bands = [
                dl.band_power(
                    stack,
                    wavelength_um,
                    float(z0),
                    p,
                    "bottom" if side == "substrate" else "top",
                    lo,
                    hi,
                )
                for _, side, lo, hi in regions
            ]
            for a in axes:
                tot[i, a, a] = r.total
                rad[i, a, a] = r.up + r.down
                for d, value in enumerate(bands):
                    det[i, d, a, a] = value
    return tot, rad, det, n_host


def response_table(
    stack: Any,
    wavelength_um: float,
    heights_um: Any,
    regions: list[tuple[str, str, float, float]],
    source: dict[str, str] | None = None,
) -> Any:
    """BindOptiX's `EmitterResponseTable` (port M) of a planar stack: grid x = y = 0 (laterally
    invariant), z = the heights (µm, sensor frame); ImportError without BindOptiX."""
    bx = _bx()
    tot, rad, det, n_host = response_tensors(stack, wavelength_um, heights_um, regions)
    z = np.asarray(heights_um, dtype=float).ravel()
    grid = bx.ports.RegularGrid((np.zeros(1), np.zeros(1), z))
    return bx.ports.EmitterResponseTable(
        grid,
        float(wavelength_um),
        n_host,
        tuple(r[0] for r in regions),
        tot.reshape(1, 1, len(z), 3, 3),
        rad.reshape(1, 1, len(z), 3, 3),
        det.reshape(1, 1, len(z), len(regions), 3, 3),
        {"solver": "rtt_fusion.couplings.dipole_layers", **(source or {})},
    )
