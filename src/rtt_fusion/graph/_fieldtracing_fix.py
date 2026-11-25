"""A fix for Fieldtracing 0.2.0, installed before the planner runs (issue #135, case 1).

`fieldtracing.selection.fit_curvature` solves its minimax linear program on a subset of the
samples and adds the violated samples until none is left. HiGHS solves to about 1e-7, the loop
asks for 1e-9: when every violated sample is in the subset already, the subset no longer grows
and the loop solves the same program forever (the camera of the FrObSy exposure head with the
head shifted by 0.25 mm: a few hundred samples, the run never ended).

The replacement is the same algorithm and stops when the subset no longer grows, with the largest
residual over all samples as the bound (≥ the optimum, so the grid step it sets is safe). It is
installed only for Fieldtracing 0.2.0 (0.2.1 carries the fix, PhotonPower/Fieldtracing#86) and
only if the module has the parts it builds on.
"""

from __future__ import annotations

from importlib import metadata
from typing import Any

import numpy as np

_installed = False


def _version() -> tuple[int, ...]:
    try:
        text = metadata.version("photonpower-fieldtracing")
    except metadata.PackageNotFoundError:
        return ()
    out: list[int] = []
    for part in text.split(".")[:3]:
        digits = "".join(c for c in part if c.isdigit())
        out.append(int(digits or 0))
    return tuple(out)


def install() -> None:
    """Installs the terminating `fit_curvature` into Fieldtracing 0.2.0 (once)."""
    global _installed
    if _installed:
        return
    _installed = True
    version = _version()
    if version[:3] != (0, 2, 0):  # fixed in Fieldtracing 0.2.1 (PhotonPower/Fieldtracing#86)
        return
    try:
        import fieldtracing.chain as chain
        import fieldtracing.selection as selection
    except ImportError:
        return
    needed = ("_curvature_lp", "_curvature_residual", "CURVATURE_SUBSET", "fit_curvature")
    if not all(hasattr(selection, n) for n in needed) or not hasattr(chain, "fit_curvature"):
        return

    def fit_curvature(points: Any, gradient: Any) -> tuple[Any, float]:
        p = np.asarray(points, dtype=np.float64)
        g = np.asarray(gradient, dtype=np.float64)
        n = len(p)
        size = int(selection.CURVATURE_SUBSET)
        if n <= size:
            c0, t0 = selection._curvature_lp(p, g)
            return c0, float(t0)
        k = size // 4
        picks = [
            np.argsort(np.hypot(p[:, 0], p[:, 1]))[-k:],
            np.argsort(np.max(np.abs(g), axis=1))[-k:],
            np.linspace(0, n - 1, 2 * k).astype(np.int64),
        ]
        subset = np.unique(np.concatenate(picks))
        while True:
            c, t = selection._curvature_lp(p[subset], g[subset])
            residual = selection._curvature_residual(p, g, c)
            tolerance = 1e-9 * max(t, float(np.max(np.abs(g))), 1e-300)
            violated = np.flatnonzero(residual > t + tolerance)
            grown = np.union1d(subset, violated[np.argsort(residual[violated])[-size:]])
            if violated.size == 0 or len(grown) == len(subset):
                return c, float(max(t, float(np.max(residual))))
            subset = grown

    fit_curvature.__doc__ = selection.fit_curvature.__doc__
    selection.fit_curvature = fit_curvature
    chain.fit_curvature = fit_curvature  # chain imported the name


__all__ = ["install"]
