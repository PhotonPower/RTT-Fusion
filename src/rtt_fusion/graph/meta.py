"""Meta blocks (stage R20, concept section 7): repetition around an acyclic working graph.

The core graph stays acyclic; a meta block repeats the part of the graph in front of one output
(the "subgraph" it encloses) with changed parameters and collects the results:

- `meta.sweep`: a parameter of a block over a list of values (wavelength, tilt, defocus, path
  difference, …),
- `meta.replicates`: the same over replicate values (seeds, realisations) with the mean of the
  collected images as result (late averaging, concept section 5),
- `meta.optimise`: one parameter within bounds for the extremum of a measure of the collected
  result (bounded scalar search; adjoint gradients of solvers come later).

The parameter is addressed as "block.param" or "block.param[i]" (an entry of a list parameter,
e.g. "mirror.point[0]"); the collected output as "block.socket". The executor's cache makes a
repetition recompute only what lies behind the changed block (concept section 6). Meta blocks are
not cached themselves: their result depends on the enclosed graph, not on links.

Measures of a collected image (`imaging.Image`): "peak", "centre" (mean of the 3×3 central
pixels), "power"; of a number: the number.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from rtt_fusion._types import FloatArray
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.imaging import Image
from rtt_fusion.graph.model import Category, Node, NodeType, ParamSpec, data_out
from rtt_fusion.graph.registry import NodeTypeRegistry

MEASURES = ("peak", "centre", "power")
_ADDRESS = re.compile(r"^(?P<node>[^.\[\]]+)\.(?P<param>[^.\[\]]+)(\[(?P<index>\d+)\])?$")


@dataclass
class MetaResult:
    """Values of the parameter, the collected results and their measure; `best` (optimisation)
    and `mean` (replicates) where they apply."""

    parameter: str
    values: list[float] = field(default_factory=list)
    results: list[Any] = field(default_factory=list)
    measures: list[float] = field(default_factory=list)
    measure: str = "peak"
    mean: Any = None
    best: float | None = None


def measure_of(result: Any, measure: str) -> float:
    if isinstance(result, Image):
        v = np.asarray(result.values)
        if measure == "peak":
            return float(v.max()) if v.size else 0.0
        if measure == "power":
            return result.power
        if measure == "centre":
            ny, nx = v.shape
            return float(v[ny // 2 - 1 : ny // 2 + 2, nx // 2 - 1 : nx // 2 + 2].mean())
        raise ExecutionError("exec.param", f"unknown measure {measure!r} ({', '.join(MEASURES)})")
    if isinstance(result, int | float | np.floating):
        return float(result)
    raise ExecutionError("exec.meta", f"no measure for a {type(result).__name__}")


class _Parameter:
    """A settable parameter "block.param" or "block.param[i]" of the graph."""

    def __init__(self, ctx: Context, address: str, where: str) -> None:
        m = _ADDRESS.match(address.strip())
        if m is None:
            raise ExecutionError("exec.param", f"parameter {address!r}: use block.param[i]", where)
        self.node_id, self.name = m["node"], m["param"]
        self.index = int(m["index"]) if m["index"] is not None else None
        if self.node_id not in ctx.graph.nodes:
            raise ExecutionError("exec.param", f"no block {self.node_id!r}", where)
        self.node = ctx.graph.nodes[self.node_id]
        self.saved = self.node.params.get(self.name)

    def set(self, value: float) -> None:
        if self.index is None:
            self.node.params[self.name] = value
            return
        current = list(self.node.params.get(self.name) or [])
        while len(current) <= self.index:
            current.append(0.0)
        current[self.index] = value
        self.node.params[self.name] = current

    def restore(self) -> None:
        if self.saved is None:
            self.node.params.pop(self.name, None)
        else:
            self.node.params[self.name] = self.saved


def _collect(ctx: Context, node: Node) -> tuple[str, str]:
    target = str(node.params.get("collect", ""))
    if "." not in target:
        raise ExecutionError("exec.param", "collect: give the output as block.socket", node.id)
    block, socket = target.split(".", 1)
    if block not in ctx.graph.nodes:
        raise ExecutionError("exec.param", f"collect: no block {block!r}", node.id)
    return block, socket


def _run(ctx: Context, node: Node, values: list[float]) -> MetaResult:
    if ctx.evaluate is None:
        raise ExecutionError("exec.meta", "meta blocks need the executor", node.id)
    block, socket = _collect(ctx, node)
    parameter = _Parameter(ctx, str(node.params.get("parameter", "")), node.id)
    measure = str(node.params.get("measure", "peak"))
    out = MetaResult(str(node.params.get("parameter", "")), measure=measure)
    try:
        for v in values:
            parameter.set(float(v))
            result = ctx.evaluate(block, socket)
            out.values.append(float(v))
            out.results.append(result)
            out.measures.append(measure_of(result, measure))
    finally:
        parameter.restore()
    return out


def sweep(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> MetaResult:
    values = [float(v) for v in node.params.get("values", [])]
    if not values:
        start, stop = float(node.params.get("start", 0.0)), float(node.params.get("stop", 1.0))
        count = int(node.params.get("count", 11))
        values = [float(v) for v in np.linspace(start, stop, max(count, 1))]
    return _run(ctx, node, values)


def replicates(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> MetaResult:
    values = [float(v) for v in node.params.get("values", [])] or [
        float(k) for k in range(int(node.params.get("count", 8)))
    ]
    out = _run(ctx, node, values)
    images = [r for r in out.results if isinstance(r, Image)]
    if images and len(images) == len(out.results):
        mean = np.mean([np.asarray(i.values) for i in images], axis=0)
        groups = [g for i in images for g in i.groups]
        out.mean = Image(
            np.asarray(mean, dtype=float),
            images[0].extent,
            float(np.mean([i.ray_power for i in images])),
            list(dict.fromkeys(groups)),
        )
    elif out.measures:
        out.mean = float(np.mean(out.measures))
    return out


def optimise(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> MetaResult:
    from scipy.optimize import minimize_scalar

    lo, hi = (float(v) for v in node.params.get("bounds", [0.0, 1.0]))
    if not lo < hi:
        raise ExecutionError("exec.param", "bounds: lower < upper", node.id)
    sign = -1.0 if str(node.params.get("goal", "maximum")) == "maximum" else 1.0
    tolerance = float(node.params.get("tolerance", 1e-6 * (hi - lo)))
    trace = MetaResult(
        str(node.params.get("parameter", "")), measure=str(node.params.get("measure", "peak"))
    )

    def merit(x: float) -> float:
        step = _run(ctx, node, [x])
        trace.values += step.values
        trace.results += step.results
        trace.measures += step.measures
        return sign * step.measures[0]

    # a coarse scan first, so that the bounded search starts in the right basin
    grid: FloatArray = np.linspace(lo, hi, int(node.params.get("scan", 9)))
    scores = [merit(float(x)) for x in grid]
    k = int(np.argmin(scores))
    a, b = float(grid[max(k - 1, 0)]), float(grid[min(k + 1, len(grid) - 1)])
    res = minimize_scalar(merit, bounds=(a, b), method="bounded", options={"xatol": tolerance})
    trace.best = float(res.x)
    return trace


def _param(name: str, type_: str, default: Any, description: str = "", **kw: Any) -> ParamSpec:
    return ParamSpec(name, type_, default, description=description, **kw)


_COMMON = (
    _param("parameter", "text", "", "block.param or block.param[i], e.g. mirror.point[0]"),
    _param("collect", "text", "", "output to collect: block.socket, e.g. camera.image"),
    _param("measure", "choice", "peak", "measure of each result", options=MEASURES),
)

TYPES = (
    NodeType(
        "meta.sweep",
        Category.META,
        (data_out("results", "result", "results"),),
        (
            *_COMMON,
            _param("values", "numbers", [], "values (empty: start, stop, count)"),
            _param("start", "number", 0.0),
            _param("stop", "number", 1.0),
            _param("count", "integer", 11, minimum=1),
        ),
        description="Sweep: a parameter over values, the output collected per value",
    ),
    NodeType(
        "meta.replicates",
        Category.META,
        (data_out("results", "result", "results"),),
        (
            *_COMMON,
            _param("values", "numbers", [], "replicate values (empty: 0 … count − 1)"),
            _param("count", "integer", 8, minimum=1),
        ),
        description="Replicates: repeated with replicate values, images averaged",
    ),
    NodeType(
        "meta.optimise",
        Category.META,
        (data_out("results", "result", "results"),),
        (
            *_COMMON,
            _param("bounds", "numbers", [0.0, 1.0], "lower and upper bound"),
            _param("goal", "choice", "maximum", options=("maximum", "minimum")),
            _param("scan", "integer", 9, minimum=3, description="coarse scan before the search"),
            _param("tolerance", "number", 0.0, description="0: 1e−6 of the bounds"),
        ),
        description="Optimisation of one parameter for the extremum of a measure",
    ),
)

EVALUATORS: dict[str, Evaluator] = {
    "meta.sweep": sweep,
    "meta.replicates": replicates,
    "meta.optimise": optimise,
}
for _ev in EVALUATORS.values():
    _ev.uncached = True  # type: ignore[attr-defined]


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    for t in TYPES:
        registry.add(t)
    return dict(EVALUATORS)
