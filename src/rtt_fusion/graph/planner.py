"""Planner, stages 1 and 2 (docs/fusionsplan.md 7.2; concept node workflow, stage R17).

The planner is a pure function of the working graph and the requirements to a plan as a data
object (`Plan`):

1. **Topology.** All fusion paths from the sources to the requested outputs, on the level of
   sockets: links between blocks, and inside a block the dependencies output -> inputs
   (`through`). Data links are dependencies, not path segments (concept 3.3).
2. **Choice of representation.** For every light link a representation, and for every block
   whose input and output representations differ a conversion operator (table 6.3; chains such
   as G -> S -> R are allowed). Inadmissible choices drop out: representations a socket does not
   accept, representations not available to the executor (`Requirements.available`) and those a
   validity hook rejects (the metrics of 7.2: Fresnel numbers, size parameter, …). A fixed
   transport is kept; "auto" is chosen by cost, but never P (no phase) while a coherent
   representation is admissible (7.3).

Costs are a provisional *ranking* (`REPRESENTATION_COST`, operator costs), not measured values;
stage 3 replaces them with the solvers' `estimate`. With only R available (stage R15) every "auto"
link becomes R and no conversion is planned, so execution gives the same results as before.

Method of stage 2: a forward dynamic programme over the links in topological order computes, per
link and representation, the cheapest way to deliver that representation (upstream cost); a
backward pass then fixes each link's representation from its upstream cost plus the conversions
into the representations already chosen downstream. On trees (every socket feeds one link and
every input is used by one output) this is exact; where links are shared (a beam splitter whose
input feeds two outputs, fan-out) the backward pass sums the demands of all consumers.

Diagnostic codes: plan.representation_unavailable, plan.no_representation, plan.no_operator
(errors); plan.incoherent_merge, plan.paths_truncated, plan.no_output (warnings).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

from rtt_fusion.graph.model import (
    LIGHT_REPRESENTATIONS,
    Category,
    Direction,
    Graph,
    Link,
    SocketKind,
    SocketSpec,
)
from rtt_fusion.graph.registry import NodeTypeRegistry
from rtt_fusion.graph.validate import Diagnostic, node_sockets, node_through

#: Order of preference for ties and the provisional cost per link (ranking, not measured; R and P
#: are rays, H is gridless, S/G/T need grids, orders or tables). Stage 3 uses `estimate`.
PREFERENCE = ("R", "H", "S", "G", "T", "P")
REPRESENTATION_COST: dict[str, float] = {"R": 1.0, "P": 1.5, "H": 3.0, "S": 4.0, "G": 6.0, "T": 8.0}


@dataclass(frozen=True)
class Operator:
    """A conversion between two representations (docs/fusionsplan.md, table 6.3)."""

    name: str
    source: str
    target: str
    cost: float
    validity: str  # the metric that decides admissibility (checked by the validity hook)
    coherent: bool = True  # False: the result carries no phase (R -> P)


#: Table 6.3. T (T-matrix, spherical waves) has no conversion yet (concept 6.2: later in M).
OPERATORS: tuple[Operator, ...] = (
    Operator("rays_to_field", "R", "H", 1.0, "ray net dense enough, no caustic at the port"),
    Operator("field_to_rays", "H", "R", 1.0, "local Fresnel number >> 1"),
    Operator("resample_residual", "H", "G", 1.0, "sampling rules (Fieldtracing ADR 0003)"),
    Operator("gridded_to_gridless", "G", "H", 1.0, "sampling rules (Fieldtracing ADR 0003)"),
    Operator("fft", "G", "S", 2.0, "table 7.1 of the paper"),
    Operator("inverse_fft", "S", "G", 2.0, "table 7.1 of the paper"),
    Operator("pft", "H", "S", 2.0, "spectrum fully captured (aperture, evanescence)"),
    Operator("orders_to_rays", "S", "R", 1.0, "orders separated, bundles large against λ"),
    Operator("drop_phase", "R", "P", 0.5, "only incoherent further computation", coherent=False),
    Operator("radiometric_rays", "P", "R", 0.5, "rays without phase", coherent=False),
)


def _operator_chains() -> dict[tuple[str, str], tuple[float, tuple[Operator, ...]]]:
    """Cheapest chain of operators between every pair of representations.

    A chain may end with an operator that drops the phase (R -> P, P -> R) but not continue after
    it: rays without phase give no coherent field again (P -> R -> H is no admissible chain)."""
    reps = sorted(LIGHT_REPRESENTATIONS)
    best: dict[tuple[str, str], tuple[float, tuple[Operator, ...]]] = {
        (a, b): ((0.0, ()) if a == b else (math.inf, ())) for a in reps for b in reps
    }

    def extend(start: str, at: str, cost: float, ops: tuple[Operator, ...]) -> None:
        for op in OPERATORS:
            if op.source != at or any(o.target == op.target or o.source == op.target for o in ops):
                continue  # simple chains only
            if op.target == start:
                continue
            c = cost + op.cost
            if c < best[(start, op.target)][0]:
                best[(start, op.target)] = (c, (*ops, op))
            if op.coherent:
                extend(start, op.target, c, (*ops, op))

    for a in reps:
        extend(a, a, 0.0, ())
    return best


CHAINS = _operator_chains()

#: validity(link, representation) -> None if admissible, else the reason (the metrics of 7.2).
Validity = Callable[[Link, str], "str | None"]


@dataclass(frozen=True)
class Requirements:
    """What the plan is for.

    Attributes:
        available: representations the executor can run (R15 and R17: only R).
        outputs: (node, socket) of the requested data outputs; None: all outputs of detectors
            and analyses.
        max_paths: limit of the enumeration of fusion paths (stage 1).
        validity: optional hook rejecting a representation on a link, with a reason.
    """

    available: frozenset[str] = LIGHT_REPRESENTATIONS
    outputs: tuple[tuple[str, str], ...] | None = None
    max_paths: int = 10_000
    validity: Validity | None = None


@dataclass(frozen=True)
class FusionPath:
    """One path from a source output to a requested output: the links in light order."""

    links: tuple[Link, ...]
    output: tuple[str, str]

    @property
    def source(self) -> tuple[str, str]:
        return self.links[0].source


@dataclass(frozen=True)
class Conversion:
    """A conversion inside a block between an input and an output socket."""

    node: str
    in_socket: str
    out_socket: str
    operators: tuple[Operator, ...]

    @property
    def chain(self) -> str:
        if not self.operators:
            return ""
        return "".join([self.operators[0].source, *(op.target for op in self.operators)])


@dataclass
class Plan:
    """Result of the planner: representation and reason per light link, conversions, paths."""

    representations: dict[str, str] = field(default_factory=dict)  # link key -> representation
    reasons: dict[str, str] = field(default_factory=dict)
    conversions: list[Conversion] = field(default_factory=list)
    paths: list[FusionPath] = field(default_factory=list)
    diagnostics: list[Diagnostic] = field(default_factory=list)
    cost: float = 0.0

    @property
    def ok(self) -> bool:
        return not any(d.severity == "error" for d in self.diagnostics)

    def representation(self, link: Link) -> str:
        return self.representations[link_key(link)]

    def to_dict(self) -> dict[str, Any]:
        """JSON form (for the node editor: what the planner made of a link or block)."""
        return {
            "links": {
                k: {"representation": r, "reason": self.reasons.get(k, "")}
                for k, r in self.representations.items()
            },
            "conversions": [
                {
                    "node": c.node,
                    "from": c.in_socket,
                    "to": c.out_socket,
                    "chain": c.chain,
                    "operators": [op.name for op in c.operators],
                }
                for c in self.conversions
            ],
            "paths": [[link_key(lk) for lk in p.links] for p in self.paths],
            "diagnostics": [d.__dict__ for d in self.diagnostics],
            "cost": self.cost,
        }


def link_key(link: Link) -> str:
    return f"{link.source[0]}.{link.source[1]} -> {link.target[0]}.{link.target[1]}"


# ---------------------------------------------------------------------------------- planner


def plan(
    graph: Graph, registry: NodeTypeRegistry, requirements: Requirements | None = None
) -> Plan:
    """Stages 1 and 2 of the planner. The graph should pass `validate` first (sockets known,
    acyclic); links with unknown sockets are left out here."""
    req = requirements or Requirements()
    result = Plan()
    specs = {n.id: {s.name: s for s in node_sockets(n, registry)} for n in graph.nodes.values()}
    through = {n.id: node_through(n, registry) for n in graph.nodes.values()}

    def spec(end: tuple[str, str]) -> SocketSpec | None:
        return specs.get(end[0], {}).get(end[1])

    light = [
        lk
        for lk in graph.links
        if (a := spec(lk.source)) is not None
        and (b := spec(lk.target)) is not None
        and a.kind is SocketKind.LIGHT
        and b.kind is SocketKind.LIGHT
    ]
    into: dict[tuple[str, str], Link] = {lk.target: lk for lk in light}
    out_of: dict[tuple[str, str], list[Link]] = {}
    for lk in light:
        out_of.setdefault(lk.source, []).append(lk)

    # ---- stage 1: topology
    outputs = _outputs(graph, registry, req, specs)
    if not outputs:
        result.diagnostics.append(
            Diagnostic(
                "plan.no_output", "warning", "no detector or analysis output to plan for", ""
            )
        )
    truncated = False

    def walk(end: tuple[str, str], tail: tuple[Link, ...], sink: tuple[str, str]) -> None:
        nonlocal truncated
        if len(result.paths) >= req.max_paths:
            truncated = True
            return
        deps = [d for d in through.get(end[0], {}).get(end[1], ()) if (end[0], d) in into]
        if not deps:
            if tail:
                result.paths.append(FusionPath(tail, sink))
            return
        for d in deps:
            lk = into[(end[0], d)]
            walk(lk.source, (lk, *tail), sink)

    for out in outputs:
        walk(out, (), out)
    if truncated:
        result.diagnostics.append(
            Diagnostic(
                "plan.paths_truncated",
                "warning",
                f"more than {req.max_paths} fusion paths; enumeration stopped",
                "",
            )
        )

    # ---- stage 2: representation per link
    order = _topological(light, through)
    domain: dict[str, list[str]] = {}
    for lk in light:
        key = link_key(lk)
        a, b = spec(lk.source), spec(lk.target)
        assert a is not None and b is not None
        common = a.accepts & b.accepts
        wanted = common if lk.transport == "auto" else common & {lk.transport}
        rejected: list[str] = []
        allowed: list[str] = []
        for r in PREFERENCE:
            if r not in wanted:
                continue
            if r not in req.available:
                rejected.append(f"{r} not available yet")
                continue
            reason = req.validity(lk, r) if req.validity is not None else None
            if reason:
                rejected.append(f"{r}: {reason}")
                continue
            allowed.append(r)
        if lk.transport == "auto" and any(r != "P" for r in allowed):
            allowed = [r for r in allowed if r != "P"]  # "auto" never drops the phase (7.3)
        domain[key] = allowed
        if not allowed:
            if lk.transport != "auto" and lk.transport not in req.available:
                result.diagnostics.append(
                    Diagnostic(
                        "plan.representation_unavailable",
                        "error",
                        f"transport {lk.transport} cannot run yet (available: "
                        f"{''.join(r for r in PREFERENCE if r in req.available)})",
                        key,
                    )
                )
            else:
                result.diagnostics.append(
                    Diagnostic(
                        "plan.no_representation",
                        "error",
                        "no admissible representation"
                        + (f" ({'; '.join(rejected)})" if rejected else ""),
                        key,
                    )
                )

    # forward: cheapest upstream cost of delivering representation r on a link
    up: dict[str, dict[str, float]] = {}
    for lk in order:
        key = link_key(lk)
        deps = [d for d in through.get(lk.source[0], {}).get(lk.source[1], ())]
        costs: dict[str, float] = {}
        for r in domain[key]:
            total = REPRESENTATION_COST[r]
            for d in deps:
                feeder = into.get((lk.source[0], d))
                if feeder is None:
                    continue
                total += min(
                    (up[link_key(feeder)][y] + CHAINS[(y, r)][0] for y in up[link_key(feeder)]),
                    default=math.inf,
                )
            costs[r] = total
        up[key] = costs

    # backward: fix each link given the choices downstream
    chosen: dict[str, str] = {}
    for lk in reversed(order):
        key = link_key(lk)
        consumers = [
            nxt
            for out_name, ins in through.get(lk.target[0], {}).items()
            if lk.target[1] in ins
            for nxt in out_of.get((lk.target[0], out_name), [])
            if link_key(nxt) in chosen
        ]
        if not domain[key]:
            continue  # reported above
        best_r, best_c = domain[key][0], math.inf
        for r in domain[key]:
            c = up[key][r] + sum(CHAINS[(r, chosen[link_key(n)])][0] for n in consumers)
            if c < best_c:
                best_r, best_c = r, c  # with no finite cost the preferred one stays; see below
        chosen[key] = best_r
        if lk.transport != "auto":
            result.reasons[key] = f"fixed by the user ({lk.transport})"
        elif len(domain[key]) == 1:
            result.reasons[key] = f"only {best_r} admissible"
        else:
            result.reasons[key] = f"lowest cost {best_c:g} among {''.join(domain[key])}"
    result.representations = {
        link_key(lk): chosen[link_key(lk)] for lk in light if link_key(lk) in chosen
    }

    # conversions inside blocks, and links whose upstream could not be reached
    for lk in order:
        key = link_key(lk)
        if key not in chosen:
            continue
        r = chosen[key]
        for d in through.get(lk.source[0], {}).get(lk.source[1], ()):
            feeder = into.get((lk.source[0], d))
            if feeder is None or link_key(feeder) not in chosen:
                continue
            y = chosen[link_key(feeder)]
            cost, ops = CHAINS[(y, r)]
            if math.isinf(cost):
                result.diagnostics.append(
                    Diagnostic(
                        "plan.no_operator",
                        "error",
                        f"no operator from {y} to {r} (table 6.3)",
                        f"{lk.source[0]}.{d} -> {lk.source[0]}.{lk.source[1]}",
                    )
                )
            elif ops:
                result.conversions.append(Conversion(lk.source[0], d, lk.source[1], ops))
    result.cost = sum(REPRESENTATION_COST[r] for r in chosen.values()) + sum(
        op.cost for c in result.conversions for op in c.operators
    )

    # coherence (7.3): contributions that meet at an output must all arrive coherently
    incoherent_links = {k for k, r in chosen.items() if r == "P"}
    incoherent_ops = {
        (c.node, c.out_socket)
        for c in result.conversions
        if any(not op.coherent for op in c.operators)
    }
    by_output: dict[tuple[str, str], list[FusionPath]] = {}
    for p in result.paths:
        by_output.setdefault(p.output, []).append(p)
    for out, paths in by_output.items():
        if len(paths) < 2:
            continue
        lossy = [
            p
            for p in paths
            if any(
                link_key(lk) in incoherent_links or lk.source in incoherent_ops for lk in p.links
            )
        ]
        if lossy:
            result.diagnostics.append(
                Diagnostic(
                    "plan.incoherent_merge",
                    "warning",
                    f"{len(lossy)} of {len(paths)} paths arrive without phase; they add as "
                    "intensities, not amplitudes",
                    f"{out[0]}.{out[1]}",
                )
            )
    return result


def _outputs(
    graph: Graph,
    registry: NodeTypeRegistry,
    req: Requirements,
    specs: dict[str, dict[str, SocketSpec]],
) -> list[tuple[str, str]]:
    if req.outputs is not None:
        return list(req.outputs)
    out: list[tuple[str, str]] = []
    for n in graph.nodes.values():
        t = registry.get(n.type_id)
        if t is None or t.category not in (Category.DETECTOR, Category.ANALYSIS):
            continue
        for s in specs[n.id].values():
            if s.direction is Direction.OUT and s.kind is SocketKind.DATA:
                out.append((n.id, s.name))
    return out


def _topological(light: list[Link], through: dict[str, dict[str, tuple[str, ...]]]) -> list[Link]:
    """Light links in light order: a link comes after every link it depends on through a block."""
    into = {lk.target: lk for lk in light}
    deps: dict[str, list[Link]] = {}
    for lk in light:
        ins = through.get(lk.source[0], {}).get(lk.source[1], ())
        deps[link_key(lk)] = [into[(lk.source[0], d)] for d in ins if (lk.source[0], d) in into]
    done: set[str] = set()
    order: list[Link] = []

    def visit(lk: Link, stack: frozenset[str]) -> None:
        key = link_key(lk)
        if key in done or key in stack:  # a cycle is reported by validate
            return
        for d in deps[key]:
            visit(d, stack | {key})
        done.add(key)
        order.append(lk)

    for lk in light:
        visit(lk, frozenset())
    return order
