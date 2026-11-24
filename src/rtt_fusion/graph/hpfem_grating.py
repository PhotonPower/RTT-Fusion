"""Grating block with hp-FEM in the loop (stage R35): `hpfem.grating`.

The block is the rigorous grating of `fg.grating` (R19a, R19d) with hp-FEM instead of
FourierGrating: a periodic sample in the place of a mirror (same sockets "in" and "out", so that
`Graph.replace_node` keeps the links) that sends every requested diffraction order on as a
channel of its own. Its grating is an hp-FEM job document (schema v1, problem "grating", as hp-FEM's
`python -m hpfem.run` reads it and the hp-FEM plugin of Raytatouille Fusion writes it: unit cell,
mesh, materials, background stack, solver), named as a file (`job`) or held by the block
(`job_data`, `graph.embedded`); the job's incidence, sweep and maps are not used – the rays set
the incidence.

Around the incident rays of one wavelength the block builds a Jones table
(`couplings.jones_table`, R6) from the hp-FEM adapter (R9, `adapters.hpfem`): one incidence
direction (a plane wave) exactly there; otherwise nodes every `sigma_step` over the σ range of the
rays plus `sigma_margin`, 4 per axis at least (cubic interpolation), at most `max_nodes`. Every
node is one `grating.jones` (2 FEM solves). The table holds the rays' wavelength only (a bundle has
one; every node is computed there, so the interpolation in λ returns it unchanged). Tables and
the prepared cells are kept per job and grid in the operator store (`Context.operators`), which
outlives a run.

Frames and phases as `fg.grating`: local x = `grating_vector` (the job's x), y the lines (the
job's z), z the surface normal; the job's x = 0 is the grating's origin at `point`, and `point`
lies on the top of the structure (the adapter's front port: the highest point of the scatterers
and the stack). σ = n·k̂⊥ with n the job's incidence medium at the rays' wavelength.

The block works in rays (R): fields and spectra arriving at it are turned into rays by the
executor (an hp-FEM solve per plane wave of a field would be too expensive for now).
"""

from __future__ import annotations

import json
from dataclasses import replace
from typing import Any

import numpy as np

from rtt_fusion._types import FloatArray
from rtt_fusion.couplings import jones_table as jt
from rtt_fusion.graph import embedded
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.grating import PLANE_WAVE, _grid, rays_through
from rtt_fusion.graph.model import (
    Category,
    Node,
    NodeType,
    ParamSpec,
    PortFrame,
    PortSpec,
    light_in,
    light_out,
)
from rtt_fusion.graph.registry import NodeTypeRegistry

#: the job file the block may name and the parameter that holds its content (`graph.embedded`)
(JOB,) = embedded.FILE_PARAMS["hpfem.grating"]

#: keys of a job the block does not use (the rays set the incidence; no sweep, no maps)
UNUSED = ("incidence", "sweep", "maps", "name")


def job_of(node: Node) -> dict[str, Any]:
    """The block's hp-FEM job without the keys it does not use (`UNUSED`)."""
    job = embedded.content(node, *JOB)
    if job is None:
        raise ExecutionError(
            "exec.param", "the hp-FEM grating needs a job (a job file or one held)", node.id
        )
    for key in UNUSED:
        job.pop(key, None)
    return job


def _store(ctx: Context) -> dict[str, Any]:
    """The operator store: kept across runs (a provider's requests, the engine's runs)."""
    store: dict[str, Any] = (
        ctx.operators if ctx.operators is not None else ctx.shared.setdefault("operators", {})
    )
    return store


def _structure(ctx: Context, node: Node, job: dict[str, Any], wavelength_um: float) -> Any:
    """The adapter's cell of the job at one wavelength (cached)."""
    from rtt_fusion.adapters import hpfem as hpa

    store = _store(ctx)
    key = "hpfem-structure|" + json.dumps([job, round(wavelength_um, 12)], sort_keys=True)
    if key not in store:
        try:
            store[key] = hpa.structure_from_job(job, wavelength_um)
        except ImportError as e:
            raise ExecutionError("exec.solver", f"hp-FEM is not installed: {e}", node.id) from e
        except (ValueError, RuntimeError, KeyError, TypeError) as e:
            raise ExecutionError("exec.param", f"hp-FEM job: {e}", node.id) from e
    return store[key]


def _one_wavelength(table: jt.JonesTable, wavelength_um: float) -> jt.JonesTable:
    """A table of one wavelength as one of 4 equal λ nodes ± 1e-7 µm (cubic stencil)."""
    offsets = np.array([-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0]) * 1e-7
    jones = np.repeat(np.asarray(table.jones)[:, :1], 4, axis=1)
    return replace(table, wavelengths_um=wavelength_um + offsets, jones=jones)


def hpfem_table(ctx: Context, node: Node, wavelength_um: float, d_loc: FloatArray) -> jt.JonesTable:
    """The block's Jones table from hp-FEM around the incident directions `d_loc` (local unit
    vectors of the live rays, z along the light), cached per job, orders and grid."""
    from rtt_fusion.adapters import hpfem as hpa

    job = job_of(node)
    side = str(node.params.get("side", "reflected"))
    orders = tuple(int(m) for m in node.params.get("orders", []) or (-1, 0, 1))
    margin = float(node.params.get("sigma_margin", 0.01))
    step = float(node.params.get("sigma_step", 0.01))
    max_nodes = int(node.params.get("max_nodes", 400))
    if not step > 0.0 or margin < 0.0:
        raise ExecutionError(
            "exec.param", "the table's σ step must be > 0 and its margin ≥ 0", node.id
        )
    structure = _structure(ctx, node, job, wavelength_um)
    sigma = float(structure.n_cover) * np.asarray(d_loc)[:, :2]
    plane = bool(np.ptp(sigma[:, 0]) <= PLANE_WAVE and np.ptp(sigma[:, 1]) <= PLANE_WAVE)
    if plane:
        s0 = (float(np.mean(sigma[:, 0])), float(np.mean(sigma[:, 1])))
        grid: list[Any] = ["plane", round(s0[0], 12), round(s0[1], 12)]
    else:
        sx = _grid(sigma[:, 0], margin, step)
        sy = _grid(sigma[:, 1], margin, step)
        if len(sx) * len(sy) > max_nodes:
            raise ExecutionError(
                "exec.param",
                f"the rays span σ_x {sx[0]:.4g}…{sx[-1]:.4g}, σ_y {sy[0]:.4g}…{sy[-1]:.4g}: "
                f"{len(sx)}×{len(sy)} table nodes (2 hp-FEM solves each) are more than "
                f"{max_nodes} – a larger σ step or a larger node limit (the block's Ray path "
                "settings)",
                node.id,
            )
        grid = [sx.tolist(), sy.tolist()]
    key = "hpfem-jones|" + json.dumps(
        [job, side, orders, round(wavelength_um, 12), grid], sort_keys=True
    )
    store = _store(ctx)
    if key not in store:
        try:
            handle = hpa.prepare(structure, (wavelength_um,))
            if plane:
                store[key] = jt.point_table(handle, wavelength_um, *s0, orders, side)
            else:
                table = jt.build_table(handle, np.array([wavelength_um]), sx, sy, orders, side)
                store[key] = _one_wavelength(table, wavelength_um)
        except ImportError as e:
            raise ExecutionError("exec.solver", f"hp-FEM is not installed: {e}", node.id) from e
        except (ValueError, RuntimeError, KeyError) as e:
            raise ExecutionError("exec.solver", f"hp-FEM: {e}", node.id) from e
        if plane:
            ctx.log.append(
                f"{node.id}: Jones matrices from hp-FEM at one incidence (σ = {s0[0]:.6g}, "
                f"{s0[1]:.6g}), orders {list(orders)}, {side}"
            )
        else:
            ctx.log.append(
                f"{node.id}: Jones table from hp-FEM, {len(sx)}×{len(sy)} nodes at λ = "
                f"{wavelength_um:.6g} µm, orders {list(orders)}, {side}"
            )
    table_out: jt.JonesTable = store[key]
    return table_out


def hpfem_grating(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
    def table_for(lam: float, d_loc: FloatArray) -> jt.JonesTable:
        return hpfem_table(ctx, node, lam, d_loc)

    return list(rays_through(ctx, node, inputs.get("in", []), table_for, True))


def _frame(n: Node) -> dict[str, PortFrame]:
    normal = np.asarray(n.params.get("normal", [0.0, 0.0, 1.0]), dtype=float)
    normal = normal / np.linalg.norm(normal)
    point = n.params.get("point", [0.0, 0.0, 0.0])
    return {
        "front": PortFrame(
            (float(point[0]), float(point[1]), float(point[2])),
            (float(normal[0]), float(normal[1]), float(normal[2])),
        )
    }


TYPE = NodeType(
    "hpfem.grating",
    Category.SPLITTER,
    (light_in("in", "RHS"), light_out("out", "RHS", label="orders")),
    (
        ParamSpec(
            "job",
            "file",
            "",
            description="hp-FEM job (JSON, problem 'grating'): cell, mesh, materials, stack and "
            "solver; Check workflow stores it in the block (job_data) and clears the path",
        ),
        ParamSpec(
            "job_data",
            "data",
            {},
            description="the hp-FEM job held by the block (before the file)",
        ),
        ParamSpec(
            "side",
            "choice",
            "reflected",
            options=("reflected", "transmitted"),
            description="the orders passed on",
        ),
        ParamSpec("orders", "numbers", [], description="orders to pass on; empty: -1, 0, 1"),
        ParamSpec(
            "sigma_step",
            "number",
            0.01,
            minimum=1e-5,
            group="Ray path",
            label="σ step of the table",
            description="spacing of the Jones table's nodes in the direction cosines σ = n·sin θ "
            "of the incident rays (cubic interpolation between them); every node is 2 hp-FEM "
            "solves",
        ),
        ParamSpec(
            "sigma_margin",
            "number",
            0.01,
            minimum=0.0,
            group="Ray path",
            label="σ margin of the table",
            description="how far the table reaches beyond the σ range of the incident rays",
        ),
        ParamSpec(
            "max_nodes",
            "integer",
            400,
            minimum=16,
            group="Ray path",
            label="Most table nodes",
            description="a wider bundle stops with a message instead of starting this many "
            "hp-FEM solves",
        ),
        ParamSpec("point", "point", [0.0, 0.0, 0.0], unit="mm", description="top of the grating"),
        ParamSpec("normal", "point", [0.0, 0.0, 1.0]),
        ParamSpec(
            "grating_vector", "point", [1.0, 0.0, 0.0], description="local x (the job's x), ⟂ lines"
        ),
    ),
    (("out", ("in",)),),
    ports=(PortSpec("front", "in", "out", "front"),),
    scatter=(("front", "front", "full", ()),),
    port_frames=_frame,
    description="Periodic sample computed by hp-FEM (finite elements): diffraction orders as "
    "channels",
)

hpfem_grating.representations = frozenset("R")  # type: ignore[attr-defined]

EVALUATORS: dict[str, Evaluator] = {"hpfem.grating": hpfem_grating}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    registry.add(TYPE)
    return dict(EVALUATORS)
