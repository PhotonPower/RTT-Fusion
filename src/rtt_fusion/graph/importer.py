"""Working graph from the explicit paths of a Raytatouille system (concept section 8.3).

Blocks: mirrors, detectors, stops and beam splitters are blocks of their own wherever they sit in
the model tree (a sample mirror inside an objective assembly stays exchangeable); the other
elements of a top-level node form sections – an objective becomes one block. A path is a
sequence of events at surfaces; consecutive events in the same
block form one *traversal*. Each traversal gets an input socket named after the block the light
comes from (`from_<block>`) and an output socket named after the next block (`to_<block>`), so the
same connection in several paths is one link, branches appear where one input feeds two outputs
(beam splitter: split) and merges where two inputs feed one output (combine).

Since R24 the import builds the model (`topology.Model`): the pair `from_X`/`to_X` of a block is
its port X (named after the neighbour; port frames follow with the poses), each pair of neighbours
one connection, and each path of the system a named branch – its runs are the steps, diffraction
orders of its events the channels. `graph_from_raytatouille` unrolls it into the socket form,
which is the same graph as before R24 (W39).

Only the structure of the system is read (root tree, surface ids, paths); no compilation, no
geometry. Accepts the JSON text / dict of a `.rtt.json` file or a `raytatouille.System`.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rtt_fusion.graph.model import Graph, Node, SocketSpec, data_out, light_in, light_out
from rtt_fusion.graph.topology import Branch, Model, Step, unroll

# representations a Raytatouille block can take and give (rays now, fields with rtt-field)
_SECTION_REPS = "RH"


@dataclass
class _Block:
    id: str
    name: str
    location: str
    element_type: str
    surfaces: list[str]
    ins: list[str] = field(default_factory=list)
    outs: list[str] = field(default_factory=list)
    sides: dict[str, str] = field(default_factory=dict)
    through: dict[str, list[str]] = field(default_factory=dict)
    kinds: set[str] = field(default_factory=set)
    paths: set[str] = field(default_factory=set)


def _system_dict(system: Any) -> dict[str, Any]:
    if isinstance(system, dict):
        return system
    if (
        isinstance(system, str | Path)
        and Path(str(system)).suffix == ".json"
        and Path(str(system)).exists()
    ):
        return dict(json.loads(Path(str(system)).read_text(encoding="utf-8")))
    if isinstance(system, str):
        return dict(json.loads(system))
    if hasattr(system, "to_json"):  # raytatouille.System
        return dict(json.loads(system.to_json()))
    raise TypeError("system: a dict, JSON text, a .rtt.json path or a raytatouille.System")


# element types that always form a block of their own (a sample mirror inside an objective
# assembly stays exchangeable; detectors end paths; splitters branch and combine)
_OWN_BLOCK = frozenset({"mirror", "detector", "beam_splitter"})
# a stop is a block of its own only directly below the root (an illumination stop); the aperture
# stop inside an objective assembly belongs to the objective's section


def _units(node: dict[str, Any], location: str) -> list[tuple[str, str, str, list[str]]]:
    """Elements with surfaces below a node, in tree order: (name, location, type, surfaces)."""
    own = [str(s["id"]) for s in node.get("surfaces", []) if "id" in s]
    out = [(str(node.get("name", "")), location, str(node.get("type", "")), own)] if own else []
    for j, child in enumerate(node.get("children", [])):
        out += _units(child, f"{location}/children/{j}")
    return out


def _safe_id(name: str, taken: set[str]) -> str:
    stem = re.sub(r"[^A-Za-z0-9_]+", "_", name).strip("_") or "block"
    candidate, k = stem, 2
    while candidate in taken:
        candidate, k = f"{stem}_{k}", k + 1
    taken.add(candidate)
    return candidate


def model_from_raytatouille(system: Any) -> Model:
    """The model of a Raytatouille system (see the module docstring): blocks with ports named
    after their neighbours, one connection per pair of neighbours, the system's paths as named
    branches."""
    d = _system_dict(system)
    root = d.get("root", {})
    blocks: list[_Block] = []
    owner: dict[str, tuple[_Block, int]] = {}  # surface id -> (block, index within the block)
    taken: set[str] = {"source"}

    def new_block(name: str, location: str, element_type: str, surfaces: list[str]) -> None:
        b = _Block(_safe_id(name, taken), name, location, element_type, surfaces)
        blocks.append(b)
        for k, sid in enumerate(surfaces):
            owner[sid] = (b, k)

    # surfaces at which some path reflects: their elements (beam splitters modelled as thin
    # elements, mirrors) are blocks of their own
    reflecting = {
        str(ev.get("surface", ""))
        for path in d.get("paths", [])
        if isinstance(path.get("events"), list)
        for ev in path["events"]
        if ev.get("kind") == "reflect"
    }

    def own_block(unit: tuple[str, str, str, list[str]]) -> bool:
        return unit[2] in _OWN_BLOCK or any(sid in reflecting for sid in unit[3])

    for i, child in enumerate(root.get("children", [])):
        location = f"/root/children/{i}"
        units = _units(child, location)  # elements with surfaces, in tree order
        if not units:
            continue
        top_name = str(child.get("name") or f"block {i}")
        # special elements are blocks of their own; consecutive ordinary ones form sections
        groups: list[list[tuple[str, str, str, list[str]]]] = []
        for unit in units:
            if own_block(unit):
                groups.append([unit])
            elif groups and not own_block(groups[-1][0]):
                groups[-1].append(unit)
            else:
                groups.append([unit])
        ordinary = [u for u in units if not own_block(u)]
        has_detector = any(u[2] == "detector" for u in units)
        for grp in groups:
            surfaces = [sid for u in grp for sid in u[3]]
            if len(grp) == 1 and (own_block(grp[0]) or len(units) == 1):
                name, loc, etype = grp[0][0] or top_name, grp[0][1], grp[0][2]
            elif len(grp) == len(ordinary) and not has_detector:
                name, loc, etype = top_name, location, "section"
            else:
                name, loc, etype = f"{top_name}: {grp[0][0]} … {grp[-1][0]}", grp[0][1], "section"
            new_block(name, loc, etype, surfaces)
    tree_order = [sid for b in blocks for sid in b.surfaces]

    source_outs: list[str] = []
    links: list[tuple[tuple[str, str], tuple[str, str]]] = []
    path_runs: dict[str, list[tuple[str, str]]] = {}
    order: list[str] = []  # first appearance of blocks, for the layout

    for p_index, path in enumerate(d.get("paths", [])):
        pname = str(path.get("name", f"path {p_index}"))
        events = path.get("events", "auto")
        if events == "auto" or not isinstance(events, list):
            events = [{"surface": sid} for sid in tree_order]
        # runs of consecutive events in the same block
        runs: list[tuple[_Block, list[int], list[str]]] = []
        run_orders: list[list[str]] = []  # the channels (diffraction orders) of each run
        for ev in events:
            sid = str(ev.get("surface", ""))
            if sid not in owner:
                raise ValueError(f"path {pname!r}: unknown surface {sid!r}")
            b, k = owner[sid]
            channel = f"{sid}:m{ev['order']}" if "order" in ev else ""
            if runs and runs[-1][0] is b:
                runs[-1][1].append(k)
                runs[-1][2].append(str(ev.get("kind", "")))
                if channel:
                    run_orders[-1].append(channel)
            else:
                runs.append((b, [k], [str(ev.get("kind", ""))]))
                run_orders.append([channel] if channel else [])
        path_runs[pname] = [(r[0].id, " ".join(o)) for r, o in zip(runs, run_orders, strict=True)]
        prev = "source"
        for j, (b, idx, kinds) in enumerate(runs):
            if b.id not in order:
                order.append(b.id)
            in_name = f"from_{prev}"
            if in_name not in b.ins:
                b.ins.append(in_name)
            b.paths.add(pname)
            b.kinds.update(k for k in kinds if k)
            forward = len(idx) < 2 or idx[-1] >= idx[0]
            b.sides.setdefault(in_name, "" if len(idx) < 2 else ("front" if forward else "back"))
            src = ("source", f"to_{b.id}") if prev == "source" else (prev, f"to_{b.id}")
            if prev == "source" and src[1] not in source_outs:
                source_outs.append(src[1])
            if (src, (b.id, in_name)) not in links:
                links.append((src, (b.id, in_name)))
            if j + 1 < len(runs):
                nxt = runs[j + 1][0]
                out_name = f"to_{nxt.id}"
                if out_name not in b.outs:
                    b.outs.append(out_name)
                b.sides.setdefault(
                    out_name, "" if len(idx) < 2 else ("back" if forward else "front")
                )
                deps = b.through.setdefault(out_name, [])
                if in_name not in deps:
                    deps.append(in_name)
            prev = b.id

    m = Model(
        meta={
            "imported_from": "raytatouille",
            "system": str(d.get("name", "")),
            "raytatouille_paths": list(path_runs),
            # the traversals the system's paths make (block, port in, port out, channel): the
            # Raytatouille backend can pass a branch through these node-wise (R25)
            "raytatouille_runs": sorted(
                {
                    (
                        bid,
                        "source" if j == 0 else runs[j - 1][0],
                        runs[j + 1][0] if j + 1 < len(runs) else "",
                        channel,
                    )
                    for runs in path_runs.values()
                    for j, (bid, channel) in enumerate(runs)
                }
            ),
        }
    )
    m.meta["raytatouille_runs"] = [list(r) for r in m.meta["raytatouille_runs"]]
    m.add(
        Node(
            "source",
            "rtt.source",
            "Source",
            sockets=tuple(light_out(s, "RHGSP", label=s[3:]) for s in source_outs),
            position=(0.0, 0.0),
        )
    )
    for b in blocks:
        if not b.ins:
            continue  # not on any path
        merges = any(len(v) > 1 for v in b.through.values())
        splits = any(sum(1 for v in b.through.values() if i in v) > 1 for i in b.ins)
        if b.element_type == "detector":
            type_id = "rtt.detector"
        elif splits or merges:
            type_id = "rtt.splitter"
        elif b.element_type == "mirror" and len(b.surfaces) == 1:
            type_id = "rtt.mirror"
        else:
            type_id = "rtt.section"
        sockets: tuple[SocketSpec, ...] = tuple(
            light_in(n, _SECTION_REPS, b.sides.get(n, ""), n[5:]) for n in b.ins
        ) + tuple(light_out(n, _SECTION_REPS, b.sides.get(n, ""), n[3:]) for n in b.outs)
        if type_id == "rtt.detector":
            sockets += (data_out("image", "result", "image"),)
        column = order.index(b.id) + 1
        m.add(
            Node(
                b.id,
                type_id,
                b.name,
                params={
                    "location": b.location,
                    "element_type": b.element_type,
                    "surfaces": b.surfaces,
                    "paths": sorted(b.paths),
                },
                position=(220.0 * column, 0.0),
                sockets=sockets,
                through=tuple((o, tuple(ins)) for o, ins in b.through.items()),
            )
        )
    # one connection per pair of neighbours: the ports are named after the neighbour
    # (sockets from_X / to_X); their frames follow with the poses (R26)
    for src, dst in links:
        a, b_end = (src[0], src[1][3:]), (dst[0], dst[1][5:])
        if not any({c.a, c.b} == {a, b_end} for c in m.connections):
            m.connect(a, b_end)
    # the system's paths are the branches of the model (content of the model, ADR P6)
    for pname, runs_of_path in path_runs.items():
        steps = [Step("source", None, runs_of_path[0][0])]
        for j, (bid, channel) in enumerate(runs_of_path):
            came = "source" if j == 0 else runs_of_path[j - 1][0]
            goes = runs_of_path[j + 1][0] if j + 1 < len(runs_of_path) else None
            steps.append(Step(bid, came, goes, channel))
        m.branches.append(Branch(pname, tuple(steps)))
    return m


def graph_from_raytatouille(system: Any) -> Graph:
    """The working graph (socket form) of a Raytatouille system: `model_from_raytatouille`
    unrolled (`topology.unroll`)."""
    from rtt_fusion.graph.registry import builtin_registry

    return unroll(model_from_raytatouille(system), builtin_registry())
