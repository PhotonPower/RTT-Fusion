"""The agent reference REFERENCE.md shipped in the package: read it, a section of it, the block
catalogue (also generated live from a registry) and the command line `python -m
rtt_fusion.reference`.

The block catalogue in REFERENCE.md is generated from the registry (`catalogue_markdown`, written
by tools/gen_llm_reference.py between `BEGIN` and `END`); tests/test_llm_reference.py checks that
it is up to date.
"""

from __future__ import annotations

import argparse
import io
import json
import re
import sys
from collections.abc import Callable
from importlib import resources
from typing import Any

FILE = "REFERENCE.md"
BEGIN = "<!-- BEGIN GENERATED: block catalogue (python tools/gen_llm_reference.py) -->"
END = "<!-- END GENERATED: block catalogue -->"

_RTT = (
    "raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in "
    "`workflow.registry()` (section 2.9)"
)
_OPTICS = "nothing (numpy/scipy); H or S on its links: fieldtracing"

#: what a block type needs beyond numpy/scipy (the registry does not say it); every registered
#: type must have an entry (tests/test_llm_reference.py)
NEEDS: dict[str, str] = {
    "analysis.camera_image": "no evaluator in this package (planned type: `exec.no_evaluator`)",
    "biooptomd.layer": "nothing; reads a run directory (`adapters.biooptomd.read_layer`, "
    "manifest `biooptomd-layer.json`); evaluator in `rtt_fusion.graph.datanodes`",
    "bindoptix.time_loop": "nothing itself (meta block; its B from `data.occupancy` needs "
    "bindoptix); evaluator in `rtt_fusion.graph.time_loop`",
    "bindoptix.signal": "bindoptix (pure Python without its C++ core; `adapters.bindoptix`); "
    "evaluator in `rtt_fusion.graph.datanodes`",
    "data.field_map": "bindoptix (reads the `raytatouille-result` file)",
    "data.occupancy": "bindoptix (reads the `raytatouille-result` file)",
    "data.response_map": "bindoptix (reads the `raytatouille-result` file)",
    "display.beam": "nothing; fields on its input: fieldtracing",
    "optics.ideal_imaging": "fieldtracing (point-source sampling of the image spectrum); "
    "evaluator in `rtt_fusion.graph.ideal_imaging`",
    "emitter.response_map": "bindoptix (builds its EmitterResponseTable); physics in "
    "`couplings.dipole_layers`",
    "emitter.dipole_layer": "fieldtracing (the spectra carry a point-source sampling); "
    "evaluator in `rtt_fusion.graph.emitter`, physics in `couplings.dipole_layers`",
    "display.image": "nothing",
    "display.plot": "nothing",
    "fg.grating": "`table`/`table_data` (Jones table): nothing for rays, the optional C++ "
    "extension `rtt_fusion._diffraction` speeds them up; `fouriergrating`/`config`: "
    "fouriergrating; H or S: fieldtracing and a configuration (a table turns fields and spectra "
    "into rays)",
    "fg.stack": "fouriergrating; H or S: also fieldtracing",
    "hpfem.grating": "hpfem (builds its Jones table from the job, 2 FEM solves per node); "
    "evaluator in `rtt_fusion.graph.hpfem_grating`; H or S arrive as rays (fieldtracing for "
    "the conversion)",
    "fusion.condensed": "made by `unroll` from a group with `condense`; what its inner blocks "
    "need; S: fieldtracing",
    "fusion.convert": "fieldtracing for every real conversion (`exec.convert_unsupported` "
    "without it)",
    "fusion.particle": "cliffordbem; S in or out: also fieldtracing",
    "fusion.sum": "nothing (structural: the executor adds what arrives; made by `unroll`)",
    "group": "nothing (structural: `unroll`/`flatten` resolve it into its blocks)",
    "meta.optimise": "nothing",
    "meta.replicates": "nothing",
    "meta.sweep": "nothing",
    "optics.aplanatic_lens": _OPTICS,
    "optics.beam_splitter": _OPTICS,
    "optics.collimated_source": _OPTICS,
    "optics.detector": _OPTICS,
    "optics.interface": _OPTICS,
    "optics.plane_mirror": _OPTICS,
    "optics.splitter": _OPTICS,
    "optics.thin_lens": _OPTICS,
    "rtt.detector": _RTT + "; H or S: also fieldtracing",
    "rtt.mirror": _RTT,
    "rtt.part.detector": _RTT + "; H or S: also fieldtracing",
    "rtt.part.lens": _RTT,
    "rtt.part.mirror": _RTT,
    "rtt.part.splitter": _RTT,
    "rtt.part.stop": _RTT,
    "rtt.section": _RTT,
    "rtt.source": _RTT,
    "rtt.splitter": _RTT,
}

_FENCE = re.compile(r"^\s*(```|~~~)")
_HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")


def text() -> str:
    """The whole reference (REFERENCE.md of the installed package)."""
    return resources.files("rtt_fusion").joinpath(FILE).read_text(encoding="utf-8")


def headings(doc: str | None = None) -> list[tuple[int, int, str]]:
    """(line index, level, title) of every Markdown heading outside code fences."""
    out: list[tuple[int, int, str]] = []
    fenced = False
    for i, line in enumerate((text() if doc is None else doc).splitlines()):
        if _FENCE.match(line):
            fenced = not fenced
            continue
        m = None if fenced else _HEADING.match(line)
        if m:
            out.append((i, len(m.group(1)), m.group(2)))
    return out


def _slug(title: str) -> str:
    return re.sub(r"[^a-z0-9._]+", "-", title.lower()).strip("-")


def section(name: str, doc: str | None = None) -> str:
    """The section whose heading matches `name` (case-insensitive): its number ("3", "3.2"), its
    title or its first word (a block type id such as "fg.grating"), else the first heading that
    contains `name`. From the heading to the next heading of the same or a higher level. Raises
    KeyError with the headings if nothing matches."""
    doc = text() if doc is None else doc
    lines = doc.splitlines()
    found = headings(doc)
    want = name.strip().lower()

    def number(title: str) -> str:
        return title.split()[0].rstrip(".").lower() if title.split() else ""

    def first_word(title: str) -> str:
        return title.split()[0].strip("`").lower() if title.split() else ""

    tests: tuple[Callable[[str], bool], ...] = (
        lambda t: number(t) == want.rstrip("."),
        lambda t: t.lower() == want or _slug(t) == _slug(want),
        lambda t: first_word(t) == want,
        lambda t: want in t.lower(),
    )
    for test in tests:
        hit = next((h for h in found if test(h[2])), None)
        if hit is not None:
            start, level, _ = hit
            end = next((i for i, lv, _ in found if i > start and lv <= level), len(lines))
            return "\n".join(lines[start:end]).rstrip() + "\n"
    listing = "\n".join("  " * (lv - 1) + t for _, lv, t in found if lv <= 3)
    raise KeyError(f"no section {name!r}; sections:\n{listing}")


def generated(doc: str | None = None) -> str:
    """The generated block catalogue of the reference (between `BEGIN` and `END`)."""
    doc = text() if doc is None else doc
    start, end = doc.index(BEGIN) + len(BEGIN), doc.index(END)
    return doc[start:end].strip("\n") + "\n"


def get(name: str | None = None) -> str:
    """The reference, or one section of it (`section`)."""
    return text() if name is None else section(name)


# ---------------------------------------------------------------------------- generation


def _value(v: Any, limit: int = 70) -> str:
    s = json.dumps(v, ensure_ascii=False, sort_keys=True)
    if len(s) > limit:
        kind = "{…}" if isinstance(v, dict) else "[…]" if isinstance(v, list) else s[:limit] + "…"
        return f"{kind} (see `registry.get(...).params`)"
    return s


def _param_line(p: Any) -> str:
    bits = [f"`{p.name}`", p.type, f"= {_value(p.default)}"]
    if p.unit:
        bits.append(f"[{p.unit}]")
    if p.minimum is not None or p.maximum is not None:
        lo = "" if p.minimum is None else f"{p.minimum:g}"
        hi = "" if p.maximum is None else f"{p.maximum:g}"
        bits.append(f"range {lo}…{hi}")
    if p.options:
        bits.append("of " + "|".join(p.options))
    text_ = " ".join(bits)
    words = "; ".join(x for x in (p.label, p.description) if x)
    return f"- param {text_}" + (f" – {words}" if words else "")


def _socket(s: Any) -> str:
    arrow = "in" if s.direction.value == "in" else "out"
    kind = "".join(sorted(s.accepts)) if s.kind.value == "light" else next(iter(s.accepts))
    side = f"@{s.side}" if s.side else ""
    return f"`{s.name}` {arrow} {s.kind.value}:{kind}{side}"


def catalogue_markdown(registry: Any = None, evaluators: dict[str, Any] | None = None) -> str:
    """The block catalogue as Markdown: per type id its category, description, what it needs,
    the representations its evaluator takes natively, sockets (name, in/out, light:accepted
    representations or data:type, @side), ports, through-paths and parameters (type, default,
    unit, range, options). Default: the workflow's registry (`workflow.registry()`) and its
    evaluators plus those of the Raytatouille backend."""
    from rtt_fusion.graph.topology import node_ports

    if registry is None:
        from rtt_fusion.graph import backend_rtt, workflow
        from rtt_fusion.graph.model import Node

        registry, extra = workflow.registry()
        evaluators = {**backend_rtt.EVALUATORS, **extra}
    from rtt_fusion.graph.model import Node

    evaluators = evaluators or {}
    types = registry.types()
    out = [
        f"{len(types)} block types. Per type: category and description; needs (optional "
        "packages); native = representations the evaluator computes in (the executor converts "
        "the rest at its ports); sockets as `name` in|out light:<accepted representations>"
        "[@side] or data:<type> (blocks of category source/domain/splitter/detector also get "
        "the data input `pose`, type drive); ports as name(in socket→out socket), the names "
        "`Model.connect` and branch steps use; through = output ← inputs it depends on; "
        "`param name type = default [unit] range of options – description`.",
        "",
    ]
    for t in types:
        out.append(f"### {t.type_id}")
        ev = evaluators.get(t.type_id)
        if ev is None:
            native = "– (no evaluator)"
        else:
            native = "".join(sorted(getattr(ev, "representations", frozenset("R")))) or (
                "– (data block, no light)"
            )
            if hasattr(ev, "produces"):
                native += " (`produces`: its output carries param `to`)"
        description = t.description.rstrip(".")
        out.append(f"{t.category.value}: {description}.")
        out.append(f"- needs: {NEEDS.get(t.type_id, '(not listed)')}; native: {native}")
        sockets = [_socket(s) for s in (*t.sockets, *t.extra_sockets)]
        if t.dynamic_sockets:
            sockets.append("+ dynamic (per node, from the import or `unroll`)")
        out.append("- sockets: " + ("; ".join(sockets) if sockets else "–"))
        ports = node_ports(Node("x", t.type_id), registry)
        if ports:
            out.append(
                "- ports: "
                + "; ".join(f"{p.name}({p.in_socket or '·'}→{p.out_socket or '·'})" for p in ports)
            )
        if t.through:
            out.append(
                "- through: " + "; ".join(f"{o} ← {', '.join(ins) or '–'}" for o, ins in t.through)
            )
        out += [_param_line(p) for p in t.params]
        out.append("")
    return "\n".join(out).rstrip() + "\n"


def replace_generated(doc: str, body: str) -> str:
    """`doc` with the generated block catalogue replaced by `body`."""
    start, end = doc.index(BEGIN) + len(BEGIN), doc.index(END)
    return doc[:start] + "\n" + body.rstrip("\n") + "\n" + doc[end:]


# ---------------------------------------------------------------------------- command line


def main(argv: list[str] | None = None) -> int:
    """`python -m rtt_fusion.reference [--list | --section NAME | --blocks [TYPE ...] [--live]]`"""
    parser = argparse.ArgumentParser(
        prog="python -m rtt_fusion.reference",
        description="Print the RTT-Fusion reference for agents (REFERENCE.md in the package).",
    )
    parser.add_argument("--list", action="store_true", help="list the section headings")
    parser.add_argument("--section", metavar="NAME", help="one section: number, title or word")
    parser.add_argument(
        "--blocks",
        nargs="*",
        metavar="TYPE",
        help="the block catalogue, or the entries of the given block type ids",
    )
    parser.add_argument(
        "--live",
        action="store_true",
        help="with --blocks: generate the catalogue from the installed registry",
    )
    args = parser.parse_args(argv)
    if isinstance(sys.stdout, io.TextIOWrapper) and not sys.stdout.isatty():
        sys.stdout.reconfigure(encoding="utf-8")  # a pipe on Windows would be the ANSI code page
    try:
        if args.list:
            print("\n".join("  " * (lv - 1) + t for _, lv, t in headings()))
        elif args.blocks is not None:
            body = catalogue_markdown() if args.live else generated()
            if args.blocks:
                print("\n".join(section(b, body) for b in args.blocks), end="")
            else:
                print(body, end="")
        elif args.section:
            print(section(args.section), end="")
        else:
            print(text(), end="")
    except KeyError as e:
        print(e.args[0], file=sys.stderr)
        return 2
    except BrokenPipeError:  # `| head`
        sys.stderr.close()
    return 0
