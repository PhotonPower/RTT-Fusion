"""Data a block holds instead of a file it names: a workflow document is complete on its own.

A workflow document (in Raytatouille Fusion the `.rttf` file) must compute on any machine without
further files – the Raytatouille blocks already carry their system data (`params["rtt"]`). Blocks
that used to name files hold the content in data parameters instead:

- `fg.grating`: `config` (the FourierGrating configuration, FG `docs/config_format.md`) instead of
  `fouriergrating`, `table_data` (the `raytatouille-result` document of a Jones table, stage R6)
  instead of `table`;
- `fg.stack`: `config` instead of `fouriergrating`;
- `rtt.source`: `rtt_catalog_data` (the glass catalogues, AGF, base64) instead of the files in
  `rtt_catalogs` (written by the Raytracer plugin); `load_catalogs` builds the material library
  from both.

The evaluators read the embedded content first and fall back to the file only for models that do
not hold it yet. `embed_nodes` / `embed_files` move a readable file's content into its block once
and clear the path (no personal paths in documents); a file that cannot be read is reported and
left as it is. Paths inside a FourierGrating configuration (`output.dir`, the STEP `file` of CAD
layers) are made absolute as `fouriergrating.load_config` does; a CAD layer's STEP file stays
outside the document and is reported (`graph.file_external`).
"""

from __future__ import annotations

import base64
import copy
import json
import os
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

from rtt_fusion.graph.execute import ExecutionError
from rtt_fusion.graph.model import Node
from rtt_fusion.graph.validate import Diagnostic

#: per block type: (file parameter, data parameter, kind of file)
FILE_PARAMS: dict[str, tuple[tuple[str, str, str], ...]] = {
    "fg.grating": (
        ("fouriergrating", "config", "fouriergrating"),
        ("table", "table_data", "result"),
    ),
    "fg.stack": (("fouriergrating", "config", "fouriergrating"),),
}

#: hint for a block that names a file it cannot read and holds no copy of it
EMBED_HINT = (
    "the block holds no copy of its own – run Check workflow once where the file can be read, "
    "it is then stored in the workflow"
)


def read_fouriergrating(path: str) -> dict[str, Any]:
    """A FourierGrating configuration file with its relative paths (`output.dir`, the `file` of
    CAD layers) made absolute, as `fouriergrating.load_config` does."""
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    if not isinstance(cfg, dict):
        raise ValueError("a FourierGrating configuration is a JSON object")
    base = os.path.dirname(os.path.abspath(path))
    out = cfg.get("output")
    if isinstance(out, dict) and isinstance(out.get("dir"), str) and not os.path.isabs(out["dir"]):
        out["dir"] = os.path.normpath(os.path.join(base, out["dir"]))
    for layer in cfg.get("layers", []):
        if not (isinstance(layer, dict) and layer.get("type") == "cad"):
            continue
        if isinstance(layer.get("file"), str) and not os.path.isabs(layer["file"]):
            layer["file"] = os.path.normpath(os.path.join(base, layer["file"]))
    return cfg


def read_result(path: str) -> dict[str, Any]:
    """A `raytatouille-result` file as its JSON document (checked when it is used)."""
    with open(path, encoding="utf-8") as f:
        doc = json.load(f)
    if not isinstance(doc, dict):
        raise ValueError("a raytatouille-result file is a JSON object")
    return doc


_READERS = {"fouriergrating": read_fouriergrating, "result": read_result}


def content(node: Node, file_param: str, data_param: str, kind: str) -> dict[str, Any] | None:
    """What the block holds in `data_param`, else the content of the file `file_param` names,
    else None (neither). A file that cannot be read raises `exec.param` with the hint to embed
    it."""
    held = node.params.get(data_param)
    if held:
        return copy.deepcopy(dict(held))
    path = str(node.params.get(file_param, "") or "")
    if not path:
        return None
    try:
        return _READERS[kind](path)
    except (OSError, ValueError) as e:
        raise ExecutionError("exec.param", f"cannot read {path}: {e}; {EMBED_HINT}", node.id) from e


def has(node: Node, file_param: str, data_param: str) -> bool:
    """True if the block holds the data or names a file for it."""
    return bool(node.params.get(data_param)) or bool(node.params.get(file_param))


@dataclass
class Embedded:
    """Outcome of embedding: what was stored (lines for the log) and what could not be."""

    stored: list[str]
    diagnostics: list[Diagnostic]


def embed_node(node: Node) -> Embedded:
    """Moves the content of the files a block names into its data parameters and clears the
    paths. Unreadable files stay named (`graph.file_unreadable`); a block that already holds the
    data loses the path only."""
    done = Embedded([], [])
    for file_param, data_param, kind in FILE_PARAMS.get(node.type_id, ()):
        path = str(node.params.get(file_param, "") or "")
        if not path:
            continue
        if node.params.get(data_param):
            node.params[file_param] = ""
            continue
        try:
            value = _READERS[kind](path)
        except (OSError, ValueError) as e:
            done.diagnostics.append(
                Diagnostic(
                    "graph.file_unreadable",
                    "warning",
                    f"{file_param}: cannot read {path} ({e}); the block holds no copy of it and "
                    "cannot compute until the file is readable here once",
                    node.id,
                )
            )
            continue
        node.params[data_param] = value
        node.params[file_param] = ""
        name = os.path.basename(path)
        done.stored.append(f"{node.id}: {name} stored in the workflow ({data_param})")
        if kind == "fouriergrating":
            for k, layer in enumerate(value.get("layers", [])):
                if isinstance(layer, dict) and layer.get("type") == "cad" and layer.get("file"):
                    done.diagnostics.append(
                        Diagnostic(
                            "graph.file_external",
                            "warning",
                            f"layer {k} of the FourierGrating configuration reads the CAD file "
                            f"{layer['file']}, which stays outside the workflow",
                            node.id,
                        )
                    )
    return done


#: what to do about a glass catalogue that cannot be read and is not in the document
CATALOG_HINT = (
    "add it under Glass catalogues in the Raytracer panel (Raytatouille Fusion) and run Check "
    "workflow – it is then stored in the workflow"
)


def catalog_name(path: str) -> str:
    """The name Raytatouille gives a catalogue file (`MaterialLibrary.add_catalog`): its stem in
    capitals."""
    return os.path.splitext(os.path.basename(path))[0].upper()


def embed_catalogs(node: Node) -> Embedded:
    """The glass catalogues of a source (`rtt_catalogs`, files) into `rtt_catalog_data`
    ({"name", "file", "base64"} per catalogue; a file read again replaces its copy). Readable
    and already held files lose their path; an unreadable file the block holds no copy of stays
    named (`graph.file_unreadable`)."""
    done = Embedded([], [])
    paths = [str(p) for p in node.params.get("rtt_catalogs") or [] if p]
    if node.type_id != "rtt.source" or not paths:
        return done
    held = {str(c.get("name")): dict(c) for c in node.params.get("rtt_catalog_data") or []}
    unread: list[tuple[str, OSError]] = []
    for path in paths:  # first what can be read: a file of the same name may stand for another
        try:
            with open(path, "rb") as f:
                data = f.read()
        except OSError as e:
            unread.append((path, e))
            continue
        file = os.path.basename(path)
        held[catalog_name(path)] = {
            "name": catalog_name(path),
            "file": file,
            "base64": base64.b64encode(data).decode("ascii"),
        }
        done.stored.append(f"{node.id}: glass catalogue {file} stored in the workflow")
    left: list[str] = []
    for path, err in unread:
        if catalog_name(path) in held:  # the document holds it: the path goes
            continue
        left.append(path)
        done.diagnostics.append(
            Diagnostic(
                "graph.file_unreadable",
                "warning",
                f"glass catalogue {path}: cannot read it ({err.strerror or err}) and the workflow "
                f"holds no copy – the glasses from it are missing; {CATALOG_HINT}",
                node.id,
            )
        )
    node.params["rtt_catalog_data"] = list(held.values())
    node.params["rtt_catalogs"] = left
    return done


def load_catalogs(library: Any, params: dict[str, Any]) -> list[str]:
    """Adds the glass catalogues of a source's parameters to a Raytatouille `MaterialLibrary`:
    those the workflow holds (`rtt_catalog_data`), then the files of `rtt_catalogs` it does not
    hold. A catalogue that cannot be loaded is left out; the reasons come back for the log."""
    problems: list[str] = []
    loaded: set[str] = set()
    for c in params.get("rtt_catalog_data") or []:
        file = str(c.get("file") or c.get("name"))
        try:
            library.add_catalog_text(
                base64.b64decode(str(c["base64"])), str(c["name"]), source=file
            )
            loaded.add(str(c["name"]))
        except Exception as e:  # a broken catalogue should not stop the others
            problems.append(f"glass catalogue {file}: {e}")
    for path in params.get("rtt_catalogs") or []:
        if not path or catalog_name(str(path)) in loaded:
            continue
        try:
            library.add_catalog(str(path))
        except Exception as e:
            problems.append(
                f"glass catalogue {path}: {e} – the glasses from it are missing; {CATALOG_HINT}"
            )
    return problems


def embed_nodes(nodes: Iterable[Node]) -> Embedded:
    """`embed_node` for every block."""
    done = Embedded([], [])
    for node in nodes:
        for one in (embed_node(node), embed_catalogs(node)):
            done.stored += one.stored
            done.diagnostics += one.diagnostics
    return done


def embed_files(model: Any) -> Embedded:
    """`embed_node` for every block of a model (`topology.Model`), also inside the variants of
    its groups."""
    done = embed_nodes(model.nodes.values())
    for group in model.groups.values():
        for variant in group.variants.values():
            inner = embed_files(variant.model)
            done.stored += inner.stored
            done.diagnostics += inner.diagnostics
    return done


__all__ = [
    "CATALOG_HINT",
    "EMBED_HINT",
    "FILE_PARAMS",
    "Embedded",
    "catalog_name",
    "content",
    "embed_catalogs",
    "embed_files",
    "embed_node",
    "embed_nodes",
    "has",
    "load_catalogs",
    "read_fouriergrating",
    "read_result",
]
