"""Port data: the values on the workflow's links as files (CAD ADR 0010, D4, stage W2a).

A block evaluated in another plugin's process gets its inputs and returns its outputs through this
format: a JSON tree (the structure: dataclasses, lists, tuples, dicts, numbers, strings) plus one
`.npy` file per NumPy array in a folder. Values come back exactly: floats travel as their hex form
(bit for bit, NaN and infinities included – the channel's JSON parser knows neither), arrays as
`.npy`. Only dataclasses and enums of the allowed packages are created when reading (`ALLOWED`),
so a plugin cannot make the reader build arbitrary objects.

    handle = write(value, folder)        # JSON-able dict; arrays in folder
    value = read(handle, folder)
"""

from __future__ import annotations

import dataclasses
import enum
import importlib
from pathlib import Path
from typing import Any

import numpy as np

FORMAT = "rtt-fusion-portdata"
VERSION = 1
#: packages whose dataclasses and enums `read` may create
ALLOWED = ("rtt_fusion.", "fieldtracing.")


class PortDataError(ValueError):
    """A handle that is not port data of this version, or names a class outside `ALLOWED`."""


def write(value: Any, folder: str | Path) -> dict[str, Any]:
    """The handle of `value`: {"format", "version", "root"}; its arrays as `a<k>.npy` in
    `folder` (created if needed)."""
    out = Path(folder)
    out.mkdir(parents=True, exist_ok=True)
    count = [0]

    def enc(v: Any) -> dict[str, Any]:
        if v is None:
            return {"t": "none"}
        if isinstance(v, bool | np.bool_):
            return {"t": "bool", "v": bool(v)}
        if isinstance(v, enum.Enum):
            cls = type(v)
            return {"t": "enum", "cls": f"{cls.__module__}.{cls.__qualname__}", "v": enc(v.value)}
        if isinstance(v, int | np.integer):
            return {"t": "int", "v": str(int(v))}
        if isinstance(v, float | np.floating):
            return {"t": "float", "v": float(v).hex()}
        if isinstance(v, complex | np.complexfloating):
            c = complex(v)
            return {"t": "complex", "re": c.real.hex(), "im": c.imag.hex()}
        if isinstance(v, str):
            return {"t": "str", "v": v}
        if isinstance(v, np.ndarray):
            name = f"a{count[0]}.npy"
            count[0] += 1
            np.save(out / name, np.ascontiguousarray(v), allow_pickle=False)
            return {"t": "array", "file": name}
        if dataclasses.is_dataclass(v) and not isinstance(v, type):
            dcls = type(v)
            return {
                "t": "dataclass",
                "cls": f"{dcls.__module__}.{dcls.__qualname__}",
                "fields": {f.name: enc(getattr(v, f.name)) for f in dataclasses.fields(v)},
            }
        if isinstance(v, tuple):
            return {"t": "tuple", "items": [enc(x) for x in v]}
        if isinstance(v, list):
            return {"t": "list", "items": [enc(x) for x in v]}
        if isinstance(v, dict):
            return {"t": "dict", "items": [[enc(k), enc(x)] for k, x in v.items()]}
        raise PortDataError(f"cannot write a {type(v).__name__} as port data")

    return {"format": FORMAT, "version": VERSION, "root": enc(value)}


def _allowed_class(name: str) -> Any:
    """The class `module.Qual.Name` of an allowed package (the longest importable module prefix)."""
    if not name.startswith(ALLOWED):
        raise PortDataError(f"port data names a class outside {ALLOWED}: {name}")
    parts = name.split(".")
    for cut in range(len(parts) - 1, 0, -1):
        try:
            obj: Any = importlib.import_module(".".join(parts[:cut]))
        except ImportError:
            continue
        for attr in parts[cut:]:
            obj = getattr(obj, attr, None)
        if isinstance(obj, type):
            return obj
        break
    raise PortDataError(f"no such class: {name}")


def read(handle: dict[str, Any], folder: str | Path, mmap: bool = False) -> Any:
    """The value of a handle written by `write`; `mmap`: arrays mapped read-only instead of read."""
    if handle.get("format") != FORMAT or handle.get("version") != VERSION:
        raise PortDataError(f"not {FORMAT} version {VERSION}")
    src = Path(folder)

    def dec(n: dict[str, Any]) -> Any:
        t = n["t"]
        if t == "none":
            return None
        if t == "bool":
            return bool(n["v"])
        if t == "int":
            return int(n["v"])
        if t == "float":
            return float.fromhex(n["v"])
        if t == "complex":
            return complex(float.fromhex(n["re"]), float.fromhex(n["im"]))
        if t == "str":
            return str(n["v"])
        if t == "array":
            name = str(n["file"])
            if "/" in name or "\\" in name or name.startswith("."):
                raise PortDataError(f"array file outside the folder: {name}")
            return np.load(src / name, mmap_mode="r" if mmap else None, allow_pickle=False)
        if t == "tuple":
            return tuple(dec(x) for x in n["items"])
        if t == "list":
            return [dec(x) for x in n["items"]]
        if t == "dict":
            return {dec(k): dec(x) for k, x in n["items"]}
        if t == "enum":
            return _allowed_class(n["cls"])(dec(n["v"]))
        if t == "dataclass":
            cls: Any = _allowed_class(n["cls"])
            if not dataclasses.is_dataclass(cls):
                raise PortDataError(f"not a dataclass: {n['cls']}")
            values = {k: dec(x) for k, x in n["fields"].items()}
            init = {f.name for f in dataclasses.fields(cls) if f.init}
            make: Any = cls  # a dataclass type (checked above)
            obj = make(**{k: v for k, v in values.items() if k in init})
            for k, v in values.items():  # fields set after __init__ (init=False)
                if k not in init:
                    object.__setattr__(obj, k, v)
            return obj
        raise PortDataError(f"unknown port data node {t!r}")

    return dec(handle["root"])


__all__ = ["ALLOWED", "FORMAT", "VERSION", "PortDataError", "read", "write"]
