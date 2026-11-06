"""Reference files in Raytatouille's result format (Raytatouille ADR 0023).

Copied from BWPM (src/bwpm/results.py, itself from Fieldtracing commit 6da9983; MIT, same
maintainer) so that all Fusion modules read one format.

JSON envelope {"format": "raytatouille-result", "schema_version", "type", "data"} with the value
encoding of Raytatouille's `raytatouille.results`:

- arrays as {"dtype", "shape", "values"} in C order; complex values as [re, im];
- NaN and ±infinity inside float arrays as "NaN", "Infinity", "-Infinity"; as scalars
  {"float": "NaN"}; complex scalars {"complex": [re, im]};
- floats in the shortest repr that gives the same double; `allow_nan=False`.

The result type `FieldReference` is a proposal for Raytatouille (docs/rtt-port/): one test case of
a building block with its parameters, inputs, the values of this Python reference implementation
(`expected`, for the port: same algorithm, tight tolerance) and, where available, an independent
value (`reference`, analytic or exact, with the physical tolerance).
"""

from __future__ import annotations

import json
import math
from typing import Any

import numpy as np

FORMAT = "raytatouille-result"
SCHEMA_VERSION = "0.1.13"  # Raytatouille at the commit in rtt-commit.txt
FIELD_REFERENCE = "FieldReference"

DTYPES = frozenset(
    {
        "bool",
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
        "float32",
        "float64",
        "complex64",
        "complex128",
    }
)
_SPECIAL = {"NaN": math.nan, "Infinity": math.inf, "-Infinity": -math.inf}


def _float_item(x: float) -> float | str:
    if math.isfinite(x):
        return x
    return "NaN" if math.isnan(x) else ("Infinity" if x > 0 else "-Infinity")


def encode(value: Any) -> Any:
    """Encode Python/NumPy values as JSON-compatible data (Raytatouille ADR 0023, point 3)."""
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [encode(v) for v in value]
    if isinstance(value, np.ndarray):
        dtype = str(value.dtype)
        if dtype not in DTYPES:
            raise ValueError(f"results: arrays of dtype {dtype} cannot be written")
        flat = value.ravel(order="C").tolist()
        if value.dtype.kind == "f":
            values: list[Any] = [_float_item(x) for x in flat]
        elif value.dtype.kind == "c":
            values = [[_float_item(z.real), _float_item(z.imag)] for z in flat]
        else:
            values = flat
        return {"dtype": dtype, "shape": list(value.shape), "values": values}
    if isinstance(value, np.generic):
        return encode(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else {"float": _float_item(value)}
    if isinstance(value, complex):
        return {"complex": [_float_item(value.real), _float_item(value.imag)]}
    return value


def _number(x: Any) -> float:
    if isinstance(x, str):
        return _SPECIAL[x]
    return float(x)


def decode(value: Any) -> Any:
    """Inverse of `encode`: arrays back to NumPy (dtype and shape), markers to floats."""
    if isinstance(value, list):
        return [decode(v) for v in value]
    if not isinstance(value, dict):
        return value
    if set(value) == {"dtype", "shape", "values"}:
        dtype = np.dtype(value["dtype"])
        if dtype.kind == "f":
            flat: Any = [_number(x) for x in value["values"]]
        elif dtype.kind == "c":
            flat = [complex(_number(re), _number(im)) for re, im in value["values"]]
        else:
            flat = value["values"]
        return np.array(flat, dtype=dtype).reshape(value["shape"])
    if set(value) == {"float"}:
        return _number(value["float"])
    if set(value) == {"complex"}:
        re, im = value["complex"]
        return complex(_number(re), _number(im))
    return {k: decode(v) for k, v in value.items()}


def to_json(type_name: str, data: dict[str, Any], indent: int | None = None) -> str:
    """JSON text of `data` under `type_name` in the format raytatouille-result."""
    envelope = {
        "format": FORMAT,
        "schema_version": SCHEMA_VERSION,
        "type": type_name,
        "data": encode(data),
    }
    return json.dumps(envelope, allow_nan=False, indent=indent, ensure_ascii=False)


def load_json(text: str) -> tuple[str, dict[str, Any]]:
    """(type, decoded data) of a raytatouille-result text; checks the format and version 0.1.x."""
    envelope = json.loads(text)
    if envelope.get("format") != FORMAT:
        raise ValueError(f"results: format must be '{FORMAT}'")
    if not str(envelope.get("schema_version", "")).startswith("0.1."):
        raise ValueError("results: only schema versions 0.1.x are read")
    data = decode(envelope["data"])
    assert isinstance(data, dict)
    return str(envelope["type"]), data


def field_reference(
    case: str,
    description: str,
    source: str,
    parameters: dict[str, Any],
    inputs: dict[str, Any],
    expected: dict[str, Any],
    tolerance: dict[str, float],
    reference: dict[str, Any] | None = None,
    reference_tolerance: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Data of one `FieldReference` (keys in a fixed order).

    Args:
        case: identifier (file stem); description: one sentence; source: test and equation.
        parameters: scalars and strings (units in the key names, e.g. `f_mm`).
        inputs: arrays the port feeds into the function.
        expected: outputs of this implementation; `tolerance` {"rel", "abs"} for the port.
        reference, reference_tolerance: independent values and the physical tolerance.
    """
    return {
        "case": case,
        "description": description,
        "source": source,
        "parameters": parameters,
        "inputs": inputs,
        "expected": expected,
        "tolerance": tolerance,
        "reference": reference or {},
        "reference_tolerance": reference_tolerance or {},
    }
