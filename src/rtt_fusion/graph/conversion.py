"""The conversion block `fusion.convert`: light from one representation into another.

R -> H (rays to a field on the port plane, `field.rays_to_field`) and H -> R (`field.field_to_rays`)
need fieldtracing; the other representations (G, S, P) come with later stages and are reported
as such. Light that already has the target representation passes unchanged. The engine converts
on links too (a link's planned representation); the block makes the conversion explicit and
visible at a chosen place.
"""

from __future__ import annotations

from typing import Any

from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.light import RayBundle
from rtt_fusion.graph.model import Node
from rtt_fusion.graph.registry import NodeTypeRegistry

#: target representations this stage can produce
SUPPORTED = frozenset("RH")


def convert(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> list[Any]:
    from rtt_fusion.graph.field import FieldBundle, FieldError, field_to_rays, rays_to_field

    to = str(node.params.get("to", "S"))
    if to not in SUPPORTED:
        raise ExecutionError(
            "exec.convert_unsupported",
            f"conversion to {to} comes with a later stage (this version converts to "
            f"{' or '.join(sorted(SUPPORTED))})",
            node.id,
        )
    out: list[Any] = []
    try:
        for item in inputs.get("in") or []:
            if to == "H" and isinstance(item, RayBundle):
                item = rays_to_field(item)
            elif to == "R" and isinstance(item, FieldBundle):
                item = field_to_rays(item)
            out.append(item)
    except FieldError as e:
        raise ExecutionError("exec.field", str(e), node.id) from e
    except ImportError as e:  # fieldtracing missing
        raise ExecutionError(
            "exec.convert_unsupported",
            f"conversion to {to} needs the package fieldtracing ({e})",
            node.id,
        ) from e
    return out


# fields arrive as fields (the engine does not turn them into rays before the block)
convert.representations = frozenset("RH")  # type: ignore[attr-defined]

EVALUATORS: dict[str, Evaluator] = {"fusion.convert": convert}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    """The type is part of the built-in registry; this adds its evaluator."""
    return dict(EVALUATORS)
