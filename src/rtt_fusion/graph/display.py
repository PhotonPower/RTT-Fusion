"""Display blocks (stage R21): views of results and of the light at any block output.

Display blocks compute nothing physical of their own; they turn what is attached to them into a
`View` – an image or curves with title and axis labels – for the GUI (Raytatouille Fusion opens
a dockable window per display block):

- `display.image`: an image result (a camera, a detector) – the image;
- `display.plot`: a meta block's result (sweep, replicates, optimisation) – the measure over the
  parameter; an image – its central row and column; numbers – over their index;
- `display.beam`: light of any block output (a probe; an output may feed several links): the
  intensity on the port plane (the plane through the centroid of the light, perpendicular to its
  mean direction) – rays as a coherent ray image (`imaging.coherent_image`), fields as the field
  image on their plane (`field.field_image`), amplitudes within, intensities across coherence
  groups. The window covers the light with a margin.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from rtt_fusion._types import FloatArray
from rtt_fusion.graph.execute import Context, Evaluator, ExecutionError
from rtt_fusion.graph.imaging import Image, coherent_image
from rtt_fusion.graph.light import RayBundle
from rtt_fusion.graph.meta import MetaResult
from rtt_fusion.graph.model import (
    Category,
    Node,
    NodeType,
    ParamSpec,
    data_in,
    data_out,
    light_in,
)
from rtt_fusion.graph.registry import NodeTypeRegistry


@dataclass
class Curve:
    name: str
    x: list[float]
    y: list[float]


@dataclass
class View:
    """What a display block shows: an image or curves."""

    kind: str  # "image" or "curves"
    title: str
    image: Image | None = None
    curves: list[Curve] = field(default_factory=list)
    x_label: str = ""
    y_label: str = ""
    unit: str = ""


def _title(node: Node, fallback: str) -> str:
    return str(node.params.get("title") or node.name or fallback)


def _single(inputs: dict[str, Any], node: Node) -> Any:
    value = inputs.get("in")
    if value is None or (isinstance(value, list) and not value):
        raise ExecutionError("exec.display", "nothing attached to the display", node.id)
    return value


def image_view(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> View:
    value = _single(inputs, node)
    if isinstance(value, View) and value.image is not None:
        value = value.image
    if not isinstance(value, Image):
        raise ExecutionError(
            "exec.display", f"an image is expected, not a {type(value).__name__}", node.id
        )
    return View(
        "image",
        _title(node, "image"),
        image=value,
        x_label="u (mm)",
        y_label="v (mm)",
        unit="W/mm² (rel.)",
    )


def plot_view(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> View:
    value = _single(inputs, node)
    title = _title(node, "plot")
    if isinstance(value, MetaResult):
        return View(
            "curves",
            title,
            curves=[Curve(value.measure, list(value.values), list(value.measures))],
            x_label=value.parameter,
            y_label=value.measure,
        )
    if isinstance(value, View) and value.image is not None:
        value = value.image
    if isinstance(value, Image):
        v = np.asarray(value.values)
        ny, nx = v.shape
        u0, u1, v0, v1 = value.extent
        u = [float(x) for x in u0 + (np.arange(nx) + 0.5) * (u1 - u0) / nx]
        w = [float(x) for x in v0 + (np.arange(ny) + 0.5) * (v1 - v0) / ny]
        return View(
            "curves",
            title,
            curves=[
                Curve("row through the centre", u, [float(x) for x in v[ny // 2, :]]),
                Curve("column through the centre", w, [float(x) for x in v[:, nx // 2]]),
            ],
            x_label="position (mm)",
            y_label="intensity",
        )
    if isinstance(value, list | tuple) and all(isinstance(x, int | float) for x in value):
        return View(
            "curves",
            title,
            curves=[
                Curve("values", [float(k) for k in range(len(value))], [float(x) for x in value])
            ],
            x_label="index",
        )
    raise ExecutionError("exec.display", f"cannot plot a {type(value).__name__}", node.id)


def _plane_of(bundles: list[RayBundle]) -> Any:
    from rtt_fusion.graph.field import port_plane

    alive = [b for b in bundles if np.any(np.asarray(b.alive))]
    if not alive:
        return None, None
    biggest = max(alive, key=lambda b: int(np.asarray(b.alive).sum()))
    plane = port_plane(biggest)
    pts: list[FloatArray] = []
    for b in alive:
        moved, uv = plane.intersect(b)
        ok = np.asarray(moved.alive)
        pts.append(uv[ok])
    xy = np.concatenate(pts) if pts else np.zeros((0, 2))
    return plane, xy


def beam_view(ctx: Context, node: Node, socket: str, inputs: dict[str, Any]) -> View:
    from rtt_fusion.graph.field import FieldBundle, FieldError, field_image

    items = list(inputs.get("in") or [])
    if not items:
        raise ExecutionError("exec.display", "no light attached to the display", node.id)
    pixels = int(node.params.get("pixels", 128))
    margin = float(node.params.get("margin", 0.1))
    title = _title(node, "beam")
    fields = [x for x in items if isinstance(x, FieldBundle)]
    rays = [x for x in items if isinstance(x, RayBundle)]
    if fields:
        plane = fields[0].plane
        pts = np.concatenate(
            [np.asarray(getattr(f.field, "points", np.zeros((0, 2)))) for f in fields]
        )
        half = float(np.max(np.abs(pts))) * (1 + margin) if len(pts) else 1.0
        try:
            image = field_image(fields, plane, half, half, pixels, pixels, ctx.log)
        except FieldError as e:
            raise ExecutionError("exec.field", str(e), node.id) from e
        return View(
            "image",
            f"{title} (field)",
            image=image,
            x_label="u (mm)",
            y_label="v (mm)",
            unit="W/mm² (rel.)",
        )
    plane, xy = _plane_of(rays)
    if plane is None or xy is None or len(xy) == 0:
        raise ExecutionError("exec.display", "no living rays at this output", node.id)
    centre = xy.mean(axis=0)
    half = float(np.max(np.abs(xy - centre))) * (1 + margin) or 1.0
    shifted = type(plane)(
        plane.origin + centre[0] * plane.u + centre[1] * plane.v, plane.normal, plane.u, plane.v
    )
    image = coherent_image(rays, shifted, half, half, pixels, pixels)
    return View(
        "image",
        f"{title} (rays)",
        image=image,
        x_label="u (mm)",
        y_label="v (mm)",
        unit="W/mm² (rel.)",
    )


beam_view.representations = frozenset("RH")  # type: ignore[attr-defined]

_TITLE = ParamSpec("title", "text", "", description="window title (default: the block name)")

TYPES = (
    NodeType(
        "display.image",
        Category.ANALYSIS,
        (data_in("in", "result", "image"), data_out("view", "result", "view")),
        (_TITLE,),
        description="Display: an image result (camera, detector) in its own window",
    ),
    NodeType(
        "display.plot",
        Category.ANALYSIS,
        (data_in("in", "result", "result"), data_out("view", "result", "view")),
        (_TITLE,),
        description="Display: curves of a sweep / optimisation, or profiles of an image",
    ),
    NodeType(
        "display.beam",
        Category.ANALYSIS,
        (light_in("in", "RH"), data_out("view", "result", "view")),
        (
            _TITLE,
            ParamSpec("pixels", "integer", 128, minimum=16, maximum=2048),
            ParamSpec("margin", "number", 0.1, minimum=0.0, description="window margin"),
        ),
        (("view", ("in",)),),
        description="Display: the light at a block output (probe), intensity on its port plane",
    ),
)

EVALUATORS: dict[str, Evaluator] = {
    "display.image": image_view,
    "display.plot": plot_view,
    "display.beam": beam_view,
}


def register(registry: NodeTypeRegistry) -> dict[str, Evaluator]:
    for t in TYPES:
        registry.add(t)
    return dict(EVALUATORS)
