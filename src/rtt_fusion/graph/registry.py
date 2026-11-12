"""Registry of block types (concept sections 3.1 and 9, A1/A2).

Built-in types cover the Raytatouille side of the graph (source, sections, splitter, mirror,
detector) and the analyses; solver adapters contribute types from their `describe()`: with the
proposed keys `sockets` (A1) and `parameters` (A2), or, for today's adapters without them, derived
from `representation` and `ports` (front/back ports, incidence at the front only).
"""

from __future__ import annotations

from typing import Any

from rtt_fusion.graph.model import (
    Category,
    NodeType,
    ParamSpec,
    PortSpec,
    SocketSpec,
    data_in,
    data_out,
    light_in,
    light_out,
)

# two-sided domain, passed in both directions (concept section 4)
TWO_SIDED = (
    light_in("front_in", side="front", label="front in"),
    light_out("back_out", side="back", label="back out"),
    light_in("back_in", side="back", label="back in"),
    light_out("front_out", side="front", label="front out"),
)
TWO_SIDED_THROUGH = (("back_out", ("front_in",)), ("front_out", ("back_in",)))


class NodeTypeRegistry:
    """Block types by id; the node library of the editor lists them by category."""

    def __init__(self) -> None:
        self._types: dict[str, NodeType] = {}

    def add(self, node_type: NodeType) -> NodeType:
        self._types[node_type.type_id] = node_type
        return node_type

    def get(self, type_id: str) -> NodeType | None:
        return self._types.get(type_id)

    def __contains__(self, type_id: str) -> bool:
        return type_id in self._types

    def types(self, category: Category | None = None) -> list[NodeType]:
        return sorted(
            (t for t in self._types.values() if category is None or t.category is category),
            key=lambda t: t.type_id,
        )

    def add_adapter(self, describe: dict[str, Any], type_id: str | None = None) -> NodeType:
        """Registers a solver adapter as a domain block (from its `describe()`)."""
        return self.add(node_type_from_describe(describe, type_id))


def node_type_from_describe(describe: dict[str, Any], type_id: str | None = None) -> NodeType:
    """A domain block type from an adapter description.

    With `sockets` (A1): taken as given. Without: from `representation` ("G"/"S") – front_in,
    front_out (reflected) and back_out (transmitted); a back_in only if a port declares
    `accepts_incidence` (A3).
    """
    name = str(describe.get("name", "adapter"))
    tid = type_id or f"{name}.domain"
    params = tuple(ParamSpec.from_dict(p) for p in describe.get("parameters", ()))
    if "sockets" in describe:
        sockets = tuple(SocketSpec.from_dict(s) for s in describe["sockets"])
        through = tuple(
            (str(o), tuple(str(i) for i in ins)) for o, ins in describe.get("through", {}).items()
        )
    else:
        rep = str(describe.get("representation", "S"))
        ports = {
            str(p.get("name", "")): p for p in describe.get("ports", ()) if isinstance(p, dict)
        }
        back_incidence = bool(ports.get("back", {}).get("accepts_incidence", False))
        sockets = (
            light_in("front_in", rep, "front", "front in"),
            light_out("front_out", rep, "front", "reflected"),
            light_out("back_out", rep, "back", "transmitted"),
        )
        through = (("front_out", ("front_in",)), ("back_out", ("front_in",)))
        if back_incidence:
            sockets += (light_in("back_in", rep, "back", "back in"),)
            through = (
                ("front_out", ("front_in", "back_in")),
                ("back_out", ("front_in", "back_in")),
            )
    return NodeType(
        tid, Category.DOMAIN, sockets, params, through, description=str(describe.get("method", ""))
    )


def builtin_registry() -> NodeTypeRegistry:
    """The built-in block types."""
    r = NodeTypeRegistry()
    r.add(
        NodeType(
            "group",
            Category.DOMAIN,
            (),
            dynamic_sockets=True,
            description="Group (R28): a sub-model with variants; in the socket form it keeps the "
            "group's pose and exported parameters for the blocks inside",
        )
    )
    r.add(
        NodeType(
            "fusion.condensed",
            Category.DOMAIN,
            (),
            dynamic_sockets=True,
            description="Condensed group (R29): plane-wave tables of its inner branches",
        )
    )
    r.add(
        NodeType(
            "fusion.sum",
            Category.SPLITTER,
            (),
            dynamic_sockets=True,
            description="Junction of unrolled branches: the sum of the light arriving at one state",
        )
    )
    r.add(
        NodeType(
            "rtt.source",
            Category.SOURCE,
            (light_out("out", "RHGSP", label="light"),),
            (
                ParamSpec("system", "file", "", label="Raytatouille system"),
                ParamSpec("fields", "text", "", description="field indices; empty: all"),
                ParamSpec("wavelengths", "text", "", description="wavelength indices; empty: all"),
            ),
            description="Field points, wavelengths, polarisation of a system (coherence groups)",
        )
    )
    r.add(
        NodeType(
            "rtt.section",
            Category.DOMAIN,
            TWO_SIDED,
            (
                ParamSpec(
                    "location", "text", "", description="node of the Raytatouille model tree"
                ),
                ParamSpec(
                    "offset",
                    "point",
                    [0.0, 0.0, 0.0],
                    unit="mm",
                    description="shift of the elements of the block (Raytatouille backend)",
                ),
                ParamSpec(
                    "drive",
                    "number",
                    0.0,
                    unit="mm",
                    description="shift along the surface normal (set by a phase shifter's drive)",
                ),
            ),
            TWO_SIDED_THROUGH,
            dynamic_sockets=True,
            extra_sockets=(data_in("drive", "drive", "drive"),),
            description="Sequential surfaces of a Raytatouille system (e.g. an objective)",
        )
    )
    r.add(
        NodeType(
            "rtt.splitter",
            Category.SPLITTER,
            (
                light_in("in", "RHGS"),
                light_out("r", "RHGS", label="reflected"),
                light_out("t", "RHGS", label="transmitted"),
                light_in("r_in", "RHGS"),
                light_in("t_in", "RHGS"),
                light_out("merged", "RHGS", label="combined"),
            ),
            (
                ParamSpec("surface", "text", ""),
                ParamSpec(
                    "offset",
                    "point",
                    [0.0, 0.0, 0.0],
                    unit="mm",
                    description="shift of the elements of the block (Raytatouille backend)",
                ),
                ParamSpec(
                    "drive",
                    "number",
                    0.0,
                    unit="mm",
                    description="shift along the surface normal (set by a phase shifter's drive)",
                ),
            ),
            (("r", ("in",)), ("t", ("in",)), ("merged", ("r_in", "t_in"))),
            dynamic_sockets=True,
            extra_sockets=(data_in("drive", "drive", "drive"),),
            description="Beam splitter: splits into channels, combines coherently per group",
        )
    )
    r.add(
        NodeType(
            "rtt.mirror",
            Category.DOMAIN,
            (light_in("in", "RHGS"), light_out("out", "RHGS")),
            (
                ParamSpec("surface", "text", ""),
                ParamSpec(
                    "offset",
                    "point",
                    [0.0, 0.0, 0.0],
                    unit="mm",
                    description="shift of the elements of the block (Raytatouille backend)",
                ),
                ParamSpec(
                    "drive",
                    "number",
                    0.0,
                    unit="mm",
                    description="shift along the surface normal (set by a phase shifter's drive)",
                ),
            ),
            (("out", ("in",)),),
            dynamic_sockets=True,
            extra_sockets=(data_in("drive", "drive", "drive"),),
            ports=(PortSpec("front", "in", "out", "front"),),
            description="Mirror (one surface, reflecting)",
        )
    )
    r.add(
        NodeType(
            "rtt.detector",
            Category.DETECTOR,
            (light_in("in", "RHGSP"), data_out("image", "result", "image")),
            (
                ParamSpec("surface", "text", ""),
                ParamSpec(
                    "pixels", "integer", 512, minimum=16, maximum=8192, label="pixels (long side)"
                ),
                ParamSpec(
                    "half_size",
                    "numbers",
                    [],
                    unit="mm",
                    description="image window [hu, hv]; empty: the detector aperture",
                ),
            ),
            dynamic_sockets=True,
            description="Detector / camera: amplitudes within, intensities across coherence groups",
        )
    )
    r.add(
        NodeType(
            "fusion.convert",
            Category.CONVERSION,
            (light_in("in"), light_out("out")),
            (ParamSpec("to", "choice", "S", options=("R", "H", "G", "S", "P")),),
            (("out", ("in",)),),
            description="Conversion between representations (R → H, H ↔ G, G ↔ S, …)",
        )
    )
    r.add(
        NodeType(
            "analysis.camera_image",
            Category.ANALYSIS,
            (light_in("in", "RHGSP"), data_out("image", "result")),
            (
                ParamSpec(
                    "rays", "integer", 200_000, minimum=1000, label="rays per field and wavelength"
                ),
                ParamSpec("pixels", "integer", 512, minimum=16, maximum=8192),
            ),
            description="Irradiance image of a detector with many rays",
        )
    )
    r.add(
        NodeType(
            "biooptomd.layer",
            Category.DATA,
            (data_out("profile", "N"), data_out("dipoles", "Q")),
            (ParamSpec("directory", "file", "", label="BioOptoMD output"),),
            description="Refractive-index profile n(z) and dipole ensembles from BioOptoMD",
        )
    )
    return r
