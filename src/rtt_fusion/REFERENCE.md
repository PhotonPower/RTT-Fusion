# RTT-Fusion reference for agents

rtt-fusion 0.3.x (PyPI `rtt-fusion`, import `rtt_fusion`, MIT). This file ships in the wheel:
`python -m rtt_fusion.reference` (`--list`, `--section 2.4`, `--blocks fg.grating`),
`rtt_fusion.reference("planner")`. English; facts checked against the code of this version.

- **What:** optical set-ups as a node workflow: blocks (sources, lenses, mirrors, gratings,
  stacks, particles, detectors, meta blocks) joined by light links. Light on a link is rays (R),
  hybrid fields (H) or plane-wave spectra (S); a planner picks the representation per link, the
  executor converts at the ports. Rigorous solvers (FourierGrating RCWA, CliffordBem BEM) are blocks.
- **Core:** numpy + scipy only. Works out of the box: `optics.*`, `meta.*`, `display.*`, rays (R),
  planner, checks, file format, `fg.grating` from a stored Jones table.
- **Optional modules** (section 0.1; fieldtracing, fouriergrating and raytatouille on PyPI, the
  others from GitHub): `fieldtracing` → H and S
  (fields, spectra, `fusion.convert`, wave-optical images); `fouriergrating` → `fg.stack`,
  `fg.grating` from a configuration, `adapters.fouriergrating`; `cliffordbem` → `fusion.particle`;
  `raytatouille` → `rtt.*` blocks; `specklebem`, `hpfem`, `bwpm` → their adapters only. Binary
  wheels contain the C++ extension `rtt_fusion._diffraction` (faster `fg.grating` with rays).
- **Check at run time:** `rtt_fusion.graph.workflow.available()` → `{'R'}` or `{'R','H','S'}`
  (fieldtracing found); `importlib.util.find_spec("cliffordbem") is not None`; section 7.
- **Units:** mm, vacuum wavelength µm, rad (`*_deg` params: degrees), k in 1/mm, irradiance W/mm².

## Contents

0.1 Installing the optional modules · 1 Mental model (1.1 graph and model, 1.2 blocks, 1.3 light
and data types, 1.4 planner, 1.5 conversions at ports, 1.6 coherence and branches, 1.7 units,
frames, signs, 1.8 polarisation) · 2 Recipes, runnable (2.1 rays, 2.2 model, file, document,
2.3 plan inspection, 2.4 fields at a focus, 2.5 converting, 2.6 stack and grating with S,
2.7 particle, 2.8 interferometer and sweep, 2.9 Raytatouille systems, 2.10 an adapter directly,
2.11 a .rttf document) · 3 Block catalogue (generated) · 4 Public API by module ·
5 Diagnostics and log lines · 6 Gotchas and performance · 7 Discovering more at run time

## 0.1 Installing the optional modules

From PyPI (binary wheels; the import names stay `fieldtracing`, `fouriergrating`, `raytatouille`):

```bash
pip install rtt-fusion photonpower-fieldtracing fouriergrating raytatouille
```

The others only from GitHub. Below are the commits the CI of this release tests against
(`requirements-modules.txt` of the repository; the first three are the alternative to PyPI). Each
needs a C++ toolchain to build (FourierGrating: Eigen 3.4, fetched if missing).

```bash
pip install "fieldtracing @ git+https://github.com/PhotonPower/Fieldtracing.git@ea6da7b2d457ee8938fc079b21ed87d1b78504cf"
pip install "fouriergrating @ git+https://github.com/PhotonPower/FourierGrating.git@174ef93d4e3b41ece7bb2f9a70055909e247ee25"
pip install "cliffordbem @ git+https://github.com/PhotonPower/CliffordBem.git@41e64e627b82f96766c7cfe4f7d663864a5e5cfe"
pip install "raytatouille @ git+https://github.com/PhotonPower/Raytatouille.git@708b6dda655a4bd0f20fb9f34b51c61d51761590"
pip install "bwpm @ git+https://github.com/PhotonPower/BWPM.git@ef2de941d226d364db000a977963b0eba24aea44"
pip install "hpfem @ git+https://github.com/PhotonPower/hp-FEM.git@dce07d5d5cc52681715044310c03a10c088d20e7" gmsh  # BLAS/LAPACK
pip install "pybind11<3" && pip install --no-build-isolation "specklebem @ git+https://github.com/PhotonPower/SpeckleBem.git@4f04efe26b71f95122ff58cc7589641ba7e8b28f#subdirectory=python"
```

## 1 Mental model

### 1.1 Graph and model

- `Graph` (socket form, file version 1): `nodes` {id: `Node`}, directed `links` (`Link(source=(node,
  socket), target=(node, socket), transport)`), `meta`. The executor and the planner work on it.
  Acyclic on socket level (inside a block only the `through` dependencies count).
- `Model` (topology form, file version 2): nodes with **ports** (a port = a reference surface with
  an *in* and an *out* channel), undirected `connections` between two ports (`Connection(a=(node,
  port), b=(node, port), transport, transport_back)`), data/observer `links`, and **branches**:
  named step lists source → detector, `Step(node, port_in, port_out, channel)`. Only the branches
  given are computed (nothing is enumerated). Groups (`Group`/`Variant`) hold sub-models.
- `unroll(model, registry) -> Graph`: one copy of a node per visit (`id`, `id@2`, …) and
  `fusion.sum` junctions; socket paths = the branches. `model_from_graph` goes back.
- Values flow on demand: `Executor.evaluate(node, socket)` computes what that output needs;
  results are cached by a content hash (type, params, upstream keys, planned representation).

### 1.2 Blocks

- `Node(id, type_id, name="", params={}, position=(0,0), muted=False, sockets=(), through=(),
  ports=())`; type from a `NodeTypeRegistry` (`NodeType`: category, sockets, params, through,
  ports, scatter, port_frames; category source, domain, splitter, detector, conversion, data,
  analysis, group or meta). Missing params take the `ParamSpec.default`.
- Light sockets accept a subset of R H G S T P; data sockets one type of `DATA_TYPES` (N, Q,
  field_map, occupancy, voxel_materials, shape, height_map, response_map, scattering_stats,
  result, drive). Every source/domain/splitter/detector block also has the data input `pose`.
- An evaluator `f(ctx, node, socket, inputs) -> value` computes one output; light values are
  lists of bundles. Attributes on it: `representations` (native, default R), `produces(node)`,
  `key_extra`, `uncached`. The node it gets has world coordinates (pose applied, section 1.7).
- Muted block: light passes unchanged, data outputs give None.

### 1.3 Light and data types

| Rep | On a link (Python type) | Module | Needs |
|---|---|---|---|
| R | `RayBundle`: one coherence group of rays (position, direction, OPL mm, complex amplitude, `prt`, `jones`, `alive`) | `graph.light` | – |
| H | `FieldBundle`: Fieldtracing `GridlessField`/`GriddedField` on a port `Plane` (u, v, n frame), phase ψ = k₀·OPL | `graph.field` | fieldtracing |
| S | `SpectrumBundle`: `ports.PlaneWaves` (kx, ky 1/mm in the plane frame; field (N, 2) transverse; `solid_angle` dΩ if continuous; `orders`) on a `Plane`, phase reference at its origin; `support` = sampling for S → H | `graph.spectrum` | fieldtracing |
| G | `ports.GriddedField` (ny, nx, 2) on a centred grid | `ports` | adapters only; not run on links |
| P | rays without phase (incoherent) | planner only | not executable |
| T | T-matrix | planner only | not executable |

- Data: `imaging.Image` (detector result: `values[j, i]` at (u_i, v_j) W/mm², `extent` (u0, u1,
  v0, v1) mm, `power` W, `ray_power`, `groups`, `note`), `meta.MetaResult`, `display.View`,
  `drive` (a number set by a meta block).
- Tables: `couplings.jones_table.JonesTable` (power-normalised Jones matrices per order over λ, σ;
  file format `raytatouille-result`), `graph.condense.Table` (plane-wave table of a condensed
  group), FourierGrating configurations (dicts, µm).

### 1.4 Planner

`plan(graph, registry, Requirements(available, outputs, max_paths, validity, native, produces))`
→ `Plan` (representation and reason per link key `"a.out -> b.in"`, conversions inside blocks,
fusion paths, diagnostics, cost). `Executor.check()` plans before every run (`Executor.plan`).

1. A fixed transport is kept (error if not available: `plan.representation_unavailable`).
2. "auto": admissible = accepted by both sockets ∩ `available` ∩ validity hook; then the cheapest
   by ranking R 1 < P 1.5 < H 3 < S 4 < G 6 < T 8 plus conversion operators (rays_to_field R→H 1,
   field_to_rays 1, pft H→S 2, orders_to_rays S→R 1, fft/inverse_fft G↔S 2, …); never P while a
   coherent one is admissible.
3. S only between two blocks that both take S natively (`Requirements.native`; decision 3).
4. **R30f:** an "auto" link between two S-native blocks is **S outright** when the light leaving
   the source block is a spectrum: the block makes S whatever comes in (`produces`:
   `fusion.convert` with `to="S"`), or it passes what it gets and *all* its inputs arrive on S
   links (fixed S, only S admissible, or this rule). Propagates along chains of native blocks
   (grating, stack, condensed group, particle, detectors). Reason text: `"spectra from <id>; S
   kept between blocks that take it natively (among RHS)"`.
5. Rays never become spectra on their own: rays → grating → stack → detector stays R. S appears
   where you ask (fixed S on a link, `fusion.convert(to="S")`, a source that emits S). Not applied
   with mixed inputs, before a block without S (lens, mirror, `rtt.*` section), when validity
   rejects S, or without `native` (abstract planning).
6. With only R available (the `Executor` default!) every "auto" link is R.

### 1.5 Conversions at ports (`Executor._convert`)

- Link planned H: rays from the source become a field on the port plane (R → H, plane through the
  centroid ⟂ mean direction); planned S: R → H → S or H → S (PFT if the local Fresnel number is
  large, else FFT on a grid).
- At the consumer: if it does not take H natively it gets rays back (H → R); a spectrum goes to
  S → H if the consumer takes H, else S → H → R; a consumer that takes S keeps it on its plane.
- A fixed transport R on the link forces rays at the consumer; a block that got a field returns
  one. Fields and spectra are not moved between port planes by the executor: the consuming block
  propagates them (`spectrum.propagate`, Fieldtracing's planner).
- `Context.used[link_key]` records what really flowed (`representation`, `chain` e.g. `"R→H→S"`,
  `reason`); `workflow.with_used` merges it into the plan.

### 1.6 Coherence and branches

- `CoherenceKey(source, field_index, wavelength_index, polarisation)`: amplitudes add within a
  key, intensities across keys. Diffraction orders keep the key of the incident bundle.
- `optics.collimated_source`: one bundle per field angle × wavelength (polarisation from the
  param). `rtt.source` (unpolarised Raytatouille source): two orthogonal polarisations = two
  groups with half the power.
- Splitters split into channels and combine coherently per group; a detector sums amplitudes of
  all bundles of one key, then intensities.
- Branches name the way (`Branch(name, steps)`); a step through a group names the inner branch as
  `channel`. `RayBundle.branch` carries (branch, position) for node-wise Raytatouille execution.

### 1.7 Units, frames, signs

- Lengths mm; vacuum wavelengths µm (`wavelengths_um`, `wavelength_um`); OPL mm vacuum-equivalent;
  k₀ = 2π/(λ·1e-3) 1/mm; angles rad except `field_angles_deg`, Raytatouille `rotation_deg`;
  FourierGrating configs µm; `fusion.particle.radius_um` µm; irradiance W/mm², powers W.
- Right-handed, main direction +z; fields ∝ exp(i(k·r − ωt)); k_z = √(n²k₀² − k⊥²) with Im ≥ 0.
- Fresnel Convention A: s = k_in × N/|…|, p = k̂ × s; at normal incidence r_p = −r_s. Plane
  mirror default `reflectance_amplitude = -1`.
- `Plane.make(origin, normal, up=(0,1,0))`: u = up × n, v = n × u. A detector with normal
  (0, 0, −1) has u = −x, v = +y (image seen from the light). Image `values[j, i]` ↔ (u_i, v_j).
- Normals: `optics.detector` and `optics.plane_mirror` `normal` points back towards the incoming
  light (default (0, 0, −1) for light along +z); `optics.interface`: `n_front` on the side `normal`
  points to; `fg.stack`: `normal` towards the front (superstrate) side; `fg.grating`: `normal`,
  `grating_vector` (local x, ⟂ grooves). Port frames: outward normal, light leaves the node
  along it; two connected ports must face each other (`graph.port_facing`).
- Ports of the S representation: phase 1 at the plane origin; `PlaneWaves.ez()` from transversality.
- Poses: `params["pose"] = {"position": [mm], "rotation": [rad], "pivot": [mm]}`, intrinsic
  X→Y→Z (R = Rx·Ry·Rz), p_parent = position + pivot + R·(p − pivot); offsets `pose_offset`
  (dx, dy, dz, rx, ry, rz) plus `pose_by_<driver>`; points (`point`, `center`) and directions
  (`normal`, `direction`, `axis`, `up`, `grating_vector`, `polarisation`) are transformed.

### 1.8 Polarisation (R31)

- Field of a ray: `a · prt · jones`; `amplitude` a (N,) complex, `prt` (N, 3, 3) complex (None =
  identity), `jones` (3,) complex source field direction. `RayBundle.field_vectors()` gives it,
  `field_direction()` gives prt·jones. Real `prt`/`jones` (older port data) are cast to complex128.
- Elliptical polarisation survives every conversion: `light.split_field(f) -> (a, ê)` with
  |a| = |f|, ê a complex unit vector, phase along the major axis of the ellipse.
- `optics.interface` and condensed groups carry the full Fresnel/table matrix in prt (mixed s/p,
  TIR phase exact); mirrors multiply prt by I − 2nnᵀ; gratings and stacks use their PRT matrices.
- `optics.collimated_source.polarisation` may be complex in Python (`[1, 1j, 0]`), not in JSON.
- The R image (`imaging.coherent_image`) is scalar (|a| and phase); H and S images are vectorial.

## 2 Recipes

Blocks marked `python run` are executed by the package's tests; `needs=` names the optional
modules they require. Each block is self-contained.

### 2.1 Rays: build, run, read

```python run
from rtt_fusion.graph import Graph, Node, validate, workflow
from rtt_fusion.graph.execute import Executor

reg, ev = workflow.registry()  # built-in + optional-module block types, and their evaluators
g = Graph()
g.add(Node("src", "optics.collimated_source",
           params={"center": [0, 0, -5], "radius": 1.0, "rings": 10, "wavelengths_um": [0.55]}))
g.add(Node("lens", "optics.thin_lens", params={"point": [0, 0, 0], "focal_length": 50.0}))
g.add(Node("cam", "optics.detector",
           params={"point": [0, 0, 40], "normal": [0, 0, -1], "half_size": [0.3, 0.3],
                   "pixels": [64, 64]}))
g.connect(("src", "out"), ("lens", "in"))           # (node, socket) -> (node, socket)
g.connect(("lens", "out"), ("cam", "in"))
assert validate(g, reg) == []                        # list of Diagnostic
ex = Executor(g, reg, ev)                            # available defaults to {"R"}
img = ex.evaluate("cam", "image")                    # imaging.Image
print(img.values.shape, img.extent, img.power, img.note)  # (64, 64) … ≈ π W (1 W/mm², r = 1 mm)
assert abs(img.power - 3.14) < 0.05
print(ex.plan.to_dict()["links"])  # {'src.out -> lens.in': {'representation': 'R', ...}, ...}
print(ex.results().keys())         # every data output of detectors and analyses
print(ex.context.log)              # notes of the run (conversions, propagation methods, ...)
```

### 2.2 Model with ports and branches; file; workflow document; .rttf container

```python run
import json
from rtt_fusion.graph import Model, Node, io, unroll, validate_model, workflow
from rtt_fusion.graph.execute import Executor

reg, ev = workflow.registry()
m = Model()
m.add(Node("src", "optics.collimated_source", params={"center": [0, 0, -5], "radius": 1.0}))
m.add(Node("lens", "optics.thin_lens", params={"point": [0, 0, 0], "focal_length": 50.0}))
m.add(Node("cam", "optics.detector",
           params={"point": [0, 0, 40], "normal": [0, 0, -1], "half_size": [0.3, 0.3]}))
m.connect(("src", "out"), ("lens", "in"))             # (node, PORT) – see "ports" in section 3
m.connect(("lens", "out"), ("cam", "in"))
m.branch("main", ("src", None, "out"), ("lens", "in", "out"), ("cam", "in", None))
assert validate_model(m, reg) == []
g = unroll(m, reg)                                    # socket graph for the executor
img = Executor(g, reg, ev, available=workflow.available()).evaluate("cam", "image")

text = json.dumps(io.model_to_dict(m), indent=1)      # rtt-fusion-graph version 2
m2 = io.model_from_dict(json.loads(text), reg)        # also reads version 1
doc = workflow.model_document(m2, reg)                # what a front end stores and shows
assert set(doc) == {"model", "types", "diagnostics", "plan", "providers"}
print(doc["plan"]["links"], doc["plan"]["connections"])  # per link / per connection+direction

# Raytatouille Fusion document (.rttf, JSON "nfcad" v2): the workflow document under "workflow"
rttf = {"format": "nfcad", "version": 2, "units": "mm", "layers": [], "groups": [],
        "entities": [], "plugins": {}, "workflow": doc}
with open("demo.rttf", "w", encoding="utf-8") as f:
    json.dump(rttf, f)
with open("demo.rttf", encoding="utf-8") as f:
    back = io.model_from_dict(json.load(f)["workflow"]["model"], reg)
assert set(back.nodes) == {"src", "lens", "cam"}
```

File format `rtt-fusion-graph` (`graph.io`):

```text
v1 (Graph):  {"format": "rtt-fusion-graph", "version": 1, "meta": {},
  "nodes": [{"id", "type", "name"?, "params"?, "position"?: [x, y], "muted"?, "sockets"?, "through"?}],
  "links": [{"from": [node, socket], "to": [node, socket], "transport"?: "auto"|"R"|"H"|"S"|...}]}
v2 (Model):  {"format": "rtt-fusion-graph", "version": 2, "meta": {},
  "nodes": [... as v1, "ports"?: [{"name", "in", "out", "side"}]],
  "connections": [{"ports": [[node, port], [node, port]], "transport"?: "R" | {"forward": "R", "back": "H"}}],
  "links": [data/observer links as v1, "branch"?: name],
  "branches": [{"name", "steps": [[node, port_in|null, port_out|null, channel?]]}]}
group node:  {"id", "type": "group", "variant": active, "variants": {name: {"nodes", "connections",
  "links", "branches", "ports": {port: [node, port]}, "params": {name: "node.param[i]"}}},
  "condense"?: {"sigma_max": 0.2, "count": 5}}
```

`.rttf` (Raytatouille Fusion CAD file, `"format": "nfcad"`): `entities` (CAD geometry), `plugins`
{plugin id: {"params": …}} (e.g. workflow: rays, pixels, gaussian_waist), `workflow` = the
document of `workflow.model_document` (`model` v2, `types` = block library, `diagnostics`,
`plan`, `providers`), optional `results` (saved display images). Old documents may hold
`workflow.graph` (v1) instead of `model`: read with `io.from_dict`. Rename `fusion.grating` →
`fg.grating` happens on load (`io.TYPE_ALIASES`).

### 2.3 Plan inspection

```python run
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import Executor

reg, ev = workflow.registry()
g = Graph()
g.add(Node("src", "optics.collimated_source", params={"center": [0, 0, -1], "radius": 0.5}))
g.add(Node("cam", "optics.detector", params={"point": [0, 0, 5], "normal": [0, 0, -1]}))
link = g.connect(("src", "out"), ("cam", "in"))
diagnostics, plan = workflow.checked(g, reg)   # graph check + plan as a front end sees it
print(plan["links"])         # {"src.out -> cam.in": {"representation": "R", "reason": "lowest cost ..."}}
print(plan["conversions"], plan["cost"], diagnostics)
g.set_transport(link, "R")   # fix a transport; "auto" lets the planner choose
ex = Executor(g, reg, ev, available=workflow.available())
ex.check()                   # validates, plans, raises ExecutionError on errors
print(ex.plan.representations, ex.plan.reasons, ex.plan.conversions)
ex.evaluate("cam", "image")
shown = workflow.with_used(ex.plan.to_dict(), ex.context.used)  # + "used", "used_chain"
print({k: v.get("used") for k, v in shown["links"].items()})
```

### 2.4 Fields (H) at a focus

A ray image at a focus is empty (all rays meet in a point, `Image.note` says so); use H on the
link into the detector.

```python run needs=fieldtracing
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import Executor

reg, ev = workflow.registry()
g = Graph()
g.add(Node("src", "optics.collimated_source",
           params={"center": [0, 0, -5], "radius": 1.0, "rings": 12, "wavelengths_um": [0.55]}))
g.add(Node("lens", "optics.aplanatic_lens",
           params={"point": [0, 0, 0], "focal_length": 10.0, "aperture_radius": 2.0}))
g.add(Node("cam", "optics.detector", params={"point": [0, 0, 10], "normal": [0, 0, -1],
                                             "half_size": [0.01, 0.01], "pixels": [32, 32]}))
g.connect(("src", "out"), ("lens", "in"))
g.connect(("lens", "out"), ("cam", "in"), "H")       # fixed transport H
ex = Executor(g, reg, ev, available=workflow.available())   # {"R", "H", "S"}
img = ex.evaluate("cam", "image")                    # Airy pattern, NA 0.1
print(ex.plan.to_dict()["conversions"])   # [{'node': 'lens', 'chain': 'RH', 'operators': ['rays_to_field'], ...}]
print(ex.context.log)                     # ['src: H propagation 9.97 mm by debye (...)']
assert img.values.max() == img.values[16, 16] and img.power > 0.9 * 3.14
```

### 2.5 Converting between representations

```python run needs=fieldtracing
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph import field as fl, spectrum as sp
from rtt_fusion.graph.execute import Executor

reg, ev = workflow.registry()
g = Graph()
g.add(Node("src", "optics.collimated_source",
           params={"center": [0, 0, -1], "radius": 0.05, "gaussian_waist": 0.02}))
g.add(Node("toS", "fusion.convert", params={"to": "S"}))   # explicit conversion block
g.add(Node("cam", "optics.detector", params={"point": [0, 0, 5], "normal": [0, 0, -1],
                                             "half_size": [0.1, 0.1], "pixels": [32, 32]}))
g.connect(("src", "out"), ("toS", "in"))
g.connect(("toS", "out"), ("cam", "in"))
ex = Executor(g, reg, ev, available=workflow.available())
img = ex.evaluate("cam", "image")
assert ex.plan.representations["toS.out -> cam.in"] == "S"   # R30f: S behind a spectrum
sb = ex.evaluate("toS", "out")[0]                  # SpectrumBundle (continuous: has dΩ)
rays = ex.evaluate("src", "out")[0]                # RayBundle
fb = fl.rays_to_field(rays)                        # R -> H (port plane from the rays)
s2 = sp.field_to_spectrum(fb)                      # H -> S (method "auto": PFT or FFT)
h2 = sp.spectrum_to_field(s2)                      # S -> H
r2 = sp.spectrum_to_rays(s2)                       # S -> H -> R
r3 = fl.field_to_rays(fb)                          # H -> R
print(fl.field_power(fb), sp.spectrum_power(s2), fl.field_power(h2), img.power)  # all ≈ 6.3e-4 W
assert abs(sp.spectrum_power(sb) / fl.field_power(fb) - 1) < 0.01
```

G ↔ S on ports (no fieldtracing): `ports.to_plane_waves(GriddedField) -> PlaneWaves`,
`ports.to_gridded(PlaneWaves, nx, ny, dx, dy)`; powers `ports.power`, `ports.spectrum_power`.
`fusion.convert` with `to` G or P: `exec.convert_unsupported`.

### 2.6 Stack and grating with spectra (S)

```python run needs=fieldtracing,fouriergrating
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.spectrum import spectrum_power

reg, ev = workflow.registry()
film = {  # FourierGrating configuration (µm); homogeneous layers, listed from the back (substrate)
    "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 0}, "period_um": {"x": 1.0},
    "materials": {"air": {"n": 1.0}, "film": {"n": 2.0}},
    "superstrate": "air", "substrate": "air",
    "layers": [{"type": "homogeneous", "thickness_um": 0.1, "material": "film"}],
}
g = Graph()
g.add(Node("src", "optics.collimated_source",
           params={"center": [0, 0, -1], "radius": 1.0, "rings": 12, "wavelengths_um": [0.55]}))
g.add(Node("lens", "optics.aplanatic_lens",
           params={"point": [0, 0, 0], "focal_length": 10.0, "aperture_radius": 2.0}))
g.add(Node("toS", "fusion.convert", params={"to": "S"}))
g.add(Node("film", "fg.stack", params={"config": film, "point": [0, 0, 5], "normal": [0, 0, -1]}))
g.add(Node("cam", "optics.detector", params={"point": [0, 0, 10], "normal": [0, 0, -1],
                                             "half_size": [0.01, 0.01], "pixels": [32, 32]}))
g.connect(("src", "out"), ("lens", "in"))
g.connect(("lens", "out"), ("toS", "in"))
g.connect(("toS", "out"), ("film", "front_in"))      # light in through the front
g.connect(("film", "back_out"), ("cam", "in"))       # transmitted through the back
ex = Executor(g, reg, ev, available=workflow.available())
img = ex.evaluate("cam", "image")
print({k: v for k, v in ex.plan.representations.items()})  # toS→film S, film→cam S (R30f)
assert ex.plan.representations["film.back_out -> cam.in"] == "S"

grating = {  # binary grating, period 2 µm, on glass
    "dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 8}, "period_um": {"x": 2.0},
    "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}}, "superstrate": "air",
    "substrate": "glass",
    "layers": [{"type": "lamellar1d", "thickness_um": 0.3, "breakpoints_um": [0, 1.0, 2.0],
                "materials": ["glass", "air"]}],
}
g.remove("film"); g.remove("cam")
g.add(Node("gr", "fg.grating", params={"config": grating, "side": "reflected",
                                       "point": [0, 0, 9.0], "normal": [0, 0, -1],
                                       "field_transfer": "interpolated"}))
g.connect(("toS", "out"), ("gr", "in"))
ex = Executor(g, reg, ev, available=workflow.available())
orders = ex.evaluate("gr", "out")                    # one SpectrumBundle per order
incident = spectrum_power(ex.evaluate("toS", "out")[0])
print({int(o.waves.orders[0][0]): spectrum_power(o) / incident for o in orders})  # m -> efficiency
assert 0 < sum(spectrum_power(o) for o in orders) < incident
```

`fg.grating` alternatives: `table` (path) / `table_data` (a `raytatouille-result` document of a
`JonesTable`) for rays without the solver; `orders` restricts the orders. With rays and a
configuration the block builds its table around the incident rays (`sigma_step`, `sigma_margin`).

### 2.7 Particle (CliffordBem) in a focus

```python run needs=cliffordbem
import numpy as np
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import Executor
from rtt_fusion.graph.particle import ray_powers

reg, ev = workflow.registry()
g = Graph()
g.add(Node("src", "optics.collimated_source",
           params={"center": [0, 0, -1], "radius": 1.0, "rings": 8, "wavelengths_um": [0.633]}))
g.add(Node("lens", "optics.aplanatic_lens",
           params={"point": [0, 0, 0], "focal_length": 2.0, "aperture_radius": 1.0}))
g.add(Node("p", "fusion.particle", params={"point": [0, 0, 2.0], "radius_um": 0.2, "index": 1.5,
                                           "refinement": 2, "n_theta": 12, "n_phi": 24}))
g.connect(("src", "out"), ("lens", "in"))
g.connect(("lens", "out"), ("p", "in"))   # rays in, backscattered rays out (Debye way)
ex = Executor(g, reg, ev, available=workflow.available())
back = ex.evaluate("p", "out")[0]         # RayBundle leaving the particle centre
inc = ex.evaluate("lens", "out")[0]
ratio = ray_powers(back, np.asarray(back.alive)).sum() / ray_powers(inc, np.asarray(inc.alive)).sum()
print(ratio, ex.context.log)              # backscattered share (≈ 0.024 here)
assert 0 < ratio < 1
```

With S on the link into the particle (`g.connect(..., "S")`, needs fieldtracing) it takes the
spectrum as it is and gives a spectrum (100 wavelengths in front of the particle).

### 2.8 Interferometer, drive and sweep

```python run
import numpy as np
from rtt_fusion.graph import Graph, Node, workflow
from rtt_fusion.graph.execute import Executor

reg, ev = workflow.registry()
g = Graph()
g.add(Node("src", "optics.collimated_source", params={"center": [0, 0, -10], "radius": 2.0}))
g.add(Node("bs", "optics.beam_splitter", params={"point": [0, 0, 0], "normal": [-1, 0, 1]}))
g.add(Node("mA", "optics.plane_mirror", params={"point": [0, 0, 20], "normal": [0, 0, -1]}))
g.add(Node("mB", "optics.plane_mirror", params={"point": [20, 0, 0], "normal": [-1, 0, 0]}))
g.add(Node("cam", "optics.detector", params={"point": [-30, 0, 0], "normal": [1, 0, 0],
                                             "half_size": [2.5, 2.5], "pixels": [32, 32]}))
for a, b in [(("src", "out"), ("bs", "in")), (("bs", "t"), ("mA", "in")),
             (("mA", "out"), ("bs", "t_in")), (("bs", "r"), ("mB", "in")),
             (("mB", "out"), ("bs", "r_in")), (("bs", "merged"), ("cam", "in"))]:
    g.connect(a, b)
lam = 0.6328e-3  # mm (default wavelength 0.6328 µm)
g.add(Node("sweep", "meta.sweep", params={"parameter": "mB.drive", "collect": "cam.image",
                                          "measure": "centre", "values": [0, lam / 8, lam / 4]}))
ex = Executor(g, reg, ev)
res = ex.evaluate("sweep", "results")      # meta.MetaResult
print(res.values, np.round(res.measures, 3))   # [0, 7.9e-05, 1.6e-04] [1. 0.5 0.]
assert np.allclose(res.measures, [1.0, 0.5, 0.0], atol=1e-6)
print(ex.calls)   # evaluator calls per output: only blocks behind mB were recomputed
```

Instead of `parameter`, link the meta block: `g.connect(("sweep", "drive"), ("mB", "drive"))` and
`g.connect(("cam", "image"), ("sweep", "observe"))`; with a `pose` input the drive moves the
offset `component` (dx … rz). `meta.replicates` averages images; `meta.optimise` searches
`bounds` for `goal` maximum/minimum.

### 2.9 Raytatouille systems (`rtt.*` blocks)

```python run needs=raytatouille
from rtt_fusion.graph import backend_rtt, new_system, unroll, workflow
from rtt_fusion.graph.execute import Executor

reg, ev = workflow.registry()
m = new_system("demo")             # rtt.source -> rtt.part.stop -> rtt.part.detector, branch "main"
g = unroll(m, reg)
evaluators = {**backend_rtt.EVALUATORS, **ev}        # rtt.* evaluators are not in registry()
shared = {"rings": 8, "pixels": [64, 64]}           # pupil rings, image size of rtt detectors
ex = Executor(g, reg, evaluators, shared, available=workflow.available())
img = ex.evaluate("detector", "image")
print(img.extent, img.power)       # window = detector aperture (±5 mm), EPD 4 mm
```

Import a system (dict, JSON text, `.rtt.json` path or `rt.System`): `model_from_raytatouille(s)`
→ `Model` (blocks named after their neighbours, system paths as branches), or
`graph_from_raytatouille(s)` → `Graph`. The blocks carry the system (`rtt_system`/`rtt_tree`/
`rtt_paths` on `rtt.source`, `rtt`/`runs` on the others); only graphs without these need
`shared["system"]`. Other `shared` keys: `"materials"` (`rt.MaterialLibrary`;
`embedded.load_catalogs(lib, source.params)` adds the stored glass catalogues), `"rings"` (30),
`"pixels"` ([256, 256]), `"gaussian_waist"` (relative to the pupil radius), `"auto_fields"`
(detector switches to H near a focus), `"rtt_execution"` ("auto"|"paths"|"nodes"),
`"field_tolerance"`. Swap a block, keeping its links: `replace_block(graph, reg, "test_mirror",
Node("gr", "fg.grating", params=...))`.

### 2.10 An adapter directly (plane waves on ports)

```python run needs=fouriergrating
import numpy as np
from rtt_fusion.adapters import fouriergrating as fga
from rtt_fusion.ports import PlaneWaves

cfg = {"dimension": "1D", "solver": "rcwa", "truncation": {"Nx": 0}, "period_um": {"x": 1.0},
       "materials": {"air": {"n": 1.0}, "glass": {"n": 1.5}},
       "superstrate": "air", "substrate": "glass", "layers": []}   # one interface air|glass
h = fga.prepare(fga.Structure(cfg), [0.55])
port = fga.front_port(h, 0.55)                     # Port(z=0, direction=+1, index 1, λ)
inc = PlaneWaves(port, np.array([0.0]), np.array([0.0]), np.array([[1.0 + 0j, 0.0]]))
out = fga.apply(h, inc)                            # {"front": reflected, "back": transmitted}
r = out["front"].field[0, 0]
print(r, out["front"].orders)                      # ≈ -0.2 (Convention A, phase at the port)
assert abs(r - (1 - 1.5) / (1 + 1.5)) < 1e-9
```

All adapters: `describe()`, `prepare(structure, wavelengths_um) -> Handle`, `apply(handle,
incident)`, `apply_many(handle, incidents)`, `estimate(handle)`; outputs `"front"` (reflected,
direction −1) and `"back"` (transmitted, +1).

### 2.11 Run a saved `.rttf` document

```python
import json
from rtt_fusion.graph import backend_rtt, embedded, io, unroll, workflow
from rtt_fusion.graph.execute import Executor

reg, ev = workflow.registry()
wf = json.load(open("set-up.rttf", encoding="utf-8"))["workflow"]
model = io.model_from_dict(wf["model"], reg)
graph = unroll(model, reg)
shared = {"rings": 30, "pixels": [256, 256]}
if any(n.type_id.startswith("rtt.") for n in graph.nodes.values()):
    import raytatouille as rt
    lib = rt.MaterialLibrary()
    for n in graph.nodes.values():
        if n.type_id == "rtt.source":
            embedded.load_catalogs(lib, n.params)     # glass catalogues stored in the document
    shared["materials"] = lib
ex = Executor(graph, reg, {**backend_rtt.EVALUATORS, **ev}, shared, available=workflow.available())
results = ex.results()          # {(node, socket): value} of every detector/analysis output
```

## 3 Block catalogue (generated)

Generated from the registry by `tools/gen_llm_reference.py`; do not edit by hand.
`python -m rtt_fusion.reference --blocks --live` prints it from the installed registry.

<!-- BEGIN GENERATED: block catalogue (python tools/gen_llm_reference.py) -->
33 block types. Per type: category and description; needs (optional packages); native = representations the evaluator computes in (the executor converts the rest at its ports); sockets as `name` in|out light:<accepted representations>[@side] or data:<type> (blocks of category source/domain/splitter/detector also get the data input `pose`, type drive); ports as name(in socket→out socket), the names `Model.connect` and branch steps use; through = output ← inputs it depends on; `param name type = default [unit] range of options – description`.

### analysis.camera_image
analysis: Irradiance image of a detector with many rays.
- needs: no evaluator in this package (planned type: `exec.no_evaluator`); native: – (no evaluator)
- sockets: `in` in light:GHPRS; `image` out data:result
- ports: in(in→·)
- param `rays` integer = 200000 range 1000… – rays per field and wavelength
- param `pixels` integer = 512 range 16…8192

### biooptomd.layer
data: Refractive-index profile n(z) and dipole ensembles from BioOptoMD.
- needs: no evaluator in this package (planned; data via `rtt_fusion.adapters.biooptomd`); native: – (no evaluator)
- sockets: `profile` out data:N; `dipoles` out data:Q
- param `directory` file = "" – BioOptoMD output

### display.beam
analysis: Display: the light at a block output (probe), intensity on its port plane.
- needs: nothing; fields on its input: fieldtracing; native: HR
- sockets: `in` in light:HR; `view` out data:result
- ports: in(in→·)
- through: view ← in
- param `title` text = "" – window title (default: the block name)
- param `pixels` integer = 128 range 16…2048
- param `margin` number = 0.1 range 0… – window margin

### display.image
analysis: Display: an image result (camera, detector) in its own window.
- needs: nothing; native: R
- sockets: `in` in data:result; `view` out data:result
- param `title` text = "" – window title (default: the block name)

### display.plot
analysis: Display: curves of a sweep / optimisation, or profiles of an image.
- needs: nothing; native: R
- sockets: `in` in data:result; `view` out data:result
- param `title` text = "" – window title (default: the block name)

### fg.grating
splitter: Periodic sample (FourierGrating table): diffraction orders as channels.
- needs: `table`/`table_data` (Jones table): nothing for rays, the optional C++ extension `rtt_fusion._diffraction` speeds them up; `fouriergrating`/`config`: fouriergrating; H or S: fieldtracing and a configuration (a table turns fields and spectra into rays); native: HRS
- sockets: `in` in light:HRS; `out` out light:HRS
- ports: front(in→out)
- through: out ← in
- param `table` file = "" – Jones table (raytatouille-result file); Check workflow stores it in the block (table_data) and clears the path
- param `fouriergrating` file = "" – FourierGrating configuration: the block builds its table around the rays (R) or computes every plane wave of a field (H); Check workflow stores it in the block (config) and clears the path
- param `config` data = {} – the FourierGrating configuration held by the block (before the file)
- param `table_data` data = {} – the Jones table held by the block (before the file)
- param `side` choice = "reflected" of reflected|transmitted – the orders passed on (with a configuration)
- param `orders` numbers = [] – orders to pass on; empty: all of the table
- param `n_in` number = 1.0 – index of the incident medium (R, built)
- param `sigma_step` number = 0.01 range 1e-05… – σ step of the table; spacing of the Jones table's nodes in the direction cosines σ = n·sin θ of the incident rays (cubic interpolation between them): smaller is more exact and needs more FourierGrating solves
- param `sigma_margin` number = 0.01 range 0… – σ margin of the table; how far the table reaches beyond the σ range of the incident rays
- param `point` point = [0.0, 0.0, 0.0] [mm]
- param `normal` point = [0.0, 0.0, 1.0]
- param `grating_vector` point = [1.0, 0.0, 0.0] – local x, ⟂ grooves
- param `field_transfer` choice = "exact" of exact|interpolated – exact: FourierGrating per plane wave; interpolated: from a verified sub-grid of the spectrum (1e-8), much faster for wide beams
- param `field_threshold` number = 1e-10 – plane waves below this share of the largest amplitude are left out (H)
- param `field_grid_max` integer = 1024 – at most this many grid points per axis for a field (H)

### fg.stack
domain: Layer stack (FourierGrating, exact per ray): coating, filter, mirror layers.
- needs: fouriergrating; H or S: also fieldtracing; native: HRS
- sockets: `front_in` in light:HRS@front; `front_out` out light:HRS@front; `back_in` in light:HRS@back; `back_out` out light:HRS@back
- ports: front(front_in→front_out); back(back_in→back_out)
- through: front_out ← front_in, back_in; back_out ← front_in, back_in
- param `fouriergrating` file = "" – FourierGrating configuration of the stack: homogeneous layers, listed from the substrate (back) upwards to the superstrate (front); Check workflow stores it in the block (config) and clears the path
- param `config` data = {} – the FourierGrating configuration held by the block (before the file)
- param `point` point = [0.0, 0.0, 0.0] [mm]
- param `normal` point = [0.0, 0.0, 1.0] – towards the front (superstrate)
- param `field_threshold` number = 1e-10 – plane waves below this share of the largest amplitude are left out (H)
- param `field_transfer` choice = "exact" of exact|interpolated – exact: FourierGrating per plane wave; interpolated: from a verified sub-grid of the spectrum (1e-8), much faster for wide beams
- param `field_grid_max` integer = 1024 – at most this many grid points per axis for a field (H)

### fusion.condensed
domain: Condensed group (R29): plane-wave tables of its inner branches.
- needs: made by `unroll` from a group with `condense`; what its inner blocks need; S: fieldtracing; native: RS
- sockets: + dynamic (per node, from the import or `unroll`)

### fusion.convert
conversion: Conversion between representations (R → H, H ↔ G, G ↔ S, …).
- needs: fieldtracing for every real conversion (`exec.convert_unsupported` without it); native: HRS (`produces`: its output carries param `to`)
- sockets: `in` in light:GHPRST; `out` out light:GHPRST
- ports: in(in→·); out(·→out)
- through: out ← in
- param `to` choice = "H" of R|H|G|S|P

### fusion.particle
domain: Particle (CliffordBem BEM): backscattered light of a focused beam.
- needs: cliffordbem; S in or out: also fieldtracing; native: RS
- sockets: `in` in light:HRS; `out` out light:HRS
- ports: front(in→out)
- through: out ← in
- param `point` point = [0.0, 0.0, 0.0] [mm] – particle centre
- param `radius_um` number = 0.25 [µm] range 0…
- param `index` number = 1.5 – refractive index of the particle
- param `n_medium` number = 1.0 – index of the medium around it
- param `refinement` integer = 4 range 1…7 – icosphere level
- param `n_theta` integer = 40 range 4… – far-field quadrature
- param `n_phi` integer = 80 range 8…

### fusion.sum
splitter: Junction of unrolled branches: the sum of the light arriving at one state.
- needs: nothing (structural: the executor adds what arrives; made by `unroll`); native: – (no evaluator)
- sockets: + dynamic (per node, from the import or `unroll`)

### group
domain: Group (R28): a sub-model with variants; in the socket form it keeps the group's pose and exported parameters for the blocks inside.
- needs: nothing (structural: `unroll`/`flatten` resolve it into its blocks); native: – (no evaluator)
- sockets: + dynamic (per node, from the import or `unroll`)

### meta.optimise
meta: Optimisation of one parameter for the extremum of a measure.
- needs: nothing; native: R
- sockets: `results` out data:result; `drive` out data:drive; `observe` in data:result
- through: results ← –
- param `parameter` text = "" – block.param or block.param[i], e.g. mirror.point[0]
- param `collect` text = "" – output to collect: block.socket, e.g. camera.image
- param `measure` choice = "peak" of peak|centre|power – measure of each result
- param `component` choice = "dz" of dx|dy|dz|rx|ry|rz – offset component when the drive goes into a block's pose input (mm, rad)
- param `bounds` numbers = [0.0, 1.0] – lower and upper bound
- param `goal` choice = "maximum" of maximum|minimum
- param `scan` integer = 9 range 3… – coarse scan before the search
- param `tolerance` number = 0.0 – 0: 1e−6 of the bounds

### meta.replicates
meta: Replicates: repeated with replicate values, images averaged.
- needs: nothing; native: R
- sockets: `results` out data:result; `drive` out data:drive; `observe` in data:result
- through: results ← –
- param `parameter` text = "" – block.param or block.param[i], e.g. mirror.point[0]
- param `collect` text = "" – output to collect: block.socket, e.g. camera.image
- param `measure` choice = "peak" of peak|centre|power – measure of each result
- param `component` choice = "dz" of dx|dy|dz|rx|ry|rz – offset component when the drive goes into a block's pose input (mm, rad)
- param `values` numbers = [] – replicate values (empty: 0 … count − 1)
- param `count` integer = 8 range 1…

### meta.sweep
meta: Sweep: a parameter over values, the output collected per value.
- needs: nothing; native: R
- sockets: `results` out data:result; `drive` out data:drive; `observe` in data:result
- through: results ← –
- param `parameter` text = "" – block.param or block.param[i], e.g. mirror.point[0]
- param `collect` text = "" – output to collect: block.socket, e.g. camera.image
- param `measure` choice = "peak" of peak|centre|power – measure of each result
- param `component` choice = "dz" of dx|dy|dz|rx|ry|rz – offset component when the drive goes into a block's pose input (mm, rad)
- param `values` numbers = [] – values (empty: start, stop, count)
- param `start` number = 0.0
- param `stop` number = 1.0
- param `count` integer = 11 range 1…

### optics.aplanatic_lens
domain: Ideal aplanatic lens (sine condition, Richards-Wolf illumination).
- needs: nothing (numpy/scipy); H or S on its links: fieldtracing; native: R
- sockets: `in` in light:HRS; `out` out light:HRS
- ports: in(in→·); out(·→out)
- through: out ← in
- param `point` point = [0.0, 0.0, 0.0] [mm]
- param `axis` point = [0.0, 0.0, 1.0]
- param `focal_length` number = 10.0 [mm]
- param `aperture_radius` number = 5.0 [mm]

### optics.beam_splitter
splitter: Plane beam splitter (Michelson / Linnik arrangement).
- needs: nothing (numpy/scipy); H or S on its links: fieldtracing; native: R
- sockets: `in` in light:HRS; `r` out light:HRS; `t` out light:HRS; `r_in` in light:HRS; `t_in` in light:HRS; `merged` out light:HRS
- ports: in(in→·); r(r_in→r); t(t_in→t); merged(·→merged)
- through: r ← in; t ← in; merged ← r_in, t_in
- param `point` point = [0.0, 0.0, 0.0] [mm]
- param `normal` point = [0.0, 0.0, -1.0]
- param `up` point = [0.0, 1.0, 0.0] – direction of the local v axis
- param `drive` number = 0.0 [mm] – shift along the normal
- param `r` number = 0.7071067811865476
- param `t` number = 0.7071067811865476

### optics.collimated_source
source: Collimated beam (plane wave per field angle and wavelength).
- needs: nothing (numpy/scipy); H or S on its links: fieldtracing; native: R
- sockets: `out` out light:HRS
- ports: out(·→out)
- param `center` point = [0.0, 0.0, 0.0] [mm]
- param `direction` point = [0.0, 0.0, 1.0]
- param `radius` number = 1.0 [mm] range 0…
- param `rings` integer = 20 range 1…
- param `wavelengths_um` numbers = [0.6328] [µm]
- param `field_angles_deg` numbers = [[0.0, 0.0]] [deg]
- param `irradiance` number = 1.0 [W/mm²]
- param `polarisation` point = [1.0, 0.0, 0.0] – field direction
- param `up` point = [0.0, 1.0, 0.0] – direction of the local v axis
- param `gaussian_waist` number = 0.0 [mm] range 0… – 0: uniform; > 0: Gaussian amplitude exp(−ρ²/w²)

### optics.detector
detector: Plane detector: coherent image (amplitudes per coherence group).
- needs: nothing (numpy/scipy); H or S on its links: fieldtracing; native: HRS
- sockets: `in` in light:HRS; `image` out data:result
- ports: in(in→·)
- through: image ← in
- param `point` point = [0.0, 0.0, 0.0] [mm]
- param `normal` point = [0.0, 0.0, -1.0]
- param `up` point = [0.0, 1.0, 0.0] – direction of the local v axis
- param `drive` number = 0.0 [mm] – shift along the normal
- param `half_size` numbers = [1.0, 1.0] [mm]
- param `pixels` numbers = [128, 128]

### optics.interface
domain: Plane interface of two media (Fresnel, s and p).
- needs: nothing (numpy/scipy); H or S on its links: fieldtracing; native: R
- sockets: `front_in` in light:HRS@front; `front_out` out light:HRS@front; `back_in` in light:HRS@back; `back_out` out light:HRS@back
- ports: front(front_in→front_out); back(back_in→back_out)
- through: front_out ← front_in, back_in; back_out ← front_in, back_in
- param `point` point = [0.0, 0.0, 0.0] [mm]
- param `normal` point = [0.0, 0.0, -1.0]
- param `up` point = [0.0, 1.0, 0.0] – direction of the local v axis
- param `drive` number = 0.0 [mm] – shift along the normal
- param `n_front` number = 1.0 – index on the side the normal points to
- param `n_back` number = 1.5 – index behind

### optics.plane_mirror
domain: Plane mirror.
- needs: nothing (numpy/scipy); H or S on its links: fieldtracing; native: R
- sockets: `in` in light:HRS; `out` out light:HRS; `drive` in data:drive
- ports: front(in→out)
- through: out ← in
- param `point` point = [0.0, 0.0, 0.0] [mm]
- param `normal` point = [0.0, 0.0, -1.0]
- param `up` point = [0.0, 1.0, 0.0] – direction of the local v axis
- param `drive` number = 0.0 [mm] – shift along the normal
- param `reflectance_amplitude` number = -1.0

### optics.splitter
splitter: Lossless splitter with four ports (Stokes relations).
- needs: nothing (numpy/scipy); H or S on its links: fieldtracing; native: R
- sockets: `front_a_in` in light:HRS@front; `front_a_out` out light:HRS@front; `front_b_in` in light:HRS@front; `front_b_out` out light:HRS@front; `back_a_in` in light:HRS@back; `back_a_out` out light:HRS@back; `back_b_in` in light:HRS@back; `back_b_out` out light:HRS@back
- ports: front_a(front_a_in→front_a_out); front_b(front_b_in→front_b_out); back_a(back_a_in→back_a_out); back_b(back_b_in→back_b_out)
- through: front_a_out ← front_b_in, back_a_in; front_b_out ← front_a_in, back_b_in; back_a_out ← back_b_in, front_a_in; back_b_out ← back_a_in, front_b_in
- param `point` point = [0.0, 0.0, 0.0] [mm]
- param `normal` point = [0.0, 0.0, -1.0]
- param `up` point = [0.0, 1.0, 0.0] – direction of the local v axis
- param `drive` number = 0.0 [mm] – shift along the normal
- param `r` number = 0.7071067811865476
- param `t` number = 0.7071067811865476

### optics.thin_lens
domain: Paraxial thin lens (matrix optics, any incidence).
- needs: nothing (numpy/scipy); H or S on its links: fieldtracing; native: R
- sockets: `in` in light:HRS; `out` out light:HRS
- ports: in(in→·); out(·→out)
- through: out ← in
- param `point` point = [0.0, 0.0, 0.0] [mm]
- param `axis` point = [0.0, 0.0, 1.0]
- param `focal_length` number = 10.0 [mm]
- param `aperture_radius` number = 1000000000.0 [mm]

### rtt.detector
detector: Detector / camera: amplitudes within, intensities across coherence groups.
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); H or S: also fieldtracing; native: HRS
- sockets: `in` in light:GHPRS; `image` out data:result; + dynamic (per node, from the import or `unroll`)
- ports: in(in→·)
- param `surface` text = ""
- param `pixels` integer = 512 range 16…8192 – pixels (long side)
- param `half_size` numbers = [] [mm] – image window [hu, hv]; empty: the detector aperture

### rtt.mirror
domain: Mirror (one surface, reflecting).
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); native: R
- sockets: `in` in light:GHRS; `out` out light:GHRS; `drive` in data:drive; + dynamic (per node, from the import or `unroll`)
- ports: front(in→out)
- through: out ← in
- param `surface` text = ""
- param `offset` point = [0.0, 0.0, 0.0] [mm] – shift of the elements of the block (Raytatouille backend)
- param `drive` number = 0.0 [mm] – shift along the surface normal (set by a phase shifter's drive)

### rtt.part.detector
detector: Raytatouille detector (a new element).
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); H or S: also fieldtracing; native: HRS
- sockets: `in_in` in light:HRS@front; `image` out data:result
- ports: in(in_in→·)
- through: image ← in_in
- param `rtt_new` data = {…} (see `registry.get(...).params`) – the Raytatouille element (pose, surfaces, glass)

### rtt.part.lens
domain: Raytatouille lens (a new element).
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); native: R
- sockets: `front_in` in light:HRS@front; `front_out` out light:HRS@front; `back_in` in light:HRS@back; `back_out` out light:HRS@back
- ports: front(front_in→front_out); back(back_in→back_out)
- param `rtt_new` data = {…} (see `registry.get(...).params`) – the Raytatouille element (pose, surfaces, glass)

### rtt.part.mirror
domain: Raytatouille mirror (a new element).
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); native: R
- sockets: `front_in` in light:HRS@front; `front_out` out light:HRS@front
- ports: front(front_in→front_out)
- param `rtt_new` data = {…} (see `registry.get(...).params`) – the Raytatouille element (pose, surfaces, glass)

### rtt.part.splitter
splitter: Raytatouille beam splitter (a new element).
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); native: R
- sockets: `front_a_in` in light:HRS@front; `front_a_out` out light:HRS@front; `front_b_in` in light:HRS@front; `front_b_out` out light:HRS@front; `back_a_in` in light:HRS@back; `back_a_out` out light:HRS@back; `back_b_in` in light:HRS@back; `back_b_out` out light:HRS@back
- ports: front_a(front_a_in→front_a_out); front_b(front_b_in→front_b_out); back_a(back_a_in→back_a_out); back_b(back_b_in→back_b_out)
- param `rtt_new` data = {…} (see `registry.get(...).params`) – the Raytatouille element (pose, surfaces, glass)

### rtt.part.stop
domain: Raytatouille stop (a new element).
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); native: R
- sockets: `front_in` in light:HRS@front; `front_out` out light:HRS@front; `back_in` in light:HRS@back; `back_out` out light:HRS@back
- ports: front(front_in→front_out); back(back_in→back_out)
- param `rtt_new` data = {…} (see `registry.get(...).params`) – the Raytatouille element (pose, surfaces, glass)

### rtt.section
domain: Sequential surfaces of a Raytatouille system (e.g. an objective).
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); native: R
- sockets: `front_in` in light:GHPRST@front; `back_out` out light:GHPRST@back; `back_in` in light:GHPRST@back; `front_out` out light:GHPRST@front; `drive` in data:drive; + dynamic (per node, from the import or `unroll`)
- ports: front(front_in→front_out); back(back_in→back_out)
- through: back_out ← front_in; front_out ← back_in
- param `location` text = "" – node of the Raytatouille model tree
- param `offset` point = [0.0, 0.0, 0.0] [mm] – shift of the elements of the block (Raytatouille backend)
- param `drive` number = 0.0 [mm] – shift along the surface normal (set by a phase shifter's drive)

### rtt.source
source: Field points, wavelengths, polarisation of a system (coherence groups).
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); native: R
- sockets: `out` out light:GHPRS
- ports: out(·→out)
- param `system` file = "" – Raytatouille system
- param `fields` text = "" – field indices; empty: all
- param `wavelengths` text = "" – wavelength indices; empty: all
- param `rtt_system` data = {} – settings of the system
- param `rtt_tree` data = {} – tree nodes no block owns
- param `rtt_paths` data = [] – order and keys of the paths
- param `rtt_catalogs` data = [] – glass catalogue files
- param `rtt_catalog_data` data = [] – glass catalogues (stored)

### rtt.splitter
splitter: Beam splitter: splits into channels, combines coherently per group.
- needs: raytatouille; evaluators in `rtt_fusion.graph.backend_rtt.EVALUATORS`, not in `workflow.registry()` (section 2.9); native: R
- sockets: `in` in light:GHRS; `r` out light:GHRS; `t` out light:GHRS; `r_in` in light:GHRS; `t_in` in light:GHRS; `merged` out light:GHRS; `drive` in data:drive; + dynamic (per node, from the import or `unroll`)
- ports: in(in→·); r(r_in→r); t(t_in→t); merged(·→merged)
- through: r ← in; t ← in; merged ← r_in, t_in
- param `surface` text = ""
- param `offset` point = [0.0, 0.0, 0.0] [mm] – shift of the elements of the block (Raytatouille backend)
- param `drive` number = 0.0 [mm] – shift along the surface normal (set by a phase shifter's drive)
<!-- END GENERATED: block catalogue -->

## 4 Public API by module

Signatures abbreviated; `reg` = `NodeTypeRegistry`; (g) = re-exported by `rtt_fusion.graph`.
`rtt_fusion`: `__version__`, `reference(section=None) -> str` (this file or one section).

### 4.1 rtt_fusion.graph.model (g)

- `Graph(nodes={}, links=[], meta={})`: `add(node)`, `remove(id)`, `connect(source, target,
  transport="auto") -> Link`, `disconnect(link)`, `set_transport(link, t) -> Link`,
  `insert_on_link(link, node, in_socket, out_socket)`, `replace_node(id, new, socket_map)`,
  `set_muted(id, muted=True)`, `node(id)`, `links_from(id, socket=None)`, `links_to(...)`,
  `unique_id(stem)`. `GraphError` for bad edits (checks are separate: `validate`).
- `Node`, `Link(source, target, transport)` (frozen), `NodeType`, `ParamSpec(name, type, default,
  unit, minimum, maximum, options, group, label, description)` (`type`: number, integer, boolean,
  choice, text, file, files, numbers, point, entity, data), `SocketSpec`, `PortSpec(name,
  in_socket, out_socket, side, label)`, `PortFrame`, `Category`, `Direction`, `SocketKind`;
  helpers `light_in(name, accepts="RHGSTP", side, label)`, `light_out`, `data_in(name, type)`,
  `data_out`; constants `LIGHT_REPRESENTATIONS`, `DATA_TYPES`, `TRANSPORTS`.

### 4.2 rtt_fusion.graph.registry (g)

- `builtin_registry() -> reg` (rtt.*, fusion.convert/sum/condensed, group, analysis/data stubs;
  no optics). `reg.add(t)`, `reg.get(id)`, `id in reg`, `reg.types(category=None)` (sorted),
  `reg.add_adapter(describe, type_id=None)`; `node_type_from_describe(describe)`.

### 4.3 rtt_fusion.graph.workflow (front-end facade)

- `registry() -> (reg, evaluators)`: built-in types plus `OPTIONAL_MODULES` (optics, condense,
  grating, stack, particle, meta, display, conversion); evaluators of those (not `rtt.*`).
- `available() -> frozenset` executable representations ({R} or {R, H, S}).
- `checked(graph, reg, evaluators=None) -> (diagnostics: list[dict], plan: dict)`.
- `native_of(graph, evaluators=None)`, `produces_of(graph, evaluators=None)`: planner callbacks.
- `catalogue(reg) -> list[dict]` (block library: type, category, description, sockets, params,
  through, ports).
- `model_document(model, reg, describe=None, keep_fields=None, providers=None, recorded=None)
  -> {"model", "types", "diagnostics", "plan", "providers"}`.
- `connection_plan(graph, model, reg, plan)`, `with_used(plan, used, graph=None, model=None,
  reg=None)`, `model_id("x@2") -> "x"`, `embed_files(model|graph)` (file contents into blocks),
  `kept_fields(document)`, `block_states(raw, model_ids)`.

### 4.4 rtt_fusion.graph.execute

- `Executor(graph, reg, evaluators, shared=None, available=frozenset("R"), cache=None,
  environment=None, operators=None)`: `check()`, `evaluate(node, socket)`, `results() ->
  {(node, socket): value}`, `status(wanted=None) -> {node: up to date}`, `needed(wanted)`,
  `native(id)`, `produces(id)`; attributes `plan`, `context`, `cache`, `operators`, `calls`.
- `Context`: `graph`, `shared` (settings, section 2.9), `log: list[str]`, `plan`, `used`,
  `evaluators`, `registry`, `evaluate(node, socket)` (for meta blocks), `progress(part, text)`,
  `on_node(id)`, `operators`; `report(part, text)`, `span(a, b)`.
- `ExecutionError(code, message, location)` (`.code`, `.location`); `ResultCache(limit_bytes=1 GiB)`
  (LRU, `warnings`); `EXECUTABLE`={R}, `WITH_FIELDS`={R,H}, `WITH_SPECTRA`={R,H,S};
  `Evaluator = Callable[[Context, Node, str, dict], Any]`; `produced_by(node, evaluators)`.

### 4.5 rtt_fusion.graph.planner (g: plan, Plan, Requirements)

- `plan(graph, reg, Requirements(...)) -> Plan`; `Plan.representations`, `.reasons`,
  `.conversions` (`Conversion(node, in_socket, out_socket, operators)`, `.chain` e.g. "RHS"),
  `.paths`, `.diagnostics`, `.cost`, `.ok`, `.representation(link)`, `.to_dict()`;
  `link_key(link)`; `OPERATORS`, `CHAINS`, `REPRESENTATION_COST`, `PREFERENCE`, `NATIVE_ONLY`.

### 4.6 rtt_fusion.graph.topology (g)

- `Model`: `add`, `connect(a, b, transport="auto", transport_back=None)`, `branch(name,
  *steps)` (steps as `Step` or tuples), `connection_at(end)`; fields `nodes`, `connections`,
  `links`, `branches`, `groups`, `observe_branch`, `meta`.
- `Connection`, `Step`, `Branch`, `Group(variants, active, condense=None).add_to(model, node)`,
  `Variant(model, ports, params)`; `ModelError(diagnostics)`.
- `unroll(model, reg) -> Graph`, `flatten(model, reg)`, `validate_model(model, reg)`,
  `connection_diagnostics(model, reg)`, `node_ports(node, reg)`, `port_frames(model, reg)`,
  `set_variant(model, group_id, name)`, `section_to_group(model, node_id)`,
  `model_from_graph(graph, reg)`.

### 4.7 rtt_fusion.graph.validate (g: validate, errors, find_cycle, replace_block, Diagnostic)

- `validate(graph, reg) -> list[Diagnostic(code, severity, message, location)]`,
  `errors(diags)`, `find_cycle(graph, reg)`, `replace_block(graph, reg, node_id, new)`,
  `node_sockets(node, reg)`, `node_through(node, reg)`.

### 4.8 rtt_fusion.graph.io (g: dumps, loads, save, load, model_to_dict, model_from_dict)

- v1 `to_dict(graph)`, `from_dict(d)`, `dumps`, `loads`, `save(graph, path)`, `load(path)`;
  v2 `model_to_dict(model)`, `model_from_dict(d, reg)` (reads v1 too); `FORMAT`, `VERSION`=1,
  `MODEL_VERSION`=2, `TYPE_ALIASES`.

### 4.9 Light: graph.light, graph.field, graph.spectrum, graph.imaging

- `light`: `RayBundle(key, wavelength_um, position, direction, opl, amplitude, source_xy, alive,
  irradiance=1, path=(), jones=None, prt=None, source_cos=1, source_frame=None, index=1,
  source_index=1, branch=None)`: `.size`, `.k0`, `.field_direction()`, `.field_vectors()`,
  `.evolve(**changes)`; `CoherenceKey`; `Plane(origin, normal, u, v)`, `Plane.make(origin,
  normal, up)`, `plane.intersect(bundle) -> (moved, uv)`; `split_field`, `split_along`,
  `physical_split`, `rotation_between`, `concat`.
- `field` (fieldtracing): `FieldBundle(key, plane, field, path, branch)`; `rays_to_field(bundle,
  plane=None)`, `field_to_rays(fb)`, `field_power(fb)`, `propagate_to(fb, plane, window=None,
  tolerance=None) -> (fb, notes)`, `field_image(fields, plane, hu, hv, nx, ny, log=None,
  progress=None, tolerance=None) -> Image`, `port_plane(bundle)`, `field_tolerance(shared)`,
  `FieldError`.
- `spectrum` (fieldtracing): `SpectrumBundle(key, plane, waves, path, branch, support)`
  (`.continuous`, `.refractive_index`, `.wavelength_um`, `.k0`); `field_to_spectrum(fb,
  method="auto", tolerance=None, notes=None, max_points=1024)`, `spectrum_to_field(sb, notes)`,
  `rays_to_spectrum(bundle, plane=None, tolerance=None, notes=None)`, `spectrum_to_rays(sb,
  notes)`, `propagate(sb, plane, notes)` (exact, in the spectrum), `spectrum_image(spectra,
  fields, plane, hu, hv, nx, ny, ...)`, `spectrum_power(sb)` (W), `overlapping_orders(spectra)`,
  `SupportError`; constants `SUM_BUDGET`, `SUM_BLOCK`, `CUT_WARNING`, `MIN_FRESNEL`.
- `imaging`: `Image(values, extent, ray_power, groups, note)` (`.power`, `.pixel_area`),
  `coherent_image(bundles, plane, hu, hv, nx, ny)` (R, scalar).

### 4.10 Block modules (`EVALUATORS`; `register(reg) -> evaluators` adds the types)

- `optics` (evaluators `source`, `mirror`, `beam_splitter`, `thin_lens`, `aplanatic_lens`,
  `interface`, `splitter`, `detector`; `hexapolar(rings)` pupil points).
- `grating`: `grating(ctx, node, socket, inputs)`, `load_table(path)`, `table_from_document(doc)`,
  `split_prt`; `stack`: `stack(...)`, `stack_config(config, from_back=False)`; `particle`:
  `particle(...)`, `ray_powers(bundle, alive)`, `rays_to_spectrum`, `spectrum_to_rays`,
  `OUT_WAVELENGTHS`; `condense`: `build_table`, `apply`, `condensed_node`; `conversion`:
  `convert`, `produces`; `meta`: `MetaResult(parameter, values, results, measures, measure,
  mean, best)`, `measure_of(result, measure)`, `MEASURES`; `display`: `View(kind, title, image,
  curves, x_label, y_label, unit)`, `Curve(name, x, y)`.
- `backend_rtt`: `EVALUATORS` (rtt.source/section/splitter/mirror/detector, rtt.part.*),
  `detector`, `foreign_blocks`.
- `importer` (g: `model_from_raytatouille`, `graph_from_raytatouille`, `new_system`),
  `system_from_model(model)`, `system_from_graph(graph)`, `NEW_SYSTEM`; `parts.PARTS`,
  `parts.prepare_parts(nodes)`.
- `pose`: `Isometry(rotation, translation)` (`from_pose`, `@`, `inverse`, `point`, `vector`),
  `rotation_xyz`, `angles_xyz`, `placement(node)`, `posed(node)`, `offset_of`, `COMPONENTS`.
- `embedded`: `FILE_PARAMS`, `embed_files(model)`, `embed_nodes(nodes)`, `load_catalogs(lib,
  params)`. `remote`: `BlockProvider(types, evaluators, prepare=None)` (`catalogue()`,
  `handle(request)`), `register_remote(reg, catalogues, call, root, local=None)`, `remote_types`.

### 4.11 Ports, adapters, couplings, data

- `rtt_fusion.ports`: `Port(z, direction ±1, index, wavelength_um)` (`.k0`, `.kz(kx, ky)`),
  `GriddedField(port, dx, dy, field (ny, nx, 2))`, `PlaneWaves(port, kx, ky, field (N, 2),
  orders=None, solid_angle=None)` (`.ez()`, `.power_density()`), `to_plane_waves(g, threshold)`,
  `to_gridded(s, nx, ny, dx, dy)`, `power(g)`, `spectrum_power(s, cone_na=None)`.
- `rtt_fusion.adapters.<fouriergrating|cliffordbem|specklebem|hpfem|bwpm>` (section 2.10);
  `Structure` per adapter (FourierGrating: `Structure(config)`; CliffordBem: `Structure(mesh,
  eps_particle, n_background, center_mm, z_front_mm, z_back_mm, tol)`); `front_port`,
  `back_port`; CliffordBem/SpeckleBem `scatter`, `scattered_plane_waves`; `adapters.biooptomd`
  (`read_profile_csv`, `profile_layers`, `stack_along_z`, `dipoles_from_tracker`).
- `couplings.jones_table`: `JonesTable(orders, side, wavelengths_um, sigma_x, sigma_y,
  period_um, n_in, n_out, jones)`, `build_table(handle, wavelengths_um, sigma_x, sigma_y,
  orders, side)`, `point_table`, `direct_jones`, `jones_many`, `interpolate(table, order, λ, σx,
  σy, method="cubic")`, `near_rayleigh_anomaly`, `prt_matrix`, `prt_many`.
  `couplings.grating_table`: `order_efficiencies`, `with_efficiencies(system, surface_id,
  table)`. `couplings.focus.richards_wolf(port, na, focus_mm, ...) -> PlaneWaves`.
- `rtt_fusion.convert`: `hemisphere(n_theta, n_phi, forward)`, `far_field_to_spectrum(...)`;
  `phasespace`: `field_to_packets`, `coupling_into_gaussian_mode`; `portdata.write(value,
  folder) -> handle`, `portdata.read(handle, folder)` (values as JSON + .npy);
  `results.to_json(type, data)`, `load_json(text)`, `encode`, `decode` (`raytatouille-result`).

## 5 Diagnostics and log lines

`Diagnostic.severity` "error" stops a run (`Executor.check` raises `ExecutionError(code, ...)`);
"warning" does not. `location`: node id, `"node.socket"` or a link key `"a.s -> b.t"`.

### 5.1 Graph checks (`validate`, `validate_model`, `connection_diagnostics`)

| Code | Meaning | Fix |
|---|---|---|
| graph.node_type_unknown | type id not in the registry | use `workflow.registry()`; register the module / plugin |
| graph.provider_missing | (document) type from a plugin that is not loaded | install/start that plugin |
| graph.socket_unknown | link names a missing socket | see section 3 for socket names |
| graph.link_direction | link not output → input | swap source/target |
| graph.kind_mismatch | light linked to data | link light to light, data to data |
| graph.data_type_mismatch | data types differ | match the data type |
| graph.no_common_representation | sockets share no representation | insert `fusion.convert` |
| graph.transport_unsupported | fixed transport not accepted by both sockets | choose one both accept or "auto" |
| graph.input_multiple | several links into one input | one link per input; in a model, branches + `fusion.sum` |
| graph.cycle | feedback loop on socket level | use branches (model) or a meta block |
| graph.param_unknown (warning) | param not in the type (typo!) | check `reg.get(type).params` |
| graph.detector_unconnected (warning) | detector without light | connect it |
| graph.no_branch | model without branches / detector no branch ends at | `model.branch(...)` |
| graph.branch_impossible | step not possible (port, connection, source→detector) | fix ports/connections of the step |
| graph.branch_duplicate | two branches with one name or the same steps | rename / drop |
| graph.port_unknown, graph.port_links | unknown port; several connections at one port | ports in section 3; split in blocks |
| graph.group_variant, _ports, _params, _branch | group inconsistent | same exported ports/params in all variants; inner branch names |
| graph.branch_backend | step through a Raytatouille block no system path makes | fix the branch or the system |
| graph.observer_ambiguous | observed node has several copies | give the link a `branch` |
| graph.port_facing, graph.port_behind, graph.port_medium | connected ports do not face / one behind the other / different media | fix `point`/`normal`/indices/pose |
| graph.port_aperture, _near_field, _tilt (warnings) | apertures miss; ports closer than 10 λ (evanescent fields lost); tilt > 1e-6 rad on an H connection (resampled) | check geometry |
| graph.file_unreadable, graph.file_external (warnings) | file param unreadable / CAD layer file stays outside | store content in `config`/`table_data` (`workflow.embed_files`) |

### 5.2 Planner (`plan.*`)

| Code | Meaning | Fix |
|---|---|---|
| plan.representation_unavailable | fixed transport not in `available` | pass `available=workflow.available()`; install fieldtracing |
| plan.no_representation | no admissible representation on a link (reasons listed) | other transport, conversion block |
| plan.no_operator | no conversion chain between two representations | change transports |
| plan.incoherent_merge (warning) | paths arrive without phase (P) | avoid P |
| plan.paths_truncated (warning) | > `max_paths` fusion paths | raise `Requirements.max_paths` |
| plan.no_output (warning) | no detector/analysis output | add a detector |

### 5.3 Execution (`ExecutionError.code`)

| Code | Meaning | Fix |
|---|---|---|
| exec.transport_unsupported | fixed transport the executor cannot run | `Executor(..., available=workflow.available())` |
| exec.no_evaluator | no evaluator for the type | add `backend_rtt.EVALUATORS` (rtt.*), the module's `register`, or the plugin |
| exec.socket_unknown, exec.cycle | bad output name; cycle (also a meta block observing itself) | names from section 3 |
| exec.param | invalid parameter (e.g. polarisation along the beam, grating without table/config) | fix params |
| exec.field | H conversion/propagation not admissible (caustic at the port, tilted planes, …) | move the port plane, other transport |
| exec.spectrum_support | continuous spectrum without sampling cannot become H or R (behind `fusion.condensed`) | keep S up to an S-native detector |
| exec.mixed_representations | rays and fields on one detector input | same transport on all links into it |
| exec.convert_unsupported | `fusion.convert` to G/P, or fieldtracing missing | to R/H/S; install fieldtracing |
| exec.solver | solver missing or failed (FourierGrating, CliffordBem) | install it; check config |
| exec.grating_range | wavelength outside the Jones table | wider table / configuration |
| exec.particle | too few rays reach the particle, or rays not towards it | aim the focus at `point` |
| exec.condense | group not laterally invariant/isotropic, probe lost, σ too large | lower `sigma_max`; don't condense lenses |
| exec.meta, exec.display | meta block without executor / no measure; display without input | link the observed output |
| exec.no_system, exec.system, exec.compile | Raytatouille: no system; generated system invalid; does not compile (missing glass catalogue) | `shared["system"]`, `shared["materials"]` with the AGF catalogue |
| exec.aiming, exec.branch, exec.split, exec.no_surface, exec.pose | Raytatouille paths cannot be aimed / followed / split; unknown surface or tree node | fix the system or the branch |
| exec.remote | a block of another process failed | see provider log |

### 5.4 Log lines (`Executor.context.log`, warnings, prefixed with the link key)

- `exec.spectrum_cut: <share> of the power left out of the spectrum (...)` – evanescent or
  backward waves dropped (> `CUT_WARNING` = 1e-6); move the port, smaller NA.
- `exec.spectrum_rays: ...` – S → R of a bundle narrower than 10 λ, or diffraction orders whose
  spectra overlap (rays lose their interference); keep S to the detector.
- Informational: `H -> S by the PFT, N waves`, `H -> S by an FFT on n×n points`, `S -> H by the
  inverse PFT/FFT`, `H propagation d mm by <method> (...)`, `image from the sum of N plane waves`,
  `the sum of N waves would alias ... – S -> H`, `field converted to rays (the table is for rays)`,
  `CliffordBem operator prepared`. `Image.note`: too few pixel centres covered (spot smaller
  than the pixels or outside the window) – adjust `half_size`/`pixels`.
- Front ends may add their own `check.*` codes (Raytatouille Fusion: `check.no_evaluator` =
  `exec.no_evaluator` reported at check time).

## 6 Gotchas and performance

### 6.1 Gotchas

1. `Executor(...)` runs only R by default: pass `available=workflow.available()` for H/S, else
   fixed H/S → `exec.transport_unsupported` and "auto" never picks H/S.
2. `workflow.registry()` has no evaluators for `rtt.*`: merge `backend_rtt.EVALUATORS`.
3. `builtin_registry()` lacks `optics.*`, gratings, …: use `workflow.registry()` or each module's
   `register(reg)` (e.g. `optics.register(reg)` returns its evaluators).
4. `from rtt_fusion.graph import validate` is the function; the module is
   `importlib.import_module("rtt_fusion.graph.validate")`. `graph.optics`, `graph.grating`,
   `graph.stack`, `graph.particle`, `graph.meta` are modules.
5. Model API takes **port** names, Graph API **socket** names (beam splitter ports `in`, `r`,
   `t`, `merged`; sockets `in`, `r`, `t`, `r_in`, `t_in`, `merged`; interface ports `front`,
   `back`). `io.load`/`loads` read only v1; a v2 file needs `model_from_dict` + `unroll`.
6. A model computes nothing without branches (`graph.no_branch`); `unroll` raises `ModelError`
   with all diagnostics.
7. Param typos are only warnings (`graph.param_unknown`) – the default is used silently.
8. R images are empty/wrong at a focus or caustic (`Image.note`); use H or S into the detector.
   Detector and field planes must be parallel for H (tilted: resampled geometrically, warning).
9. Units mix: geometry mm, wavelengths µm, FourierGrating configs µm, particle radius µm,
   `field_angles_deg` degrees, poses rad (Raytatouille JSON `rotation_deg` degrees).
10. `fg.stack` layers are listed from the substrate (back) upwards; the light from the back sees
    them reversed. `fg.grating` with only a table converts fields/spectra to rays (log line).
11. `rtt.detector`: the Python backend takes the image size from `shared["pixels"]`, the window
    from `half_size` or the aperture; its own `pixels` param is not read there.
12. The cache key covers node params and plain `shared` values; objects in `shared` (a compiled
    system, a material library) are not hashed: pass `environment=` or a fresh `cache`.
    Mutating `node.params` between runs is fine (keys change).
13. Complex `polarisation` only from Python; JSON documents hold real vectors.

### 6.2 Performance knobs

- Rays: `rings` of `optics.collimated_source` (3k² + 3k + 1 rays); `shared["rings"]` for
  `rtt.*`; detector `pixels`.
- Gratings/stacks with H or S: `field_transfer="interpolated"` (transfer from a verified sub-grid
  of the spectrum, tolerance 1e-8; much faster for wide beams; falls back to exact with a log line);
  `field_threshold` (drop weak plane waves, default 1e-10), `field_grid_max` (grid points per
  axis, 1024). Rays with a configuration: `sigma_step`, `sigma_margin` (table nodes), `orders`.
- Spectra: S images sum waves × pixels directly up to `spectrum.SUM_BUDGET` (5e8 terms) in blocks
  of `SUM_BLOCK`; beyond, or if the sum would alias, S → H and Fieldtracing propagation. Narrow
  windows near the focus are cheap; large windows far away go through H.
- `shared["field_tolerance"]` = ε of Fieldtracing's selection rules (default 0.02).
- Particle: `refinement` (icosphere level, 4× the triangles per level), `n_theta`, `n_phi`; BEM
  operators are kept in `Executor(operators=ResultCache())` across runs and poses.
- Condensed groups (`Group(..., condense={"sigma_max": 0.2, "count": 5})`): one table per inner
  branch and wavelength, reused across runs (operator store).
- Caching: reuse `cache=ResultCache(limit_bytes)` across executors; `Executor.status()` tells what
  is up to date; meta blocks recompute only what lies behind the changed block.
- C++ extension `rtt_fusion._diffraction` (binary wheels): per-ray loop of `fg.grating`.

## 7 Discovering more at run time

```python run
import rtt_fusion
from rtt_fusion.graph import workflow
from rtt_fusion.graph.execute import Executor

print(rtt_fusion.__version__, sorted(workflow.available()))
reg, ev = workflow.registry()
print([t.type_id for t in reg.types()])                      # every block type
t = reg.get("optics.detector")
print([s.to_dict() for s in t.sockets], [p.to_dict() for p in t.params])
print(t.to_dict()["through"], t.description)
print({k: sorted(getattr(e, "representations", "R")) for k, e in ev.items()})  # native reps
print(workflow.catalogue(reg)[0]["type"])                    # the editor's block library
print(rtt_fusion.reference("1.4")[:200])                     # one section of this file
help(Executor.evaluate)
```

- Command line: `python -m rtt_fusion.reference --list | --section NAME | --blocks [TYPE ...]
  [--live]`.
- Module docstrings are detailed (stages, formulas, sources): `help(rtt_fusion.graph.spectrum)`,
  `help(rtt_fusion.graph.planner)`, `help(rtt_fusion.graph.backend_rtt)`.
- Optional modules: `importlib.util.find_spec(name) is not None` for fieldtracing,
  fouriergrating, cliffordbem, raytatouille, specklebem, hpfem, bwpm; C++ extension:
  `importlib.util.find_spec("rtt_fusion._diffraction")`.
- Adapters: `rtt_fusion.adapters.<name>.describe()` (method, representation, ports, validity).
- Repository (docs in German, tests as usage examples): https://github.com/PhotonPower/RTT-Fusion.
