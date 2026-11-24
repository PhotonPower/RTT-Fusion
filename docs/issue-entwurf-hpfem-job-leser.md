# Issue-Entwurf für hp-FEM: öffentlicher Leser für Gitter-Jobs

> Entwurf nach AGENTS.md, Regel 4 (Löser-Repositories nur lesen). Anlass: Stufe R35.

## Befund

Der Block `hpfem.grating` (R35) baut die Gitterzelle aus einem hp-FEM-Job. Damit ein Job dort
genauso rechnet wie mit `python -m hpfem.run`, nutzt `adapters.hpfem.structure_from_job` die
Funktionen, mit denen `hpfem.run` den Job liest: `_unit_cell`, `_material`, `_stack` und
`_build_mesh`. Diese Funktionen sind privat. Ändert hp-FEM sie, bricht der Adapter, und das
fällt erst in der CI von RTT-Fusion auf, wenn der Commit in `requirements-modules.txt` erhöht
wird.

## Vorschlag

Eine öffentliche Funktion in `hpfem.run`, etwa

```python
def grating_cell(job, wavelength, base=None) -> GratingCell:
    """mesh, materials (by tag), stack at the wavelength, order, pml, options – what
    run_job hands to grating.solve for one point of a scattering job."""
```

Dazu gehört `structure_extent(mesh, materials, stack, omega)` (heute privat in
`grating._structure_extent`): Ober- und Unterkante der Streuer, auf die der Adapter seine Ports
legt.

Bis dahin prüft `tests/test_hpfem_grating.py` in der CI (mit `RTT_FUSION_REQUIRE_MODULES=1`), ob
die privaten Funktionen noch passen.
