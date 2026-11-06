# Schnittstelle der Adapter

Nach dem Fusionsplan (Fieldtracing `docs/rtt-fusion/plan.md`, Abschnitt 6.4). Ein Adapter ist
ein Python-Modul unter `rtt_fusion.adapters` mit vier Funktionen; Daten an der Grenze sind Werte
in den Konventionen dieses Repositorys (`docs/konventionen.md`, Typen in `rtt_fusion.ports`).

| Funktion | Inhalt |
| --- | --- |
| `describe() -> dict` | `name`, `version`, `method`, `representation` ("G" oder "S"), `ports` (Name, Ebene, Richtungen), `linear`, `reciprocal`, `operator` (ob eine Streumatrix angeboten wird), `deterministic`, `conventions_of_solver` (was der Adapter umrechnet), `validity` (gemessene Grenzen) |
| `prepare(structure, wavelengths_um) -> Handle` | prüft die Struktur, bildet einen Inhalts-Hash (Struktur, Wellenlängen, Löser-Version, Adapter-Version); rechnet nichts |
| `apply(handle, incident) -> dict[str, Ports]` | einfallendes Feld am vorderen Port (Richtung +1) → `"front"` (reflektiert, −1) und `"back"` (durchgelassen, +1) in derselben Darstellung |
| `estimate(handle, ...) -> dict` | `seconds`, `memory_bytes` aus Messungen |

Darstellungen: Löser mit Gitterfeldern (BWPM) nehmen `GriddedField`, Löser mit Beugungsordnungen
(FourierGrating) nehmen `PlaneWaves`; die Umwandlung G ↔ S liegt in `rtt_fusion.ports`, nicht im
Adapter.

**Konformitätstests** (Fusionsplan 6.5) je Adapter, gegen unabhängige Werte: Platte gegen Airy
und Transfer-Matrix (s, p, schräg, Azimut), Energieerhaltung, Reziprozität, Determinismus,
Phasenbezug der Ports; verfahrensspezifisch (Gittergleichung, Mie, …); und, wo zwei Löser
dasselbe rechnen können, die Querprobe zwischen ihnen.
