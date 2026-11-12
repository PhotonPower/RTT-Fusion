# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/), Versionen nach SemVer.

## [Unreleased]

### Hinzugefügt

- R21: Darstellungs-Blöcke (`rtt_fusion.graph.display`: `display.image`, `display.plot`,
  `display.beam` als Sonde an jedem Lichtausgang) mit dem Ergebnis `View`; Fortschritt der
  Ausführung (`Context.progress`, `Context.span`, `Context.report`; Raytatouille-Detektor je
  Pfad/Feld/λ, Feldausbreitung, Meta-Blöcke je Wert) und laufender Block (`Context.on_node`).
  Tests W32–W34.

### Geändert

- FourierGrating eingestellt auf 0.18.0 (Commit `174ef93`, statt 0.16.0 `959c93b`): baut auch mit
  Apple clang (FourierGrating PR #2), so dass FourierGrating, Raytatouille und CliffordBem in
  einem Prozess dieselbe C++-Laufzeit teilen; Konventionsdateien unverändert.

### Hinzugefügt

- R19c: Partikel als Probe – Block `fusion.particle` (`rtt_fusion.graph.particle`): fokussierte
  Strahlen als Winkelspektrum (`rays_to_spectrum`), Streuung mit CliffordBem (Adapter R4),
  Rückstreuung als Strahlen (`spectrum_to_rays`); kollimiertes Licht als eine ebene Welle.
  Tests W29–W31. CliffordBem lokal mit Apple clang gebaut (eine Zeile in den Python-Bindings,
  Vorschlag für CliffordBem).

- R20b: `offset` an Raytatouille-Blöcken verschiebt deren Elemente (Kopie des Systems,
  übersetzt je Verschiebung zwischengespeichert); Sweeps und Optimierung wirken damit auf
  importierte Systeme. Test W28.

- R19d: Der Gitter-Block baut seine Jones-Tabelle mit FourierGrating selbst (`fouriergrating`:
  Konfigurationsdatei, `side`: reflektiert/transmittiert), um die einfallenden Strahlen herum
  und zwischengespeichert. Tests W26, W27 (laufen mit FourierGrating, 1D lokal in Sekunden).

- R20: Meta-Blöcke (`rtt_fusion.graph.meta`) – Sweep, Replikate (Mittel der Bilder) und
  Optimierung eines Parameters um den azyklischen Graphen; Parameter als „block.param[i]“,
  gesammelter Ausgang „block.socket“, Messgrößen Spitze/Mitte/Leistung; die Ausführung bietet
  verschachtelte Auswertung (`Context.evaluate`) und ungecachte Evaluatoren. Tests W24, W25.

- R19b: Ausgetauschte Blöcke im Raytatouille-Backend – der Pfad eines importierten Systems
  wird am fremden Block geteilt (Ereignisse den Blöcken zugeordnet, Teilpfade A und B in einer
  Kopie des Systems), der Block wirkt auf die Strahlen, Raytatouille rechnet weiter; Felder
  (skalar entlang einer reellen Polarisationsrichtung, Rest in der PRT), OPL und PRT
  verkettet, zwei Polarisationen; für Strahlen- und Feldbilder. `replace_block` (Austausch mit
  Zuordnung der Sockets nach Richtung und Art). Evaluatoren und Registry im Kontext.
  Testsystem `tests/data/linnik_mini.rtt.json`. Tests W22, W23.

- R19a: Gitter als Probe – Block `fusion.grating` (`rtt_fusion.graph.grating`) mit den Sockets
  eines Spiegels: Beugungsordnungen als Kanäle aus einer Jones-Tabelle (R6, FourierGrating;
  `load_table` liest `raytatouille-result`-Dateien), Gittergleichung, PRT der Ordnung, Phase der
  Phasenschicht m·2πx/Λ, Ordnungen in derselben Kohärenzgruppe; `light.rotation_between`.
  Tests W19–W21.

- R18b: Felder im Raytatouille-Backend – ein Band in einen Detektor eines importierten Systems
  kann H tragen: Felder aus den aufgezeichneten Strahlen nach der letzten Fläche vor dem
  Detektor (PRT und Gewicht von Raytatouille, Brechzahlen aus dem OPL, zwei Polarisationen der
  unpolarisierten Quelle), Ausbreitung mit Fieldtracings Planer. Die Quellachsen werden
  geometrisch mitgeführt (`RayBundle.source_frame`), damit Drehung und Spiegelung des
  Strahlenbilds durch die Geometrie den Maslov-Index nicht verfälschen; Brechzahl des Mediums
  (`index`, `source_index`); Gauß-Beleuchtung (`gaussian_waist`), Bildfenster `half_size`,
  Plan im Kontext der Ausführung. Tests W17, W18.

- R18: Field Tracing auf Bändern – Darstellung H (`rtt_fusion.graph.field`: `FieldBundle`,
  Strahlen → Feld nach Fieldtracing ADR 0004, Feld → Strahlen, Bild aus Feldern mit Ausbreitung
  durch Fieldtracings Planer), Polarisation der Strahlbündel (Jones der Quelle, geometrische PRT,
  Ausgleich der Spiegel-Orientierung im Maslov-Index), `optics.aplanatic_lens`, Gauß-Quelle
  (`gaussian_waist`), Ausführung mit wählbaren Darstellungen (`Executor(available=…)`,
  `WITH_FIELDS`), Umwandlung an den Ports. fieldtracing bleibt optional (Import erst bei H).
  Tests W13–W16. Der Fall in W8 mit festem Transport H (abgelehnt) nutzt jetzt G, weil die
  Blöcke des analytischen Backends H annehmen.

- R17: Planer, Stufen 1–2 (`rtt_fusion.graph.planner`) – Fusionspfade auf Socket-Ebene,
  Darstellung je Licht-Band und Umwandlungen in Blöcken nach Tabelle 6.3 (auch Ketten; nach
  Phasenverlust keine kohärente Darstellung mehr), feste Transporte, verfügbare Darstellungen,
  Gültigkeits-Hook, „auto“ ohne P solange kohärent möglich, Warnung `plan.incoherent_merge`;
  Plan als Datenobjekt mit JSON-Form für den Node-Editor; die Ausführung plant vor dem Rechnen
  (`Executor.plan`). Tests W10–W12.

- R15: Ausführung des Arbeitsgraphen ohne Planer – `graph.light` (Strahlenbündel R mit
  Kohärenz-Kennung), `graph.execute` (Auswertung je Ausgang, Cache nach Inhalts-Hash,
  Stummschalten, Fehlercodes `exec.*`), `graph.imaging` (kohärentes Kamerabild: Amplituden je
  Kohärenzgruppe, Intensitäten über Gruppen), `graph.optics` (analytisches Backend mit ebenen
  Elementen), `graph.backend_rtt` (importierte Raytatouille-Graphen, auch nicht zielbare Pfade);
  Pilot W9 (Tests W4–W9).

- R14: Graph-Modell des Node-Workflows `rtt_fusion.graph` – Blöcke, Sockets (Licht: R/H/G/S/T/P,
  Daten: N, Q und benannte Typen), Bänder mit Transport, Bearbeitung (Verbinden, Einfügen auf ein
  Band, Austauschen mit erhaltenen Bändern, Stummschalten), Prüfung mit Diagnosecodes und
  Kreissuche auf Socket-Ebene, Format `rtt-fusion-graph`, Block-Typen aus `describe()`, Import aus
  den Pfaden eines Raytatouille-Systems; Tests W1–W3.

- L11: Fusionsplan aus Fieldtracing (Zweig `docs-fusion-plan`, Commit `3e0d5f8`) als
  `docs/fusionsplan.md` übernommen; Verweise angepasst.

- R13: Konzept Node-Workflow (freigegeben, E3) `docs/konzept-node-workflow.md` – Arbeitsgraph aus
  Blöcken, Sockets und Bändern als bearbeitbare Form des Planers, Kohärenz-Kennungen,
  Meta-Blöcke für Schleifen, Ort von Logik (RTT-Fusion / `rtt-fusion`) und GUI (Raytatouille
  Fusion), Vorschläge A1–A6 zur Adapter-Schnittstelle, Inventar aller Module, Stufen R14–R20.
  Vorschläge L1–L12 zu den offenen Punkten: ein G-Typ mit optionaler Trägerphase, Umrechnung der
  FFT-Normierung, Port-Typen R/H/P/M, `apply_many`/`operator`, Einfall am hinteren Port, Namen der
  Daten-Sockets, Streustatistik, maschinenlesbare Gültigkeit, Beistelldateien für große Daten,
  Kopplung GUI ↔ Bibliothek mit Formatskizze `rtt-fusion-graph`.

- R0: Gerüst (Regeln, Plan mit D1, D2, D10, D11), Konventionen, Schnittstelle der Adapter,
  Port-Typen `GriddedField` und `PlaneWaves` mit Umwandlung G ↔ S und Poynting-Leistung, CI mit
  den Löser-Modulen auf festen Commits.
- R1: FourierGrating-Adapter `rtt_fusion.adapters.fouriergrating` (Drehung der Koordinaten,
  kartesische Rekonstruktion mit FourierGratings Basen, Einfall als Überlagerung von δ = 0°/90°,
  Ports an oberer und unterer Grenzfläche) mit Konformitätstests F1–F5, darunter die Querprobe
  gegen den BWPM-Adapter.
- R2: BWPM-Adapter `rtt_fusion.adapters.bwpm` (Gitterfeld-Ports), bitgleich zu `bwpm.adapter`;
  Energie über `ports.power`, Querprobe mit FourierGrating über die Umwandlung G → S.
- R3: Pilot K1 – `rtt_fusion.couplings.grating_table` (Effizienzen je Ordnung aus FourierGrating,
  als `diffraction_efficiency` in ein Raytatouille-System) mit Test gegen Gittergleichung,
  Pfadtransmission und Energie; Raytatouille als Modul (CI baut es mit Eigen, nlohmann-json, TBB).
  Befunde für den Fusion-ADR in `docs/testfaelle.md`.
- R4: CliffordBem-Adapter `rtt_fusion.adapters.cliffordbem` (Einfall als Winkelspektrum, Fernfeld →
  Port-Spektren über Weyl), `couplings.focus.richards_wolf`, `ports.spectrum_power` und
  `PlaneWaves.solid_angle`; Tests C1–C3 und Pilot K3 (Teilchen im Fokus, Rückstreuung ins
  Objektiv). Befunde für den Fusion-ADR in `docs/testfaelle.md`.
- R5: ADR-Entwurf „Raytatouille Fusion“ für Raytatouille (`docs/adr-entwurf-fusion.md`): Scope,
  Gebiete und Ports, Darstellungen R/H/G/S/T, tabellierte Jones-Interaktion als Erweiterung von
  ADR 0025, Adapter mit Konformitätssatz, Modell und Pfade, Planung, Kohärenz, Cache, GPU;
  begründet mit den Befunden B1–B7 aus R1–R4 und BWPM.
- Fusion-ADR freigegeben (E1): Punkt 9 „`rtt-fusion` in C++“ mit drei Wegen für Adapter;
  Vorschläge zu Nummer, Zeitpunkt und eigenem ADR für die tabellierte Interaktion; Stufen R6, R7.
- R6: Prototyp der tabellierten Beugungsinteraktion `rtt_fusion.couplings.jones_table`
  (leistungsnormierte Jones-Matrix je Ordnung in kartesischer Basis, PRT-Matrix, kubische
  Interpolation, Erkennung von Rayleigh-Anomalien) mit Tests J1–J5 und Konvergenzmessung;
  ADR-Entwurf `docs/adr-entwurf-beugungstabelle.md`.
- R7: Sollwertdateien für den C++-Port (`tools/export_reference.py`, `reference/`, sechs Fälle mit
  unabhängigen Werten); Umwandlung Fernfeld → Spektrum und Halbkugel-Quadratur als Modul
  `rtt_fusion.convert` (der CliffordBem-Adapter nutzt sie).
- R8: SpeckleBem-Adapter `rtt_fusion.adapters.specklebem` für geschlossene Körper (Konjugation
  von ε, Einfall und Fernfeld wegen exp(+jωt), SI-Einheiten) mit Tests X1–X4: Mie für eine
  absorbierende Kugel, Vorzeichen der Absorption, Querprobe gegen CliffordBem, Überlagerung.
  CI baut SpeckleBem mit `pybind11<3` (`constraints-modules.txt`).
- CI: jedes Löser-Modul in einem eigenen Schritt; SpeckleBem ohne Isolation gegen pybind11 < 3;
  Issue-Entwurf für SpeckleBem `docs/vorschlaege/specklebem-pybind11.md`.
- R9: hp-FEM-Adapter für Gitter `rtt_fusion.adapters.hpfem` (`grating.jones` → kartesische
  Felder in den Port-Konventionen) mit Tests H1–H3: Fresnel (3e−15), Querprobe gegen
  FourierGrating mit allen Ordnungen (starke < 0,6 %, schwache < 3 %), Energie.
- R10: BioOptoMD-Datenadapter `rtt_fusion.adapters.biooptomd` (n(z)-CSV, Schichten nach
  BioOptoMDs Regel, optische Überschussdicke, Dipole aus dem Tracker, Beleuchtungsseite als
  Spiegelung) und Pilot K4 (SPR-Biosensor): de Feijter formunabhängig auf 0,07 %, Vergröberung
  auf 1 nm < 5e−4; Nachtrag zum Fusion-ADR (Datenports, weitere Querproben).
- CI: Löser-Module als Wheels in einem Cache (actions/cache); ohne Änderung der Pins entfällt der
  Bau der Module.
- E2: ADR-Entwurf tabellierte Beugungsinteraktion freigegeben.
- R11: C++-Prototyp von `rtt-fusion` (`cpp/`): Port-Typen, Umwandlungen, Ergebnisleser,
  `SolverAdapter`/`AdapterRegistry`; BWPM-Adapter in C++ gegen BWPMs Kern (FetchContent); Tests
  gegen die Sollwertdateien (C++1, C++2); CI-Job `cpp`.
- CI: Cache-Schlüssel ohne den Workflow (Änderungen an anderen Jobs bauen die Module nicht mehr
  neu); Version `v1` im Schlüssel für Änderungen der Bauschritte.
- R12: tabellierte Beugungsinteraktion in C++ (`DiffractionTable`: kubische Interpolation,
  Rayleigh-Anomalien, PRT-Matrix, Prüfung) mit Test C++3 gegen die Sollwertdatei; Leser für
  boolesche Werte.
- Phasenraum-Übergang `rtt_fusion.phasespace` (Feld → Husimi-Pakete, Kopplung in eine Gauß-Mode)
  für Scattuccino S9b; Tests PS1, PS2.
- `phasespace.field_to_packets(periodic=True)`: wrapping windows for periodic fields (BWPM ohne
  Absorberrand), Fenstergröße geprüft.
