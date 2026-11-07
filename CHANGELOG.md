# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/), Versionen nach SemVer.

## [Unreleased]

### Hinzugefügt

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
