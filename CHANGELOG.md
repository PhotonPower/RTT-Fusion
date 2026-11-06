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
