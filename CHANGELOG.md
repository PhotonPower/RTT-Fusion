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
