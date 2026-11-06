# Plan und Stufen

> Maßgebliche Fassung für dieses Repository. Gesamtplan: Fieldtracing `docs/rtt-fusion/plan.md`
> (Phasen P0–P5); dieses Repository ist der Prototyp aus P2.

## Stufen

| Stufe | Inhalt | Abnahme |
| --- | --- | --- |
| **R0 Gerüst** | Regeln, Konventionen, Schnittstelle, Port-Typen G und S mit Umwandlung und Leistung, CI mit den Modulen | Tests der Port-Typen gegen analytische Werte ✅ (P1) |
| **R1 FourierGrating-Adapter** | Adapter mit Umrechnung der Konventionen, Konformität: Fresnel/Airy, Transfer-Matrix, Gittergleichung, Energie, Querprobe gegen den BWPM-Adapter | alle Konformitätstests grün ✅ (F1–F5) |
| **R2 BWPM-Adapter** | `bwpm.adapter` in die gemeinsamen Port-Typen eingebunden | Konformität wie R1, Querprobe ✅ (B1) |
| **R3 Pilot K1** | Gitterspektrometer: FourierGrating-Effizienzen als Tabelle im Raytracer (Python-API von Raytatouille) | Gittergleichung, Energie über die Ordnungen ✅ (K1, Befunde für R5 in `docs/testfaelle.md`) |
| **R4 Pilot K3** | Partikel im Fokus mit CliffordBem (zweiter Löser) | Mie |
| **R5 Fusion-ADR** | Entwurf für Raytatouille aus den Erfahrungen R1–R4, vor dem M9-ADR | Freigabe |

## Entscheidungen

| Nr. | Frage | Entscheidung |
| --- | --- | --- |
| D1 | Scope von Raytatouille | Kopplung rigoroser Löser über Ports, Löser bleiben außerhalb des Kerns; ADR in Raytatouille später (Maintainer, 2026-11-06, grundsätzlich) |
| D2 | Ort des Prototyps | dieses Repository `PhotonPower/RTT-Fusion` (2026-11-06) |
| D10 | Reihenfolge der Löser | FourierGrating zuerst, CliffordBem als zweiter Löser (2026-11-06) |
| D11 | BWPM | eigenes Repository `PhotonPower/BWPM` (2026-11-05) |
