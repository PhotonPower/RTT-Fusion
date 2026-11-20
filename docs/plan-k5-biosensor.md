# Plan K5: Biosensor – Emitter aus BioOptoMD an der Sensoroberfläche, Abbildung durch ein Objektiv

> Fusionsplan Abschnitt 8 (K5), D13 Schritt 4 (RTT-Fusion#106). Stand 2026-11-20.

## Befund

- Q (Dipol-Reihe je Rezeptorzustand im Sensorsystem) und N (n(z)) liefert `biooptomd.layer`
  (R33); BindOptiX rechnet aus F, M, Q und B die Raten (R34).
- Ein Fluorophor wenige nm über Glas oder Gold koppelt sein Nahfeld ins Substrat: überkritisches
  Licht (SAF), über einen Goldfilm SPCE, und eine veränderte Zerfallsrate. Der S-Transport wirft
  evaneszente Wellen weg (`spectrum.py`); ein Fernfeld, das erst danach durch `fg.stack` läuft,
  verfehlt beides.

## Entscheidung des Maintainers (2026-11-20)

Analytischer Schichtblock: Dipol in bzw. über einem planaren Stapel (Weyl/Sommerfeld mit
Transfermatrix aus N), Spektrum S im Substrat einschließlich überkritischem Licht und die Raten
A_tot, A_rad, A_det je NA – damit auch die Antwortkarte M für BindOptiX ohne hp-FEM. Rigorose
Löser (CliffordBem, hp-FEM) später als Gegenprobe (DP3, K7).

## Stufen

| Stufe | Inhalt | Abnahme |
| --- | --- | --- |
| **K5a** | `couplings.dipole_layers`: Stapel je Wellenzahl mit dem Dipol als Primärfeld, `rates`, `angular_spectrum` | freier Dipol, Spiegeldipol, Nahfeld-Rate = abgestrahlte Leistung, Konvergenz mit Gold, Fluss des Spektrums ✅ (W94) |
| **K5b** | Block `emitter.dipole_layer`: Q (Zustand, Stichprobe) + N bzw. Schichtparameter → S im Substrat auf einer Port-Ebene, je Emitter eine Kohärenzgruppe (inkohärent über Emitter); Raten als Ergebnis; Pose Sensor ↔ RTT mit Spiegelregel (DP8, `in_rtt_frame`) | S gegen `angular_spectrum`; Leistung = `rates`; Spiegelung ✅ (W95) |
| **K5c** | Abbildung durch ein ideales aplanatisches Mikroskop für S (`optics.ideal_imaging`) auf die Kamera; ein Raytatouille-Objektiv später | Leistung in der NA, Bild bei kleiner NA = Airy, Verschiebung −M·Δ, Ring des z-Dipols ✅ (W96) |
| **K5d** | Antwortkarte M (`EmitterResponseTable`) aus `rates` über einem Gitter von Höhen und Detektorbereichen (NA-Kegel) für BindOptiX | DP4 (Chance–Prock–Silbey, Lukosz–Kunz), DP5 (Anteil in der NA), DP6 (Tensoren) |

Die Methode und ihre Konventionen stehen im Modulkopf von `couplings/dipole_layers.py`.
