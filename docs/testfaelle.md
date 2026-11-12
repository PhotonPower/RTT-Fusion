# Testfälle

Status: ✅ umgesetzt und grün, ◐ teilweise, ○ offen.

| Nr. | Prüfung | Referenz | Abnahme | Stufe | Status | Test |
| --- | --- | --- | --- | --- | --- | --- |
| P1 | Port-Typen: ebene Welle G ↔ S, Leistung, E_z | analytisch (Dispersionsrelation, Parseval) | < 1e−13 | R0 | ✅ | `test_ports.py` |
| F1 | FourierGrating-Adapter: einzelne Grenzfläche Luft/Glas, senkrecht und schräg (bis 69°), beliebiger Azimut, komplexe Polarisation | vektorieller Fresnel-Operator in Convention A (unabhängig aufgebaut) | < 1e−10 | R1 | ✅ | `test_adapter_fouriergrating.py` |
| F2 | Schichtstapel 2,1 / 1,8 + 0,05i / 1,45 auf Glas | Transfer-Matrix der tangentialen Felder (s, p getrennt), Phasenbezug an oberer und unterer Grenzfläche | < 1e−9 | R1 | ✅ | `test_adapter_fouriergrating.py` |
| F3 | Lamellengitter (Λ = 0,7 µm, n = 2,1/Luft), schräg mit Azimut | Gittergleichung k_x = k_x,ein + n·2π/Λ; Energieerhaltung über alle Ordnungen | < 1e−6 /mm bzw. < 1e−9 | R1 | ✅ | `test_adapter_fouriergrating.py` |
| F4 | Reziprozität (asymmetrischer, absorbierender Stapel, Vakuum beidseitig) | Transmission von beiden Seiten gleich | < 1e−10 | R1 | ✅ | `test_adapter_fouriergrating.py` |
| F5 | Querprobe FourierGrating gegen BWPM, Schichtstapel, 45° Azimut, 30° | BWPM-Adapter (zweiter Löser) und Transfer-Matrix | < 1e−9 | R1 | ✅ | `test_adapter_fouriergrating.py` |

| B1 | BWPM-Adapter in den gemeinsamen Port-Typen | `bwpm.adapter` direkt (bitgleich); Energieerhaltung über `ports.power`; Querprobe mit FourierGrating über G → S | bitgleich bzw. < 1e−12 / 1e−9 | R2 | ✅ | `test_adapter_bwpm.py` |

| K1 | Pilot: FourierGrating-Effizienzen als Tabelle im sequenziellen Tracer von Raytatouille (Lamellengitter Λ = 0,7 µm in Vakuum, λ = 0,6 µm, Ordnungen −1, 0, +1) | Gittergleichung; FourierGrating-Richtungen; Effizienzen; Energieerhaltung | Richtungen < 1e−12 (gegen Gittergleichung) bzw. 1e−9 (gegen FourierGrating), Pfadtransmission = Effizienz < 1e−12, Σ = 1 < 1e−9 | R3 | ✅ | `test_pilot_k1.py` |

| C1 | CliffordBem-Adapter: ebene Welle über den Adapter = Lösung von CliffordBem selbst (Phase am Port, Polarisation) | `ScatteringProblem.solve_plane_wave` | < 1e−7 | R4 | ✅ | `test_adapter_cliffordbem.py` |
| C2 | Kugel R = 0,25 µm, n = 1,5, λ = 0,6 µm (720 Dreiecke) | Mie (eigene Reihe, Bohren & Huffman); σ_sca = σ_ext | Mie < 4 % (Genauigkeit flacher Elemente), Energie < 2 % | R4 | ✅ | `test_adapter_cliffordbem.py` |
| C3 | Verschiebung des Teilchens: Port-Amplituden ändern sich um e^(i(k_ein − k_aus)·Δ) | analytisch | < 1e−5 | R4 | ✅ | `test_adapter_cliffordbem.py` |
| K3 | Pilot: Teilchen im Fokus eines Objektivs (NA 0,8), Rückstreuung ins Objektiv | Fokusfeld gegen `BeamField.focused`; Leistungsanteil gegen die Rechnung allein in CliffordBem | < 1e−10 bzw. < 1e−6; Anteil 3,70 % | R4 | ✅ | `test_pilot_k3.py` |

| J1 | Tabellierte Beugungsinteraktion: PRT-Matrix je Ordnung, P k_ein = k_m, Transversalität, Energie über alle Ordnungen bei beliebiger Polarisation | Energieerhaltung (verlustfreies Gitter) | < 1e−9 | R6 | ✅ | `test_jones_table.py` |
| J2 | Senkrechter Einfall ohne Basis-Entartung (TE ≠ TM) | FourierGrating direkt | anisotrop, keine Kopplung < 1e−12 | R6 | ✅ | `test_jones_table.py` |
| J3 | Leistungsnormierung | Effizienzen aus K1 (s, p, 0° und 4,6°) | < 1e−12 | R6 | ✅ | `test_jones_table.py` |
| J4 | Phasenbezug am Gitterursprung | verschobenes Gitter: Faktor e^(−i m 2π x₀/Λ) | < 1e−10 | R6 | ✅ | `test_jones_table.py` |
| J5 | Interpolation (kubisch) und Rayleigh-Anomalien | FourierGrating direkt an Kontrollpunkten | Tabelle 4 × 5 × 4: < 5e−3; Anomaliezellen erkannt | R6 | ✅ | `test_jones_table.py` |

| S1 | Sollwertdateien für den C++-Port (`reference/*.result.json`): G → S, Halbkugel-Quadratur, Richards-Wolf-Strahl, Fernfeld → Spektrum, Beugungstabelle, FourierGrating-Stapel | unabhängig: Einzelwelle, Σw = 2π, ∫cos θ dΩ = π, Strahlleistung (2π)²/k·π NA², Dipolleistung k·4π/3, Tabelle gegen FourierGrating direkt | Dateien aktuell, unabhängige Werte < 1e−12 (Tabelle 5e−3), Rundreise | R7 | ✅ | `test_reference_files.py` |

| X1 | SpeckleBem-Adapter (exp(+jωt)): absorbierende Kugel ε = 2,25 + 0,3i, Fernfeld | Mie | < 3 % | R8 | ✅ | `test_adapter_specklebem.py` |
| X2 | Vorzeichen der Absorption (falsche Konjugation wäre Verstärkung) | Mie: σ_ext, σ_sca, σ_abs | σ_abs > 0,2 σ_ext; σ_ext < 3 % | R8 | ✅ | `test_adapter_specklebem.py` |
| X3 | Querprobe SpeckleBem gegen CliffordBem, schräger Einfall, elliptische Polarisation, Port-Spektren vorn und hinten | zweiter BEM-Löser mit anderen Konventionen | < 4 % | R8 | ✅ | `test_adapter_specklebem.py` |
| X4 | Überlagerung zweier Einfallswellen (Operator wiederverwendet) | Summe der Einzellösungen | < 1e−6 | R8 | ✅ | `test_adapter_specklebem.py` |

| H1 | hp-FEM-Adapter: ebene Grenzfläche Luft/Glas, senkrecht, schräg, konisch | vektorieller Fresnel (Convention A) | < 1e−12 (gemessen 3e−15) | R9 | ✅ | `test_adapter_hpfem.py` |
| H2 | Querprobe hp-FEM gegen FourierGrating: Lamellengitter aus K1, 0°, 4,6°, 4,6° konisch (φ = 40°), komplexe Vektoramplituden aller Ordnungen | zweites rigoroses Verfahren | starke Ordnungen < 0,6 %, R±1 < 3 % | R9 | ✅ | `test_adapter_hpfem.py` |
| H3 | Energie über alle Ordnungen (verlustfrei) | Energieerhaltung | < 5e−3 | R9 | ✅ | `test_adapter_hpfem.py` |

| K4 | Pilot: SPR-Biosensor (Prisma 1,515 / Gold 50 nm / Adlayer n(z) aus BioOptoMD-Format / Wasser, 633 nm); Datenport N | FourierGrating-Adapter gegen Transfer-Matrix; de Feijter (Verschiebung ∝ ∫(n − n_w) dz); Linearität; Vergröberung der Schichten | < 1e−9; 0,5 %; 6 %; 1 % | R10 | ✅ | `test_pilot_k4.py` |

| C++1 | C++-Prototyp von `rtt-fusion`: G → S, E_z, Leistung, Halbkugel-Quadratur, Richards-Wolf-Strahl, Fernfeld → Spektrum | Sollwertdateien `reference/*.result.json` (R7) samt unabhängigen Werten | < 1e−13 bzw. Dateitoleranz | R11 | ✅ | `cpp/tests/test_reference.cpp` |
| C++2 | BWPM-Adapter in C++ (Weg 2: außerhalb gegen die Header gebaut, BWPM-Kern per FetchContent) über die Registry | BWPMs Sollwertdatei `bwpm3d_tilted`; Energie über `power` | 3e−16 bzw. 6e−15 | R11 | ✅ | `cpp/tests/test_bwpm_adapter.cpp` |

| C++3 | Tabellierte Beugungsinteraktion in C++ (`DiffractionTable`): kubische Interpolation, Anomalie-Erkennung, PRT-Matrix, Prüfung der Tabelle | Sollwertdatei `jones_table_lamellar` (Interpolation aus der gespeicherten Tabelle, PRT aus FourierGratings direkten Matrizen); Anomalie-Zellen wie J5 | Interpolation bitgleich (0), PRT 1e−16, P k = k_m 1e−16 | R12 | ✅ | `cpp/tests/test_diffraction_table.cpp` |
| W1 | Import des Linnik aus Raytatouille-Pfaden (Struktur `tests/data/linnik_structure.rtt.json`) | Graph aus dem Konzept, Abschnitt 4: Objektiv als ein Block, Probe als eigener Block, Strahlteiler teilt und vereinigt, Detektoren; Hin- und Rückweg je Objektiv | Blockarten, Sockets, Seiten, interne Abhängigkeiten, 13 Bänder; keine Fehler, kein Kreis | R14 | ✅ | `test_graph.py` |
| W2 | Prüfung des Graphen | Definitionen der Diagnosecodes (Darstellungen, Datentypen, Richtung, Mehrfacheingang, Kreis auf Socket-Ebene, unbekannte Typen und Parameter) | je Fall genau die erwarteten Codes; Resonator = Kreis, Hin- und Rückweg ≠ Kreis | R14 | ✅ | `test_graph.py` |
| W3 | Bearbeitung und Format | Austausch Spiegel → Gitter mit erhaltenen Bändern, Einfügen auf ein Band, Stummschalten; `rtt-fusion-graph` | Bänder erhalten; JSON-Rundreise bitgleich, Version 2 und fremde Formate abgelehnt | R14 | ✅ | `test_graph.py` |
| W4 | Pilot W9: Michelson mit ebenen Elementen, Referenzspiegel um θ = 1 mrad gekippt (eigenes analytisches Backend `graph.optics`) | Zwei-Strahl-Interferenz ebener Wellen: Streifenabstand λ/sin 2θ auf einem Detektor senkrecht zu einem Strahl, Kontrast 1 | Abstand < 1e−6 rel., Kontrast 1 ± 2e−3 | R15 | ✅ | `test_execute.py` |
| W5 | Pilot W9: Wegunterschied δ = 0, λ/8, λ/4 | I = 2 I_Arm (1 + cos 2k₀δ), I_Arm = |r t|² = 1/4 | < 1e−3 absolut | R15 | ✅ | `test_execute.py` |
| W6 | Energie (δ = λ/8, Interferenzterm 0) | ∫I dA = 2 I_Arm π R² | < 5e−3 rel. | R15 | ✅ | `test_execute.py` |
| W7 | Kohärenzgruppen | zwei Wellenlängen bzw. zwei Feldpunkte: Bild = Summe der Einzelbilder (inkohärent) | < 1e−12 | R15 | ✅ | `test_execute.py` |
| W8 | Ausführung: Cache, Stummschalten, Fehler | Definitionen (nur Nachfolger neu, stummer Spiegel lässt durch → nur Arm A, I = I_Arm; Codes) | exakt | R15 | ✅ | `test_execute.py` |
| W10 | Planer, Entscheidungstabellen (R17) | Regeln aus Fusionsplan 7.2/7.3 und Tabelle 6.3 mit der vorläufigen Kosten-Rangfolge; je Fall ausgerechnet im Test: auto → R; nur-H-Detektor → Umwandlung R→H im Block davor; fester Transport bleibt; nicht verfügbar → `plan.representation_unavailable`; Kette G→H→R; P→H unzulässig (`plan.no_operator`); Gültigkeits-Hook; geteilter Eingang summiert die Anforderungen; inkohärente Vereinigung warnt | exakt (Wahl, Operatoren, Kosten, Codes) | R17 | ✅ | `test_planner.py` |
| W11 | Planer, Topologie am Linnik | Konzept Abschnitt 4: zwei Arme (Probe, Referenz) treffen sich an der Kamera | 2 Pfade zur Kamera, je einer über Probe und Referenz, alle ab der Quelle; Grenze der Pfadzahl meldet `plan.paths_truncated` | R17 | ✅ | `test_planner.py` |
| W12 | „auto“ wie R15 | Abnahme R17: bei „auto“ alle Bänder R, keine Umwandlung (Linnik mit R bzw. allen Darstellungen verfügbar; Michelson der Ausführung); W4–W9 unverändert | exakt | R17 | ✅ | `test_planner.py`, `test_execute.py` |
| W13 | Fokus einer aplanatischen Linse (NA 0,5, f = 2 mm) über ein H-Band | Richards-Wolf (Bessel-Integrale, Fieldtracing `focus.richards_wolf`), x-polarisiert | Spitze 2e−3 rel., Bild < 5e−3 der Spitze (gemessen 1,2e−3) | R18 | ✅ | `test_field.py` |
| W14 | Gauß-Strahl w₀ = 50 µm nach z = 0,5 / 1 / 2 z_R über ein H-Band | I = I₀ (w₀/w)² exp(−2r²/w²); Leistung π w₀² I₀/2 | < 2e−2 der Spitze (gemessen ≤ 1e−2); Leistung 1e−2 | R18 | ✅ | `test_field.py` |
| W15 | Michelson mit Gauß-Strahl, Band Strahlteiler → Kamera als H; alle Bänder H; Arm B um 2 mrad gekippt | Mitte H/R = 1/(1 + (L/z_R)²) mit L = 6 mm; δ = λ/8 halbiert die Intensität (W5); Hin- und Rückwandlung an den Spiegeln ändert nichts; gekippter Arm wie R | 3e−3; 2e−3; 1e−6 rel.; 3e−2 der Spitze | R18 | ✅ | `test_field.py` |
| W16 | Strahlen → Feld → Strahlen ohne Ausbreitung | Strahlenbild unverändert (Definition) | < 2e−3 der Spitze | R18 | ✅ | `test_field.py` |
| W17 | Paraboloid f = 100 mm, D = 2 mm (Raytatouille m6/paraboloid_f50), Band in den Fokus-Detektor als H | Airy-Muster [2 J₁(v)/v]², Spitze P π NA²/λ² mit P = π; Strahlenbild (R) im Fokus dunkel | normiert < 2e−3 (gemessen 5e−4); Spitze 1e−2 | R18b | ✅ | `test_field_rtt.py` |
| W18 | Michelson (Raytatouille m0, Testspiegel 2 mrad, EPD 1 mm, Gauß-Beleuchtung w = 0,5) mit H in die Kamera | Strahlenbild W9 bis auf Beugung über 50 mm ≪ z_R; gleiche Leistung | beleuchteter Bereich < 0,12 der Spitze (gemessen 0,088); Leistung 2e−2 | R18b | ✅ | `test_field_rtt.py` |
| W19 | Gitter-Block mit der Jones-Tabelle aus R6 (Lamellengitter Λ = 0,7 µm, λ = 0,6 µm), Einfall 0°, 3°, −2° | Gittergleichung sin θ_m = sin α + mλ/Λ; bei 0° gleiche Leistung in ±1 (symmetrisches Profil) | 1e−9; 1e−9 rel. | R19a | ✅ | `test_grating.py` |
| W20 | Interferenz der Ordnungen 2 µm hinter dem Gitter | Intensität periodisch mit Λ: Anteile bei 1/Λ (0 mit ±1) und 2/Λ (−1 mit +1), dazwischen nichts | Spitzen 1e−2 rel.; Rest < 5 % der Spitzen | R19a | ✅ | `test_grating.py` |
| W21 | Arbeitsgraph: Spiegel → Gitter mit `replace_node`; H-Band ins Gitter | Bänder bleiben, Graph rechnet; H wird am Port zu Strahlen, gleiche Ordnungen | Bänder gleich; Richtungen 1e−9 | R19a | ✅ | `test_grating.py` |
| W22 | Linnik-Arm mit idealen Linsen (`tests/data/linnik_mini.rtt.json`), Probe ersetzt durch eine Tabelle „Spiegel“ (Ordnung 0, J = −I) | Bild des ganzen Pfads mit dem Spiegel (Definition) – Strahlen und Felder | 1e−9 bzw. 1e−3 der Spitze | R19b | ✅ | `test_grating_rtt.py` |
| W23 | Derselbe Linnik, Probe dünnes sinusförmiges Phasengitter (φ₀ = 1, Λ = 10 µm, Ordnungen −1, 0, 1) | Raman–Nath: I/I_Spiegel = J₀² + 4J₁² sin²(K′x); Mittel J₀² + 2J₁², Extreme J₀² und J₀² + 4J₁², Periode MΛ/2, kein Anteil bei 1/(MΛ) | 5e−3 absolut; Periode 1e−2; Anteil < 2 % | R19b | ✅ | `test_grating_rtt.py` |
| W24 | Sweep über den Wegunterschied δ (Spiegel B, 9 Werte über λ/2), Messgröße Bildmitte | I(δ) = 2 I_Arm (1 + cos 2k₀δ) (W5); Quelle nur einmal gerechnet (Cache) | 2e−3 absolut; exakt | R20 | ✅ | `test_meta.py` |
| W25 | Replikate über δ gleichverteilt über eine Periode; Optimierung von δ in [λ/4, 3λ/4] | Mittel 2 I_Arm; Maximum bei δ = λ/2 | 2e−3; 1e−3 λ | R20 | ✅ | `test_meta.py` |
| W26 | Gitter-Block mit FourierGrating (Lamellengitter der Sollwertdatei, λ = 0,6 µm, TE, Einfall 0° und 3°), Tabelle vom Block gebaut | FourierGrating direkt (`direct_jones`, ohne Tabelle und Interpolation) | 1e−4 rel. (gemessen ~1e−6) | R19d | ✅ (ohne FourierGrating übersprungen) | `test_grating_fg.py` |
| W27 | Dasselbe Gitter, reflektierte und transmittierte Ordnungen | verlustfrei: Σ R_m + Σ T_m = 1 | 1e−5 (gemessen 3e−7) | R19d | ✅ (ohne FourierGrating übersprungen) | `test_grating_fg.py` |
| W28 | Sweep der Verschiebung (`offset`) des Testspiegels im Michelson aus Raytatouilles Referenzen, gleiche Arme | Zwei-Strahl-Interferenz: Periode λ/2 in der Verschiebung, Kontrast 1 | Periode 1e−3 rel., Kontrast 1e−2, Rest 1e−3 | R20b | ✅ (ohne Raytatouille übersprungen) | `test_meta.py` |
| W9 | Michelson aus Raytatouilles Referenzen über Import, Graph und Raytatouille-Backend, Testspiegel um 1 mrad gekippt | λ/sin 2θ; eine Kohärenzgruppe für beide Arme | Abstand < 1e−2 rel. (Strahlraster), Kontrast > 0,99 | R15 | ✅ (ohne Raytatouille übersprungen) | `test_execute.py` |

| PS1 | Phasenraum-Übergang Feld → Pakete (Husimi): Leistung, zweite Momente (Wigner + Fenster), Stokes, Kopplung in eine Gauß-Mode | analytische Gauß-Integrale | Leistung 1e−12, Momente ≤ 0,5 %, Kopplung 2 % | S9b-1 (Scattuccino D8/D10) | ✅ | `test_phasespace.py` |
| PS2 | Fokussierter Strahl: Pakete geradlinig weitergeführt gegen Wellenausbreitung (Winkelspektrum) | Wellenmoment + Fensterverbreiterung | < 2 % (gemessen ≤ 0,8 %) | S9b-1 | ✅ | `test_phasespace.py` |

## F1–F5: FourierGrating-Adapter (R1)

FourierGrating 0.16.0 (Commit `959c93b`), RCWA, CPU-Kern. Alle Vergleiche über die
Adapter-Schnittstelle in den Konventionen von Raytatouille.

**Befund zur Basis bei senkrechtem Einfall:** FourierGrating nimmt für die Ordnungen bei κ = 0
s = (0, 0, 1), für die Einfallswelle aber s_ein = (sin ψ, 0, −cos ψ). Bei φ = 0, ψ = 0 haben
dadurch amp_s und amp_p der Ordnungen das umgekehrte Vorzeichen gegenüber Convention A
(gemessen: r = +0,2 statt r_s = −0,2 an Luft/Glas). Physikalisch ist das richtig – es ist nur
eine andere Basis –, aber wer amp_s direkt als r_s liest, liegt bei senkrechtem Einfall falsch.
Der Adapter rekonstruiert deshalb die kartesischen Felder mit FourierGratings eigener Basis und
gibt nur diese weiter; F1 prüft gerade den senkrechten Einfall.

**Querprobe (F5):** Derselbe Schichtstapel (2,1 über 1,45 auf Glas, Dicken Vielfache von
λ/40) mit dem BWPM-Adapter (Gitter 80², 60 Iterationen) und dem FourierGrating-Adapter: Reflexion
und Transmission stimmen untereinander und mit der Transfer-Matrix auf < 10⁻⁹ überein (Transmission
nach Verschiebung um eine Scheibe, weil der hintere Port der BWPM eine Scheibe hinter der letzten
Grenzfläche liegt). Damit ist die offene Querprobe aus BWPM K16 erledigt.

## K1: Pilot Gitter im Raytracer (R3)

Raytatouille 0.9.0 (Commit `708b6dd`), FourierGrating 0.16.0. Effizienzen je Ordnung (Leistungsanteil,
Transmission und Reflexion) des Lamellengitters:

| Einfall | Pol. | T(−1) | T(0) | T(+1) | R(−1) | R(0) | R(+1) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0° | s | 0,3162 | 0,0768 | 0,3162 | 0,0323 | 0,2263 | 0,0323 |
| 0° | p | 0,3677 | 0,1307 | 0,3677 | 0,0518 | 0,0302 | 0,0518 |
| 4,6° | s | 0,2942 | 0,0762 | 0,3074 | 0,0143 | 0,2592 | 0,0486 |
| 4,6° | p | 0,4107 | 0,1349 | 0,3183 | 0,0481 | 0,0207 | 0,0672 |

Ergebnisse: Raytatouille lenkt Ordnung m in dieselbe Richtung wie FourierGrating die Ordnung n = m
(gleiches Vorzeichen bei Gittervektor = lokale x-Achse); die Pfadtransmission jeder Ordnung ist
genau die eingetragene Effizienz (keine doppelte Zählung von Grenzflächenverlusten, weil das Gitter
zwischen gleichen Medien liegt); die Energie über alle Ordnungen ist erhalten.

**Befunde für den Fusion-ADR (R5):**

1. `diffraction_efficiency` ist in Raytatouille eine Zahl je Ordnung und Fläche, unabhängig von
   Wellenlänge, Einfallswinkel und Polarisation (ADR 0025, Punkt 5). Schon bei 4,6° Einfall
   ändern sich die Effizienzen um bis zu 12 % (p, −1: 0,368 → 0,411), und s und p unterscheiden
   sich um bis zu 70 % (0. Ordnung). Für eine echte Kopplung braucht Raytatouille eine
   **Tabelle über (λ, Einfallsrichtung)** mit **Jones-Matrix je Ordnung** statt einer Zahl;
   der Pilot behilft sich mit einem System je Wellenlänge, Einfall und Polarisation.
2. Die Effizienz wirkt nach der eigenen Interaktion der Fläche (Fresnel der Ordnung 0). Liegt ein
   rigoros gerechnetes Gitter auf einer Grenzfläche zwischen verschiedenen Medien, muss die
   Kopplung durch diesen Faktor teilen – oder die Fläche bekommt eine ideale Interaktion.
3. Reflektierte Ordnungen eines Transmissionsgitters (hier bis 26 %) sind in einem sequenziellen
   Pfad nur über einen zweiten Pfad mit `reflect` erreichbar; die Tabelle muss beide Seiten
   liefern (der Adapter tut es).

## C1–C3, K3: CliffordBem und Teilchen im Fokus (R4)

CliffordBem 0.63 (Commit `41e64e6`), Kugel aus 720 flachen Dreiecken, R = 0,25 µm, n = 1,5 in
Vakuum, λ = 0,6 µm.

- **C1:** Eine ebene Welle über den Adapter (Phase 1 am vorderen Port, 2 µm vor dem Mittelpunkt)
  ergibt dieselbe Spur wie `solve_plane_wave` mal e^(ik·2 µm), auf 10⁻⁷.
- **C2:** σ_ext 3,0 % unter Mie, Fernfeld bei 0,3 und 1,0 rad 2,8 % – die bekannte Genauigkeit
  flacher Elemente bei dieser Auflösung (CliffordBem selbst: 1,1 % bei 2 000 Dreiecken).
  σ_sca aus den Port-Spektren gegen σ_ext 1,4 %.
- **C3:** Verschiebung um (0,3; −0,1; 0,4) µm: Amplituden auf 10⁻⁵ wie vorhergesagt.
- **K3:** Das Richards-Wolf-Spektrum (`couplings.focus`) gleicht `BeamField.focused` auf 10⁻¹³
  nach Normierung (reeller Faktor: gleiche Phasenkonvention). Vom Strahl kehren **3,70 %** der
  Leistung in die NA 0,8 zurück, über Ports und Spektren wie direkt in CliffordBem (< 10⁻⁶);
  insgesamt streut die Kugel das 2,16-fache der Strahlleistung (möglich bis 4×, weil die
  Extinktion die Interferenz mit dem Strahl enthält).

**Befunde für den Fusion-ADR (R5):**

1. **Fernfeld → Port-Spektrum:** Die Weyl-Darstellung E_s = (ik/2π)∫E∞ e^(ik k̂·r) dΩ ist nur
   mit evaneszenten Anteilen exakt. Lässt man sie weg, fehlt in einer Port-Ebene in endlichem
   Abstand ein Beitrag, der nur wie 1/r abfällt (Randterm streifender Richtungen): gemessen 6 %
   auf der Achse in 2–6 µm Abstand. Für alles, was ein Objektiv mit NA < n erreicht, ist das
   Spektrum dagegen exakt. Ports zu Streukörpern sollten deshalb **Spektren** (S) sein, nicht
   Gitterfelder in einer nahen Ebene.
2. **Kontinuierliche Spektren brauchen Gewichte:** Ein auf Quadraturrichtungen diskretisiertes
   Spektrum (Strahl, Fernfeld) hat eine Leistung (2π)²/k Σ|A|²dΩ; die Port-Typen tragen dafür
   `solid_angle`. Diskrete Wellen (Gitterordnungen) und Spektren lassen sich nicht ohne diese
   Unterscheidung addieren.
3. **Numerische Transversalität:** Das Fernfeld flacher Elemente ist nur auf etwa 1 % transversal;
   die Rekonstruktion von E_z aus (E_x, E_y) verstärkt das nahe streifender Richtungen
   (σ_sca vorher 6 % zu groß). Der Adapter projiziert transversal.

## J1–J5: Tabellierte Beugungsinteraktion (R6)

Prototyp `rtt_fusion.couplings.jones_table` für den ADR-Entwurf
`docs/adr-entwurf-beugungstabelle.md`; Gitter aus K1.

**Konvergenz der Interpolation** (Messung außerhalb der Tests, 30 Zufallspunkte, λ 0,585–0,615 µm,
σ_x ±0,08, Tabellen λ × σ_x mit n × n Stützstellen, σ_y 4 Stützstellen; größter Fehler eines
Eintrags von J_m gegen FourierGrating direkt):

| n | Schritt σ_x | linear | kubisch |
| --- | --- | --- | --- |
| 5 | 0,05 | 9,0·10⁻³ | 2,5·10⁻³ |
| 9 | 0,025 | 1,6·10⁻³ | 2,5·10⁻⁴ |
| 17 | 0,0125 | 5,9·10⁻⁴ | 1,4·10⁻⁵ |

**Nahe der Rayleigh-Anomalie** (Ordnung ±1 streifend bei σ_x ≈ ±0,143, λ = 0,6 µm; Tabelle
σ_x 0,10–0,20): Zellen, die die Anomalie enthalten, haben Fehler 0,11 / 0,087 / 0,031 bei
Schritt 0,025 / 0,0125 / 0,0063 (≈ √h), Zellen daneben 5,9·10⁻⁴ / 1,4·10⁻⁴. Die Prüfung
`near_rayleigh_anomaly` (Vorzeichenwechsel von |σ + m′λ/Λ| − n über der Zelle, alle Ordnungen,
beide Medien) erkennt genau die schlechten Zellen.

## X1–X4: SpeckleBem (R8)

SpeckleBem 0.1.0 (Commit `4f04efe`), Kugel R = 0,25 µm aus `make_sphere` mit Kanten ≤ λ/10
(1 280 Dreiecke), ε = 2,25 + 0,3i (in SpeckleBem 2,25 − 0,3j), λ = 0,6 µm, dicht, GMRES 10⁻⁸.

| Größe | Adapter | Mie | Abweichung |
| --- | --- | --- | --- |
| σ_ext (µm²) | 0,5353 | 0,5402 | −0,9 % |
| σ_sca aus den Port-Spektren | 0,3716 | 0,3751 | −0,9 % |
| σ_abs = σ_ext − σ_sca | 0,1637 | 0,1650 | −0,8 % |

Querprobe gegen CliffordBem (720-Dreieck-Kugel aus `make_icosphere(8)`, 17,5° Einfall,
Polarisation (0,8; 0,2 + 0,4i)): Port-Spektren vorn 1,7 %, hinten 0,7 % (größte Abweichung
relativ zum Maximum), Leistungen 2,7 % bzw. 1,2 %. Die Unterschiede liegen im Bereich der
Diskretisierung beider Löser; ein Konventionsfehler (Konjugation, Einheiten, Phasenbezug) gäbe
Abweichungen der Größenordnung 1 (X2: Verstärkung statt Absorption).

**Befund für SpeckleBem (Vorschlag, dort nur lesend):** Der Bau scheitert mit pybind11 3.x
(`bind_excitation.cpp:151`: Lambda mit zwei verschiedenen `pybind11::typing::Tuple`-Rückgabetypen,
3er-Vektor und N×3-Matrix). Mit pybind11 < 3 baut es. Eine Beschränkung über `PIP_CONSTRAINT`
reicht in der CI nicht: Ein früher installiertes Modul bringt pybind11 3.x in die Umgebung, und
CMake findet es auch bei isoliertem Bau. Die CI installiert deshalb pybind11 < 3 und baut
SpeckleBem ohne Isolation (`.github/workflows/ci.yml`). Issue-Entwurf für SpeckleBem:
`docs/vorschlaege/specklebem-pybind11.md`.

## H1–H3: hp-FEM-Gitter und Querprobe mit FourierGrating (R9)

hp-FEM 0.4.0 (Commit `dce07d5`), `grating.jones` auf einer Einheitszelle aus `hpfem.meshing`
(gmsh, `element_sizes(p=2)`, Ordnung 5; 822 Zellen), Lamellengitter wie K1. Größte Abweichung
der komplexen Vektoramplituden gegen FourierGrating (Nx = 15, selbst auf 10⁻⁴ konvergiert),
relativ zum Betrag der Ordnung:

| Einfall | R(−1) | R(0) | R(+1) | T(−1) | T(0) | T(+1) |
| --- | --- | --- | --- | --- | --- | --- |
| 0° | 1,3 % | 0,21 % | 0,82 % | 0,10 % | 0,22 % | 0,21 % |
| 4,6° | 2,2 % | 0,23 % | 1,1 % | 0,45 % | 0,23 % | 0,45 % |
| 4,6°, φ = 40° | 2,0 % | 0,23 % | 0,96 % | 0,31 % | 0,15 % | 0,40 % |

Ordnung 6 statt 5 auf demselben Netz ändert nichts: Die Ecken der Lamelle begrenzen die
p-Konvergenz, verfeinert werden muss das Netz (h). Auf dem feinsten Netz (`p=1`, Ordnung 6,
2 426 Zellen) liegen die Effizienzen bei s und 0° auf 0,1–0,3 % an FourierGrating (R(±1) 0,03221
/ 0,03223 gegen 0,0323; R(0) 0,22611 gegen 0,2263; T(0) 0,07676 gegen 0,0768; T(±1) 0,3165 gegen
0,3162). Die Querprobe bestätigt beide Adapter mit allen Konventionsumrechnungen: Rahmen, s/p-
Vorzeichen (hp-FEM: ŝ = k̂ × ŷ, p̂ = ŝ × k̂), Leistungsnormierung, Bezug auf Zellrand und
Stapelgrenzen; ein Fehler darin gäbe Abweichungen der Größenordnung 1.

## K4: Biosensor mit Brechzahlprofil aus der Molekulardynamik (R10)

BioOptoMD braucht CUDA und lief hier nicht; der Adapter liest seine Ausgabeformate
(`analysis.export_refractive_index`, Tracker-Zeitreihen; Commit `4a4c492`), die Testprofile sind
synthetisch in diesem Format (2 Replicas, Δz = 0,05 nm, 400 Bins). Kretschmann-Anordnung,
p-Polarisation, Resonanz ohne Schicht bei 72,214°:

| Profil | ∫(n − n_w) dz | Verschiebung | 1-nm-Schichten statt 0,05 nm |
| --- | --- | --- | --- |
| exponentiell, 0,12·e^(−z/3 nm) | 0,3559 nm | 0,5403° | +1,4·10⁻⁴ |
| Kasten, 0,06 auf 6 nm | 0,3564 nm | 0,5407° | 0 |
| exponentiell, doppelte Amplitude | 0,7119 nm | 1,0581° | +4,7·10⁻⁴ |

Verschiebung je optischer Überschussdicke: 1,5182 bzw. 1,5172 °/nm – formunabhängig auf
0,07 % (de Feijter); bei doppelter Amplitude 1,958-fach (2 % nichtlinear). Die Vergröberung auf
1-nm-Schichten ändert die Verschiebung um < 5·10⁻⁴: Für Plasmonensensoren genügt ein grobes
Profil, die MD-Auflösung von 0,05 nm wird nicht gebraucht.

**Befund für den Fusion-ADR:** Ein Datenport (N, Q) unterscheidet sich von einem Lichtport: Er
hat keine Richtung und keine Leistung, wohl aber eine Wellenlänge (n(z) gilt nur bei der λ des
MD-Laufs), eine Streuung über Replicas und eine Lage relativ zur Sensoroberfläche. Die Wahl der
Beleuchtungsseite ist eine Spiegelung (z → −z), keine Drehung – für Orientierungen von Dipolen
wichtig.

## PS1, PS2: Phasenraum-Übergang (Scattuccino S9b-1)

`rtt_fusion.phasespace`: Husimi-Verteilung des Feldes (Gauß-Fenster σ, Kurzzeit-FFT nur in den
gezogenen Fenstern, speicherarm), Pakete mit Ort, Richtung, Gewicht, Stokes; Kopplung in eine
Gauß-Mode über die Wigner-Dichte der Mode (Entscheidungen D8, D10 in Scattuccino #22).

- Kollimierter Gauß-Strahl (w₀ = 30 µm, σ = 10 µm, 2·10⁵ Pakete): Varianzen in Ort und Wellenzahl
  = Wigner + Fenster + Zellenzittern auf 0,1–0,5 %; Leistung exakt; Kopplung in dieselbe Mode
  0,385 = Π[(v_x,m + v_x,Q)(v_k,m + v_k,Q)]^(−1/2) (die Glättung verringert die Kopplung kohärenten
  Lichts – ballistisches Licht wird deshalb über das Feld gekoppelt, S9b-3).
- Fokussierter Strahl (w₀ = 60 µm, f = 1 mm, NA ≈ 0,06): Varianz der geradlinig weitergeführten
  Pakete in 0,5/1,0/1,5 mm gegen die Wellenausbreitung plus Fensterverbreiterung σ² + z²/(4σ²k₀²):
  0,1 / 0,8 / 0,6 %.
- **Befund für D9:** Im Fokus (1 mm hinter der Übergabe) ist der Paketfleck rms 12 µm breit, der
  Wellenfokus 2 µm – die Fensterverbreiterung dominiert. Die Übergabe gehört an
  oder hinter den Fokus, oder der Fokus in den wellenoptischen Bereich; das kleinste Verschmieren
  bei σ_opt = sqrt(z/(2nk₀)).
