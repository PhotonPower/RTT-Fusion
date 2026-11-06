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
3er-Vektor und N×3-Matrix). Mit pybind11 < 3 baut es; die CI von RTT-Fusion setzt deshalb
`PIP_CONSTRAINT` mit `pybind11<3`. Abhilfe in SpeckleBem: expliziter Rückgabetyp der Lambda oder
zwei getrennte Überladungen.
