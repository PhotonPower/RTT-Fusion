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
