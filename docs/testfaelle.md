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
