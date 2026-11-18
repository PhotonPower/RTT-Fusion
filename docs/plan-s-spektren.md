# Plan R30: Spektren ebener Wellen (S) auf Verbindungen

Stand: freigegeben (2026-11-18, Entscheidungen unten). Grundlage: Darstellungen und Umwandlungen in
[fusionsplan.md](fusionsplan.md) (6.2, 6.3), Stufe F5 im
[Entwurf Gitter-Feldoperator](adr-entwurf-gitter-feldoperator.md) („S und G als
Verbindungsdarstellungen, wenn ein zweiter Block Spektren tauscht“) und eine Bestandsaufnahme des
Codes vom 2026-11-18. Anlass: Der Block `fusion.convert` bietet „S“ an, die Ausführung kann es aber
nicht rechnen (`exec.convert_unsupported`).

## Ziel

Licht kann zwischen Blöcken als **Spektrum ebener Wellen** fließen, nicht nur als Strahlen (R) oder
hybrides Feld (H). Blöcke, die ohnehin im Spektrum rechnen – Gitter, Schichtstapel, verdichtete
Gruppen, Streukörper –, tauschen es direkt aus, statt an jeder Verbindung zurück nach H oder R zu
gehen.

**Nutzen:**

- **Genauigkeit:** Ein Gitter oder Stapel wirkt auf jede ebene Welle exakt (Jones bzw. S-Matrix je
  Richtung). Heute tastet jeder Block das Feld neu ab, rechnet im Spektrum und interpoliert zurück auf
  die Stützpunkte von H (`graph/spectral.py`). Zwei Gitter hintereinander kosten zwei solche
  Rundreisen und deren Fehler.
- **Streukörper:** CliffordBem und SpeckleBem nehmen ein Winkelspektrum als Einfall und liefern ein
  Fernfeld, also wieder ein Spektrum. `fusion.particle` wandelt dafür heute intern R → S → R
  (`particle.rays_to_spectrum`, `spectrum_to_rays`). Mit S auf Verbindungen werden K3, K5 und K6 aus
  dem Fusionsplan als Ketten von Blöcken möglich, ohne Strahlen dazwischen.
- **Freie Ausbreitung** zwischen parallelen Ebenen ist im Spektrum exakt: nur die Phase
  e^(i k_z Δz), ohne Fresnel-Näherung und ohne Kaustikproblem.

## Was es schon gibt

| Baustein | Ort | Stand |
|---|---|---|
| Datentyp S: `PlaneWaves` (k_x, k_y, E_x, E_y, Ordnungen, Gewichte dΩ), E_z aus Transversalität, Leistung | `ports.py` | fertig, getestet (P1) |
| G ↔ S per FFT (`to_plane_waves`, `to_gridded`) | `ports.py` | fertig (P1) |
| H ↔ H im Frequenzbereich (PFT, inverse PFT) | Fieldtracing `GridlessField.pft`, `inverse_pft` | fertig; RTT-Fusion nutzt es noch nicht |
| R → S (Debye/Richards–Wolf je Strahl) und S → R (Strahlen vom Fokus aus) | `graph/particle.py` | fertig, aber nur blockintern |
| Feld → Spektrum → Gitter/Stapel je Welle → Feld | `graph/spectral.py`, `couplings/grating_field.py` | fertig, aber nur blockintern (F1–F4b) |
| Tabellen ebener Wellen einer Gruppe | `graph/condense.py` | fertig, wendet sie auf Strahlen an |
| Fernfeld → S | `convert.far_field_to_spectrum` (CliffordBem-, SpeckleBem-Adapter) | fertig |
| S als zulässige Darstellung im Modell, Operatoren und Kosten im Planer (`pft`, `fft`, `orders_to_rays`) | `graph/model.py`, `graph/planner.py` | nur geplant: Kosten ja, Ausführung nein |
| Konvertierungsblock `fusion.convert` | `graph/conversion.py` | R ↔ H; S meldet `exec.convert_unsupported` |

## Was fehlt

1. **Ein Lichttyp S in der Ausführung.** Die Ausführung kennt nur `RayBundle` und `FieldBundle`
   (`graph/light.py`, `graph/field.py`); `EXECUTABLE`/`WITH_FIELDS` sind {R} bzw. {R, H}
   (`graph/execute.py`).
2. **Umwandlungen als Operatoren der Ausführung** (heute nur blockintern): R → S, H → S, S → H,
   S → R, jeweils mit der Gültigkeitsprüfung aus Tabelle 6.3.
3. **Blöcke, die S annehmen und abgeben.** Heute nimmt kein Block S: Detektor `RH`, Gitter `RH`,
   Partikel R.
4. **CAD:** Transport „S“ an Verbindungen wählen, Anzeige R→S/H→S nach dem Lauf, Beschriftung der
   Ports, Handbuch.

## Entscheidungen (zur Freigabe)

**S1 Datentyp.** `SpectrumBundle(key, plane, waves, path, branch)` in `graph/spectrum.py`, analog zu
`FieldBundle`: eine Kohärenzgruppe, Port-Ebene `plane` (Bezugspunkt der Phase = Ursprung der Ebene,
Rahmen u, v, n), `waves` vom Typ `ports.PlaneWaves`. Zwei Arten:

- **diskret:** einzelne Wellen bzw. Beugungsordnungen (`orders` gesetzt). Ebene Wellen ohne
  Ausdehnung; Leistung je Welle als Fluss durch eine Bezugsfläche.
- **kontinuierlich:** ein abgetastetes Winkelspektrum (`solid_angle` gesetzt). Es beschreibt einen
  begrenzten Strahl, die Ausdehnung steckt im Spektrum.

*Empfehlung:* beide Arten in einem Typ (wie `PlaneWaves` schon), die Art aus den gesetzten Feldern.

**S2 Wie ein begrenzter Strahl als S weitergeht.** Ein Gitter, das von einem begrenzten Strahl
getroffen wird, gibt je Ordnung wieder einen begrenzten Strahl ab. Als diskrete Ordnungen ginge die
Ausdehnung verloren.

*Empfehlung:* Gitter und Stapel geben **kontinuierliche** Spektren ab: das Einfallsspektrum je
Ordnung um n·K verschoben, Gewichte unverändert. Diskrete Wellen nur für unbegrenzt beleuchtete
periodische Strukturen und für Probewellen (`condense`, Tabellen).

**S3 Umwandlungen.** Als Operatoren in der Ausführung, mit Diagnose bei Verletzung der Gültigkeit:

| Umwandlung | Verfahren | Gültigkeit (Diagnose) |
|---|---|---|
| H → S | PFT von Fieldtracing, dann Stützstellen mit dΩ | Spektrum vollständig erfasst: Apertur, evaneszenter Anteil < Schwelle (`exec.spectrum_cut`, Warnung mit Anteil) |
| S → H | inverse PFT auf die Port-Ebene | Abtastung nach Fieldtracing ADR 0003 |
| R → S | über H (R → H → S); im Fokus Debye je Strahl wie `particle.rays_to_spectrum` | wie R → H; Debye nur bei Fresnel-Zahl ≫ 1 |
| S → R | eine Strahlfamilie je Spektralbereich bzw. Ordnung, Start auf der Port-Ebene über der Ausdehnung aus der inversen PFT | Ordnungen getrennt, Bündel groß gegen λ (`exec.spectrum_rays`, Warnung) |
| G ↔ S | `ports.to_plane_waves` / `to_gridded` | Tabelle 7.1 des Papers |

*Empfehlung:* R → S immer über H; kein eigener Weg außer im Fokus (dort übernimmt der Partikelblock
seinen Debye-Weg in den gemeinsamen Operator).

**S4 Freie Ausbreitung.** Eine Verbindung mit S zwischen zwei Port-Ebenen mit gleicher Normale
verschiebt nur die Phase je Welle (e^(i k·Δr)). Bei gekippten Ebenen wird jede Welle in den neuen
Rahmen gedreht; evaneszente Anteile fallen mit Diagnose weg. Raytatouille-Abschnitte (Linsen)
nehmen S nicht: Davor wandelt die Ausführung S → H → R.

**S5 Planer.** Mit fieldtracing ist die verfügbare Menge {R, H, S}. „auto“ wählt S nur zwischen zwei
Blöcken, die beide S nativ können (Gitter → Stapel, Gitter → Partikel); sonst bleibt es wie heute.
Kosten aus `REPRESENTATION_COST` bleiben vorläufig, bis `estimate` (Stufe 3) sie ersetzt.

**S6 Standard von `fusion.convert`.** Bis R30a fertig ist, steht `to` auf H; danach kann S wieder
der Standard werden.

## Stufen

| Stufe | Inhalt | Abnahme (unabhängige Sollwerte) |
|---|---|---|
| **R30a Lichttyp und Umwandlungen** | `SpectrumBundle`; Operatoren H ↔ S, R → S, S → R, G ↔ S in `Executor._convert`; `WITH_SPECTRA` = {R, H, S}; Planer lässt S auf Verbindungen zu; `fusion.convert` mit `to = S`; Diagnosen `exec.spectrum_cut`, `exec.spectrum_rays` | Gauß-Strahl: H → S gegen das analytische Winkelspektrum, S → H → S und R → S → R Rundreisen, Leistung erhalten (1e−3); Fokus gegen Richards–Wolf wie W29 (W76–W78) |
| **R30b Detektor und freie Ausbreitung** | Detektor nimmt S (Bild aus der Summe der Wellen auf dem Pixelgitter bzw. über inverse PFT); Ausbreitung zwischen Port-Ebenen im Spektrum (S4) | zwei ebene Wellen: Streifen analytisch; Gauß-Fokus über S = über H (1e−3); Ausbreitung = Fresnel-Integral im Gültigkeitsbereich (W79, W80) |
| **R30c Gitter, Stapel, Verdichten mit S** | `fg.grating`, `fg.stack`, `fusion.condensed` nehmen S und geben S ab (je Welle exakt, Ordnungen als verschobene Spektren, S2); bei H am Eingang wie heute | zwei Gitter in Folge: S-Kette = FourierGrating mit beiden Gittern im Stapel bzw. = H-Kette bei großen Strahlen; Σ Ordnungen = R + T; Stapel gegen Transfermatrix (W81–W83) |
| **R30d Streukörper mit S** | `fusion.particle` nimmt S am Eingang (ohne inneres R → S) und gibt das Fernfeld als S ab; SpeckleBem-Adapter als Block, sobald L4 ihn mit vielen Einfällen hat | W29–W31 unverändert mit S-Verbindung; K3 ohne Strahlen zwischen Fokus und Teilchen (W84, W85) |
| **R30e CAD** | Transport „S“ im Verbindungsmenü des Node-Editors, Anzeige der genutzten Darstellung (H→S, S→R), Hinweise der Diagnosen, Handbuch | Linnik mit Si-Gitter (`auflicht_si_gitter`): Bild mit S am Gitter = mit H (Toleranz nach GF5) |

R30a und R30b sind die Voraussetzung; R30c und R30d sind unabhängig voneinander. R30e folgt R30a
(der Transport kann erscheinen, sobald die Ausführung S kann) und wächst mit jeder Stufe.

## Was R30 nicht löst

- **Linsen und Spiegel aus Raytatouille** bleiben Strahlen; S wird davor nach H bzw. R gewandelt.
- **Darstellung G** auf Verbindungen (Gitterfeld) bleibt intern; sie kommt erst, wenn ein Block sie
  braucht (BWPM).
- **Radiometrisch (P)** und **Kugelwellen (T-Matrix)** sind eigene Stufen (M9 bzw. später in M).
- **Leistung:** Viele Spektralpunkte × FourierGrating je Richtung bleiben teuer; L4 (S-Matrix je
  Löser, `apply_many`) ist die Antwort, nicht R30.

## Entscheidungen des Maintainers (2026-11-18)

| Nr. | Frage | Entscheidung |
|---|---|---|
| 1 | S2: begrenzte Strahlen hinter Gittern | **kontinuierliche Spektren** (Empfehlung übernommen) |
| 2 | S3: Weg R → S | **über H**, Debye direkt nur im Fokus (Empfehlung übernommen) |
| 3 | S5: wählt „auto“ S? | **Ja:** der Planer wählt S selbst zwischen zwei Blöcken, die S nativ können |
| 4 | Reihenfolge nach R30b | **zuerst Gitter und Stapel (R30c)**, dann Streukörper (R30d) |

## Nachtrag R30a (Umsetzung, 2026-11-18)

**Gebaut.** `graph/spectrum.py` mit `SpectrumBundle(key, plane, waves, path, branch, support)` und
den Operatoren `field_to_spectrum` (H → S, G → S), `spectrum_to_field` (S → H), `rays_to_spectrum`
(R → H → S) und `spectrum_to_rays` (S → H → R); `Executor._convert` wandelt an den Ports nach dem
Plan der Verbindung und dem, was der Block nativ nimmt (S → H, wenn er H nimmt, sonst S → R);
`execute.WITH_SPECTRA` = {R, H, S}, `workflow.available()` gibt es mit fieldtracing zurück;
`fusion.convert` mit `to = S`, Standard H (S6). Abnahme W76–W78 in `tests/test_spectrum_links.py`
(Messwerte in [testfaelle.md](testfaelle.md)): PFT gegen das Winkelspektrum des Gauß-Strahls mit
genau dem erwarteten Fehler z_R/z der stationären Phase, FFT 2·10⁻⁵, Leistung ≤ 6·10⁻⁵,
Rundreisen über die PFT auf Rundung, Fokus NA 0,8 aus dem Spektrum gegen Richards–Wolf 9·10⁻⁴.

**Abweichungen und Präzisierungen.**

1. **`support` im Lichttyp.** Aus den ebenen Wellen allein (k, Amplitude) lässt sich S → H nicht
   rechnen: Die inverse PFT braucht die Spektralphase ψ̃ mit Gradient (−ρ) und Hesse-Matrix, die
   inverse FFT das Gitter. Der Typ hat deshalb ein sechstes Feld `support`: das Spektrum als
   Fieldtracing-Feld im Frequenzbereich (eine Stützstelle je Welle) bzw. das Gitter mit der Lage
   jeder Welle. S → H nimmt die Amplituden aus `waves`, nur die Geometrie aus `support`; ein Block,
   der die Wellen ändert (R30c), lässt die Geometrie stehen. Diskrete Wellen ohne `support` werden
   noch kein Feld (`FieldError`, kommt mit R30b).
2. **H → S nicht nur über die PFT.** Die PFT braucht eine reguläre Hesse-Matrix und eine große
   lokale Fresnel-Zahl F (Fehler O(1/F)): ein kollimierter Strahl oder ein Gauß-Strahl an der Taille
   hat keine. „auto“ nimmt die PFT, wo 1/F ≤ ε (ε = `field_tolerance` des Laufs, Standard 0,02 wie
   Fieldtracings Auswahlregeln), sonst eine FFT des Felds auf dem Gitter von `graph.spectral` (Träger
   k_c, Residuum kubisch interpoliert – linear verlöre 2·10⁻³ der Leistung); passt kein Gitter, die
   PFT bis F = 10 mit einem Vermerk ihres Fehlers.
3. **Gewichte dΩ.** Δ²k je Welle aus den Delaunay-Dreiecken der Wellenvektoren (dΩ = Δ²k/(k k_z)),
   nicht über die Jacobi-Abbildung der Flächen in der Port-Ebene: Ein fokussierter Strahl hat seine
   Stützstellen gleichmäßig in der Pupille, also in k∥, und dort gilt die Trapezregel. Die
   Leistung des Spektrums ist dann die der Strahlen (2·10⁻¹⁶ bei NA 0,8); `field.field_power` der
   Felder liegt bei konvergierenden Bündeln am Rand etwa 1 % darunter (Quadratur in der Port-Ebene).
4. **G → S, S → G.** Gitterfelder von Fieldtracing (nach einer Ausbreitung in einen Fokus) gehen über
   `ports.to_plane_waves`; der Rückweg ist die inverse FFT mit derselben Konvention wie
   `ports.to_gridded` (das je Welle über das Gitter summiert und für große Gitter zu langsam ist).
   Wellen unter 10⁻⁸ der größten Amplitude fallen weg.
5. **Sockets.** Damit eine Verbindung S tragen kann, müssen beide Sockets S zulassen: Die Blöcke
   `optics.*` nehmen jetzt R, H und S (die Ausführung wandelt an ihren Ports). Die `rtt.*`-Typen
   ließen S schon seit R14 zu; Gitter, Stapel, Partikel und Sonde bleiben bei R und H bis R30c/R30d.
6. **Planer.** `Requirements.native` (Block → native Darstellungen) und `planner.NATIVE_ONLY` = {S}:
   „auto“ lässt S nur zwischen zwei Blöcken zu, die es beide nativ nehmen; ein fester Transport S
   gilt immer. Ohne `native` (abstrakte Planung, Planer-Tests) gilt die Regel nicht. Die Ausführung
   nimmt die Attribute `representations` ihrer Auswerter, das Workflow-Dokument die der optionalen
   Module (`workflow.native_of`). Da noch kein Block S nativ nimmt, ändert „auto“ keinen Plan.
7. **Diagnosen** stehen wie die übrigen Hinweise der Ausführung als Zeilen in `Context.log`
   (`<Verbindung>: exec.spectrum_cut: …`), nicht als `Diagnostic` des Plans. `exec.spectrum_rays`
   prüft in R30a nur die Größe des Bündels (< 10 λ/n); „Ordnungen getrennt“ folgt mit den Ordnungen
   in R30c.
8. **Test angepasst:** `test_convert_block.py` prüfte `exec.convert_unsupported` mit `to = S`; jetzt
   mit G und P, dazu ein Test für `to = S`.

**Offen.** `SpectrumBundle` in den Port-Daten (`portdata`, Blöcke anderer Prozesse) fehlt noch; es
wird gebraucht, sobald ein Plugin-Block S nativ nimmt (R30c). Ein Raytatouille-Detektor
(`rtt.detector`) mit S auf der Eingangsverbindung rechnet im Modus ganzer Pfade mit Strahlen wie bei
R (er wertet nur H aus).
