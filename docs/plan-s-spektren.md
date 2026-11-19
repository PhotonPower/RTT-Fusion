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

## Nachtrag R30b (Umsetzung, 2026-11-18)

**Gebaut.** In `graph/spectrum.py`: `spectrum_image` (Bild aus Spektren und Feldern, Amplituden
je Kohärenzgruppe), `wave_sum` (Σ_j E_j e^(i k_j·(r − o)) an Punkten im Raum), `grid_sum` (dieselbe
Summe für Gitterspektren als zwei Matrixprodukte), `sum_phase_step` (Prüfung der Abtastung) und
`propagate` (S4). `optics.detector` und `rtt.detector` nehmen S nativ (`representations` = RHS);
in `graph/field.py` ist `field_image` in `field_sums` und `image_of` geteilt, damit Felder und
Spektren auf einem Eingang interferieren. Abnahme W79–W80 in `tests/test_spectrum_detector.py`
(Messwerte in [testfaelle.md](testfaelle.md)): Streifen zweier ebener Wellen auf Rundung
(2,5·10⁻¹², auch auf einem gekippten Detektor), Gauß-Fokus NA 0,5 über S gegen H 7,5·10⁻⁴ und
gegen Richards–Wolf mit Gauß-Apodisation 6,8·10⁻⁵ (H: 7,8·10⁻⁴), Gauß-Strahl im Spektrum um
0,5–2 z_R verschoben gegen das Fresnel-Integral 1,6·10⁻⁵, verschobenes PFT-Spektrum über die
inverse PFT = geometrische Ausbreitung von Fieldtracing auf Rundung, Paraboloid (W17) mit S gegen
Airy 5,6·10⁻⁴, knotenweise = ganze Pfade.

**Entscheidungen und Abweichungen.**

1. **Bild als Summe der Wellen im Raum.** Der Detektor wertet die Wellen an den Pixelmitten im
   Raum aus, ohne das Spektrum erst auf die Detektorebene zu bringen: Eine ebene Welle ist in
   jedem Rahmen dieselbe, so verliert auch ein gekippter Detektor nichts (S4 ließe dort Wellen
   fallen). Für diskrete Wellen ist das exakt und braucht keine Abtastung (`support`) – damit sind
   diskrete Spektren am Detektor nutzbar. Für kontinuierliche Spektren ist es die Quadratur des
   Winkelspektrums: exakte freie Ausbreitung ohne Fresnel- oder Debye-Näherung, solange die Wellen
   für das Fenster dicht genug liegen. Geprüft wird der Phasenschritt |r − ρ_j|·δk_j zwischen
   Nachbarwellen über die Ecken des Fensters, ρ_j der Auftreffpunkt der Welle (aus −∇ψ̃ der PFT
   bzw. der Ausdehnung des Felds auf dem FFT-Gitter, entlang k_j verschoben): Grenze 2π für
   Gitterspektren (Poisson: Kopien im Abstand des Fensters), π für die verstreuten Wellen der PFT.
   Im Fokus ist das immer erfüllt (dort ist die Summe das Debye-Integral ohne dessen Näherung und
   genauer als der H-Weg), weit vom Fokus nicht; dann wandelt der Detektor S → H auf der Ebene des
   Spektrums (die inverse Transformation, exakt) und rechnet weiter wie mit H (Fieldtracings
   Planer), mit Vermerk im Protokoll – das Bild ist dann das von H. Statt „inverse PFT auf die
   Detektorebene“: diese bräuchte im Fokus eine reguläre Hesse-Matrix, die es dort nicht gibt.
2. **Kosten.** Eine Exponentialfunktion je Welle und Pixel (PFT-Spektrum mit 2800 Wellen auf 256²
   Pixeln: 1,8·10⁸ Terme, wenige Sekunden). Gitterspektren (FFT, oft 65 000 Wellen) auf einem zur
   Spektrumsebene parallelen Detektor mit gleichen Achsen (bis auf Vorzeichen und Reihenfolge)
   zerfallen in zwei Matrixprodukte über Zeilen und Spalten des Gitters (`grid_sum`, gleich der
   direkten Summe auf 10⁻¹²). Über `SUM_BUDGET` = 5·10⁸ Terme geht ein Spektrum mit Abtastung über
   H. Beim Paraboloid (W17) ist S rund 150-mal schneller als H (Debye auf großem Gitter).
3. **Freie Ausbreitung (S4) als Helfer, nicht an der Verbindung.** `propagate(sb, plane)`
   verschiebt die Phase je Welle (e^(i k·Δr)), dreht bei gekippter Ebene jede Welle in den neuen
   Rahmen (Normale entlang des Lichts) und lässt evaneszente und nicht vorwärts laufende Wellen mit
   `exec.spectrum_cut` (Anteil der Leistung) fallen; auf parallelen Ebenen klingen evaneszente
   Wellen ab und fallen nur hinter der Ebene weg. dΩ gehört zur Richtung und bleibt. Die Abtastung
   eines PFT-Spektrums wandert mit (ψ̃ + k·Δr, Gradient −ρ um die Strahlversetzung −(k∥/k_n)Δn,
   Hesse-Matrix um −Δn ∂²k_n/∂k∥²), sodass S → H auf der neuen Ebene gilt; Gitterspektren und
   gekippte Ebenen verlieren sie (`support` None, Vermerk). Die Ausführung bewegt Spektren wie
   Felder nicht zwischen Port-Ebenen (`_convert` kennt die Port-Ebene des Verbrauchers nicht, und
   H wird ebenfalls erst im Block ausgebreitet); ein Block, der S nimmt, ruft `propagate` selbst
   auf (R30c: Gitter und Stapel auf ihrer Ebene).
4. **`rtt.detector`** (offener Punkt aus R30a): Mit S auf der Eingangsverbindung nimmt er im
   Modus ganzer Pfade die Felder hinter der letzten Fläche wie bei H und wandelt sie in Spektren
   (H → S, meist PFT); knotenweise wandelt die Ausführung an der Verbindung. Damit eine Verbindung
   in einen Raytatouille-Detektor S tragen kann, nehmen die Licht-Sockets importierter Systeme
   (`importer`, `topology.unroll`, `parts`) jetzt R, H und S; vor den übrigen Raytatouille-Blöcken
   (Strahlen) wandelt die Ausführung S → H → R, und auf einer S-Verbindung hinter einem Block
   bekommt das Bündel wie bei H das Medium hinter seiner letzten Fläche.
5. **Planer.** Detektoren nehmen S nativ; „auto“ wählt S aber nur, wenn auch der Block davor S
   nativ nimmt. Vor Detektoren stehen bisher nur Blöcke, die Strahlen rechnen, oder
   `fusion.convert` (RHS). Zwischen `fusion.convert` und einem Detektor ist S nun zulässig, die
   Kostenrangfolge (R 1, S 4) lässt „auto“ aber bei R; ein Spektrum aus `fusion.convert` mit
   `to = S` erreicht den Detektor dann unverändert (wie Felder auf einer R-Verbindung seit R30a),
   statt an seinem Port nach H gewandelt zu werden. Bestehende Pläne mit „auto“ bleiben unverändert
   (Linnik, Fokus: `test_existing_plans_are_unchanged`); S wählt „auto“ nur, wo R gar nicht
   zulässig ist (Quelle, die nur S abgibt, W79).
6. **Leistung im Bild** (`ray_power`): `spectrum_power` der kontinuierlichen Spektren, diskrete
   Wellen tragen keine endliche Leistung (0).

**Offen.** `SpectrumBundle` in den Port-Daten fehlt weiter: Ein Gitterspektrum trägt in
`GridSupport` eine Funktion (`Sampled.back`), die nicht über eine Prozessgrenze geht – nötig,
sobald ein Raytatouille-Detektor in einem anderen Prozess knotenweise S bekommt oder ein
Plugin-Block S nativ nimmt (R30c). Die Anzeige „geplant R, genutzt S“ für `fusion.convert` →
Detektor (`Context.used`) gehört zu R30e. Ein gekippter Detektor hinter einem großen
Gitterspektrum fällt auf die direkte Summe zurück und über `SUM_BUDGET` auf H.

## Nachtrag R30c (Umsetzung, 2026-11-18)

**Gebaut.** `fg.grating`, `fg.stack` und `fusion.condensed` nehmen Spektren nativ und geben
Spektren ab; mit H oder R am Eingang rechnen sie wie bisher. Gitter und Stapel
(`couplings.grating_field.spectrum_through`): das Spektrum auf die Blockebene (S4), jede Welle
für sich durch FourierGrating – derselbe Weg wie der Feldoperator (`_transfer_at`,
`apply_many`, derselbe Cache je Richtung), nur ohne Gitter –, Ordnung (n, m) bei
k∥ + n·K_x + m·K_y, je Ordnung ein eigenes Spektrum (Pfad `block[n]` wie bei Strahlen und
Feldern). Kontinuierliche Spektren bleiben kontinuierlich (Entscheidung 1: Δ²k unverändert,
dΩ = Δ²k/(k k_z) im Medium dahinter), diskrete bleiben diskret; die Wellen tragen ihre Ordnung
(`orders`). Ausgangsebene: Blockebene mit den Achsen des Eingangs, reflektiert (u, −v, −n).
Verdichtete Gruppe (`condense.apply_spectrum`): je Welle die Tabelle bei ihrem σ (bilinear, die
Toleranz der Tabelle), Phase k₀·OPL + k_out·(o_x − o_e − Versatz), das ist bei einer Platte die
Übertragung e^(ik_z′d) der ebenen Welle. Diagnose `exec.spectrum_rays` „Ordnungen nicht
getrennt“ (`spectrum.overlapping_orders`) bei S → R in `Executor._convert`: Spektren einer
Kohärenzgruppe mit verschiedenen Ordnungen, deren wesentliche Wellen sich in k∥ überlappen.
Abnahme W81–W83 in `tests/test_spectrum_blocks.py` (Messwerte in [testfaelle.md](testfaelle.md)):
zwei Gitter als Kette gegen FourierGrating mit beiden Gittern in einer Struktur 5,9·10⁻⁵ (die
Mehrfachreflexionen zwischen den Gittern, ∝ Δn²: 5,7·10⁻⁶ bei einem Drittel des Kontrasts),
Gauß-Strahl durch zwei Gitter gegen die H-Kette 5,8·10⁻⁴, Σ Ordnungen = R + T auf 4·10⁻¹⁵,
Stapel gegen Airy je Welle 1,9·10⁻¹³, verdichtete Platte gegen die Airy-Teilsumme 2·10⁻¹².

**Entscheidungen und Abweichungen.**

1. **Abtastung (`support`) als Daten.** `GridSupport` hält jetzt Gitter (Mitte, Schritt, Größe,
   Träger), Bins und die Stützstellen, auf die S → H zurückgeht (das `GridlessField` mit ψ,
   Gradient und Hesse-Matrix bzw. das `GriddedField` mit Krümmung) statt `Sampled.back`;
   `GridSupport.to_field` baut das Feld daraus. Damit gehen Spektren als Port-Daten über die
   Prozessgrenze (offener Punkt aus R30a/R30b, Rundreise bitgenau, Gitter im FourierGrating-
   Plugin = im Prozess). Die Eigenschaft `sampled` bleibt als Sicht für bestehenden Code.
2. **Abtastung hinter einer Ordnung** (`spectrum.order_support`): PFT-Spektren behalten ψ̃, −ρ
   und Hesse-Matrix auf den verschobenen Punkten; Gitterspektren behalten ihre Bins, das Gitter
   bekommt den Träger k_c + n·K und die Stützstellen die Phase ψ + n·K·x – wie beim Feldoperator
   (Ordnungen auf denselben Bins, kein mit der Periode kommensurables Fenster nötig). Spiegelung
   bei Reflexion: v, Gradient, Hesse-Matrix und Bins (κ_v → −κ_v).
3. **Ausbreitung eines Gitterspektrums** (`propagate`, Erweiterung von S4): zwischen parallelen
   Ebenen behält ein Gitterspektrum seine Bins; das Gitter wandert mit dem Licht, die Stützstellen
   kommen geometrisch von der alten Ebene (H → R → H wie `onto_plane`), solange der Strahl im
   Fenster bleibt (Ausdehnung · `PAD` ≤ Fenster), sonst entfällt die Abtastung mit Vermerk. S → H
   bleibt dabei die exakte inverse FFT der verschobenen Wellen.
4. **Gekippte Ebenen.** Die Port-Ebene einer Verbindung steht senkrecht zum Licht, die Ebene eines
   Gitters oder Stapels oft nicht. Gedrehte Wellen lägen nicht mehr auf dem Gitter ihrer
   Abtastung; ein Spektrum mit Abtastung wird deshalb vor dem Block über H auf dessen Ebene
   gebracht (S → H → S, wie Felder), mit Vermerk; ohne Abtastung (diskrete Wellen) exakt gedreht.
   Exakt im Spektrum ist die Kette also zwischen parallelen Ebenen – Gitter → Stapel, Gitter →
   Gitter.
5. **Befund W81 (b):** mit 20 µm zwischen den Gittern weicht die H-Kette um 2,7·10⁻³ ab, genau
   die Beugung z/z_R (2,9·10⁻³), die ihr geometrischer Schritt zwischen den Blöcken auslässt
   (`onto_plane` wählt Wellenausbreitung nur an Kaustiken); S breitet exakt aus. Bei einer
   senkrecht einfallenden Gauß-Strahl 1 mm vor dem Gitter waren es 14 %: dort ist S die genauere
   Kette.
6. **Interpolierte Übertragung auch für S.** `field_transfer = interpolated` gilt für
   Gitterspektren (die verifizierte Interpolation auf ihrem Gitter, `interpolated_transfer`, im
   Rahmen des Spektrums ausgewertet und in den des Gitters gedreht); sonst und wenn die Prüfung
   scheitert exakt je Welle. Ohne sie kostet ein Strahl mit 64² Wellen je Gitter rund 10⁴
   FourierGrating-Lösungen.
7. **Verdichtete Gruppe:** nativ R und S (keine H-Wege; Felder kommen als Strahlen). Hinter ihr hat
   das Spektrum keine Abtastung: Der seitliche Versatz in der Gruppe steckt in den Phasen, nicht
   in einer Abtastung – S → H dahinter ist offen; der Detektor summiert die Wellen.
8. **Planer und bestehende Pläne.** Sockets von Gitter, Stapel und verdichteter Gruppe nehmen
   RHS, die Auswerter S nativ: „auto“ lässt S zwischen ihnen (und zum Detektor) zu. Die
   Kostenrangfolge (R 1 < H 3 < S 4, Umwandlung S → R 1) behält R überall, wo R zulässig ist, auch
   hinter einem Block, der S bekommt; kommt ein Spektrum an einen nativen Block, bleibt es auf der
   Verbindung mit „auto“ ein Spektrum (wie seit R30b vor dem Detektor). Geprüft: die Pläne des
   Linnik mit Gitterprobe (`linnik_mini`, mit und ohne H am Gitter) und der Beispiele
   `auflicht_si_gitter` und `auflicht_si_gitter_100x` sind auf `main` und hier gleich (nur R
   bzw. dasselbe H mit denselben Umwandlungen) – keine Planänderung.

**Offen.** S → H hinter einer verdichteten Gruppe; die Phase des Übertrags eines dicken Stapels
steckt im Residuum, nicht in ψ̃ (für die inverse PFT eine Näherung wie beim Feldoperator); die
Anzeige „geplant R, genutzt S“ (R30e).

## Nachtrag R30d (Umsetzung, 2026-11-19)

**Gebaut.** `fusion.particle` nimmt Spektren nativ (Auswerter R und S, Sockets RHS) und gibt das
rückgestreute Fernfeld als Spektrum ab; Strahlen am Eingang rechnen wie bisher (Debye je Strahl,
Strahlen vom Mittelpunkt). Eingang (`particle.incident_waves`): das Spektrum im Spektrum auf den
vorderen Port von CliffordBem verschoben (S4, exakte Phase e^(ik·Δr) je Welle, ohne inneres
R → S), Achse des Teilchens = Normale des Spektrums, diskrete Wellen bleiben diskret.
Ausgang (`particle.scattered_spectrum`): die Wellen des Ports (Richtung −1) auf den Mittelpunkt
bezogen, im gespiegelten Rahmen (u, −v, −n) wie hinter einem Gitter (R30c), mit der Abtastung
einer Punktquelle im Mittelpunkt (ψ̃ = 0, −ρ = 0, Hesse-Matrix 0; das Residuum ist das
Fernfeld) und mit `spectrum.propagate` auf eine Ebene `OUT_WAVELENGTHS` = 100 Wellenlängen vor
dem Teilchen verschoben (ψ̃ = k·Δr, −ρ = Strahlversatz, Hesse-Matrix regulär). S → H dort ist die
inverse PFT der Kugelwelle E∞ e^(ikr)/r, S → R gibt Strahlen, die vom Teilchen weglaufen.
Abnahme W84–W85 in `tests/test_particle_spectra.py` (Messwerte in [testfaelle.md](testfaelle.md)):

| Prüfung | Referenz | Toleranz | gemessen |
|---|---|---|---|
| W84 (a) Einfallsspektrum einer S-Verbindung (NA 0,8): \|Faktor\| / Phase / Feld | Richards–Wolf (`couplings.focus`) wie W29 | 1e−3 / 1e−6 / 2e−3 | 6,5e−5 / 4e−12 / 9,3e−4 |
| W84 (a) dasselbe gegen den Debye-Weg je Strahl; Leistung gegen die Strahlen | W29 (`particle.rays_to_spectrum`) | 1e−9; 1e−3 | 5e−13; 0 |
| W84 (a) diskrete Welle auf den Port verschoben | Phase der ebenen Welle | 1e−12 | 0 |
| W84 (b) ebene Welle (diskretes S) auf die Kugel, Netzstufe 6 / 8 / Richardson | Mie | – / 5e−2 / 5e−3 | −6,1e−2 / −3,6e−2 / −3,3e−3 (Strahlweg W31: 7e−3) |
| W84 (b) dieselbe Rückstreuung über den Strahlweg (Fernfeld-Quadratur) | W31 | 1e−9 | 5e−15 |
| W85 (a) Dipol-Fernfeld als S, S → H in 0,1 und 1 mm | Kugelwelle E∞ e^(ikr)/r | 1e−9 | 3,4e−12 |
| W85 (a) dasselbe gegen das exakte Dipolfeld, Abweichung / (1/(kD)) | Nahfeldterm | 1 ± 0,02 | 0,9991 |
| W85 (a) Leistung; Richtung der Stützstellen; Strahlen durch den Mittelpunkt | (4π/3) k⁴\|p\|²; Wellenvektoren | 1e−12; 1e−9; 1e−9 | 2e−16; 3e−16; 2e−13 |
| W85 (a) Leistung der Strahlen in der NA (S → H → R) | Leistung der Wellen | 1e−2 | 3,9e−3 (Delaunay-Flächen) |
| W85 (b) K3 mit S auf beiden Verbindungen: Anteil in der NA (3,92 %) | CliffordBem allein (`BeamField.focused`) wie W30 | 5e−3 | 4,3e−4 (Strahlweg: 2,2e−3) |
| W85 (b) derselbe Anteil gegen den Strahlweg | W30 | 5e−3 | 1,8e−3 |
| W85 (b) Bild 1 mm vor dem Teilchen je Einfallsleistung | \|E∞\|²/r² von CliffordBem allein | 1e−2 der Spitze | 5,5e−3 |

**Entscheidungen und Abweichungen.**

1. **Ausgangsebene 100 λ vor dem Teilchen, nicht die gespiegelte Ebene des Einfalls.** Die
   inverse PFT der Punktquellen-Abtastung ist für jede Entfernung D > 0 genau die Kugelwelle
   (W85: 3·10⁻¹²); was ihr zum vollen Feld fehlt, der Nahfeldterm (1/(kD)) und die evaneszenten
   Wellen, erreicht ohnehin kein Objektiv. D muss also nur positiv sein und vor dem nächsten
   Block liegen. Die Ebene, aus der das Einfallsspektrum kam (durch den Schwerpunkt der Strahlen
   hinter der letzten Fläche), kann im Pfeil einer gekrümmten Fläche liegen, und ein Detektor
   zwischen ihr und dem Teilchen bräuchte eine Ausbreitung von H rückwärts, die der Planer von
   Fieldtracing nicht kann (`math domain error` im ersten Versuch). 100 λ sind 60 µm bei 0,6 µm.
2. **Strahlen bleiben Strahlen, Spektren bleiben Spektren.** Mit R am Eingang bleibt der
   Debye-Weg (Entscheidung 2: Debye direkt nur im Fokus) und die Ausgabe Strahlen vom
   Mittelpunkt, mit S am Eingang wird die Ausgabe S. H auf der Verbindung in das Teilchen wird wie
   bisher zu Strahlen (der Auswerter nimmt R und S; H → S im Block wäre möglich, änderte aber die
   Ergebnisse bestehender Graphen mit festem H).
3. **Eingang.** Die Abtastung des Einfalls wird nicht gebraucht und vor dem Verschieben
   weggelassen (`propagate` rechnete sonst für ein Gitterspektrum eine geometrische Abtastung auf
   der neuen Ebene). Am Port evaneszente Wellen fallen mit `exec.spectrum_cut` weg (CliffordBem
   nimmt nur laufende). Ein Spektrum in einem anderen Medium als `n_medium` lehnt der Block ab
   (`exec.particle`): die Grenzfläche dazwischen ist nicht modelliert.
4. **W29 mit S ist der Debye-Weg** (5·10⁻¹³): die PFT des geometrischen Felds eines fokussierten
   Bündels hat ihren stationären Punkt in der Richtung des Strahls, mit dessen Amplitude.
5. **W31 mit S:** Das Spektrum zählt die Leistung mit der Fernfeld-Quadratur (Gauß–Legendre),
   `ray_powers` des Strahlwegs mit den Delaunay-Flächen der Richtungskosinus (Netzstufe 6:
   1 % weniger). Mit derselben Quadratur sind beide Wege gleich (5·10⁻¹⁵); extrapoliert liegt S
   näher an Mie (3,3·10⁻³ statt 7·10⁻³). Zusätzlich geprüft: extrapoliert 5·10⁻³.
6. **Detektor hinter dem Teilchen.** Weit vom Teilchen faltet die Summe der Wellen immer (die
   fast streifenden Richtungen landen bei D·tan θ, bis etwa 550·D), das Bild geht über H. Das
   Feld der ganzen Halbkugel ließ den Planer von Fieldtracing „mixed“ wählen (F_min = 0,55 an den
   äußersten Stützstellen bei 68 mm) und Gigabytes anlegen. Deshalb fallen vor S → H im Detektor
   die Wellen eines PFT-Spektrums weg, deren Strahlen mehr als `REACH` = 10 halbe Diagonalen an
   der Mitte des Fensters vorbeigehen (`spectrum.towards`, mit Vermerk und Leistungsanteil; Dipol,
   Detektor 1 mm davor, Fenster ±0,6 mm: 720 von 3200 Wellen, 9 % der Leistung); danach geometrisch. Bestehende Bilder ändern sich nicht (alle
   Wellen von W79/W80 erreichen das Fenster; `test_spectrum_detector.py` unverändert grün).
   Ebenso lässt `field._onto` (gekippter Detektor) Stützstellen weg, deren Strahlen die Ebene
   nicht vorwärts kreuzen, statt abzubrechen.
7. **Offener Punkt aus R30c (S → H ohne Abtastung).** Es gab keinen Absturz, aber eine falsche
   Meldung („diskrete Wellen“) als `exec.field`. Jetzt `spectrum.SupportError` und in der
   Ausführung `exec.spectrum_support` mit dem Block, aus dem das Spektrum kommt, und dem Ausweg (H
   oder R auf die Verbindung in diesen Block; der Detektor nimmt das Spektrum weiter als Summe).
   Ein Rückfall auf ein Feld bräuchte eine Abtastung aus den Wellen selbst (Spektralphase mit
   Gradient und Hesse-Matrix aus den Nachbarn, also Phasenentfaltung auf verstreuten k-Punkten) –
   nicht billig, deshalb nur die klare Meldung.
8. **Planer und bestehende Pläne.** Die Sockets nehmen RHS, der Auswerter S nativ: „auto“ lässt S
   zwischen Teilchen und nativen Blöcken (Detektor, Gitter, Stapel, verdichtete Gruppe) zu, die
   Kostenrangfolge behält R. Geprüft auf `main` und hier: K3 (Quelle → Linse → Teilchen →
   Detektor), Gitter → Teilchen → Stapel → Detektor und K3 mit festem H haben dieselben
   Darstellungen und Umwandlungen; nur die Begründung der Verbindungen zwischen Teilchen und
   nativen Blöcken lautet „among RHS“ statt „among RH“. **Keine Planänderung.** Ein Spektrum, das
   auf einer „auto“-Verbindung ankommt, bleibt ein Spektrum (wie seit R30b/R30c).
9. **SpeckleBem-Block zurückgestellt.** Die Bedingung der Stufe („sobald L4 ihn mit vielen
   Einfällen hat“) ist nicht erfüllt: `specklebem.apply_many` ist eine Schleife über `apply`, und
   `scatter` löst je einfallender ebener Welle einmal (`set_excitation` + GMRES, 0,2–22 s je
   Lösung nach plan-l4-operatoren.md). Ein fokussierter Strahl mit rund 2800 Wellen wären 2800
   Lösungen. Nötig ist ein Einfall aus einem beliebigen Spektrum als eine Anregung (D12,
   `AngularSpectrumBeam` kann bisher nur Gauß) oder Mehrfach-RHS mit behaltener LU (O3,
   optional). SpeckleBem ist zudem nur in der CI gebaut (Commit `4f04efe`), nicht lokal.

**Offen.**

- **H → R und elliptische Polarisation:** `field.field_to_rays` behält von einem elliptisch
  polarisierten Feld nur den linearen Anteil (Projektion auf die reelle Polarisationsrichtung).
  Die Rückstreuung der Kugel ist neben den Achsen elliptisch (S1 und S2 mit verschiedener
  Phase); ein gegen die Spektrumsebene um 1,5 mrad gekippter Detektor (H → R → H in `_onto`)
  verliert deshalb an den Diagonalen bis 13 %. Der Strahlweg des Teilchens behält die Leistung,
  macht die Polarisation aber linear. Beides bestand vor R30d; eine komplexe Polarisation in
  `prt` würde es beheben (Frage an den Maintainer, betrifft R18).
- Mit Strahlen am Eingang und festem S (oder H) auf der Ausgangsverbindung scheitert R → H an
  den Strahlen, die alle im Mittelpunkt starten (wie vor R30d mit H).
- SpeckleBem-Block (Punkt 9); die Anzeige der genutzten Darstellung (R30e).

## Nachtrag R30e (Umsetzung, 2026-11-19)

**Gebaut.** RTT-Fusion: `Executor._convert` trägt nach jeder Umwandlung an einem Port in
`Context.used` ein, was die Lichtverbindung getragen hat (`Executor._record_used`):
`representation` – die Darstellung *auf der Verbindung*, also nach der Umwandlung am Ausgang der
Quelle (R, H → S auf einer S-Verbindung; R, S → H auf einer H-Verbindung) und vor der am Port des
Verbrauchers; auf einer mit R geplanten Verbindung kommen Felder und Spektren unverändert an –,
`chain` – die Schritte vom Ausgang der Quelle bis zu dem, womit der Verbraucher rechnet, R und S
immer über H („R→H→S“, „H→S“, „S→H→R“, „S“) – und `reason` („lens gave rays; spectra on the link;
cam takes them as they are“ bzw. „… converted to rays at the port of lens“). Eingetragen wird,
wo nicht nur Strahlen über eine R-Verbindung gehen, und auf jeder S-Verbindung. Der Detektor
nahe einem Fokus überschreibt den Eintrag seiner Verbindung danach wie bisher (R→H).
`workflow.with_used` gibt die Kette als `links[key]["used_chain"]` und je Verbindung und Richtung
eines Modells als `connections_used_reason` („Kette – Grund“) weiter. Tests in
`tests/test_spectrum_links.py`: Fokus mit S (R→H→S, nativ am Detektor), `fusion.convert` →
Detektor auf „auto“ (geplant R, getragen S), S in einen Strahlenblock (R→H→S→H→R),
`with_used`.

CAD (Raytatouille-Fusion #138, Branch `feat/r30e-spectra`):

- **Transportmenü** einer Verbindung bzw. Richtung: jede Darstellung mit Beschreibung („S –
  plane-wave spectra“), die aktuelle abgehakt, nicht rechenbare ausgegraut mit Grund im Tooltip.
  Was rechenbar ist, schreibt das Workflow-Plugin beim Prüfen und nach dem Lauf als `available`
  (`workflow.available()`) ins Dokument, der Raytracer beim Anlegen eines Modells; ohne den
  Schlüssel (ältere Dokumente) bleibt alles wählbar.
- **Anzeige:** Auf einer Verbindung auf „auto“ steht nach dem Lauf geplant→getragen, wo es
  abweicht (R→H wie bisher, neu R→S); der Tooltip nennt in der Zeile „last run“ Darstellung, Kette
  und Grund, auch bei festem Transport (dort bleibt das Schild beim Stern, z. B. „S*“). Die Zeile
  erscheint, wo die Darstellung vom Plan abweicht, wo umgewandelt wurde oder wo ein Spektrum floss
  (auch nativ). Gilt solange das Ergebnis des Blocks am Ende aktuell ist (wie R→H).
- **Diagnosen:** Das Workflow-Plugin liest die Zeilen `<Ort>: exec.spectrum_*: …` aus
  `Context.log` (Ort = Verbindung oder Block; eine Blockmeldung hinter einer Verbindung zählt für
  die Verbindung) und `exec.spectrum_support` aus dem `ExecutionError`, hängt sie an den Plan
  (`links[key]["notes"]`, `node_notes`, `connections_notes` je Richtung) und nennt die erste
  Warnung in der Statuszeile. Der Editor setzt ein ⚠ hinter die Darstellung bzw. in den Kopf des
  Blocks, die Meldungen in den Tooltip.
- Die Prüfung eines Socket-Graphen im Raytracer plant über `workflow.checked` (S verfügbar,
  native Darstellungen der Blöcke) statt mit {R, H} ohne `native`.
- Handbuch (Kapitel Node-Workflow): Darstellungen R/H/S/G/T/P, wann S lohnt und welche Blöcke es
  nativ nehmen, „auto“ wählt S nur, wo R nicht zulässig ist, Umwandlungen an den Ports, Anzeige
  R→S und „last run“, Transportmenü, Meldungen `exec.spectrum_cut`/`_rays`/`_support`; Anhang
  Workflow-Blöcke: Sockets RHS, `fusion.convert` mit `to = S` (Standard H).

**Entscheidungen und Abweichungen.**

1. **„Genutzt“ heißt: auf der Verbindung**, nicht „womit der Verbraucher rechnet“. Die
   Umwandlung am Port eines Blocks, der S nicht nimmt (S → R vor einer Linse), steht schon als
   Umwandlung im Plan; wichtig ist die Abweichung vom Plan selbst – auf „auto“ mit R geplant, aber
   ein Spektrum geflossen. Womit der Verbraucher rechnet, steht im letzten Glied der Kette.
2. **Diagnosen bleiben Zeilen im Protokoll** (Entscheidung 7 aus R30a); die CAD-Seite liest
   sie aus dem Text. Ein `Diagnostic` je Umwandlung mit Ort wäre sauberer, ändert aber alle
   Aufrufer von `notes`; zurückgestellt, bis eine zweite GUI sie braucht.
3. Das Schild zeigt geplant→getragen nur auf „auto“; ein fester Transport behält sein Schild mit
   Stern (der Planer hält ihn), die Kette steht im Tooltip.

**Offen.**

- **Abnahme am Linnik mit Si-Gitter** (`auflicht_si_gitter`: Bild mit S am Gitter = mit H, Toleranz
  nach GF5): nicht gerechnet – braucht FourierGrating mit vielen Lösungen je Welle und gehört auf
  den Server. Die Planprüfung aus R30c (Pläne mit „auto“ unverändert) gilt weiter.
- Die Meldungen von Blöcken in anderen Prozessen kommen als Protokollzeilen herüber; ihr Ort ist
  der Block in diesem Prozess, wie ihn die Zeile nennt.
