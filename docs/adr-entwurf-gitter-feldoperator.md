# ADR: Gitter als Feldoperator und Schichtstapel (`fg.*`)

Stand: **angenommen** 2026-11-13 – der Maintainer folgt den Empfehlungen zu allen sechs Fragen
(Tabelle am Ende). Anlass: CAD ADR 0010, Stufe W4, nennt als Blöcke des FourierGrating-Plugins
„Gitter als Jones-Tabelle, als S-Matrix-Operator, Stapel“. Die Jones-Tabelle gibt es (`fg.grating`,
R19, W72/W73); dieses ADR legt die beiden anderen fest. Reihenfolge: F4 (`fg.stack`) zuerst, dann
F1–F3.

Bezüge: Fusionsplan 6.2 (Darstellungen S, M), 6.3 (Umwandlungen), 6.4 (`operator`);
[plan-l4-operatoren.md](plan-l4-operatoren.md) (O3: `Simulation.smatrix` in der s/p-Basis der
Ordnungen); [adr-entwurf-beugungstabelle.md](adr-entwurf-beugungstabelle.md) (Tabelle, Phasenbezug,
Rayleigh-Anomalien); [adr-entwurf-ports-posen-aeste.md](adr-entwurf-ports-posen-aeste.md) (Ports,
Streutabelle); [konventionen.md](konventionen.md) (Zeitfaktor e^(−iωt), Port-Konventionen); Pilot K2
(Metafläche bzw. DOE im Feldpfad).

## Kontext

**Heute (`fg.grating`, Darstellung R).** Der Block behandelt jeden Strahl als lokale ebene Welle:
Richtung aus der Gittergleichung, Feld aus der Jones-Matrix der Ordnung, interpoliert über (λ, σ);
die Phase der Ordnung m am Auftreffort x ist m·2πx/Λ. Die Tabelle baut FourierGrating über
`apply_many` mit seiner S-Matrix. Das ist richtig, solange das Licht auf der Gitterfläche lokal eine
ebene Welle ist – die beleuchtete Fläche groß gegen Λ, die Wellenfront über viele Perioden glatt.

**Was damit nicht geht.**

1. *Felder mit Struktur auf der Skala weniger Perioden*: ein Fokus auf dem Gitter (Linnik mit
   hoher NA, Spot ≈ Λ), Strahlen mit wenigen beleuchteten Perioden. Dann überlagern sich die
   Ordnungen im Nahfeld kohärent, und die Strahlnäherung verliert die Interferenz der benachbarten
   Ordnungen innerhalb eines Strahlbündels.
2. *Kohärente Ausbreitung hinter dem Gitter*: Talbot-Bilder, Selbstabbildung, die Überlagerung der
   Ordnungen in kurzem Abstand (Gitter vor einem zweiten Element) – dafür braucht der nächste Block
   ein Feld (H), nicht Strahlfamilien je Ordnung.
3. *Schichtstapel ohne Gitter* (Beschichtungen, Filter, Spiegelschichten einer Probe): heute nur über
   `fusion.condensed` aus Raytatouille-Flächen oder gar nicht; FourierGrating rechnet sie exakt
   (N = 0), aber kein Block bietet das an.

**Was vorhanden ist.** Die Darstellungen S (`ports.PlaneWaves`, ebene Wellen mit Ordnungen) und G
(`ports.GriddedField`) mit G ↔ S (FFT), H (`graph.field.FieldBundle`) mit R ↔ H; der Planer kennt die
Operatoren `fft`, `inverse_fft`, `pft`, `orders_to_rays`, ausführbar sind auf Verbindungen bisher nur
R und H (`execute.WITH_FIELDS`). FourierGrating liefert je Einfallsrichtung die S-Matrix aller
Ordnungen in der s/p-Basis (`Simulation.smatrix`, Plan eq:Ssp); die Wiederverwendung über eine
Bloch-Familie ist wegen der um die Einfallsordnung zentrierten Abschneidung nur genähert (≈ 7·10⁻³
bei N = 6, L4/O3).

## Entscheidung

### 1. Ein Gitterblock, zwei Darstellungen

`fg.grating` bleibt der eine Gitterblock (gleiche Sockets, gleiche Ports, Austausch gegen einen
Spiegel unverändert). Sein Evaluator nimmt künftig **R oder H** an (`representations = {"R", "H"}`):

- **R**: wie heute, Jones-Tabelle je Ordnung (unverändert, bitgleich zu W72/W73).
- **H**: der Feldoperator nach Abschnitt 2.

Welche Darstellung die Verbindung trägt, entscheidet wie überall der Planer bzw. der Nutzer über den
Transport der Verbindung (Rechtsklick → Transport → H). „auto“ bleibt R, bis ein Gültigkeitskriterium
(Abschnitt 5) die Wahl automatisch treffen kann.

Begründung: Der Nutzer tauscht eine Probe, nicht eine Rechenmethode; die Methode gehört an die
Verbindung (Konzept: „Die Darstellung auf einem Band wählt der Planer bzw. der Nutzer“). Ein zweiter
Blocktyp (`fg.grating_field`) würde Workflows bei jedem Methodenwechsel umbauen.

### 2. Der Feldoperator: kommensurables Fenster, S-Matrix je Spektralpunkt

Für jede Kohärenzgruppe des einfallenden Feldes (H, in der Port-Ebene des Gitters, Rahmen:
x = Gittervektor, y = Furchen, z = Flächennormale wie in `fg.grating`):

1. **H → G** auf ein Fenster der Breite W_x = q·Λ_x (q ganzzahlig) bzw. W_y = r·Λ_y bei
   Kreuzgittern, das die beleuchtete Fläche samt Rand enthält; Abtastung nach den Regeln von
   Fieldtracing ADR 0003 (Restfeld). W ist **kommensurabel** mit der Periode: Dann ist der
   Gittervektor K = 2π/Λ ein ganzes Vielfaches (q) des Spektralabstands 2π/W, und jede Ordnung eines
   Spektralpunkts fällt wieder auf einen Spektralpunkt – **keine Interpolation im Spektrum**.
2. **G → S** (`ports.to_plane_waves`): Spektralpunkte k∥ = (k_x, k_y) mit ihren Amplituden; Punkte
   unterhalb einer Schwelle (relativ, Parameter `spectrum_threshold`, Standard 10⁻⁸) entfallen, ebenso
   evaneszente des Einfalls.
3. **S → S mit M**: je verbleibendem Spektralpunkt die S-Matrix von FourierGrating für die Richtung
   dieses Punkts (λ, φ, ψ aus k∥), angewandt auf (E_s, E_p) des Punkts; die Ordnung n landet auf
   k∥ + n·K. Ausgänge desselben k∥ aus verschiedenen Einfällen addieren sich kohärent. Das ist
   `adapters.fouriergrating.apply_many` mit dem Spektrum als Einfallsliste – derselbe Code wie für die
   Tabelle.
4. **S → G → H** je Ausgangsport (reflektiert vorn, transmittiert hinten): inverse FFT auf dasselbe
   Fenster, dann das Restfeld als `FieldBundle` mit der Phase des Gitterbezugs (x = 0, wie
   Abschnitt 4 der Beugungstabelle).

Die Darstellungen S und G bleiben **innerhalb des Blocks**; auf Verbindungen kommen weiter nur R und
H vor. S und G als Verbindungsdarstellungen sind eine spätere Stufe (F5), wenn mehrere Blöcke
Spektren austauschen (CliffordBem, SpeckleBem).

**Genauigkeit und Kosten.** Je Spektralpunkt eine Lösung von FourierGrating (die S-Matrix beider
Polarisationen auf einmal). Standard ist die **exakte** Lösung je Richtung. Opt-in `bloch = true`
nutzt eine Lösung je Bloch-Familie (k_x modulo K): q statt aller Punkte je Zeile, mit dem
Abschneidefehler der verschobenen Ordnung (L4/O3); der Block meldet ihn als Warnung mit einer
Schätzung aus zwei Abschneidungen. Die Operatoren landen im Operator-Speicher (Schlüssel: Struktur,
λ, k∥), so dass Optimierung und wiederholte Läufe sie wiederverwenden. Die Kostenschätzung des
Plugins (*Estimate*) zählt Spektralpunkte × Lösungen.

### 3. Schichtstapel: `fg.stack`

Ein neuer Block **`fg.stack`** für lateral invariante Schichtfolgen (homogene Schichten, auch
absorbierend und anisotrop, wie FourierGrating sie kennt) zwischen zwei Medien – Beschichtung,
Filter, Spiegelschicht. Sockets und Ports wie ein Spiegel bzw. eine Platte (`front` mit Reflexion,
`back` mit Transmission, beides wählbar). FourierGrating rechnet ihn mit N = 0 exakt; der Operator
ist diagonal in k∥, also

- **R**: je Strahl die Koeffizienten (r_s, r_p, t_s, t_p) für seine Richtung und λ, interpoliert aus
  einer Tabelle über (λ, σ) wie in `fusion.condensed` (die Tabelle aus derselben Rechnung);
- **H**: je Spektralpunkt dieselben Koeffizienten, ohne Fenster-Bedingung (keine Ordnungen).

Ein Stapel mit Gittern **gleicher oder kommensurabler Periode** ist kein eigener Block: Er ist eine
FourierGrating-Konfiguration mit mehreren Schichten und gehört in `fg.grating` (die Konfiguration
darf beliebig viele Gitterschichten und Abstandsschichten enthalten, die evaneszente Kopplung ist dann
exakt). Gitter mit **inkommensurablen** Perioden lassen sich nicht in eine Einheitszelle legen; sie
bleiben zwei Blöcke mit freier Ausbreitung (H) dazwischen.

### 4. Ports, Rahmen, Konventionen

Unverändert gegenüber `fg.grating`: Port `front` an `point` mit `normal`, `grating_vector` als lokale
x-Achse; Zeitfaktor e^(−iωt); s ⟂ Einfallsebene, p = s × k (Beugungstabelle, Abschnitt 2);
Phasenbezug am Ursprung der Port-Ebene. Für `fg.stack` entfällt `grating_vector`. Rückwärtseinfall
(Licht von hinten) bleibt wie bei allen FourierGrating-Blöcken offen (Konzept A3).

### 5. Gültigkeit und Diagnosen

| Code | Schwere | Wann |
| --- | --- | --- |
| `fg.window_small` | Fehler | das Feld reicht über das Fenster hinaus (Leistung im Rand > 10⁻⁶) – der periodische Umlauf der FFT würde das Ergebnis verfälschen |
| `fg.spectrum_truncated` | Warnung | Spektralanteile außerhalb der Abschneidung von FourierGrating (|n| > N_x) tragen mehr als 10⁻⁶ der Leistung |
| `fg.rayleigh` | Warnung | ein Spektralpunkt liegt näher als die Toleranz an einer Rayleigh-Anomalie (wie Beugungstabelle 6) |
| `fg.bloch_error` | Warnung | mit `bloch = true`: geschätzter Abschneidefehler über der Toleranz |
| `fg.energy` | Warnung | verlustfreie Struktur: Σ P_aus ≠ P_ein über 10⁻⁸ |
| `fg.field_advised` | Warnung | R auf der Verbindung, aber weniger als etwa 20 Perioden beleuchtet oder merklich gekrümmte Wellenfront über eine Periode: H wählen |

Ein späteres Kriterium für „auto“ (R oder H): die beleuchtete Breite in Perioden und die
Krümmung der Wellenfront über eine Periode; bis dahin wählt der Nutzer.

### 6. Ort im Code

- RTT-Fusion: `couplings/grating_field.py` (Schritte 1–4, ohne Graph), `graph/grating.py` (Evaluator
  für H), `graph/stack.py` (`fg.stack`).
- FourierGrating: keine Änderung (`smatrix` vorhanden).
- CAD: keine Änderung am Kern; das FG-Plugin liefert die Typen und Evaluatoren aus RTT-Fusion wie
  heute `fg.grating`, `fg.stack` kommt automatisch hinzu.

## Verworfen

- **Eigener Blocktyp für den Feldoperator.** Siehe 1.
- **Spektrum interpolieren statt kommensurablem Fenster.** Die Ordnungen landen sonst zwischen den
  Spektralpunkten; Interpolation komplexer Amplituden mit schneller Phase ist genau das, was
  Fieldtracing ADR 0003 vermeidet.
- **S und G sofort als Verbindungsdarstellungen.** Erst wenn ein zweiter Block Spektren nimmt; bis
  dahin bleibt der Planer bei R und H, und die Abnahme ist kleiner.
- **Bloch-Familie als Standard.** Schneller, aber mit einem Fehler, der von der Abschneidung abhängt
  und nicht von der Rundung; die Regel „Bitgleichheit oder erklärte Toleranz“ (D8) spricht für exakt
  als Standard.

## Abnahme (Sollwerte unabhängig)

| Nr. | Fall | Sollwert | Toleranz |
| --- | --- | --- | --- |
| GF1 | eine ebene Welle (ein Spektralpunkt) durch den Feldoperator | Effizienzen und Amplituden von `fouriergrating.solve` direkt | 1e−12 |
| GF2 | Gauß-Strahl auf ein verlustfreies Lamellengitter, Taille 5 … 50 Λ | Leistung je Ordnung → Effizienz der ebenen Welle für große Taille (Konvergenz), Σ P = P_ein | Konvergenzordnung; 1e−10 für die Bilanz |
| GF3 | Talbot-Selbstabbildung hinter einem binären Amplitudengitter, paraxial | Bild in z_T = 2Λ²/λ gleich dem Gitter, in z_T/2 um Λ/2 verschoben (Rayleigh 1881) | Kontrast und Lage auf 1e−3 |
| GF4 | breiter kollimierter Strahl: Feldpfad (H) gegen Strahlpfad (R, Jones-Tabelle) bis zur Kamera | gleiches Bild | Tabellentoleranz |
| GF5 | Linnik mit Phasengitter, kleine NA: H gegen R (W23); große NA: Unterschied ausgewiesen | W23 bei kleiner NA | 1e−3 |
| GF6 | Feldoperator im Plugin-Prozess (H über Port-Daten) | gleich im Prozess der Engine | 1e−12 |
| ST1 | λ/4-Antireflexschicht n = √n_s bei der Auslegungswellenlänge | R = 0 (analytisch) | 1e−12 |
| ST2 | Bragg-Spiegel aus N Paaren | Transfermatrix (analytisch, eigener Code im Test) für s und p über Winkel und λ | 1e−10 |
| ST3 | `fg.stack` mit R und mit H | gleiche Koeffizienten je Richtung | 1e−12 |

## Stufen

| Stufe | Inhalt | Abnahme |
| --- | --- | --- |
| F1 | `couplings/grating_field.py`: Fenster, H → G → S → M → S → G → H, ohne Graph | GF1, GF2, GF3 |
| F2 | `fg.grating` nimmt H; Planer und Ausführung; Diagnosen | GF4, GF6 |
| F3 | Linnik im Feldpfad | GF5 |
| F4 | `fg.stack` mit R und H (empfohlen als erste Stufe) | ST1–ST3, GF6 sinngemäß |
| F5 (später) | S und G als Verbindungsdarstellungen, wenn ein zweiter Block Spektren tauscht | – |

## Entscheidungen des Maintainers (2026-11-13: Empfehlungen übernommen)

| Nr. | Frage | Entscheidung |
| --- | --- | --- |
| 1 | Ein Block, zwei Darstellungen (Abschnitt 1) statt eines eigenen Feld-Gitterblocks? | **Ja.** Die Methode gehört an die Verbindung; ein zweiter Typ baut Workflows bei jedem Wechsel um. Die Feldparameter (Fenster, Schwelle, Bloch) bilden eine eigene Gruppe „Feldpfad“, die bei R unbeachtet bleibt. |
| 2 | Fenster automatisch oder vom Nutzer? | **Automatisch mit Überschreiben:** Träger der Feldleistung (ohne die letzten 10⁻⁶) plus Rand – mindestens 10 Perioden und der seitliche Versatz durch die Gitterdicke –, aufgerundet auf ein Vielfaches von Λ; `window_um` (0 = automatisch) überschreibt, zu klein ergibt `fg.window_small`. |
| 3 | Exakt als Standard, Bloch-Familie als Opt-in (Abschnitt 2)? | **Ja** (Regel D8). Größenordnung Linnik-Fokus NA 0,8, λ 0,6 µm, Fenster 20 µm: etwa 2 800 Spektralpunkte = Lösungen (Sekunden bis Minuten, im Operator-Speicher wiederverwendet); Bloch etwa halb so viele bei Λ = 0,7 µm, mit Fehler ≈ 10⁻³. *Estimate* zeigt die Zahl vorher; bei vielen Lösungen weist der Block auf Bloch hin. |
| 4 | Was heißt „Stapel“? | **`fg.stack`** für lateral invariante Schichtfolgen; Gitter gleicher oder kommensurabler Periode in einer Konfiguration von `fg.grating` (exakt, mit evaneszenter Kopplung); inkommensurable Gitter bleiben zwei Blöcke mit H dazwischen. |
| 5 | „auto“ bleibt R? | **Ja, mit Hinweis:** Diagnose `fg.field_advised` (Warnung), wenn weniger als etwa 20 Perioden beleuchtet sind oder die Wellenfront über eine Periode merklich krümmt. Automatisch umschalten erst, wenn GF5 zeigt, wo R und H tatsächlich auseinanderlaufen. |
| 6 | Reihenfolge? | **F4 zuerst, dann F1–F3:** `fg.stack` ist klein und sofort nützlich (beschichteter Referenzspiegel, AR-Schichten) und erprobt den Mechanismus „ein Block, Evaluator für R und H“ an einem diagonalen Operator ohne Fenster und Ordnungen; GF6 (Feld über die Prozessgrenze) wird schon dort geprüft. |
