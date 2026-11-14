# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/), Versionen nach SemVer.

## [Unreleased]

### Behoben

- Leistung eines Felds hängt nicht mehr von der Ringzahl ab, wenn eine Blende die abgetastete
  Pupille beschneidet (Nebenbefund in #69): Ein Feld reicht nur bis zu seinen äußersten
  Stützstellen, mit dem Ringraster allein also bis zu einem Ringabstand vor der wahren Kante. Für
  Felder (H) bekommen die Startstrahlen jetzt Punkte an der Kante (Bisektion zwischen benachbarten
  Strahlen, von denen einer durchkommt und einer verloren geht, `PointListPupil`). Beispiel: eine
  Blende, die 0,7² der Pupille freigibt, gab vorher 4/9 (Raster), jetzt 0,490 bei 6 und 12 Ringen.

- Meta-Blöcke wirken auf Blöcke anderer Prozesse (#70): Ein `meta.sweep`, `meta.replicates` oder
  `meta.optimise`, der einen Parameter eines `rtt.*`-Blocks ändert (Antrieb, Pose, Parameter),
  blieb wirkungslos, wenn die Blöcke ein anderer Prozess rechnet (in Raytatouille Fusion der
  Normalfall: der Raytracer-Plugin-Prozess) – der Lieferant rechnete mit dem einmal vorbereiteten
  Graphen. Jede Auswertung nennt jetzt die Blöcke, deren Parameter vom vorbereiteten Graphen
  abweichen; der Lieferant setzt seine Kopie vor dem Auswerten auf vorbereitete Parameter plus
  diese Änderungen. Abnahme W74.

- Felder (H) außeraxialer Feldpunkte fokussieren richtig (#69): Raytatouille startet die Strahlen
  eines Feldpunkts im Unendlichen mit OPL 0 auf ihrer Startebene, die nur auf der Achse ihre
  Wellenfront ist. Der Raytatouille-Backend gibt den parallelen Startstrahlen jetzt die Phase der
  ebenen Welle mit, n·d̂·(p − p_ref) mit dem Index des Objektraums; vorher widersprach ψ = k₀·OPL
  den Strahlrichtungen um etwa 90 rad über die Pupille, und Kamerabilder im Feldpfad lagen bis zu
  0,2 mm daneben mit 5- bis 25-fach zu niedriger Spitze. Die Quelle im knotenweisen Lauf trägt
  zudem den Brechungsindex des Objektraums.

### Geändert

- Felder werden bikubisch mit FITPACK (`spectral.bicubic`) auf ihre Stützstellen zurückgelegt;
  `RegularGridInterpolator(method="cubic")` konvergiert auf solchen Gittern nicht.

- Eine Verbindung mit Transport R liefert Strahlen, auch wenn der Block davor ein Feld ausgibt
  (vorher kam das Feld bei einem Empfänger, der Felder nimmt, als Feld an).

- Der Gitterblock heißt jetzt `fg.grating` (vorher `fusion.grating`; CAD ADR 0010, Entscheidung 3:
  Typkennungen tragen das Präfix ihres Plugins). In Raytatouille Fusion liefert ihn das
  FourierGrating-Plugin (W4). Ältere Dateien laden weiter: `io.TYPE_ALIASES` bildet alte
  Kennungen beim Laden auf die neuen ab, gespeichert wird die neue. Abnahme W72.

### Hinzugefügt

- Linnik im Feldpfad (ADR Gitter/Stapel, Stufe F3, Abnahme GF5): Strahlen werden vor einer
  Gitterprobe zum Feld und danach wieder zu Strahlen, im knotenweisen Lauf; `FieldBundle` trägt
  dafür den Ast der Strahlen (`branch`). Der Feldoperator rechnet die Übertragung je Richtung
  einmal und hält sie für den Lauf (die zweite Polarisation kostet nichts); der interpolierte
  Übertrag gibt bei einem Knick im Spektrum (Rayleigh-/Wood-Anomalie) früh auf.

- `fg.grating` nimmt Felder (ADR Gitter/Stapel, Stufe F2): mit einer FourierGrating-Konfiguration
  rechnet der Block jede ebene Welle des Feldes und gibt jede Ordnung der gewählten Seite als
  eigenes Feld aus (Kanal `[m]` wie bei Strahlen); mit nur einer Jones-Tabelle wird ein Feld wie
  bisher an seinem Port zu Strahlen. Der Typ nennt jetzt auch `fouriergrating`, `side`, `n_in` und
  die Feldparameter. Abnahme GF3, GF4, GF6.

- Gitter als Feldoperator, ohne Graph (ADR Gitter/Stapel, Stufe F1): `couplings.grating_field`
  bringt ein Feld in die Ebene eines periodischen Elements, zerlegt es in ebene Wellen, rechnet
  jede mit FourierGrating und gibt je Seite und Ordnung ein eigenes Feld mit dem Träger k_c + nK
  zurück (Nachtrag 1 zum ADR). `fg.stack` nutzt denselben Weg. Neue Wahl `field_transfer =
  "interpolated"`: Übertrag von einem geprüften Teilgitter des Spektrums, für breite Strahlen
  etwa zehnmal schneller. Abnahme GF1, GF2.

- `fg.stack` nimmt Felder (H, Stufe F4b): das Feld wird in die Ebene des Stapels gebracht
  (geometrisch, bei einem Fokus auf dem Stapel durch Wellenausbreitung), ohne seinen Träger auf
  ein Gitter gelegt, in ebene Wellen zerlegt, jede signifikante Welle exakt von FourierGrating
  durchgerechnet und auf die eigenen Stützstellen des Feldes zurückgelegt (ψ bleibt).
  `graph.spectral` hält das für den Feldoperator des Gitters (F1) bereit. Abnahme ST3, GF6.

- Block `fg.stack` (ADR Gitter als Feldoperator und Schichtstapel, Stufe F4a): lateral invarianter
  Stapel homogener Schichten (Beschichtung, Filter, Spiegelschichten) zwischen zwei Medien,
  exakt je Strahlrichtung von FourierGrating gerechnet (N = 0, keine Tabelle); Sockets wie
  `optics.interface`, Licht von beiden Seiten. Die Schichten stehen wie bei FourierGrating vom
  Substrat (hinten) nach oben. Abnahme ST0–ST2 gegen Transfermatrix, Airy-Rekursion und
  `optics.interface`.

- Abnahme W73 (CAD ADR 0010, W4): der Linnik mit einem Phasengitter als Probe, das Gitter bei
  einem Lieferanten hinter einem JSON-Kanal ausgewertet, gibt das Kamerabild des Gitters im
  Prozess der Engine (1e−12).

- Vorbereitete Läufe für fremde Blöcke (CAD ADR 0010, D6, Stufe W3): Blöcke, die mehr als ihre
  Eingänge brauchen – Raytatouille-Blöcke bauen ein System aus den Daten aller `rtt.*`-Blöcke und
  folgen den Ästen des ganzen Graphen –, bekommen vor ihrer ersten Auswertung den abgerollten Graphen
  des Laufs (`prepare`, ein Rückruf des Lieferanten ergänzt z. B. eine Materialbibliothek); die
  Auswertungen rechnen in diesem Kontext. `BlockProvider.handle` nimmt beide Anfragen. Abnahme W71.

- Lieferanten im Workflow-Dokument (CAD ADR 0010, D7): `model_document(…, providers, recorded)`
  vermerkt unter `providers`, welches Plugin (mit Version) einen im Modell benutzten Blocktyp
  liefert, und markiert die Bibliothekseinträge mit ihrem Plugin. Fehlt ein vermerktes Plugin,
  meldet die Prüfung `graph.provider_missing` mit Plugin und Version statt eines unbekannten Typs;
  der Block bleibt mit seinen Daten. Abnahme W70.

- `graph.remote` (CAD ADR 0010, D2/D3, Stufe W2b): Blöcke in einem anderen Prozess. Lieferant:
  `BlockProvider` beschreibt seine Blocktypen (`catalogue`, mit den nativen Darstellungen) und
  rechnet eine Anfrage mit Port-Daten (`evaluate`). Engine: `register_remote` übernimmt die Kataloge
  in die Registry, `RemoteEvaluator` rechnet über eine Anfrage. `NodeType.to_dict`/`from_dict` für
  den statischen Teil (ohne Port-Rahmen). Ein fremder Block bekommt seine Eingänge, Parameter und die
  JSON-fähigen Einstellungen des Laufs. Ein Fehler des Lieferanten oder des Kanals (ein Lieferant,
  der abbricht) ist der `ExecutionError` „exec.remote“ des Blocks; der Lauf geht anderswo weiter.
  Abnahme W69.
- Die Blockbibliothek der Workflow-Engine (`graph.workflow.registry`) enthält auch die allgemeinen
  Optik-Blöcke (`optics.*`: Quellen, Detektoren, Grenzflächen, Teiler, dünne Linse) und den
  Evaluator verdichteter Gruppen (R29), die der Editor schon einstellen kann.

- `rtt_fusion.portdata` (CAD ADR 0010, D4, Stufe W2a): die Werte auf den Links des Workflows als
  Dateien – ein JSON-Baum (Datenklassen, Listen, Zahlen; Gleitkommazahlen bitgenau als Hex, auch NaN
  und ±∞) und je Array eine `.npy`-Datei. Damit rechnet ein Block im Prozess eines anderen Plugins.
  Beim Lesen entstehen nur Klassen aus `rtt_fusion` und `fieldtracing`. Abnahme W68.
- L4, Stufe O3 für FourierGrating: `apply_many` nutzt die S-Matrix in der s/p-Basis der Ordnungen
  (FourierGrating `Simulation.smatrix`, FG #4): ein Eigenproblem je Richtung für beide
  Polarisationen statt zweier Lösungen, gleich `apply`. Mit `bloch=True` ein Eigenproblem je
  Bloch-Familie; dann weicht eine Welle außerhalb der Mitte um den Trunkierungsfehler ab (das
  Ordnungsfenster liegt um die Einfallsordnung, etwa 7e−3 bei N = 6). Ohne `smatrix` (älteres
  FourierGrating) wie bisher. Die CI baut FourierGrating auf dem Stand b3e2c31. Abnahme W67.
- `graph.workflow` (CAD ADR 0010, Stufe W1a): das Workflow-Dokument, das der Node-Editor zeigt –
  Modell mit Sockets und Ports, Blockbibliothek, Diagnosen von Modell und abgerolltem Graphen, Plan
  je Link und je Verbindung – samt Blockzuständen je Block des Modells. Bisher im Raytracer-Plugin
  des CAD; jetzt in der Bibliothek, damit Workflow-Engine und Plugins es gemeinsam nutzen. Die
  bearbeitbaren Felder liefert der Lieferant über einen Rückruf; die Engine übernimmt sie. Abnahme
  W66.

- Bilder warnen, wenn das Licht zu wenige Pixelmitten bedeckt (`Image.note`): Das Feld wird an den
  Pixelmitten abgetastet, ein Fleck kleiner als wenige Pixel fehlt dann ganz oder teilweise, ohne
  dass etwas falsch ist. Die Warnung nennt Pixelgröße, Fleckgröße und ein passendes Fenster. Sie
  prüft geometrisch, nicht über die Leistung, damit Dunkelheit durch Interferenz nicht warnt.
  Abnahme W65.

- Ein Raytatouille-System lässt sich ohne Import im Workflow anlegen (`new_system`): eine Quelle mit
  Standard-Einstellungen (mm/µm, d-Linie, Eintrittspupille 4 mm, ein Feldpunkt auf der Achse),
  bearbeitbar wie die eines importierten Systems, ein leerer Baum, eine Blende und ein Detektor aus
  der Bibliothek und der Ast „main“. Neues Bauteil `rtt.part.stop` (Raytatouille braucht eine
  Blende im Pfad). Abnahme W64.
- L4, Stufe O1 (Plan `docs/plan-l4-operatoren.md`): `apply_many(handle, incidents)` in allen fünf
  Adaptern. CliffordBem nutzt dabei Krylov-Recycling (`use_recycling`, je Wellenlänge), SpeckleBem
  den einmal assemblierten Operator; FourierGrating, hp-FEM und BWPM schleifen, bis die Löser
  Operatoren anbieten (O3). Abnahme W62.
- L4, Stufe O2: Die Aufrufer nutzen `apply_many`. Der Partikel-Block gibt alle Bündel mit demselben
  CliffordBem-Operator in einem Aufruf weiter (Recycling ab zwei Bündeln), die Jones-Tabellen alle
  Punkte einer Wellenlänge (bitgleich, bis FourierGrating einen Operator hat). Abnahme W62b.
- L5: CliffordBem und SpeckleBem nehmen Licht auch am hinteren Port (Richtung −1) an; Ausgänge
  `"back"` (reflektiert), `"front"` (durchgelassen) und `"front_scattered"`. Abnahme W63.

- `section_to_group` (ADR Nachtrag 1): ein Raytatouille-Abschnitt (mehrere Elemente einer Baugruppe
  in einem Block) wird eine Gruppe mit einem Block je Element. Jeder Block hält die Daten seines
  Elements und seinen Teil jedes Durchgangs; die Durchgänge des Abschnitts werden innere Äste, die
  Baugruppe selbst kommt in den Rest des Baums an der Quelle. Das erzeugte System bleibt dasselbe.
  Abgelehnt: ein Abschnitt aus einem Element, ein verschobener Abschnitt, ein Geist durch ihn.
  Abnahme W61.

- S4 (ADR Nachtrag 1): neue Raytatouille-Bauteile aus der Bibliothek (`graph/parts.py`:
  `rtt.part.mirror`, `rtt.part.splitter`, `rtt.part.lens`, `rtt.part.detector`). Das Element ist
  der Parameter `rtt_new` (bearbeitbar wie die Daten importierter Blöcke), die Ports sind fest, die
  Ereignisse eines Durchgangs folgen aus den Ports. Vor dem Abrollen bekommt jedes Bauteil seine
  Stelle im System, eindeutige Flächen-IDs und seine Durchgänge (`prepare_parts`) und rechnet dann
  wie ein importierter Block. Knotenweise Durchgänge werden über die Ports des Ast-Schritts
  gesucht. Ein im Workflow gelöschter Block hinterlässt im erzeugten System eine leere Baugruppe.
  Abnahme W59.

- S2 (ADR Nachtrag 1): Das Raytatouille-Backend rechnet mit dem System, das der Workflow ergibt
  (`importer.system_from_graph`: Baum aus den Daten der Blöcke, Pfade aus den Ästen; kompiliert je
  Zustand). Ohne gehaltene Daten (alte Graphen) weiter mit `shared["system"]`. Knotenweise Durchgänge
  kommen aus den Daten der Blöcke (`runs`), nicht mehr aus Systempfaden; ein eigener Ast aus
  bekannten Durchgängen wird ein Pfad. Ein ersetzter Block behält die Raytatouille-Daten des
  Elements (`replace_block`), das Gitter wirkt an dessen Stelle (R19b). Der Import hält auch die
  Durchgänge von „auto“-Pfaden (ausdrückliche Pfade haben Vorrang, Spiegel reflektieren). Ein Ast
  mit einem Schritt ohne Ereignisse (Geist) gibt einen Zielpfad „<Ast> (aim)“. Abnahme W58.
- S1 (ADR Nachtrag 1, der Workflow ist die Wahrheit): Der Import hält das ganze System im Modell –
  jeder Baumknoten (Element, Baugruppe) mit seiner Stelle im Block, dem alle seine Flächen gehören
  (`rtt`), die Ereignisse jedes Durchgangs im Block (`runs`), die Systemangaben, die übrigen
  Baumknoten und die Reihenfolge der Pfade in der Quelle (`rtt_system`, `rtt_tree`, `rtt_paths`);
  `importer.system_from_model` erzeugt das System daraus zurück. Kanäle umfassen auch
  Ereignisarten, die die Ports nicht unterscheiden (`importer.event_channel`, z. B.
  `P.S1:extraordinary`). Abnahme W57 mit drei Referenzsystemen von Raytatouille. W39 vergleicht
  ohne die neu gehaltenen Daten.
- `Executor.needed(wanted)` und `Executor.status(wanted=…)`: Marken aktuell/veraltet nur für die
  Blöcke, die die gewünschten Ausgänge (die Darstellungen) brauchen, je Block mit den benötigten
  Ausgängen – ein Teiler ist aktuell, auch wenn sein Weg zu einer nicht gezeigten Kamera fehlt
  (W56).
- R29 Verdichten (`graph/condense.py`): `Group.condense` macht eine Gruppe zu einem Block
  `fusion.condensed`, der je innerem Ast eine Tabelle ebener Wellen anwendet (s- und
  p-Koeffizient, Ausfallsrichtung, optischer Weg und Versatz zwischen Ein- und Austrittsebene je
  Einfallsrichtung σ, bilinear zwischen den Stützstellen). Die Tabellen entstehen aus Probewellen
  durch die Blöcke der Gruppe, prüfen seitliche Invarianz und Isotropie (`exec.condense`) und
  liegen im Operator-Speicher im Rahmen der Gruppe. Äste durch die Gruppe nennen weiter ihren
  inneren Ast; der Block wendet genau diese an. Format Version 2 mit `condense`. Abnahme W55.
- R28 Gruppen und Varianten: `Group` und `Variant` im Modell (exportierte Ports und Parameter,
  innere Äste, aktive Variante, `set_variant`), ein Schritt durch eine Gruppe nennt ihren inneren
  Ast als Kanal; `flatten` löst Gruppen auf (verschachtelt, innere Blöcke als `gruppe/block`,
  Prüfungen `graph.group_variant`, `graph.group_ports`, `graph.group_params`,
  `graph.group_branch`); der Gruppenknoten bleibt in der Socket-Form, seine Pose, sein Eingang
  `pose` und seine exportierten Parameter wirken zur Laufzeit auf die Blöcke darin (`_parent`,
  `_bind`, `pose.bind_exported`; der Schlüssel des Ergebnisspeichers enthält die Gruppen).
  Format Version 2 mit Varianten. `optics.thin_lens` (paraxiale dünne Linse für jede Einfallsrichtung).
  Abnahme W53–W54.
- R27 allgemeiner Teiler und Mehrfachumläufe: `optics.interface` (Grenzfläche zweier Medien,
  Fresnel s und p, leistungsnormierte Amplitude, Totalreflexion) und `optics.splitter`
  (verlustfreier Teiler mit vier Ports, Stokes-Beziehungen); `Plane.intersect` rechnet die OPL im
  Medium des Bündels. Umläufe trägt der Modellierer als Äste ein (Abnahme: Airy-Teilsummen).
  Geister in Raytatouille-Abschnitten: Streueintrag `ghost`, Kanal `ghost:<Fläche>,…`, das
  Backend setzt die Ereignisse aus der Flächenfolge zusammen; importierte Abschnitte und Teiler
  haben Ports mit beiden Kanälen, fehlende Sockets ergänzt `unroll`; „auto“-Pfade ohne Ast
  werden für die Runs ausgeschrieben. Testsystem `ghost_plates` (Raytatouille M4), Abnahme
  W50–W52.
- R26 Posen (`graph/pose.py`): Isometrien in Raytatouilles Konvention (Translation, Rotation um
  den Pivot, intrinsisch X → Y → Z); Pose je Knoten ohne Geometrie (Parameter `pose`), Versatz
  `pose_offset` und je Treiber `pose_by_<id>` (sie addieren), Eingang `pose` an jedem Lichtblock,
  Meta-Blöcke treiben eine Komponente (`component`: dx … rz). Die Ausführung übergibt Auswertern
  die geometrischen Parameter in Weltkoordinaten (`posed`); Raytatouille-Blöcke bewegen ihren
  Knoten im Modellbaum starr (`exec.pose` für Teilabschnitte). Port-Lagen (`PortFrame`,
  `NodeType.port_frames` für Quelle, Spiegel, Linse, Detektor, Gitter) und Prüfung der
  Verbindungen (`connection_diagnostics`: `graph.port_facing`, `graph.port_behind`,
  `graph.port_medium`, `graph.port_aperture`, `graph.port_near_field`, `graph.port_tilt`).
  Operator-Speicher `Executor.operators` / `Context.operators`: der CliffordBem-Operator des
  Partikels ohne Lage im Schlüssel. Abnahme W45–W49.
- R25 knotenweise Ausführung im Raytatouille-Backend: `rtt.section`, `rtt.mirror`, `rtt.splitter`
  geben Strahlen je Ast und Besuch weiter (`RayBundle.branch`), die Quelle startet je Ast, der
  Detektor bildet ab bzw. nimmt Felder auf einem H-Band; Runs aus allen Systempfaden; ausgetauschte
  Blöcke (R19b) dazwischen, an der Stelle der ersetzten Fläche; Modus `shared["rtt_execution"]`
  ("nodes", "paths", "auto": knotenweise, sobald ein Ast kein Systempfad ist). Der Import merkt
  sich die Durchgänge der Systempfade (`meta["raytatouille_runs"]`), `graph.branch_backend` prüft
  Äste dagegen statt gegen die Pfadnamen. Abnahme W43–W44: knotenweise = ganze Pfade (≤ 3e−10).
- R24 Ports, Äste im Modell, Format 2 (`graph/topology.py`): Licht-Ports mit zwei Kanälen
  (`PortSpec`, `NodeType.ports`, abgeleitet aus Socket-Namen), Streutabelle (`NodeType.scatter`,
  sonst aus `through`), ungerichtete Verbindungen (`Connection`, Transport je Richtung), Äste nur
  vom Modellierer (`Branch`, `Step`) mit Prüfung `validate_model` (`graph.no_branch`,
  `graph.branch_impossible`, `graph.branch_duplicate`, `graph.port_unknown`, `graph.port_links`,
  `graph.branch_backend`, `graph.observer_ambiguous`); `unroll` erzeugt die Socket-Form als
  minimalen Automaten der Äste (Kopien `id@k` je Besuch, Sammelblock `fusion.sum`), Planer,
  Ausführung und Speicher unverändert; Format Version 2 (`io.model_to_dict`,
  `io.model_from_dict`, Version 1 über `model_from_graph`); `model_from_raytatouille`, der Import
  läuft über das Modell (Pfade als Äste, Beugungsordnungen als Kanäle). Abnahme W39–W42; Ausgabe
  des alten Imports als Referenz in `tests/data/import_v1`.
- ADR-Entwurf „Bidirektionale Ports, Posen, Äste und Gruppen“ (`docs/adr-entwurf-ports-posen-aeste.md`)
  mit der Festlegung P6 des Maintainers (Äste nur vom Modellierer, keine Schleifen im Framework)
  und den Stufen R24–R29; W6 im Konzept als ersetzt markiert.

- R23: Ergebnisse über Läufe – `ResultCache` (Grenze, LRU, Warnung), Schlüssel ohne Rechnen mit
  Umgebung, Treffer sparen alles davor, Meta-Blöcke gespeichert, `Executor.status()`. Tests W37,
  W38.

- R22: Phasenschieber am Spiegel – Meta-Blöcke mit `drive` (aus) und `observe` (ein); Eingang
  `drive` (Datentyp `drive`) an Spiegeln, Abschnitten und Teilern (`NodeType.extra_sockets` für
  importierte Blöcke) sowie am ebenen Spiegel des analytischen Backends: Verschiebung entlang der
  Flächennormalen. Tests W35, W36.

- R21: Darstellungs-Blöcke (`rtt_fusion.graph.display`: `display.image`, `display.plot`,
  `display.beam` als Sonde an jedem Lichtausgang) mit dem Ergebnis `View`; Fortschritt der
  Ausführung (`Context.progress`, `Context.span`, `Context.report`; Raytatouille-Detektor je
  Pfad/Feld/λ, Feldausbreitung, Meta-Blöcke je Wert) und laufender Block (`Context.on_node`).
  Tests W32–W34.

### Geändert

- FourierGrating eingestellt auf 0.18.0 (Commit `174ef93`, statt 0.16.0 `959c93b`): baut auch mit
  Apple clang (FourierGrating PR #2), so dass FourierGrating, Raytatouille und CliffordBem in
  einem Prozess dieselbe C++-Laufzeit teilen; Konventionsdateien unverändert.

### Hinzugefügt

- R19c: Partikel als Probe – Block `fusion.particle` (`rtt_fusion.graph.particle`): fokussierte
  Strahlen als Winkelspektrum (`rays_to_spectrum`), Streuung mit CliffordBem (Adapter R4),
  Rückstreuung als Strahlen (`spectrum_to_rays`); kollimiertes Licht als eine ebene Welle.
  Tests W29–W31. CliffordBem lokal mit Apple clang gebaut (eine Zeile in den Python-Bindings,
  Vorschlag für CliffordBem).

- R20b: `offset` an Raytatouille-Blöcken verschiebt deren Elemente (Kopie des Systems,
  übersetzt je Verschiebung zwischengespeichert); Sweeps und Optimierung wirken damit auf
  importierte Systeme. Test W28.

- R19d: Der Gitter-Block baut seine Jones-Tabelle mit FourierGrating selbst (`fouriergrating`:
  Konfigurationsdatei, `side`: reflektiert/transmittiert), um die einfallenden Strahlen herum
  und zwischengespeichert. Tests W26, W27 (laufen mit FourierGrating, 1D lokal in Sekunden).

- R20: Meta-Blöcke (`rtt_fusion.graph.meta`) – Sweep, Replikate (Mittel der Bilder) und
  Optimierung eines Parameters um den azyklischen Graphen; Parameter als „block.param[i]“,
  gesammelter Ausgang „block.socket“, Messgrößen Spitze/Mitte/Leistung; die Ausführung bietet
  verschachtelte Auswertung (`Context.evaluate`) und ungecachte Evaluatoren. Tests W24, W25.

- R19b: Ausgetauschte Blöcke im Raytatouille-Backend – der Pfad eines importierten Systems
  wird am fremden Block geteilt (Ereignisse den Blöcken zugeordnet, Teilpfade A und B in einer
  Kopie des Systems), der Block wirkt auf die Strahlen, Raytatouille rechnet weiter; Felder
  (skalar entlang einer reellen Polarisationsrichtung, Rest in der PRT), OPL und PRT
  verkettet, zwei Polarisationen; für Strahlen- und Feldbilder. `replace_block` (Austausch mit
  Zuordnung der Sockets nach Richtung und Art). Evaluatoren und Registry im Kontext.
  Testsystem `tests/data/linnik_mini.rtt.json`. Tests W22, W23.

- R19a: Gitter als Probe – Block `fusion.grating` (`rtt_fusion.graph.grating`) mit den Sockets
  eines Spiegels: Beugungsordnungen als Kanäle aus einer Jones-Tabelle (R6, FourierGrating;
  `load_table` liest `raytatouille-result`-Dateien), Gittergleichung, PRT der Ordnung, Phase der
  Phasenschicht m·2πx/Λ, Ordnungen in derselben Kohärenzgruppe; `light.rotation_between`.
  Tests W19–W21.

- R18b: Felder im Raytatouille-Backend – ein Band in einen Detektor eines importierten Systems
  kann H tragen: Felder aus den aufgezeichneten Strahlen nach der letzten Fläche vor dem
  Detektor (PRT und Gewicht von Raytatouille, Brechzahlen aus dem OPL, zwei Polarisationen der
  unpolarisierten Quelle), Ausbreitung mit Fieldtracings Planer. Die Quellachsen werden
  geometrisch mitgeführt (`RayBundle.source_frame`), damit Drehung und Spiegelung des
  Strahlenbilds durch die Geometrie den Maslov-Index nicht verfälschen; Brechzahl des Mediums
  (`index`, `source_index`); Gauß-Beleuchtung (`gaussian_waist`), Bildfenster `half_size`,
  Plan im Kontext der Ausführung. Tests W17, W18.

- R18: Field Tracing auf Bändern – Darstellung H (`rtt_fusion.graph.field`: `FieldBundle`,
  Strahlen → Feld nach Fieldtracing ADR 0004, Feld → Strahlen, Bild aus Feldern mit Ausbreitung
  durch Fieldtracings Planer), Polarisation der Strahlbündel (Jones der Quelle, geometrische PRT,
  Ausgleich der Spiegel-Orientierung im Maslov-Index), `optics.aplanatic_lens`, Gauß-Quelle
  (`gaussian_waist`), Ausführung mit wählbaren Darstellungen (`Executor(available=…)`,
  `WITH_FIELDS`), Umwandlung an den Ports. fieldtracing bleibt optional (Import erst bei H).
  Tests W13–W16. Der Fall in W8 mit festem Transport H (abgelehnt) nutzt jetzt G, weil die
  Blöcke des analytischen Backends H annehmen.

- R17: Planer, Stufen 1–2 (`rtt_fusion.graph.planner`) – Fusionspfade auf Socket-Ebene,
  Darstellung je Licht-Band und Umwandlungen in Blöcken nach Tabelle 6.3 (auch Ketten; nach
  Phasenverlust keine kohärente Darstellung mehr), feste Transporte, verfügbare Darstellungen,
  Gültigkeits-Hook, „auto“ ohne P solange kohärent möglich, Warnung `plan.incoherent_merge`;
  Plan als Datenobjekt mit JSON-Form für den Node-Editor; die Ausführung plant vor dem Rechnen
  (`Executor.plan`). Tests W10–W12.

- R15: Ausführung des Arbeitsgraphen ohne Planer – `graph.light` (Strahlenbündel R mit
  Kohärenz-Kennung), `graph.execute` (Auswertung je Ausgang, Cache nach Inhalts-Hash,
  Stummschalten, Fehlercodes `exec.*`), `graph.imaging` (kohärentes Kamerabild: Amplituden je
  Kohärenzgruppe, Intensitäten über Gruppen), `graph.optics` (analytisches Backend mit ebenen
  Elementen), `graph.backend_rtt` (importierte Raytatouille-Graphen, auch nicht zielbare Pfade);
  Pilot W9 (Tests W4–W9).

- R14: Graph-Modell des Node-Workflows `rtt_fusion.graph` – Blöcke, Sockets (Licht: R/H/G/S/T/P,
  Daten: N, Q und benannte Typen), Bänder mit Transport, Bearbeitung (Verbinden, Einfügen auf ein
  Band, Austauschen mit erhaltenen Bändern, Stummschalten), Prüfung mit Diagnosecodes und
  Kreissuche auf Socket-Ebene, Format `rtt-fusion-graph`, Block-Typen aus `describe()`, Import aus
  den Pfaden eines Raytatouille-Systems; Tests W1–W3.

- L11: Fusionsplan aus Fieldtracing (Zweig `docs-fusion-plan`, Commit `3e0d5f8`) als
  `docs/fusionsplan.md` übernommen; Verweise angepasst.

- R13: Konzept Node-Workflow (freigegeben, E3) `docs/konzept-node-workflow.md` – Arbeitsgraph aus
  Blöcken, Sockets und Bändern als bearbeitbare Form des Planers, Kohärenz-Kennungen,
  Meta-Blöcke für Schleifen, Ort von Logik (RTT-Fusion / `rtt-fusion`) und GUI (Raytatouille
  Fusion), Vorschläge A1–A6 zur Adapter-Schnittstelle, Inventar aller Module, Stufen R14–R20.
  Vorschläge L1–L12 zu den offenen Punkten: ein G-Typ mit optionaler Trägerphase, Umrechnung der
  FFT-Normierung, Port-Typen R/H/P/M, `apply_many`/`operator`, Einfall am hinteren Port, Namen der
  Daten-Sockets, Streustatistik, maschinenlesbare Gültigkeit, Beistelldateien für große Daten,
  Kopplung GUI ↔ Bibliothek mit Formatskizze `rtt-fusion-graph`.

- R0: Gerüst (Regeln, Plan mit D1, D2, D10, D11), Konventionen, Schnittstelle der Adapter,
  Port-Typen `GriddedField` und `PlaneWaves` mit Umwandlung G ↔ S und Poynting-Leistung, CI mit
  den Löser-Modulen auf festen Commits.
- R1: FourierGrating-Adapter `rtt_fusion.adapters.fouriergrating` (Drehung der Koordinaten,
  kartesische Rekonstruktion mit FourierGratings Basen, Einfall als Überlagerung von δ = 0°/90°,
  Ports an oberer und unterer Grenzfläche) mit Konformitätstests F1–F5, darunter die Querprobe
  gegen den BWPM-Adapter.
- R2: BWPM-Adapter `rtt_fusion.adapters.bwpm` (Gitterfeld-Ports), bitgleich zu `bwpm.adapter`;
  Energie über `ports.power`, Querprobe mit FourierGrating über die Umwandlung G → S.
- R3: Pilot K1 – `rtt_fusion.couplings.grating_table` (Effizienzen je Ordnung aus FourierGrating,
  als `diffraction_efficiency` in ein Raytatouille-System) mit Test gegen Gittergleichung,
  Pfadtransmission und Energie; Raytatouille als Modul (CI baut es mit Eigen, nlohmann-json, TBB).
  Befunde für den Fusion-ADR in `docs/testfaelle.md`.
- R4: CliffordBem-Adapter `rtt_fusion.adapters.cliffordbem` (Einfall als Winkelspektrum, Fernfeld →
  Port-Spektren über Weyl), `couplings.focus.richards_wolf`, `ports.spectrum_power` und
  `PlaneWaves.solid_angle`; Tests C1–C3 und Pilot K3 (Teilchen im Fokus, Rückstreuung ins
  Objektiv). Befunde für den Fusion-ADR in `docs/testfaelle.md`.
- R5: ADR-Entwurf „Raytatouille Fusion“ für Raytatouille (`docs/adr-entwurf-fusion.md`): Scope,
  Gebiete und Ports, Darstellungen R/H/G/S/T, tabellierte Jones-Interaktion als Erweiterung von
  ADR 0025, Adapter mit Konformitätssatz, Modell und Pfade, Planung, Kohärenz, Cache, GPU;
  begründet mit den Befunden B1–B7 aus R1–R4 und BWPM.
- Fusion-ADR freigegeben (E1): Punkt 9 „`rtt-fusion` in C++“ mit drei Wegen für Adapter;
  Vorschläge zu Nummer, Zeitpunkt und eigenem ADR für die tabellierte Interaktion; Stufen R6, R7.
- R6: Prototyp der tabellierten Beugungsinteraktion `rtt_fusion.couplings.jones_table`
  (leistungsnormierte Jones-Matrix je Ordnung in kartesischer Basis, PRT-Matrix, kubische
  Interpolation, Erkennung von Rayleigh-Anomalien) mit Tests J1–J5 und Konvergenzmessung;
  ADR-Entwurf `docs/adr-entwurf-beugungstabelle.md`.
- R7: Sollwertdateien für den C++-Port (`tools/export_reference.py`, `reference/`, sechs Fälle mit
  unabhängigen Werten); Umwandlung Fernfeld → Spektrum und Halbkugel-Quadratur als Modul
  `rtt_fusion.convert` (der CliffordBem-Adapter nutzt sie).
- R8: SpeckleBem-Adapter `rtt_fusion.adapters.specklebem` für geschlossene Körper (Konjugation
  von ε, Einfall und Fernfeld wegen exp(+jωt), SI-Einheiten) mit Tests X1–X4: Mie für eine
  absorbierende Kugel, Vorzeichen der Absorption, Querprobe gegen CliffordBem, Überlagerung.
  CI baut SpeckleBem mit `pybind11<3` (`constraints-modules.txt`).
- CI: jedes Löser-Modul in einem eigenen Schritt; SpeckleBem ohne Isolation gegen pybind11 < 3;
  Issue-Entwurf für SpeckleBem `docs/vorschlaege/specklebem-pybind11.md`.
- R9: hp-FEM-Adapter für Gitter `rtt_fusion.adapters.hpfem` (`grating.jones` → kartesische
  Felder in den Port-Konventionen) mit Tests H1–H3: Fresnel (3e−15), Querprobe gegen
  FourierGrating mit allen Ordnungen (starke < 0,6 %, schwache < 3 %), Energie.
- R10: BioOptoMD-Datenadapter `rtt_fusion.adapters.biooptomd` (n(z)-CSV, Schichten nach
  BioOptoMDs Regel, optische Überschussdicke, Dipole aus dem Tracker, Beleuchtungsseite als
  Spiegelung) und Pilot K4 (SPR-Biosensor): de Feijter formunabhängig auf 0,07 %, Vergröberung
  auf 1 nm < 5e−4; Nachtrag zum Fusion-ADR (Datenports, weitere Querproben).
- CI: Löser-Module als Wheels in einem Cache (actions/cache); ohne Änderung der Pins entfällt der
  Bau der Module.
- E2: ADR-Entwurf tabellierte Beugungsinteraktion freigegeben.
- R11: C++-Prototyp von `rtt-fusion` (`cpp/`): Port-Typen, Umwandlungen, Ergebnisleser,
  `SolverAdapter`/`AdapterRegistry`; BWPM-Adapter in C++ gegen BWPMs Kern (FetchContent); Tests
  gegen die Sollwertdateien (C++1, C++2); CI-Job `cpp`.
- CI: Cache-Schlüssel ohne den Workflow (Änderungen an anderen Jobs bauen die Module nicht mehr
  neu); Version `v1` im Schlüssel für Änderungen der Bauschritte.
- R12: tabellierte Beugungsinteraktion in C++ (`DiffractionTable`: kubische Interpolation,
  Rayleigh-Anomalien, PRT-Matrix, Prüfung) mit Test C++3 gegen die Sollwertdatei; Leser für
  boolesche Werte.
- Phasenraum-Übergang `rtt_fusion.phasespace` (Feld → Husimi-Pakete, Kopplung in eine Gauß-Mode)
  für Scattuccino S9b; Tests PS1, PS2.
- `phasespace.field_to_packets(periodic=True)`: wrapping windows for periodic fields (BWPM ohne
  Absorberrand), Fenstergröße geprüft.

### Behoben

- Die Prüfung warnte bei jedem Bauteil aus der Bibliothek vor den Parametern, die `prepare_parts`
  ihm gibt (`rtt`, `runs`, `surfaces`, `location`, `element_type`; `graph.param_unknown`). Diese
  abgeleiteten Parameter gelten jetzt als bekannt; ein wirklich unbekannter warnt weiter (W64).
- Raytatouille-Backend: Ein Strahl galt als angekommen, wenn er alle Plätze des aufgezeichneten
  Pfads füllte. Raytatouille gibt aber auch einem Strahl, der am letzten Ereignis verloren geht,
  alle Plätze (`count = j + 2`). Knotenweise liefen so z. B. an der letzten Fläche eines Laufs
  abgeschattete Strahlen weiter (Linnik mit kleiner Objektivapertur: 289 statt 79). Jetzt zählt
  zusätzlich `lost_at == -1`, an allen vier Stellen. Die Sonderregel für Strahlen, die auf oder
  hinter der ersten Fläche eines Laufs starten, entfällt: Raytatouille meldet sie als verfehlt;
  sie war kein Unterschied zwischen ganzem Pfad und Teilpfad, sondern dieser Lesefehler (W43b).
- `system_from_model` las nur die oberste Ebene: Raytatouille-Blöcke in Gruppen fehlten im erzeugten
  System (das Zeichnen im CAD verlor sie). Gruppen zählen jetzt mit ihrer aktiven Variante, aufgelöst
  wie zum Rechnen (W61).
- `graph.branch_backend` prüfte die Durchgänge eines Raytatouille-Blocks gegen die Liste des Imports
  nach Block-ID; in einer Gruppe (`<Gruppe>/<Block>`) fand sie nichts. Jetzt zählen die Durchgänge,
  die der Block selbst hält (S2), bei älteren Modellen die Liste ohne den Pfad der Gruppen (W60).
- Geister durch gekrümmte Flächen (R27): durchlaufene Flächen wurden als `transmit` (ungebrochen)
  gerechnet, jetzt mit der Art des Vorwärtsdurchgangs (gebrochen); Test W52b mit Raytatouilles
  `ghost_singlet`.

