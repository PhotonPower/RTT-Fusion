# ADR-Entwurf: Bidirektionale Ports, Posen, Äste und Gruppen

> **Entwurf für RTT-Fusion**, Nachtrag zum Konzept Node-Workflow
> ([konzept-node-workflow.md](konzept-node-workflow.md), E3) und Präzisierung von Fusionsplan 6.1
> („Richtung ein, aus, beides“) und 7.2 Stufe 1. Entstanden am 2026-11-12 aus der Durchsicht von
> RTT-Fusion (Stand R22, `310fe71`) und Raytatouille-Fusion (`1b54aae`); ergänzt am selben Tag um die
> Festlegungen des Maintainers zu Ästen (Abschnitt 6) und die Durchsicht gegen R23 (Abschnitte 4,
> 5, 8). Ändert die Abschnitte 3.2, 3.3, 4 und 7 des Konzepts sowie W6; alles Übrige
> (Darstellungen, Kohärenz, Planer Stufe 2, Adapter) bleibt.

**Status:** vorgeschlagen; P6 vom Maintainer festgelegt (2026-11-12), P1–P5 und P7–P9 offen.

## Kontext

Der Arbeitsgraph hat heute **gerichtete Sockets** (`front_in`, `back_out`, …) und je Blocktyp eine
Tabelle `through` (Ausgang ← Eingänge). Die Verdrahtung ist damit zugleich die Auswahl der
Lichtwege. Das trägt den Linnik, hat aber Folgen:

- Ein Doppeldurchgang braucht vier Bänder statt zwei; der Nutzer zeichnet Wege, nicht den Aufbau.
- Der Strahlteiler (`in, r, t, r_in, t_in, merged`) ist auf den Linnik zugeschnitten. Ein
  allgemeiner Teiler hat vier Ports mit je zwei Richtungen; Licht zurück zur Quelle oder von der
  Kameraseite ist nicht darstellbar.
- Ein zusätzlicher Reflex (Ghost) verlangt neue Sockets am Blocktyp.
- `Category.GROUP` existiert nur als Enum-Wert; Modell, Format, Prüfung, Planer, Ausführung und
  Editor sind flach. Ein Objektiv ist ein Block oder viele, nicht wahlweise beides.
- Blöcke haben keine eigene Lage. Ebenen stehen als absolute Parameter im Block bzw. im
  `.rtt.json`; `drive` verschiebt nur entlang einer Normalen. In der CAD werden Knoten aber
  verschoben und gedreht.
- Im Raytatouille-Backend geben Blöcke kein Licht weiter (`_structural`), der Detektor rechnet
  ganze Pfade. Der Graph ist dort Anzeige, kein Datenfluss.

Ziel des Maintainers (2026-11-12): Knoten mit Ports, die zugleich Ein- und Ausgang sind; zwischen
Knoten wird propagiert; der Modellierer gibt über Äste vor, welche Kopplungen zählen; mehrere
Knoten lassen sich zu einer Black-Box zusammenfassen; Knoten sind in der CAD starr beweglich.

## Entscheidung

### 1. Drei Ebenen

| Ebene | Inhalt | Wer bearbeitet |
| --- | --- | --- |
| **Topologie** | Knoten mit Pose, Licht-Ports, ungerichtete Verbindungen, Gruppen: was physikalisch da ist | Nutzer (Node-Editor, 3D-Ansicht), Import |
| **Äste** | gerichtete Wege über die Topologie: was gerechnet wird | **nur** der Modellierer (Editor, Import, Skript) |
| **Ausführungsplan** | abgerollte Äste als DAG gemeinsamer Präfixe, Darstellung je Abschnitt, Operatoren | Planer (Stufe 2 wie bisher) |

Die bisherige Ebene „Arbeitsgraph“ teilt sich in Topologie und Äste. Der heutige Socket-Graph
bleibt als **abgeleitete** Zwischenform bestehen (Abschnitt 8), Planer Stufe 2 und Ausführung
arbeiten weiter auf ihm.

### 2. Licht-Port

Ein Licht-Port ist eine Bezugsfläche am Knoten mit

- einem **Rahmen** im lokalen System des Knotens: Ursprung, Außennormale n, Achsen u, v in der
  Ebene (volle Lage, sechs Freiheitsgrade; die Drehung um n legt Gitterachsen und
  Polarisationsbasis fest);
- Ausdehnung (Apertur) und Medium ñ außen;
- den zulässigen Darstellungen (R, H, G, S, T, P) wie bisher;
- **zwei Kanälen**: *ein* (Welle läuft gegen n in den Knoten) und *aus* (Welle läuft entlang n aus
  dem Knoten). Das ersetzt „Richtung ±1“ und die Paare `*_in` / `*_out`.

Der Knoten ist ein **Streuoperator**: aus = S(ein) über alle Ports. Für kohärente lineare Knoten
ist S eine Blockmatrix, für Strahlen und radiometrische Pakete eine Abbildung; die Schnittstelle
ist dieselbe. Quellen und Detektoren sind Ein-Port-Knoten (Rückreflex in die Quelle, Reflex am
Sensor sind dann Einträge von S, im Normalfall null).

**Daten-Ports** (N, Q, `drive`, `pose`, `result`, …) bleiben gerichtet.

### 3. Streutabelle des Knotentyps

`describe()["scatter"]` nennt je Paar (Port ein → Port aus), was der Knoten liefern kann:

| Wert | Bedeutung | Beispiel |
| --- | --- | --- |
| `none` | physikalisch null oder nicht modelliert | Spiegel: hinten |
| `oneway` | geliefert, ohne innere Mehrfachreflexe | Tracer-Abschnitt: vorn → hinten, hinten → vorn |
| `full` | geliefert einschließlich aller inneren Reflexe | FourierGrating, BWPM, hp-FEM, CliffordBem |
| `ghost` | lieferbar (Flächenreflexe eines Abschnitts) | Tracer-Abschnitt: vorn → vorn |

Je Eintrag zusätzlich die **Kanäle** (Ordnung, r/t) wie bisher. Die Tabelle ersetzt `through`;
sie sagt, **was ein Ast an diesem Knoten nutzen darf**, nicht, was gerechnet wird (Abschnitt 6).
Ein Ast über einen Eintrag `none` ist ein Fehler (`graph.branch_impossible`). Adapter melden mit
`accepts_incidence` je Port, wo Einfall möglich ist (A3/L5 bleibt Voraussetzung für `full` von
beiden Seiten). Für reziproke Knoten wird S = Sᵀ in symmetrischer Basis zum Konformitätstest.

### 4. Pose des Knotens

- Jeder Knoten hat eine **Pose** gegenüber seinem Elternknoten (Gruppe oder Welt): starre
  Bewegung in den Konventionen von Raytatouille (`Pose`: Position in mm, Drehungen intrinsisch
  X → Y → Z in rad, Pivot; `docs/architecture.md`, „Transformationen“).
- Port-Rahmen in der Welt = Pose aller Vorfahren ∘ Pose des Knotens ∘ lokaler Port-Rahmen. Die
  Port-Rahmen gehören zum Knotentyp bzw. zur Geometrie und bewegen sich starr mit.
- **Eine Wahrheit:** Hat der Knoten Geometrie, ist seine Pose die des CAD-Objekts bzw. des Knotens
  im Raytatouille-Modellbaum; der Graph verweist darauf (stabile Kennung, z. B. `rtt.location`)
  und speichert sie nicht doppelt. Nur Knoten ohne Geometrie (ideale Linse, ebene Welle, Sonde)
  tragen die Pose selbst.
- **Intern** rechnet die Bibliothek mit Isometrien (R, t); das Sechser-Tupel ist nur die Form an
  Schnittstelle und GUI, weil Euler-Winkel von Reihenfolge und Pivot abhängen.
- **Daten-Eingang `pose`** an jedem Knoten mit Pose. Die **Nulllage gibt die Geometrie**; der
  Eingang trägt nur den zusätzlichen Versatz: sechs Werte (dx, dy, dz in mm; rx, ry, rz in rad).
  - Bezug: lokales System des Knotens in der Nulllage, Drehung um dessen Pivot, Reihenfolge wie
    bei `Pose` (erst verschieben, dann X → Y → Z). Lage = Nulllage ∘ Versatz.
  - Unbeschaltet gilt null. Mehrere Treiber an einem Eingang (z. B. Sweep und Toleranz) sind
    erlaubt und werden addiert; das ist nur für kleine Drehungen von der Reihenfolge unabhängig,
    sonst Diagnose `graph.pose_order`.
  - Der Versatz ändert das CAD-Modell nicht. Die 3D-Ansicht zeigt wahlweise Nulllage oder
    verschobene Lage; „Versatz übernehmen“ schreibt ihn als einen Undo-Schritt in die Geometrie.
  - An einer Gruppe bewegt er alle Kinder starr; ein Knoten, der mehrfach durchlaufen wird,
    sieht in jedem Ast denselben Versatz (wie Raytatouilles Toleranzen über `SurfaceId`).
  - Treiber: Sweep (Michelson-Spiegel), Replikate (Toleranzanalyse kippt Linsen), Optimierung
    (Kompensatoren). `drive` (R22) bleibt als Kurzform für dz entlang der Normalen.
- **Spiegelung** ist keine Pose. Ein gespiegelt eingebauter Knoten ist ein Parameter des Knotens
  (wie „Seite“ bei den Daten-Sockets).
- **Zwei Speicher (Ergänzung, R23):** Der Ergebnisspeicher (R23) hält Licht an Ports in
  Weltkoordinaten; deren Schlüssel enthalten die Pose des Knotens und aller Vorgänger, sonst
  lieferte ein verschobener Knoten alte Werte. Davon getrennt hält ein **Operator-Speicher** die
  lageunabhängigen Teile eines Löser-Knotens im lokalen System (`prepare`, Faktorisierung,
  Streumatrix) ohne Pose. Ein starr bewegter Löser-Knoten behält Vorbereitung und Operator; neu
  gerechnet werden seine Ausgaben, die Verbindungen und alles danach (PA9).

### 5. Verbindung

Eine Verbindung koppelt **zwei Ports** (ungerichtet) und trägt beide Laufrichtungen. Im homogenen
Medium koppeln die Richtungen nicht: Die Verbindung ist zwei unabhängige Einweg-Ausbreitungen.
Verfahren mit gekoppelter Vor- und Rückwärtsrichtung sind immer Knoten, nie Verbindungen.

- Die Ausbreitung folgt aus der **relativen Lage** der beiden Port-Rahmen: Abstand entlang der
  Normalen, Querversatz, Kippung, Drehung um die Normale. Decken sich die Ebenen, entfällt sie.
- Transport je Verbindung `auto` oder fest, auf Wunsch je Richtung verschieden.
- Prüfung mit Diagnosecodes: `graph.port_facing` (Normalen zeigen nicht zueinander),
  `graph.port_behind` (Zielport liegt hinter dem Startport), `graph.port_medium` (Medien
  verschieden), `graph.port_aperture` (Warnung: Aperturen überdecken sich nicht),
  `graph.port_near_field` (Warnung: Abstand unter einer Schwelle in Wellenlängen bei Ports, die nur
  propagierende Wellen tragen), **`graph.port_tilt`** (Warnung bei Transport H: Felder zwischen
  gekippten Ebenen werden heute geometrisch umgetastet, R18; gültig bei großen lokalen
  Fresnel-Zahlen, bis Fieldtracing gekippte Ebenen direkt ausbreitet).
- Ein Port hat höchstens eine **Verbindung**. Aufteilen und Vereinigen geschieht in Knoten.
- **Beobachter (Ergänzung):** Darstellungs-Sonden (`display.beam`) und die Beobachtung eines
  Meta-Blocks (`observe`) sind keine Verbindungen. Sie heften sich an einen Port eines Asts (Ast,
  Schritt, Kanal), ändern weder Topologie noch Äste und lesen das Licht, das dort im Plan ohnehin
  entsteht.

### 6. Äste (festgelegt vom Maintainer, 2026-11-12)

**Der Modellierer entscheidet alles.** Ein **Ast** ist eine endliche Liste von Schritten (Knoten,
Port ein, Port aus, Kanal) von einer Quelle zu einem Detektor. Gerechnet werden **genau die Äste,
die im Modell stehen.**

- Das Framework **zählt keine Wege auf, ergänzt keine und lässt keine weg.** Es gibt keine
  Grenzen wie Tiefe, Zahl der Umkehrungen oder Leistungsschwelle und keine Voreinstellung.
- Ein Graph ohne Ast zu einem Detektor ist ein **Fehler** (`graph.no_branch`), kein Anlass für
  eine Voreinstellung. Ein Ast über einen Streueintrag `none`, über einen Port ohne Verbindung
  oder mit einem nicht vorhandenen Kanal ist ebenfalls ein Fehler (`graph.branch_impossible`).
- Ein Ast darf einen Knoten mehrfach besuchen (Objektiv hin und zurück); ein Kreis in der
  Topologie ist deshalb kein Fehler. Höhere Umläufe zwischen zwei Knoten (Doppelreflex, Etalon)
  trägt der Nutzer **als eigene Äste** ein.
- **Keine Schleifen im Framework.** Äste sind endliche Listen, der Plan ist damit immer ein DAG mit
  vorab bekannten Kosten. Das Framework **summiert nie eine Reihe auf und löst keine Rückkopplung.**
  Ersetzt W6 „Kern azyklisch, Schleifen über Meta-Blöcke“ durch: *Äste sind endlich und vom
  Modellierer gegeben.* Der Meta-Block „Iteration“ (Konzept, Abschnitt 7) entfällt; Sweep,
  Replikate, Optimierung und die Zeitschleife (BindOptiX) bleiben, sie wiederholen den Graphen,
  koppeln aber kein Licht zurück.
- **Woher Äste kommen:** aus dem Editor (Liste, Hervorheben im Editor und in 3D, Schritte per
  Klick entlang der Verbindungen), aus Skripten und aus dem **Import**: Der Import aus
  Raytatouille legt dessen Pfade als benannte Äste an – das ist Inhalt des Modells, keine
  Voreinstellung des Frameworks.

Kohärenz bleibt wie in Konzept 5: Äste derselben Kohärenz-Kennung addieren Amplituden am
gemeinsamen Port.

### 7. Gruppen und Varianten

Eine **Gruppe** ist ein Knoten, dessen Inhalt ein Untergraph ist (Kompositum):

- eigene Pose; Kinder liegen relativ dazu (Baum wie System → Baugruppe → Element in
  Raytatouille). Hierarchie (Enthaltensein) und Lichtweg (Verbindungen auf einer Ebene) sind zwei
  verschiedene Strukturen;
- **exportierte Ports**: benannte Ports der Gruppe, je einer auf einen inneren Port abgebildet.
  Ihr Rahmen ist im System der Gruppe festgelegt und in allen Varianten gleich, damit ein
  Wechsel der Variante die äußeren Verbindungen nicht ändert;
- **exportierte Parameter**; Adressen sind Pfade (`objA/L3.radius`), auch für Sweep, Replikate
  und Optimierung;
- **innere Äste**: Ein äußerer Ast, der die Gruppe über (Port ein, Port aus) durchläuft, verweist
  auf einen benannten inneren Ast dieser Gruppe; auch die inneren Äste stehen im Modell. Außen
  erscheint die Streutabelle der Gruppe, abgeleitet aus ihren inneren Ästen.

Ein Knoten kann mehrere **Varianten** mit gleichen Ports haben, eine ist aktiv: etwa *ideale
Linse*, *Tracer-Abschnitt als Ganzes*, *Linse für Linse mit wählbarer Ausbreitung*, *rigoros*.
`replace_block` wird zum Wechsel der Variante.

Eine Gruppe lässt sich **verdichten**: Ihr Operator für genau ihre inneren Äste wird einmal je
(λ, Basis) berechnet und unter dem Inhalts-Hash im Operator-Speicher abgelegt (Tabelle T oder
Streumatrix, L4). Das ist die rechnerische Black-Box; das Zuklappen im Editor ist nur die
Anzeige. Verdichten summiert nichts auf, was nicht als Ast im Modell steht.

**Grenze:** Ein exportierter Port muss dort liegen, wo evaneszente Felder abgeklungen sind, oder
die Darstellung trägt sie mit. Nahfeldgekoppelte Strukturen bleiben ein Knoten.

### 8. Format und Übergang

`rtt-fusion-graph` Version 2 (Skizze):

```json
{"format": "rtt-fusion-graph", "version": 2,
 "nodes": [
   {"id": "bs", "type": "rtt.splitter", "pose": {"ref": "cad:BS"}},
   {"id": "objA", "type": "group", "variant": "section",
    "pose": {"position": [0, 0, 40], "rotation": [0, 0, 0]},
    "ports": {"front": {"inner": ["L1", "front"]}, "back": {"inner": ["L5", "back"]}},
    "variants": {"section": {"nodes": [], "links": [], "branches": []},
                 "ideal": {"nodes": [], "links": [], "branches": []}}},
   {"id": "probe", "type": "fusion.grating", "pose": {"position": [0, 0, 62.5]}}
 ],
 "links": [
   {"ports": [["bs", "t"], ["objA", "front"]], "transport": "auto"},
   {"ports": [["objA", "back"], ["probe", "front"]], "transport": {"forward": "R", "back": "H"}}
 ],
 "branches": [
   {"name": "messarm",
    "steps": [["src", null, "out", 0], ["bs", "in", "t", "t"], ["objA", "front", "back", 0],
              ["probe", "front", "front", "m0"], ["objA", "back", "front", 0],
              ["bs", "t", "cam", "r"], ["cam", "in", null, 0]]}
 ]}
```

**Abgeleitete Form = abgerollte Äste.** Jeder Besuch eines Knotens in einem Ast wird ein eigener
Knoten der Socket-Form (Knoten@Besuch) mit genau einem Ein- und einem Aus-Socket; gemeinsame
Präfixe verschiedener Äste werden zusammengelegt. Damit ist die Socket-Form auch bei
Mehrfachbesuchen (Ghost, Doppelreflex, zweimal dieselbe Richtung) ein DAG, und Planer Stufe 2,
Ausführung, Cache (R23) und Kohärenz arbeiten unverändert auf ihr. Die bloße Umsetzung „jeder Port
→ zwei Sockets“ reicht dafür nicht, sobald ein Kanal je Ast mehrfach durchlaufen wird.

**Übergang:** Version 1 wird beim Lesen umgesetzt; die gerichteten Bänder und `through` einer
Datei von Version 1 werden zu Ports und zu den Ästen, die diese Bänder beschreiben. Der Import aus
Raytatouille wird neu geschrieben: Ports mit Seiten (vorn/hinten aus der Lage der Flächen) statt
`from_X`/`to_Y`, Raytatouilles Pfade als benannte Äste.

## Abnahme (Sollwerte unabhängig, Regel 2)

| Nr. | Fall | Sollwert |
| --- | --- | --- |
| PA1 | Rundreise Version 1 → 2 → abgerollte Socket-Form | Graphen und Ergebnisse von W4–W38 unverändert |
| PA2 | Linnik mit zwei Verbindungen je Arm statt vier Bändern, Äste aus dem Import | wie W9 |
| PA3 | Platte zwischen zwei Halbräumen aus zwei Grenzflächen-Knoten, Äste bis Umlauf N vom Nutzer eingetragen | Partialsumme der Airy-Reihe bis Glied N (unabhängig ausgerechnet); ohne den N-ten Ast genau das Glied N weniger |
| PA4 | allgemeiner Strahlteiler, alle vier Ports beleuchtet | Energie, Reziprozität S = Sᵀ, Phasenbeziehung r/t des verlustfreien Teilers |
| PA5 | Gruppe aufgeklappt gegen zugeklappt gegen verdichtet | gleich auf 10⁻¹² (R) bzw. in der Toleranz der Tabelle |
| PA6 | Variante wechseln (ideale Linse ↔ Abschnitt) | äußere Verbindungen und Port-Rahmen unverändert; paraxial gleiche Brennweite |
| PA7 | starre Bewegung des ganzen Aufbaus (zufällige Pose, Seed) | alle Ergebnisse invariant |
| PA8 | Michelson, Referenzspiegel über `pose` um α gekippt | Streifenabstand λ / (2 sin α); dz wie W35 |
| PA9 | Löser-Knoten verschoben | Operator-Speicher: `prepare` und Operator nicht neu; Ergebnisspeicher: Ausgaben, Verbindungen und Folgende neu |
| PA10 | Modell ohne Ast zu einem Detektor, Ast über `none`, Ast über einen Port ohne Verbindung | `graph.no_branch` bzw. `graph.branch_impossible`; nichts wird gerechnet |

## Folgen

- **GUI (Raytatouille-Fusion):** Ports als ein Punkt je Seite statt zwei; Verbindungen ohne
  Pfeil; **Äste als eigene Liste** (anlegen, benennen, Schritt für Schritt entlang der
  Verbindungen klicken, hervorheben im Editor und in 3D, prüfen mit Diagnosen); Gruppen auf- und
  zuklappen, Variante wählen; Verschieben und Drehen eines Knotens in der 3D-Ansicht schreibt
  seine Pose, Port-Flächen bewegen sich mit, Diagnosen der Verbindung erscheinen sofort; Marken
  aktuell/veraltet (R23) je Knoten über alle seine Besuche.
- **Raytatouille-Backend:** Abschnitte müssen Licht je Knoten weitergeben (heute ganze Pfade am
  Detektor), sonst sind Äste, Varianten und Gruppen dort nur Anzeige. R19b (Pfad am
  ausgetauschten Block teilen) ist die Vorlage; deshalb kommt diese Stufe vor den Ästen.
- **Planer Stufe 1** (R17, Aufzählung der Fusionspfade aus dem Socket-Graphen) entfällt; die Äste
  des Modells sind die Topologie des Plans.
- **Adapter:** `scatter` in `describe()`; Einfall am hinteren Port (L5) wird zur Voraussetzung für
  `full` beidseitig, Operatoren (L4) für das Verdichten.
- **C++-Port (W2):** erst nach dieser Entscheidung, damit das Graph-Modell nur einmal portiert wird.
- **Nicht Teil:** nichtlineare und zeitabhängige Module (Voxelspaetzle, BindOptiX) bleiben
  Daten-Knoten und Meta-Blöcke; geometrisch gefundene Verbindungen (nicht-sequenziell, M9) später.

## Stufen (Vorschlag)

| Stufe | Inhalt | Abnahme |
| --- | --- | --- |
| **R24 Ports, Äste im Modell, Format 2** | Port mit zwei Kanälen, Streutabelle, Äste als Liste im Modell mit Prüfung, abgerollte Socket-Form, Leser für Version 1, neuer Import aus Raytatouille | PA1, PA2, PA10 |
| **R25 Knotenweise Ausführung im Raytatouille-Backend** | Abschnitte geben Strahlen und Felder je Besuch weiter | Linnik knotenweise = ganze Pfade (W9, W17, W18) |
| **R26 Posen** | Pose je Knoten, Port-Rahmen, relative Lage je Verbindung, Diagnosen, Eingang `pose`, Operator-Speicher | PA7, PA8, PA9 |
| **R27 Allgemeiner Teiler und Mehrfachumläufe** | Teiler mit vier Ports, Ghost-Einträge, vom Nutzer eingetragene Umläufe | PA3, PA4 |
| **R28 Gruppen und Varianten** | Kompositum, exportierte Ports und Parameter, innere Äste, Pfad-Adressen, Auflösen | PA5 (auf/zu), PA6 |
| **R29 Verdichten** | Operator einer Gruppe für ihre inneren Äste | PA5 (verdichtet) |

## Entscheidungen

| Nr. | Frage | Vorschlag |
| --- | --- | --- |
| P1 | Drei Ebenen (Topologie, Äste, Ausführungsplan) | ja |
| P2 | Licht-Ports bidirektional, Daten-Ports gerichtet | ja |
| P3 | Streutabelle `none / oneway / full / ghost` statt `through`, als Erlaubnis für Äste | ja |
| P4 | Pose je Knoten in Raytatouilles Konvention; Geometrie ist die Wahrheit; Eingang `pose`; Ergebnis- und Operator-Speicher getrennt | ja |
| P5 | Verbindung = zwei Einweg-Ausbreitungen aus der relativen Lage der Port-Rahmen; Beobachter sind keine Verbindungen | ja |
| P6 | Äste nur vom Modellierer; keine Aufzählung, keine Grenzen, keine Voreinstellung; `graph.no_branch`; keine Schleifen, keine Reihen, keine Rückkopplung im Framework; ersetzt W6 | **festgelegt** (Maintainer, 2026-11-12) |
| P7 | Gruppen als Kompositum mit Varianten, festen Port-Rahmen und inneren Ästen | ja |
| P8 | Verdichten einer Gruppe zum Operator ihrer inneren Äste | ja, nach L4 |
| P9 | Reihenfolge R24–R29; knotenweise Ausführung (R25) vor Posen und Teiler; vor dem C++-Port des Graphen; R16 (Editor) folgt je Stufe | ja |
