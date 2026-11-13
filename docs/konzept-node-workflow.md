# Konzept: Node-Workflow der Raytatouille Fusion

> **Freigegeben** vom Maintainer am 2026-11-11 (Entscheidung E3 in [roadmap.md](roadmap.md)):
> W1–W9, A1–A6, L1–L12 und die Stufen R14–R20; L11: Plan übernommen als [fusionsplan.md](fusionsplan.md).
>
> Entstanden am 2026-11-11. Anlass: Wunsch des Maintainers, Systeme
> der Fusion wie in Blender als Graph aus Blöcken („Nodes“) mit Ein- und Ausgängen zu entwerfen,
> in einem Fenster „Node-Workflow“ der GUI Raytatouille Fusion (CAD), Repository
> `PhotonPower/Raytatouille-Fusion`. Die Logik entsteht in diesem Repository, beim Port nach C++
> (Entscheidung E1). Grundlagen: Fusion-ADR-Entwurf ([adr-entwurf-fusion.md](adr-entwurf-fusion.md),
> E1, E2), Fusionsplan ([fusionsplan.md](fusionsplan.md), aus Fieldtracing übernommen,
> Abschnitte 6 und 7), Inventar der Module (Anhang A dieses Entwurfs, gelesen am 2026-11-11).
> Nichts davon ist in den Modul-Repositories angelegt (AGENTS.md Regel 4).

## 1. Idee

Ein System der Fusion wird als **Graph** beschrieben: Blöcke mit Namen und Parametern, an denen
Anschlüsse („Sockets“) für Ein- und Ausgänge sitzen, und Verbindungen („Bänder“) zwischen ihnen.
Blöcke lassen sich verschieben (die Verbindungen bleiben), einfügen, austauschen, gruppieren und
stummschalten. Beispiel Linnik-Interferometer:

```
                    ┌──────────────┐ t  ┌────────────┐ hinten ┌──────────┐
 ┌────────┐  ┌────┐ │              ├───►│ Objektiv A ├───────►│  Probe   │
 │ Quelle ├─►│Bel.├►│ Strahlteiler │    │ (Messarm)  │◄───────┤ (Spiegel)│
 └────────┘  └────┘ │      BS      │◄───┤            │        └──────────┘
                    │              │    └────────────┘
                    │              │ r  ┌────────────┐        ┌──────────┐
                    │              ├───►│ Objektiv B ├───────►│ Referenz │
                    │              │◄───┤ (Ref.-Arm) │◄───────┤ (Spiegel)│
                    └──────┬───────┘    └────────────┘        └──────────┘
                           │ (vereinigt)
                    ┌──────▼──────┐   ┌─────────┐
                    │ Tubuslinse  ├──►│ Kamera  │  → Bild, Interferogramm
                    └─────────────┘   └─────────┘
```

- **Linnik-Mikroskop statt Interferometer:** Der Block *Probe (Spiegel)* wird durch einen Block
  *FourierGrating* (periodische Probe, RCWA/FMM) ersetzt; die Verbindungen bleiben.
- **Field Tracing statt Strahlen:** Ein Band zwischen zwei Blöcken trägt standardmäßig, was der
  Planer wählt; der Nutzer kann es auf *Field Tracing* (H) oder ein Spektrum (S) festlegen.
- **Interferenz:** Am Strahlteiler laufen Mess- und Referenzarm wieder zusammen. Gehören beide
  Beiträge zur selben Kohärenzgruppe (Abschnitt 5), addieren sich ihre Amplituden, und die Kamera
  zeigt Streifen; sonst addieren sich die Intensitäten.

Der Graph ist damit die sichtbare, bearbeitbare Form dessen, was der Fusionsplan als **Planer**
vorsieht (Plan Abschnitt 7.2: „Plan als Datenobjekt …, der Nutzer kann jeden Punkt
überschreiben“). Er ersetzt den Planer nicht, sondern gibt ihm Struktur und Nutzerwünsche vor.

## 2. Zwei Ebenen

| Ebene | Inhalt | Wer bearbeitet | Entsteht aus |
| --- | --- | --- | --- |
| **Arbeitsgraph** | Blöcke (Quellen, Gebiete, Strahlteiler, Detektoren, Daten, Auswertungen), Sockets, Bänder mit gewünschtem Transport („auto“ oder fest), Parameter | Nutzer im Node-Editor, Skripte über die Python-API | von Hand; oder automatisch aus den Pfaden eines Raytatouille-Systems (Abschnitt 8.3) |
| **Ausführungsplan** | Fusionspfade mit Kanälen, Darstellung und Abtastung je Port, Operator je Übergang, Fehlerbudget, Kohärenzgruppen, Begründungen; DAG der Löser-Aufrufe | Planer (reine Funktion); Nutzer überschreibt einzelne Punkte im Arbeitsgraphen | Arbeitsgraph + Anforderungen (Ausgaben, Toleranz, Ressourcen), Plan Abschnitt 7.2 |

Der Arbeitsgraph bleibt klein und verständlich (ein Objektiv ist **ein** Block, nicht 17 Flächen);
der Ausführungsplan darf groß werden (Ordnungen, Mehrfachdurchgänge, Wellenlängen, Feldpunkte).
Ein ausgewählter Block oder ein Band zeigt im Editor, was der Planer daraus gemacht hat
(Darstellung, Abtastung, Kosten, Fehlerbeitrag, Begründung).

## 3. Elemente des Arbeitsgraphen

### 3.1 Blöcke

| Art | Beispiele | Sockets | Bemerkung |
| --- | --- | --- | --- |
| **Quelle** | Punkt/Feldpunkte, ebene Welle, Gauß-Strahl, Köhler-Beleuchtung, Dipol-Ensemble (aus Q) | aus: Licht | legt Wellenlängen, Feldpunkte, Polarisation und damit die Kohärenzgruppen fest |
| **Gebiet** | Tracer-Abschnitt (sequenzielle Flächen, z. B. ein Objektiv), freie Ausbreitung, Löser-Gebiet (FourierGrating, BWPM, CliffordBem, SpeckleBem, hp-FEM) | ein/aus je Seite (vorn, hinten) | ein Gebiet = genau ein Verfahren (Fusion-ADR Punkt 2); zweiseitig, damit Licht zurücklaufen kann |
| **Teiler / Vereiniger** | Strahlteiler, Gitter mit mehreren Ordnungen, polarisierender Teiler | mehrere aus bzw. ein; je Ausgang ein **Kanal** (Ordnung, r/t) | der Kanal entspricht `order` am Ereignis (ADR 0025) bzw. dem Kanal eines Fusionspfads (Plan 7.1) |
| **Detektor** | Kamera, Leistungsmesser, Lochblende + Detektor, Fernfeld | ein: Licht; aus: Bild/Werte | addiert je Kohärenzgruppe Amplituden, über Gruppen Intensitäten |
| **Umwandlung** | R → H, H ↔ G, G ↔ S, Fernfeld → S, S → R, G → P | ein/aus verschiedener Darstellung | meist vom Planer eingefügt; sichtbar nur, wenn der Nutzer sie festlegt |
| **Daten** | BioOptoMD (N, Q), Voxelspaetzle (Struktur, Höhenkarte), Material, Phasenfunktion/Mueller-Tabelle | Daten-Sockets | liefern Strukturen und Parameter für Gebiete; tragen keine Leistung |
| **Auswertung** | Kamerabild, Interferogramm, PSF/MTF, Spot, Effizienz je Ordnung, Merit | ein: Licht oder Daten; aus: Ergebnis | Ergebnisse gehen an die GUI (Bild, Kurve, Feld auf Fläche) |
| **Gruppe** | „Linnik-Objektarm“ = Objektiv + Probe | Sockets der Gruppe nach außen | Untergraph; wiederverwendbar und austauschbar wie ein Block |
| **Meta** | Parameter-Sweep, Replikate (Monte Carlo, Prozessstreuung), Optimierung, Zeitschleife, Iteration | umschließt einen Untergraphen | Abschnitt 7 |

Jeder Block hat **Parameter** als Daten (keine Python-Callbacks im Format, Entscheidung O2/F3 von
Raytatouille) mit Typ, Einheit, Grenzen und Gruppe – dieselbe Beschreibung, die die GUI heute für
Plugin-Parameter verwendet (`nfcad.Param`: Zahl mit Einheit, Auswahl, Liste, Datei, Punkt, Objekt
aus dem Modell). Blöcke, die Geometrie haben (Tracer-Abschnitte, Löser-Gebiete), verweisen auf
Objekte bzw. Gruppen des CAD-Modells oder auf Knoten des Raytatouille-Modellbaums.

### 3.2 Sockets

Ein Socket hat einen **Typ** (im Editor eine Farbe wie in Blender) und eine Seite (vorn/hinten,
ein/aus). Zwei Familien:

**Licht-Sockets** sitzen an einem **Port** (Fusion-ADR Punkt 2): Ebene z = const im lokalen
System des Blocks mit Pose, Richtung (+1/−1), Medium und den Wellenlängen. Sie tragen Licht in
einer der Darstellungen der abgeschlossenen Menge (Fusion-ADR Punkt 3):

| Kürzel | Inhalt | Stand |
| --- | --- | --- |
| R | `RayBatch` (Strahlen, mit OPL) | Raytatouille |
| H | hybrides Feld U e^(iψ) | Fieldtracing, `rtt-field` geplant |
| G | Gitterfeld (E_x, E_y) | RTT-Fusion `GriddedField` |
| S | ebene Wellen, Ordnungen, Gewichte dΩ | RTT-Fusion `PlaneWaves` |
| T | Jones-Tabelle je Ordnung über (λ, Einfall) | RTT-Fusion R6/R12 |
| P | radiometrische Pakete (Ort, Richtung, Gewicht, Weg, Stokes) | Scattuccino; mit M9 |

Ein Licht-Socket legt die Darstellung nicht fest, sondern nur, **welche** er annehmen kann (z. B.
FourierGrating: S; BWPM: G). Die Darstellung auf einem Band wählt der Planer bzw. der Nutzer.

**Daten-Sockets** tragen keine Leistung und keine Richtung, aber Wellenlänge (bzw. Tabelle über λ),
Bezugsfläche/Pose, Seite (Spiegelung statt Drehung), Herkunft (Löser, Version, Inhalts-Hash,
Konvergenz) und Replikate (Fusion-ADR, Nachtrag R8–R10):

| Typ | Inhalt | Liefert | Nutzt |
| --- | --- | --- | --- |
| N | Brechzahlprofil n(z) | BioOptoMD | Schichtlöser, FourierGrating, hp-FEM |
| N3 | Voxelgitter aus Material-Kennungen + Materialtabelle je λ (optional Δñ je Voxel) | Voxelspaetzle | BWPM, FourierGrating, Raytatouille |
| Q | Dipol-Ensemble je Rezeptorzustand | BioOptoMD | Quellen (Emitter), CliffordBem, hp-FEM |
| F | komplexe Feldkarte E ∈ ℂ³ auf geschachtelten Gittern, je 1 mW | hp-FEM, CliffordBem, `rtt-field` | BindOptiX, Voxelspaetzle (Fokus), BioOptoMD (optische Kraft) |
| Form | benannte Körper (STEP) | Voxelspaetzle, CAD | Löser-Gebiete, Tracer |
| Höhenkarte | z(x, y) | Voxelspaetzle | SpeckleBem, Gitter |
| M | Antwortkarte / Tabelle (z. B. Emitter-Antworttensoren) | Kette K5 | BindOptiX |
| Statistik | Phasenfunktion, Mueller-Tabelle, μ_s, μ_a, g | CliffordBem, Mie | Scattuccino |
| Belegung | f_s(x, y, t) + Ereignisliste | BindOptiX | Daten-Nodes, label-freie Sensoren |
| Ergebnis | Bild, Kurve, Kennzahl mit Standardfehler | Auswertungen | GUI, Optimierung |

Namenskonflikt: „B“ bezeichnet in BindOptiX die Belegung, im Vorschlag von Voxelspaetzle
(RTT-Fusion#20) die STEP-Form. Vorschlag: Kürzel nur für Licht-Darstellungen (R, H, G, S, T, P)
und die eingeführten Datenports N, Q; alle weiteren Daten-Sockets mit Namen statt Buchstaben.

### 3.3 Bänder

Ein Band verbindet einen Aus-Socket mit einem Ein-Socket und trägt einen **Transport**:

- `auto` (Standard): Der Planer wählt Darstellung und Operator (Plan 7.2, Stufe 2: kürzester Weg
  im Graphen der Zustände mit Kosten aus `estimate` und Gültigkeitskennzahlen).
- fest: `R` (Strahlen), `H` (Field Tracing), `G`, `S`, `T`, `P`; gegebenenfalls mit Parametern
  (Strahlzahl, Abtastung, Ordnungen).

Liegen die beiden Port-Ebenen nicht aufeinander, enthält das Band die **freie Ausbreitung** in
einem homogenen Medium dazwischen (Winkelspektrum, Field Tracing oder Strahlen). Alles mit
Flächen ist ein Block, kein Band. Die Pose des Bands ergibt sich aus den Posen der Ports; der
Editor zeigt an, wenn zwei verbundene Ports nicht zueinander passen (Medium, Richtung, Abstand).

Daten-Bänder tragen keinen Transport; sie sind Abhängigkeiten im Ausführungsplan, keine
Pfadabschnitte (BindOptiX, ADR-Entwurf Datenports).

## 4. Beispiel Linnik im Einzelnen

| Block | Art | Inhalt | Sockets |
| --- | --- | --- | --- |
| Quelle | Quelle | Feldpunkte und Wellenlängen des Systems (`fields`, `wavelengths`) | aus |
| Beleuchtung | Gebiet (Tracer) | Beleuchtungsblende | ein, aus |
| BS | Teiler/Vereiniger | Strahlteiler, Kanäle `r` und `t`; vereinigt die Rückwege | ein (Quelle), aus r, aus t, ein r, ein t, aus (vereinigt) |
| Objektiv A, B | Gebiet (Tracer), zweiseitig | Flächen A.G5…A.L1 bzw. B.… aus dem `.rtt.json` | vorn ein/aus, hinten ein/aus |
| Probe, Referenz | Gebiet | Spiegel (Tracer) – austauschbar gegen FourierGrating, CliffordBem, hp-FEM | ein, aus |
| Tubuslinse | Gebiet (Tracer) | T1, T2 | ein, aus |
| Kamera | Detektor | CAM | ein; aus: Bild |

- **Hin und zurück:** Ein Objektiv wird zweimal durchlaufen. Im Arbeitsgraphen ist es ein Block
  mit vier Sockets; die Bänder BS → vorn-ein, hinten-aus → Probe, Probe → hinten-ein, vorn-aus → BS
  bilden auf der Ebene der Sockets einen azyklischen Graphen. Erst Mehrfachreflexionen erzeugen
  Kreise (Abschnitt 7).
- **Herkunft:** Raytatouille beschreibt den Linnik heute mit vier expliziten Pfaden (`messarm`,
  `referenzarm` und zwei Pupillenpfade). Der Import (Abschnitt 8.3) legt die gemeinsamen
  Abschnitte als einen Block an und die Verzweigung am BS als Teiler.
- **Kamerabild mit Interferenz:** Für Strahlen-Bänder genügt die OPL je Strahl: an der Kamera aus
  beiden Armen R → H bzw. der Zwei-Strahl-Ansatz von `rt.analysis.interferogram` (Raytatouille
  #210, erste Ordnung im Strahlversatz), je Kohärenzgruppe kohärent, über Gruppen inkohärent.

## 5. Kohärenz

Regeln nach Fusion-ADR Punkt 7 und Plan 7.3, im Graphen so umgesetzt:

- Jede Quelle vergibt **Kohärenz-Kennungen**: (Quelle, Feldpunkt bzw. Mode, Wellenlänge,
  Polarisationszustand der Quelle, Replikat). Sie wandern über die Bänder mit.
- Ein **Vereiniger** (Strahlteiler, Detektor) addiert Amplituden innerhalb einer Kennung und
  Intensitäten bzw. Stokes-Vektoren über Kennungen. Wegunterschiede über der Kohärenzlänge (aus
  der spektralen Breite der Quelle) machen Beiträge derselben Kennung inkohärent.
- Kommt ein Beitrag in einer inkohärenten Darstellung an (P, Intensitätsbild), warnt der Planer,
  statt still Intensitäten zu addieren.
- Ensembles (Dipole aus BioOptoMD, Speckle-Realisierungen aus SpeckleBem, Prozessstreuung aus
  Voxelspaetzle) werden **spät** gemittelt, je Probe und nicht über gemittelte Tensoren.

## 6. Auswertung

- **Inkrementell wie in Blender:** Ändert sich ein Parameter, werden nur die nachfolgenden Blöcke
  neu gerechnet. Jedes Zwischenergebnis wird unter einem Inhalts-Hash zwischengespeichert
  (Struktur, Wellenlängen, Einfall, Löser- und Adapter-Version; Fusion-ADR Punkt 7).
- **Operatoren wiederverwenden:** Alle kohärenten Löser sind linear und haben intern
  Wiederverwendung (hp-FEM: eine Faktorisierung je Frequenz; SpeckleBem: `set_excitation` behält
  den Operator; CliffordBem: Krylov-Recycling; FourierGrating: Schichtmodelle einmal je Struktur).
  Keiner der Adapter bietet sie bisher an (`operator: False`). Vorschlag: optionales
  `operator()` bzw. ein langlebiges Handle mit mehreren `apply`, damit der Ausführungsplan bei
  vielen Feldpunkten, Wellenlängen und Optimierungsschritten einmal tabelliert (Plan 7.2, Stufe 4).
- **Vorschau je Band:** Ein Klick auf ein Band zeigt Licht oder Daten an diesem Port (Strahlen in
  3D, Feld als Bild, Spektrum) – das wichtigste Werkzeug bei der Fehlersuche einer Kopplung.
- **Stummschalten/Überbrücken** eines Blocks: Licht läuft ungeändert durch bzw. der Block fällt
  weg. Das entspricht für Tracer-Elemente dem Schalter `active` aus Raytatouille #380.
- **Fortschritt, Abbruch, Ressourcen** je Block; Löser mit GPU oder viel Speicher laufen in einem
  eigenen Prozess (Fusion-ADR Punkt 5, Weg 3).

## 7. Schleifen und Meta-Blöcke

Der Kern bleibt ein azyklischer Graph. Alles, was wiederholt, läuft über Meta-Blöcke, die einen
Untergraphen umschließen:

| Meta-Block | Zweck | Beispiel |
| --- | --- | --- |
| Sweep | Parameter variieren | Wellenlänge, Einfallswinkel, Defokus |
| Replikate | Streuung | Prozessstreuung (Voxelspaetzle), Speckle-Realisierungen, MD-Replikate |
| Optimierung | Parameter gegen eine Merit | Probe-Gitter, Vorkompensation des Drucks (Voxelspaetzle V8); adjungierte Gradienten, wo ein Löser sie hat (FourierGrating, hp-FEM) |
| Zeitschleife | quasistatische Optik in einem zeitabhängigen Prozess | BindOptiX: Optik einmal je Design, Biologie über s bis h; Neuberechnung nur bei label-freien Sensoren |
| Iteration | Rückkopplung bis zur Konvergenz | Mehrfachreflexion, Resonator; optische Kraft ↔ MD (BioOptoMD, Zukunft) |

Die Zeitschleife ist eine eigene Ebene **über** dem optischen Graphen: Sie ruft den gecachten
optischen Untergraphen auf und hat eigene Datentypen (Belegung, Zeitreihen, Kamerabilder über der
Zeit).

## 8. Wo was liegt

### 8.1 Bibliothek (RTT-Fusion, später `rtt-fusion` in C++)

- **Graph-Modell:** Block, Socket, Band, Typen; Prüfung (Typen, Ports, Kreise, Kohärenz);
  Serialisierung als JSON-Format `rtt-fusion-graph` mit Schema-Version.
- **Planer und Ausführung:** Arbeitsgraph → Ausführungsplan (Plan 7.2); Ausführung mit Cache,
  Parallelität und Fortschritt; Ergebnisse im Format `raytatouille-result`.
- **Block-Typen** registrieren sich wie Adapter (`AdapterRegistry`); jeder Adapter beschreibt in
  `describe()` zusätzlich seine Sockets und Parameter (Abschnitt 9).
- **Python-Anbindung** für Skripte, Tests und die GUI.

So bleibt der Graph ohne GUI nutzbar, skriptbar und mit Sollwerten testbar – nach der Regel
„Erst Test, dann Code“.

### 8.2 GUI (Raytatouille Fusion)

- **Fenster „Node-Workflow“:** Blöcke als Kästen mit Namen und Parametern, Sockets als farbige
  Punkte („Magnete“), Bänder als Kurven; Verschieben mit erhaltenen Verbindungen, Einfügen eines
  Blocks auf ein Band (teilt es), Austauschen, Gruppieren, Stummschalten; Bibliothek der
  Block-Typen aus den angemeldeten Adaptern.
- **Kopplung mit dem 3D-Modell:** Ports erscheinen als kleine Flächen in den Ansichten; die Auswahl
  eines Blocks markiert seine Geometrie und umgekehrt.
- **Parameter** im Plugin-Panel, wie heute (`nfcad.Param`).
- **Ergebnisse** in den vorhandenen Ansichten: Strahlen, Felder auf Flächen, Bildfenster,
  Plot-Fenster.
- Die GUI rechnet nicht selbst; sie spricht mit der Bibliothek über deren Python-Anbindung, wie
  heute mit dem Raytracer-Plugin (eigener Prozess).

### 8.3 Raytatouille

- Gebiete und Fusionspfade im Systemmodell (Fusion-ADR Punkt 6, mit M9).
- **Import:** Aus den expliziten Pfaden eines Systems entsteht ein Arbeitsgraph – gemeinsame
  Pfadabschnitte als ein Block, Verzweigungen als Teiler. Ein einfaches Objektiv wird ein einziger
  Block; der Graph ist für sequenzielle Systeme nicht lästiger als heute.
- **Export:** Der Ausführungsplan schreibt Fusionspfade mit Kanälen zurück, soweit sie
  Raytatouille rechnen kann.

## 9. Erweiterungen der Adapter-Schnittstelle (Vorschläge)

| Nr. | Erweiterung | Grund |
| --- | --- | --- |
| A1 | `describe()["sockets"]`: Name, Seite, Licht oder Daten, zulässige Darstellungen bzw. Datentyp, Port-Ebene | Block im Editor ohne Sonderwissen |
| A2 | `describe()["parameters"]`: Name, Typ, Einheit, Grenzen, Standard, Gruppe | Parameterformular, Prüfung, Sweeps |
| A3 | `apply(handle, incident, port=…)`: Einfall auch am hinteren Port | Rückwege (Linnik, Reflexe); FourierGrating und BWPM nehmen heute nur vorn |
| A4 | optional `operator()` oder langlebiges Handle mit mehreren `apply` | Wiederverwendung (Abschnitt 6) |
| A5 | `describe()["validity"]` maschinenlesbar (Kennzahlen und Grenzen) | Darstellungswahl des Planers; heute freier Text, nur Fieldtracing hat `selection.Situation` |
| A6 | binäre Nutzdaten neben dem JSON-Umschlag | Voxelgitter mit ~10⁸ Zellen (Voxelspaetzle), große Feldkarten |

## 10. Lücken aus dem Inventar

Vorschläge dazu in Abschnitt 10a.

- **Zwei Gitterfeld-Typen:** RTT-Fusion `GriddedField` (ny, nx, 2), ohne Trägerphase; Fieldtracing
  `GriddedField` quadratisch N × N × 3 mit quadratischer Trägerphase C und Gitter L. Vorschlag: ein
  G-Typ mit optionaler Trägerphase, Umwandlung mit Sollwerten.
- **FFT-Normierung:** RTT-Fusion G → S diskret (1/N)Σ, Fieldtracing kontinuierlich (2π)^(−d)∫; die
  Umrechnung gehört in die Umwandlung, mit Test.
- **Fehlende Typen:** R, H, P, M und die neuen Daten-Sockets haben noch keine Port-Typen in
  RTT-Fusion; Leistungsnormierung Strahl ↔ Feld nach Fieldtracing ADR 0004 folgt mit dem
  Strahlen-Port.
- **Adapter:** Keiner bietet einen Operator; FourierGrating und BWPM nur Einfall vorn; hp-FEM nur
  für Gitter angebunden (3D, axialsymmetrisch und Wellenleiter-Ports nicht); SpeckleBem nur
  geschlossene Körper (raue Fläche fehlt, Einfall mit beliebigem Spektrum fehlt, D12).
- **Scattuccino** braucht P und eine Statistik-Schnittstelle (Phasenfunktion, Mueller), die noch
  kein Modul liefert.

## 10a. Vorschläge zu den offenen Punkten

Je Lücke ein Vorschlag mit Abnahme. Was ein Modul-Repository ändern müsste, wird dort als
Issue-Entwurf vorgeschlagen (Regel 4), nicht hier umgesetzt.

### L1 Ein Gitterfeld-Typ G

- **Kanonische Form** (wie heute in RTT-Fusion): transversale Komponenten (E_x, E_y) auf einem
  zentrierten, rechteckigen Gitter (n_x, n_y, dx, dy) in der Port-Ebene; E_z folgt aus der
  Transversalität (Port in homogenem Medium, `docs/konventionen.md`).
- **Optionale Trägerphase** wie in Fieldtracing: physikalisches Feld = Werte · exp(i/2·ρᵀCρ) mit
  C (2 × 2, 1/mm²) und optional ein linearer Kipp k₀⊥ (exp(i k₀⊥·ρ)). Ohne Träger gilt C = 0,
  k₀⊥ = 0 – alle heutigen Ports bleiben gültig.
- **Gitter (Lattice) L** und der Frequenzbereich bleiben intern in Fieldtracing bzw. `rtt-field`;
  an einem Port der Fusion liegt G immer im Ortsraum auf einem achsparallelen Gitter. Die
  Umwandlung „Träger abziehen / neu abtasten“ ist ein eigener geprüfter Operator.
- **Abnahme:** Rundreise Fieldtracing-Gitterfeld (quadratisch, 3 Komponenten, mit C) → G mit
  Träger → G ohne Träger (abgetastet) → S gegen das analytische Spektrum eines Gauß-Strahls mit
  Krümmung; Leistung nach Parseval vor und nach jeder Stufe gleich.

### L2 FFT-Normierung

Die beiden Normierungen beschreiben dasselbe Feld E(ρ) = ∫∫ A(k⊥) e^(i k⊥·ρ) d²k⊥ (kontinuierliche
Winkelspektrum-Dichte A, Fieldtracing) bzw. E(ρ) = Σ Ẽ_mn e^(i k_mn·ρ) (diskrete ebene Wellen Ẽ,
RTT-Fusion S). Auf einem periodischen Fenster L_x × L_y mit Δk_x = 2π/L_x, Δk_y = 2π/L_y ist die
Summe die Riemann-Summe des Integrals, also (Herleitung hier, keine fremde Gleichung)

    Ẽ_mn = A(k_mn) · Δk_x Δk_y = A(k_mn) · (2π)² / (L_x L_y).

- **Festlegung:** S trägt Amplituden Ẽ (Phase 1 am Port-Ursprung, wie heute); kontinuierliche
  Spektren tragen zusätzlich das Gewicht dΩ (B5). Die Dichte A gibt es nur innerhalb von
  Fieldtracing; der Operator H/G → S rechnet um.
- **Abnahme:** ebene Welle (ein Ẽ), Gauß-Strahl (analytisches A) und Fieldtracing-Gitterfeld liefern
  über beide Wege dieselbe Leistung (Parseval) und dieselben Ẽ auf 10⁻¹².

### L3 Fehlende Port-Typen

Jeder Typ kommt mit dem ersten Nutzer (Fusion-ADR Punkt 3) und Sollwertdateien:

| Typ | Felder | Erster Nutzer | Abnahme |
| --- | --- | --- | --- |
| R | Port + `RayBatch` (Ort in der Port-Ebene, Richtung, OPL, PRT-Matrix P, Gewicht, Wellenlängen-Index, Feld-Index) + Kohärenz-Kennung | Linnik (R15) | Strahlen eines Pfads, an der Port-Ebene geschnitten, gegen Raytatouilles eigene Spur; Leistung = Σ Gewicht |
| H | Fieldtracing `GridlessField` (Punkte, ψ, ∇ψ, Hesse, U) | R18 | R → H nach Fieldtracing ADR 0004 (Leistung n|U|² cos γ dA = Σ Gewicht); Fokus gegen Richards-Wolf |
| P | Pakete (Ort, Richtung, Gewicht, Weg, Stokes, Seed) | Scattuccino (mit M9) | G → P aus `phasespace` mit dokumentierter Unschärfe; Energie |
| M | Tabelle über Orte/Zustände (Antworttensoren) | BindOptiX (K8) | wie dort DP1–DP8 |

**Leistung zwischen R und Feldern:** Strahlen sind leistungsnormiert (P, Gewicht), Felder sind
Amplituden; die Umrechnung an jedem Übergang folgt Fieldtracing ADR 0004 (U_aus aus P_T, s, cos γ,
Flächenverhältnis und Maslov-Index). Ein R-Port trägt deshalb die Quellleistung je Kohärenz-Kennung
als absolute Größe (W), damit Strahl- und Feldbeiträge an einem Vereiniger dieselbe Einheit haben.

### L4 Wiederverwendbare Operatoren

Plan mit Bestandsaufnahme der Löser und Stufen O1–O4: [plan-l4-operatoren.md](plan-l4-operatoren.md).

- **Schnittstelle:** `apply_many(handle, incidents) -> list[PortFields]` als Pflicht (darf intern
  einzeln rechnen), dazu optional `operator(handle, basis_in, basis_out) -> SMatrix`. `SMatrix` ist
  eine Blockmatrix zwischen den Basen der Ports (bei S: Liste der ebenen Wellen bzw. Ordnungen und
  ihre zwei Polarisationen; bei G: Gitterpunkte × 2), mit Inhalts-Hash.
- **Planer:** nimmt den Operator, sobald die Zahl der Einfälle (Feldpunkte × Wellenlängen ×
  Polarisationen × Optimierungsschritte) die geschätzte Schwelle aus `estimate` übersteigt.
- **Bei den Lösern** (Issue-Entwürfe): FourierGrating – die S-Matrix in Python freigeben; hp-FEM –
  `ScatteringOperator`/`solve_many` über `grating` hinaus; SpeckleBem – `set_excitation` im Adapter
  nutzen (ohne Änderung im Löser möglich); CliffordBem – Recycling im Adapter nutzen.
- **Abnahme:** `operator` und `apply` liefern für 20 zufällige Einfälle (Seed) dasselbe auf 10⁻¹²;
  Reziprozität Sᵀ = S in einer symmetrischen Basis, wo der Löser reziprok ist.

### L5 Einfall am hinteren Port

- `apply(handle, incident, port="back")`; `describe()["ports"]` erklärt je Port
  `accepts_incidence`.
- **FourierGrating:** im Adapter die Struktur an z spiegeln und den Lauf als Einfall von vorn
  rechnen (Spiegelung, keine Drehung; Ports und Phasenbezüge entsprechend vertauscht) – ohne
  Änderung im Löser.
- **BWPM:** Die Spezifikation erlaubt Rückwärts-Einfall, `apply` setzt ihn nicht um → Issue-Entwurf
  für BWPM.
- **Abnahme:** symmetrische Struktur (Platte, symmetrisches Gitter) liefert vorn und hinten dasselbe;
  unsymmetrischer Stapel gegen Transfer-Matrix von beiden Seiten; Reziprozität t₁₂ = t₂₁.

### L6 Namen der Daten-Sockets

Kürzel bleiben R, H, G, S, T, P (Licht) und N, Q (eingeführt). Alle weiteren Daten-Sockets heißen
mit englischem Typnamen im Format und deutschem Namen in der Doku: `field_map` (BindOptiX F),
`occupancy` (BindOptiX B), `voxel_materials` (Voxelspaetzle N3), `shape` (Voxelspaetzle B),
`height_map`, `response_map` (M), `scattering_stats`. Abnahme: Prüfung des Graphen lehnt Bänder
zwischen unverträglichen Typen mit Diagnose ab.

### L7 Streustatistik für Scattuccino

- **Daten-Socket `scattering_stats`** je Wellenlänge: μ_s, μ_a (1/mm), g und entweder eine
  Phasenfunktion p(cos θ) als Tabelle oder eine Mueller-Tabelle (S11, S12, S33, S34 über cos θ,
  Bohren-Huffman-Normierung) für kugelsymmetrische bzw. orientierungsgemittelte Streuer.
- **Erzeuger:** Mie (unabhängige Referenz, `tests/mie_reference.py` liegt vor), dann CliffordBem
  (Fernfeld über Orientierungen gemittelt → Mueller; Issue-Entwurf für einen Export, falls die
  Mittelung im Löser bleiben soll).
- **Abnahme:** Mie-Kugel über beide Erzeuger auf < 1 %; Scattuccino mit diesen Daten gegen
  Adding-Doubling (dort vorhanden).

### L8 Raue Flächen (SpeckleBem) und weitere hp-FEM-Löser

Nicht für die ersten Stufen. SpeckleBem braucht Einfall mit beliebigem Spektrum (D12) und den Weg
für offene Flächen; hp-FEM 3D/axialsymmetrisch/Wellenleiter-Ports bekommen Adapter, sobald ein
Pilot sie braucht (K7 als Brute-Force-Referenz des Planers ist der Anlass).

### L9 Gültigkeit maschinenlesbar (A5)

`describe()["validity"]` als Liste von Regeln mit Kennzahl, Bereich, Folge und Beleg:

```json
{"validity": [
  {"metric": "interface_tilt_rad", "max": 0.0, "effect": "reflection_quantitative",
   "else": "reflection_direction_and_phase_only", "evidence": "BWPM E5"},
  {"metric": "absorber_width_wavelengths", "min": 20, "effect": "required", "evidence": "BWPM K9"},
  {"metric": "fresnel_min", "min": "<aus Fieldtracings Auswahlregeln>", "effect": "geometric_propagation",
   "evidence": "<Testfall in Fieldtracing>"}
]}
```

Nur zur Veranschaulichung des Formats: Die BWPM-Regeln geben die gemessenen Grenzen aus dem
Fusion-ADR (B6) wieder; Schwellen und Belege der übrigen Einträge liefern die Adapter selbst.

Das Vokabular der Kennzahlen beginnt mit Fieldtracings `selection.Situation` (lokale
Fresnel-Zahlen, Kantenverhältnis, Abtastzahlen, Krümmungsradius) und den Kennzahlen aus Plan 7.2
(Größenparameter, Winkelbreite, Zahl der Ordnungen); jede neue Kennzahl mit Definition in
`docs/konventionen.md`. Abnahme: Der Planer wählt in den Fieldtracing-Fällen F2/F3 dieselben
Operatoren wie `plan_propagation`.

### L10 Große Nutzdaten (A6)

JSON-Umschlag im Format `raytatouille-result` mit Verweisen auf Beistelldateien statt Inline-Arrays:
`{"array": "voxels.bin", "dtype": "<u2", "shape": [nz, ny, nx], "sha256": "…"}`, rohe Daten in
Little-Endian, zeilenweise (C-Reihenfolge), per Memory-Map lesbar. Keine HDF5-Abhängigkeit im
Kern. Abnahme: Rundreise eines 10⁸-Voxel-Gitters ohne Kopie; Hash geprüft.

### L11 Fusionsplan auf `main`

Der Fusionsplan liegt nur auf dem Fieldtracing-Zweig `docs-fusion-plan`. Vorschlag: dort einen PR
auf `main` (Entscheidung des Maintainers in Fieldtracing), oder den Plan als `docs/fusionsplan.md`
in dieses Repository übernehmen, weil RTT-Fusion ihn heute am häufigsten zitiert.

### L12 Kopplung GUI ↔ Bibliothek

- Die GUI hält den Arbeitsgraphen im Dokument (`.rttf`) als `rtt-fusion-graph`-JSON und schickt ihn
  bei jeder Änderung an einen Plugin-Prozess mit `rtt_fusion`; der Prozess antwortet mit Prüfung,
  Plan (je Block und Band: Darstellung, Kosten, Fehlerbeitrag, Begründung) und auf Anforderung mit
  Ergebnissen über die vorhandene Ergebnis-API (Strahlen, Felder auf Flächen, Bilder, Kurven).
- Der Graph ist die einzige Wahrheit über die Kopplung; die Geometrie bleibt im CAD-Modell bzw.
  im `.rtt.json`, Blöcke verweisen darauf.
- Skizze des Formats für den Linnik:

```json
{"format": "rtt-fusion-graph", "version": 1,
 "nodes": [
   {"id": "src", "type": "rtt.source", "params": {"from_system": "linnik.rtt.json"}},
   {"id": "bs", "type": "rtt.splitter", "params": {"surface": "BS"}},
   {"id": "objA", "type": "rtt.section", "params": {"surfaces": ["A.AP", "A.L1.S1"]}},
   {"id": "probe", "type": "rtt.mirror", "params": {"surface": "SPM"}},
   {"id": "cam", "type": "rtt.detector", "params": {"surface": "CAM", "pixels": [512, 384]}}
 ],
 "links": [
   {"from": ["src", "out"], "to": ["bs", "in"]},
   {"from": ["bs", "t"], "to": ["objA", "front_in"], "transport": "auto"},
   {"from": ["objA", "back_out"], "to": ["probe", "in"]},
   {"from": ["probe", "out"], "to": ["objA", "back_in"]},
   {"from": ["objA", "front_out"], "to": ["bs", "t_in"]},
   {"from": ["bs", "merged"], "to": ["cam", "in"], "transport": "R"}
 ]}
```

## 11. Entscheidungen

| Nr. | Frage | Vorschlag |
| --- | --- | --- |
| W1 | Zwei Ebenen (Arbeitsgraph, Ausführungsplan) | ja (Abschnitt 2) |
| W2 | Ort der Logik | RTT-Fusion, beim Port in `rtt-fusion` (C++20, E1); Python-Prototyp als Referenz mit Sollwertdateien wie bisher |
| W3 | Ort der GUI | Raytatouille Fusion (CAD), Fenster „Node-Workflow“; Zugriff über die Python-Anbindung im Plugin-Prozess |
| W4 | Socket-Typen | Licht: die abgeschlossene Menge des Fusion-ADR (+ P mit M9); Daten: benannte Typen, je neuer Typ mit dem ersten Modul, das ihn braucht |
| W5 | Kohärenz | Kennungen von der Quelle, Regeln in Vereinigern (Abschnitt 5) |
| W6 | Schleifen | nur über Meta-Blöcke; Kern azyklisch – **ersetzt** durch den ADR-Entwurf [Ports, Posen, Äste](adr-entwurf-ports-posen-aeste.md), P6: Äste nur vom Modellierer, endlich, keine Schleifen und keine Reihen im Framework |
| W7 | Format | `rtt-fusion-graph` (JSON, versioniert), im Dokument der GUI (`.rttf`) eingebettet oder verwiesen; Parameter sind Daten |
| W8 | Namen der Daten-Sockets | Kürzel nur für R, H, G, S, T, P, N, Q; sonst Namen (löst „B“) |
| W9 | Erster Pilot | Linnik-Interferometer (Abschnitt 4) mit ebenem Spiegel und Kippung/Versatz der Referenz; Sollwert: Zwei-Strahl-Interferenz mit analytischem Wegunterschied (Streifenabstand, Kontrast), Energie |

## 12. Stufen (Vorschlag)

| Stufe | Inhalt | Abnahme |
| --- | --- | --- |
| **R13 Konzept** | dieses Dokument, Entscheidungen W1–W9 | Freigabe |
| **R14 Graph-Modell** | Block, Socket, Band, Typen, Prüfung, JSON `rtt-fusion-graph`, Python-API; Import aus Raytatouille-Pfaden | Rundreise JSON; Linnik-Import ergibt den Graphen aus Abschnitt 4 |
| **R15 Ausführung ohne Planer** | feste Transporte (R), Kohärenz-Kennungen, Vereiniger, Detektor; Cache | Pilot W9 gegen den analytischen Sollwert |
| **R16 Node-Editor** | Fenster in Raytatouille Fusion (dort eigene Stufe) | Bedienung am Linnik; Vorschau je Band |
| **R17 Planer, Stufen 1–2** | Topologie, Darstellungswahl (Plan 7.2) | Entscheidungstabellen; gleiche Ergebnisse wie R15 bei `auto` |
| **R18 Field Tracing auf Bändern** | H-Transport (Fieldtracing, später `rtt-field`) | Fokus gegen Richards-Wolf; Linnik mit H gegen R |
| **R19 Rigorose Blöcke austauschen** | Probe → FourierGrating (Linnik-Mikroskop), Partikel (CliffordBem) | Gittergleichung bzw. K3 im Graphen |
| **R20 Meta-Blöcke** | Sweep, Replikate, Optimierung; später Zeitschleife (BindOptiX) und Druck-Schleife (Voxelspaetzle) | je Block eigene Abnahme |

## Anhang A: Inventar der Module (gelesen am 2026-11-11)

| Modul | Rolle im Graphen | Licht-/Daten-Sockets | Sprache, Reife | Adapter in RTT-Fusion | Operator |
| --- | --- | --- | --- | --- | --- |
| Raytatouille | Tracer-Gebiete, Quellen, Detektoren | R; erzeugt Pfade, Felder über `rtt-field` | C++, v0.9 | (Bezug der Konventionen) | – |
| Fieldtracing | Transport auf Bändern (H), Umwandlungen R → H, H ↔ G, G ↔ S | H, G | Python, F0–F5, Port F6 geplant | keiner (eigene Kette `chain.Plan`) | linear, reine Funktionen |
| FourierGrating | Gebiet: periodische Schicht | S (Ordnungen), T | C++20, v0.18, CI, CUDA | ✅ R1 (nur vorn, nur propagierende Ordnungen) | intern ja, nicht angeboten |
| BWPM | Gebiet: Brechzahlvolumen | G | Python + C++/CUDA, v0.1, B0–B4 | ✅ R2, C++ R11 | nein |
| CliffordBem | Gebiet: Partikel (auch chiral, beschichtet) | S (Fernfeld → S), Q als Quelle | C++17, v0.63, ohne CI | ✅ R4 | intern (Recycling) |
| SpeckleBem | Gebiet: geschlossener Körper; raue Fläche (fehlt) | S | C++20, CI, CUDA experimentell | ✅ R8 (geschlossene Körper) | intern (`set_excitation`) |
| hp-FEM | Gebiet: Gitter (angebunden); 3D, axialsymmetrisch, Wellenleiter (nicht angebunden) | S/T (Jones je Ordnung), F | C++20, v0.4, CI, cuDSS | ✅ R9 (Gitter) | intern (Faktorisierung) |
| Scattuccino | Gebiet: Gewebe/Volumen (radiometrisch) | P; Statistik ein | Python + C++/CUDA, jung | – (G → P in `phasespace`) | Umgewichtung |
| BioOptoMD | Daten-Block | N, Q aus; F ein (optische Kraft, Zukunft) | C++/CUDA, Python | ✅ R10 (Dateien) | – |
| Voxelspaetzle | Daten-Block (Druckprozess) | Fokus/F ein; N3, Form, Höhenkarte aus | Python, CI | – (Vorschlag RTT-Fusion#20) | nichtlinear (Dosis ∝ I²) |
| BindOptiX | Zeitschleife über dem optischen Graphen | F, M, Q ein; Belegung, Bilder, Kennzahlen aus | Python + C++/CUDA | – (ADR-Entwurf Datenports, D13) | nichtlinear, stochastisch |
| Mikroptim | Entwurf sequenzieller Optik; Tracer-Gebiet über Dateiadapter | Systemdatei; Fokusfeld als möglicher Ausgang | C++20/CUDA, groß | – (Dateiadapter zu Raytatouille geplant) | Optimierer |
| FMM | kein Block: Theorie-Grundlage von FourierGrating | – | nur Dokumente | – | – |
