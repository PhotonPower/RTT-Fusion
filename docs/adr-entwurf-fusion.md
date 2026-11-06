# ADR 00xx (Entwurf): Raytatouille Fusion – Kopplung wellenoptischer Löser über Ports

> **Entwurf aus RTT-Fusion (Stufe R5) für Raytatouille.** Nicht dort angelegt (RTT-Fusion,
> AGENTS.md Regel 4); Nummer und Annahme vergibt der Maintainer von Raytatouille nach dessen
> Prozess. Stand Raytatouille: Commit `708b6dd` (0.9.0, M8 offen). Belege: RTT-Fusion
> `docs/testfaelle.md` (Tests P1, F1–F5, B1, C1–C3, Pilote K1, K3), BWPM `docs/testfaelle.md`
> (K1–K16, Pilot K6), Fieldtracing ADR-Entwurf `rtt-field` (angenommen als Raytatouille #284) und
> Fusionsplan `docs/rtt-fusion/plan.md`. Grundsatzentscheidungen D1, D2, D10, D11 des Maintainers
> vom 2026-11-05/06 in RTT-Fusion `docs/roadmap.md`.

**Status:** vorgeschlagen (für Raytatouille; in RTT-Fusion freigegeben am 2026-11-06, Entscheidung E1)

## Kontext

**Scope heute.** `docs/architecture.md`: „Nicht geplant: rigorose EM-Solver (RCWA, FDTD) …
Dafür sind Import-Schnittstellen (Jones-Tabellen, Coating-Dateien) vorgesehen.“ ADR 0013: keine
GPU in v1. Gitter wirken als dünne Phasenschicht mit einer Effizienz je Ordnung (ADR 0025,
Punkt 5: eine Zahl, unabhängig von Wellenlänge, Winkel und Polarisation). Mit `rtt-field` (ADR
aus #284, Start nach v1.0) kommen Felder dazu: Strahlen → hybrides Feld (H) → Gitterfeld (G).

**Bedarf.** Mikro- und Nanooptik im System (Gitter, Metaflächen, 3D-gedruckte Mikrooptik,
Faserfacetten, Partikel im Fokus, Biosensoren) braucht rigorose oder halbrigorose Löser in
einzelnen Gebieten, deren Ergebnis in den Strahl- bzw. Feldpfad zurückgeht. Solche Löser liegen
beim Maintainer als eigene Programme mit Python-Anbindung vor (MIT): FourierGrating (RCWA/FMM),
CliffordBem (BEM für Partikel), hp-FEM, SpeckleBem (raue Flächen), BWPM (bidirektionale
Wellenpropagation) und BioOptoMD (Strukturdaten für Sensorschichten).

**Zeitpunkt.** M9 (nicht-sequenziell: BVH, Strahlteilung, Streumodelle) beginnt nach v1.0. Die
Pfad- und Teilungsmodelle von M9 sollen Gebiete mit Ports nicht ausschließen; deshalb kommt
dieser Entwurf vor dem M9-ADR (wie der GRIN-Entwurf X1 vor dem Tracer-Ausbau).

**Prototyp und Befunde.** In RTT-Fusion (Python) sind drei Löser über eine gemeinsame
Port-Schnittstelle gekoppelt und gegen unabhängige Referenzen geprüft:

| Nr. | Befund | Beleg |
| --- | --- | --- |
| B1 | Die Kopplung über Ports trägt: zwei völlig verschiedene Verfahren (FourierGrating, BWPM) rechnen denselben Schichtstapel über dieselbe Schnittstelle auf 10⁻⁹ gleich; ein Teilchen im Fokus liefert über Ports und Spektren dieselbe rückgestreute Leistung wie die Rechnung allein im Löser (10⁻⁶). | F5, B1, K3 |
| B2 | Konventionen der Löser sind die teuerste Fehlerquelle und nur im Adapter zu behandeln: FourierGrating schichtet in y, bezieht Phasen lokal und nimmt bei senkrechtem Einfall für die Ordnungen eine andere s-Basis als für den Einfall – amp_s hat dort das umgekehrte Vorzeichen gegenüber Convention A. | F1 |
| B3 | Eine Zahl je Ordnung reicht für rigoros gerechnete Gitter nicht: Schon bei 4,6° Einfall ändern sich die Effizienzen um bis zu 12 %, s und p unterscheiden sich um bis zu 70 %; die Effizienz von ADR 0025 wirkt zusätzlich zum Fresnel-Faktor der Fläche. | K1 |
| B4 | Ports zu Streukörpern müssen Spektren sein: Ohne evaneszente Anteile fehlt einer Gitterfeld-Ebene in endlichem Abstand ein Beitrag, der nur wie 1/r abfällt (6 % in 2–6 µm); für alles, was eine Optik mit NA < n erreicht, ist das Spektrum exakt. | K3 |
| B5 | Kontinuierliche Spektren (Strahl, Fernfeld) brauchen Quadraturgewichte dΩ, sonst ist ihre Leistung nicht bestimmt und sie lassen sich nicht mit diskreten Ordnungen zusammenführen. | K3 |
| B6 | Jeder Löser hat gemessene Gültigkeitsgrenzen, die der Planer kennen muss (BWPM: Reflexe nur an z-senkrechten Flächen quantitativ, gefangene Wellen konvergieren nicht; flache BEM-Elemente: 1–3 % in Wirkungsquerschnitten). | BWPM E5, K4, K9; C2 |
| B7 | Rechenzeit ist kein Grund für einen Kern-Umbau: Die Löser laufen schnell genug in eigenen Prozessen (BWPM-CUDA: Faserfacette 1142² × 1011 in 4 min, 130–150-mal schneller als NumPy). | BWPM K15 |

## Entscheidung

### 1. Scope

Raytatouille koppelt wellenoptische Löser über **Ports**; die Löser bleiben **außerhalb des
Kerns**. Kein Löser wird gelinkt, keiner wird Abhängigkeit einer Kernbibliothek. Der Satz in
`docs/architecture.md` wird ersetzt durch: „Rigorose und halbrigorose Löser (RCWA, BEM, FEM,
Wellenpropagation) sind nicht Teil des Kerns; sie werden über Ports und Adapter gekoppelt
(ADR 00xx). Import von Jones-Tabellen und Coating-Dateien bleibt.“

### 2. Gebiete und Ports

- Ein **Gebiet** ist ein Bereich, in dem genau ein Verfahren rechnet: sequenzielle Flächen
  (Tracer), freie Ausbreitung (`rtt-field`) oder ein externer Löser.
- Ein **Port** ist eine Ebene z = const im lokalen System des Gebiets mit Richtung (+1/−1),
  Medium (homogen in der Ebene) und Vakuumwellenlänge. Geschlossene Flächen sind keine Ports;
  ein Streukörper hat Ports in Ebenen vor und hinter ihm (B4).
- Gebiete tauschen Licht nur über Ports. Ein Pfadabschnitt in einem Löser-Gebiet trägt einen
  **Kanal** (Ausgangsport, bei periodischen Strukturen die Ordnung), analog zu `order` am
  Ereignis (ADR 0025).

### 3. Darstellungen an Ports (abgeschlossene Menge)

| Kürzel | Inhalt | Herkunft |
| --- | --- | --- |
| R | `RayBatch` (Strahlen) | vorhanden |
| H | hybrides Feld U e^(iψ) | `rtt-field` |
| G | Gitterfeld (E_x, E_y) auf zentriertem Gitter, Spektrum auf (0, 0) bezogen | `rtt-field` |
| S | ebene Wellen (k_x, k_y, E_x, E_y), Phase 1 am Port-Ursprung; optional Ordnungen (n, m) und Gewichte dΩ für kontinuierliche Spektren | neu (B4, B5) |
| T | tabellierte Interaktion: Jones-Matrix je Ordnung über (λ, Einfallsrichtung) | neu, Punkt 4 |

E_z folgt aus der Transversalität. Leistung je Port als Poynting-Fluss: für G über Parseval,
für S je Welle Re(k_z)|E|², für kontinuierliche S (2π)²/k Σ |A|² dΩ. Umwandlungen sind geprüfte
Operatoren mit Sollwertdateien (`raytatouille-result`): R → H und H ↔ G (`rtt-field`), G ↔ S
(FFT mit NumPy-Rundung der Frequenzen; im Prototyp P1), Fernfeld → S (Weyl, transversal
projiziert; C2, K3), S → R (eine Strahlfamilie je Welle bzw. Ordnung; folgt mit P3).
Radiometrische Ports (Stokes/Strahldichte) kommen mit M9, Datenports (Strukturen, Emitter) mit
dem ersten Modul, das sie braucht; jede weitere Darstellung braucht ein ADR.

### 4. Tabellierte Interaktion an Flächen (Erweiterung von ADR 0025)

- Neues optionales Flächenfeld (Arbeitstitel `diffraction_table`), nur an Flächen mit
  Phasenschicht: Verweis auf eine Tabellendatei im Format `raytatouille-result`.
- Inhalt je Ordnung m und Seite (durchgelassen/reflektiert): **Jones-Matrix 2 × 2** von der
  s/p-Basis der einfallenden Welle in die s/p-Basis der Ordnung (Convention A, s = k × N), als
  Funktion von (λ, tangentialen Richtungskosinus der einfallenden Welle); Interpolation glatt in
  den Stützstellen (Regel und Prüfung in der Umsetzung).
- Die Tabelle **ersetzt** die Interaktion der Fläche und die Effizienz: P_event ist die
  leistungsnormierte PRT-Matrix aus der Jones-Matrix (ADR 0021), sodass Grenzflächenverluste
  nicht doppelt zählen (B3). Fehlt die Tabelle, gilt ADR 0025 unverändert.
- Erzeugt wird die Tabelle von einem Adapter (z. B. FourierGrating); Raytatouille rechnet sie
  nicht selbst.
- Präzisierung als eigenes ADR möglich, weil sie unabhängig vom Rest umgesetzt werden kann.

### 5. Adapter

- Ein Adapter je Löser, in einem eigenen Paket (Lizenzen, Abhängigkeiten, GPU-Treiber bleiben
  draußen), mit `describe()`, `prepare(struktur, wellenlängen)`, `apply(handle, einfall)`,
  `estimate(handle)`; optional `operator()` (Streumatrix) für lineare Löser.
- `describe()` erklärt Darstellung, Ports, Linearität, Reziprozität, Determinismus,
  umgerechnete Konventionen und **gemessene Gültigkeitsgrenzen** (B6).
- **Konventionen werden nur im Adapter umgerechnet** (B2). Ein Adapter wird aufgenommen, wenn er
  den Konformitätssatz besteht: Platte gegen Airy/Transfer-Matrix (s, p, schräg, Azimut, auch
  senkrecht), Energieerhaltung, Reziprozität, Determinismus, Phasenbezug am Port; dazu
  verfahrensspezifisch (Gittergleichung, Mie) und, wo möglich, die Querprobe gegen einen zweiten
  Löser.
- Ausführung: zunächst im Python-Prozess über die Python-API (Schicht Workflows); Unterprozess
  als Standard, sobald ein Löser GPU oder viel Speicher braucht (B7). Daten an der Grenze sind
  Werte im Format `raytatouille-result`.

### 6. Modell und Pfade

- Neuer Abschnitt `regions` im Systemmodell: Gebiet mit Löser-Kennung, Verweis auf eine
  Strukturdatei, Pose, Ports (Ebenen, Medien). Löser-Parameter sind Daten, kein Code (wie die
  Entscheidung zu Python-Callbacks, O2/F3).
- Neues Ereignis in Pfaden: Durchgang durch ein Gebiet mit Kanal. Im sequenziellen Tracer
  zunächst nur über T (Punkt 4); Felder und Spektren im Workflow über `rtt-field`.
- M9: Ein Löser-Gebiet ist für den nicht-sequenziellen Tracer ein Körper mit Ports, an dem
  Strahlen enden und neue Strahlfamilien entstehen; eine tabellierte BSDF eines Lösers und ein
  Streumodell von M9 sind zwei Formen derselben Interaktion. Der M9-ADR nimmt das auf.

### 7. Planung, Kohärenz, Cache

- **Planer** als reine Funktion; Ergebnis ist ein Plan als Datenobjekt (Abschnitte,
  Darstellungen, Operatoren, Abtastung, Fehlerbudget, Begründung), vom Nutzer überschreibbar –
  wie die Auswahlregeln von `rtt-field`. Einzelheiten in einem eigenen ADR mit P4 des
  Fusionsplans.
- **Kohärenz:** Amplituden addieren sich innerhalb einer Kohärenzgruppe (gleiche Wellenlänge,
  gleicher Quellpunkt/Modus, gleicher Polarisationszustand), Intensitäten über Gruppen. Diskrete
  Wellen und gewichtete Spektren werden nur getrennt oder nach Umrechnung addiert (B5).
- **Cache und Determinismus:** Inhalts-Hash aus Struktur, Wellenlängen, Löser- und
  Adapter-Version. Der Kern bleibt bitgleich deterministisch; Löser erklären ihre Toleranz (z. B.
  GPU), sie steht im Ergebnis.

### 8. GPU

ADR 0013 bleibt für den Kern. Löser dürfen die GPU in ihrem eigenen Prozess nutzen (B7).

### 9. `rtt-fusion` in C++ (Maintainer, 2026-11-06)

- Neue Bibliothek `rtt-fusion` in C++20 in der Schicht Workflows, wie die übrigen Bibliotheken
  (ADR 0001, 0003, 0004): Port-Typen (Port, GriddedField, PlaneWaves mit Ordnungen und dΩ),
  Umwandlungen (G ↔ S, Fernfeld → S), Leistung, Kohärenzregeln, Inhalts-Hash, die abstrakte
  Adapter-Schnittstelle mit Registry und später der Planer. Python-Anbindung über `rtt-py`
  (nanobind, ADR 0002).
- **Adapter bleiben außerhalb.** Die Schnittstelle ist eine abstrakte C++-Klasse
  (`describe`, `prepare`, `apply`, `estimate`); Umsetzungen gibt es auf drei Wegen:
  1. Adapter in Python (die vorhandenen Löser haben Python-Anbindungen): Sie implementieren die
     Schnittstelle über `rtt-py` und werden nur aus Python-Workflows aufgerufen – im Einklang
     mit der Regel, dass Python-Callbacks nur über die API kommen (O2/F3), nie aus dem
     Dateiformat;
  2. Adapter in C++ als eigene Bibliothek außerhalb von Raytatouille, die gegen die öffentlichen
     Header von `rtt-fusion` baut (für Löser mit C++-Kern, z. B. BWPM);
  3. Unterprozess-Adapter: ein Protokoll mit Werten im Format `raytatouille-result` über
     Dateien oder Pipes, für Löser mit GPU, eigenem MPI oder fremder Lizenz.
- Der Python-Prototyp in RTT-Fusion bleibt die Referenz: Seine Sollwertdateien sind die Abnahme
  des Ports, wie bei `rtt-field` (Fieldtracing → Raytatouille).

## Verworfen

- **Löser in den Kern:** Lizenzen, Build-Abhängigkeiten (Eigen, CUDA, MPI) und Wartung von fünf
  Verfahren im Kern; widerspricht ADR 0013 und dem Scope.
- **Nur Datei-Import von Jones-Tabellen:** reicht für Gitter (Punkt 4), nicht für Felder,
  Partikel, Mehrfachreflexe und Rückwege durch das System.
- **Ein gemeinsames Gitterfeld als einziger Port:** falsch bei Streukörpern (B4) und unhandlich
  bei Gittern (diskrete Ordnungen sind natürlicher als S).
- **Python-Callbacks im Dateiformat:** bereits verworfen (O2/F3); Löser-Parameter sind Daten.
- **Erweiterung von ADR 0025 um eine Effizienz je Polarisation:** löst B3 nur halb (keine
  Phase, keine Winkelabhängigkeit, Fresnel weiter doppelt).

## Folgen

- **`docs/architecture.md`:** Scope-Satz (Punkt 1), Schicht Workflows um `rtt-fusion`
  (Ports, Darstellungen, Umwandlungen, Adapter-Registry, später Planer).
- **`rtt-model`/`rtt-compile`/`rtt-trace`:** Feld `diffraction_table` (Punkt 4) mit
  Diagnosecodes (ADR 0022) und Schemasprung; später Abschnitt `regions` (Punkt 6).
- **`rtt-field`:** S als zusätzliche Darstellung neben H und G; G ↔ S.
- **Tests:** Sollwertdateien aus RTT-Fusion (Ports, FourierGrating-Stapel, Gitter-Pilot,
  Teilchen-Pilot) als Abnahme; Konformitätssatz als Testvorlage für jeden Adapter.
- **Reihenfolge** (Fusionsplan P3–P5): (1) `diffraction_table` und Adapter FourierGrating
  (direkter Nutzen im sequenziellen Tracer), (2) S in `rtt-field` mit Port-Typen, (3)
  `regions` und Pfade zusammen mit M9, (4) Planer.
- **Vorschläge zu den offenen Punkten** (RTT-Fusion, 2026-11-06, zur Entscheidung in Raytatouille):
  1. **Nummer und Zeitpunkt:** Dieser Grundsatz-ADR wird das erste ADR nach dem Release v1.0,
     mit der dann nächsten freien Nummer (heute 0043), **vor** dem ADR zu M9. M8 bleibt frei von
     Fusion-Themen; der M9-ADR verweist auf Punkt 6.
  2. **Tabellierte Interaktion als eigenes ADR:** Ja – Punkt 4 wird ein eigener, genauer ADR
     („Tabellierte Beugungsinteraktion“, Nummer direkt nach diesem) mit Tabellenformat, Basis,
     Interpolationsregel, Gültigkeitsbereich, Diagnosecodes und Abnahme aus Pilot K1. Er hängt
     nicht von `regions` ab und kann in einem Release 1.x vor M9 umgesetzt werden; er bringt den
     ersten Nutzen der Fusion im sequenziellen Tracer.
  3. **Reihenfolge nach v1.0:** (a) Port von `rtt-field` (Fieldtracing, angenommen), (b)
     tabellierte Beugungsinteraktion mit FourierGrating-Adapter in Python, (c) `rtt-fusion` in
     C++ mit Port-Typen und Umwandlungen, Abnahme gegen die Sollwertdateien von RTT-Fusion,
     (d) `regions` und Pfade zusammen mit M9, (e) Planer.
