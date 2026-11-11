> **Übernommen** aus Fieldtracing, Zweig `docs-fusion-plan` (Commit `3e0d5f8`, 2026-11-09),
> `docs/rtt-fusion/plan.md`, nach Entscheidung E3 / L11 (Maintainer, 2026-11-11). Ab jetzt wird der
> Plan hier gepflegt; Verweise auf Abschnitte des Plans gelten unverändert.

# Plan: Raytatouille nach v1.0 – Port von Field Tracing und „Raytatouille Fusion“

> **Entwurf** (2026-11-05) zur Abstimmung mit dem Maintainer. Nichts davon ist in Raytatouille
> angelegt (Fieldtracing, AGENTS.md Regel 8); Entscheidungen fallen dort nach dessen Prozess.
> Stand Raytatouille: Commit `6078e55` (nach 0.9.0, M8 offen), gelesen am 2026-11-05.
> Inventar der Module (Abschnitt 2.1) aus den Repositories des Maintainers, gelesen am 2026-11-05;
> BindOptiX ergänzt am 2026-11-09 (K8, D13, Anhang A.1).

## 1. Ziel

Raytatouille soll nach v1.0 optische Systeme rechnen können, in denen geometrische Optik,
Field Tracing und rigorose Wellenoptik zusammenwirken. Vorgesehene Module: sequenzielles und
später nicht-sequenzielles Raytracing (vorhanden bzw. M9), Field Tracing (`rtt-field`), eine
schnelle bidirektionale vektorielle Wellenpropagation (BWPM, geplant), Fourier-Modal-Methode und
RCWA für periodische Strukturen (FourierGrating), BEM für Partikel und raue Flächen (CliffordBem,
SpeckleBem), hp-FEM für allgemeine Nahfeldprobleme und die Molekulardynamik BioOptoMD, die
optische Parameter biologischer Sensorschichten liefert. Bis auf die BWPM liegen alle als headless
Programme mit Python-Anbindung unter MIT-Lizenz vor (Abschnitt 2.1).

Daraus folgen drei Aufgaben, die dieser Plan nacheinander und mit Abhängigkeiten beschreibt:

1. **Port von Field Tracing** (F6, Bibliothek `rtt-field`) nach den freigegebenen Entwürfen in
   `docs/rtt-port/`. Er liefert die Brücke Strahl ↔ Feld, ohne die keine Fusion möglich ist.
2. **Wellenoptische Blackbox:** eine einheitliche Schnittstelle, über die externe Löser als
   Gebiete mit Ein- und Ausgangsports in eine Rechnung eingebunden werden, ohne dass ihr Code in
   den Kern von Raytatouille wandert.
3. **Pfadplanung:** eine Verallgemeinerung der expliziten Pfade von Raytatouille („Pfade sind
   Daten“) auf Pfade durch Gebiete verschiedener Verfahren, mit einem Planer, der je Abschnitt
   Verfahren, Darstellung, Abtastung und Fehlerbudget wählt, wie es Kapitel 7 des Papers für die
   Feldbausteine tut.

## 2. Ausgangslage

**Raytatouille** (C++20, Python-API, MIT, headless und deterministisch):

- Strahlen als Structure-of-Arrays mit leistungsnormierter PRT-Matrix P und Gewicht (ADR 0005,
  0021); Pfade als explizite Ereignislisten, mehrere benannte Pfade, Ghost-Generator (ADR 0027),
  Gitter als Phasenschicht mit Ordnung am Ereignis (ADR 0025).
- Ergebnisformat `raytatouille-result` (ADR 0023), Diagnosecodes (ADR 0022), Änderungen am Modell
  über JSON Patch (ADR 0024), Python-Callbacks nur über die API und nie im Dateiformat
  (Entscheidung zu O2/F3 in `docs/roadmap-erweiterungen.md`).
- Scope laut `docs/architecture.md`: **„Nicht geplant: rigorose EM-Solver (RCWA, FDTD) …
  Dafür sind Import-Schnittstellen (Jones-Tabellen, Coating-Dateien) vorgesehen.“** Keine GPU in
  v1 (ADR 0013). Nach v1.0 als erste Priorität M9 (nicht-sequenziell mit Streuung), dann M10
  (Python-GUI); Beamlets, GPU-Backend und globale Optimierung unter „v2+“.
- Field Tracing: ADR `rtt-field` angenommen (Raytatouille #284), Start nach v1.0 (M8); #285 (C++-
  Ersatz für SciPy-Bausteine) offen.

**Fieldtracing** (dieses Repository): alle Bausteine geprüft, Sollwertdateien unter
`reference/field/`, ADR-Entwurf, API-Skizze und Issue-Entwürfe I1–I3 unter `docs/rtt-port/`.

**Folgerung für die Fusion:** Der bestehende Scope schließt rigorose Löser *im Kern* aus, nicht
ihre *Kopplung*. Der Plan hält daran fest: Löser bleiben außerhalb, der Kern bekommt
Darstellungen, Umwandlungen, Schnittstelle und Planer. Das ändert trotzdem die Scope-Aussage in
`architecture.md` und braucht dort ein ADR (Abschnitt 9, Entscheidung D1).

### 2.1 Inventar der Module

Gelesen am 2026-11-05 (Repositories unter `PhotonPower`, alle MIT, nur lesend).

| Modul | Verfahren | Stand | Sprache, Anbindung | GPU |
| --- | --- | --- | --- | --- |
| **FourierGrating** | Fourier-Modal-Methode: Differentialmethode (implizites RK, S-Matrix) und RCWA, schnelle Fourierfaktorisierung; 1D/2D, isotrop/anisotrop, Nahfeld | v0.12, M0–M13 umgesetzt | C++, Python `fouriergrating` (pybind11, Wheels), CLI `fgsolve`, JSON-Konfiguration | CUDA (cuBLAS/cuSOLVER), Abnahme auf dem Server offen |
| **hp-FEM** | Maxwell zeitharmonisch, Nédélec-Elemente hoher Ordnung, a-posteriori-Fehler, hp-Adaptivität, PML; Streuung, Eigenmoden, Wellenleiter, periodisch (Bloch), geschichteter Hintergrund, Dipole, Formableitungen, Optimierung/UQ | 0.2.0, M0–M9 | C++20, Python `hpfem` (pybind11), JSON-Projekte | cuDSS als optionaler Direktlöser |
| **CliffordBem** | BEM in Cl₃(ℂ) (Dirac-Form, resonanzfreie Transmissionsgleichung), H-Matrix/ACA; Partikel, mehrere Körper, beschichtete und chirale Medien, Dipole, Nahfeld, Kräfte | 0.63 | C++, Python `cliffordbem` (pybind11) | – |
| **SpeckleBem** | Oberflächen-Integralgleichungen (PMCHWT u. a.), RWG, MLFMM; raue Flächen und beliebige Körper, Speckle-Statistik | Phase 4 (MLFMM) weit fortgeschritten, Python-API für Simulation und Felder | C++20, Python `specklebem` (pybind11), `.npy`/HDF5 | geplant (Phase 7) |
| **BioOptoMD** | Molekulardynamik von Biorezeptoren an Sensoroberflächen, Optik-Observer: Übergangsdipole, FRET-Geometrie, Brechzahlprofil n(z) (Lorentz-Lorenz) | Phase 9.3 abgenommen | C++/CUDA, Python `biooptomd` (nanobind, DLPack) | nur GPU |
| **BWPM** (geplant) | schnelle bidirektionale vektorielle Wellenpropagation nach Wende et al. 2023 auf Basis der FPWPM | – | – | vorgesehen |
| **BindOptiX** (geplant) | kein Lichtlöser: Zeit-, Transport- und Biologieschicht für Biosensoren (Langevin-Partikel für Zellen, Konvektion–Diffusion, Bindungskinetik, Signal- und Kameramodell aus optischen Karten) | Z0 (Spezifikation) | Python-Referenz, später C++/CUDA | vorgesehen |

**Konventionen** – der wichtigste Befund, denn hier entstehen die teuersten Kopplungsfehler:

| Modul | Zeitfaktor | Einheiten | Geometrie und Polarisation | Anpassung im Adapter |
| --- | --- | --- | --- | --- |
| Raytatouille / `rtt-field` | e^(−iωt) (Wellen exp(i(k·r − ωt))) | mm, µm, rad | +z Ausbreitung; s = k × N, Fresnel Convention A | Referenz |
| FourierGrating | e^(−iωt) | intern k₀-normiert, API µm und Grad | Schichtung in **y**, Licht von oben (−y); δ = 0 s/TE; Amplituden tangentiales E mit lokalem Phasenbezug je Grenzfläche | Koordinatendrehung y ↔ z, Phasenbezug auf die Port-Ebene, s/p-Basis |
| hp-FEM | e^(−iωt), Im ε > 0 | SI (m), Python mit `units` | frei; Streufeld- oder Gesamtfeldformulierung | nur Einheiten |
| CliffordBem | e^(−iωt) | freie Längeneinheit L, ε₀ = μ₀ = 1 | frei; Strahlen auf Leistung 1 normiert | Einheiten, absolute Leistung |
| SpeckleBem | **e^(+jωt)**, Im ε ≤ 0 | SI (m) | Licht aus z < 0 in +z; p in xz-Ebene, s entlang y | **komplex konjugieren** (Felder, Materialien), Basis |
| BioOptoMD | – (keine Wellenoptik) | nm, ps, kJ/mol, e·nm | Sensoroberfläche bei z, n(z)-Profil | Einheiten; nur Datenkopplung |
| BindOptiX | – (nur Daten); Feldkarten mit e^(−iωt), je mW | µm, s, K, mW; Teilchen/µm³ | Sensoroberfläche z = 0, Lösung z > 0, Schwerkraft als Vektor | Pose und Spiegelregel z → −z zum RTT-System |

**Ein- und Ausgänge**, soweit für Ports relevant:

- FourierGrating: Einfall ebene Welle (λ, Polar- und Azimutwinkel, Polarisation); Ausgabe je
  Ordnung Effizienz und Amplituden s/p, Energiebilanz, Nahfeld. `Simulation.solve` ändert nur die
  Einfallsrichtung ohne Neuaufbau, gibt den GIL frei – Tabellen über k∥ sind damit billig.
- hp-FEM: Quellen ebene Welle, Dipol, frei definierte Felder (Callback, langsam); Ausgabe
  Feld, Beugungsordnungen periodischer Strukturen mit abgezogenem Einfall, Leistungsbilanz.
- CliffordBem: `BeamField` ist bereits ein **Winkelspektrum ebener Wellen** (Richtung, Amplitude,
  Raumwinkelgewicht), aber nur über `focused` (Richards-Wolf) und `gaussian` erzeugbar und auf
  Leistung 1 normiert; Ausgabe Fernfeld, Nahfeld, Extinktion, Kräfte, Emitterraten.
- SpeckleBem: `AngularSpectrumBeam` (rigoros, bisher nur Gaußspektrum), Gauß, ebene Welle;
  Ausgabe Nah- und Fernfeld, RCS/DRC, Speckle-Statistik.
- BioOptoMD: keine Lichtports, sondern **Strukturdaten**: n(z) mit Varianz über Frames,
  Dipolorte und -richtungen, FRET-Größen.

**Folgerungen für den Plan:**

1. Das **Ebene-Wellen-Spektrum** ist die gemeinsame Sprache: CliffordBem und SpeckleBem nehmen
   Winkelspektren, FourierGrating und hp-FEM ebene Wellen je Einfallsrichtung, und `rtt-field`
   erzeugt Spektren mit PFT/FFT. Kugelwellen (T-Matrix) braucht zunächst niemand.
2. Kleine Erweiterungen in zwei Modulen sind nötig: ein Konstruktor `BeamField` aus beliebigen
   Stützstellen mit absoluter Leistung (CliffordBem) und ein `AngularSpectrumBeam` mit beliebigem
   Spektrum (SpeckleBem). Diese Repositories sind für dieses Projekt nur lesbar; die Erweiterungen
   werden dort als Issues vorgeschlagen (Entscheidung D12).
3. BioOptoMD ist kein Lichtlöser, sondern eine **Struktur- und Quellenquelle**: n(z) wird zur
   Schicht in Raytatouilles Coatings, in FourierGrating oder als Dünnschicht in CliffordBem; die
   Dipole werden zu Emittern in CliffordBem und hp-FEM. Der Plan bekommt dafür eine eigene,
   einseitige Datenkopplung (6.2). Ein inkohärentes Strahlungstransport-Modul gibt es nicht; die
   radiometrische Darstellung bleibt für M9 (Streuung, SpeckleBem-BSDF).
4. Alle Löser außer Raytatouille rechnen auf der GPU oder können es. Das stützt die Ausführung als
   Unterprozess (D3) und lässt den Kern von Raytatouille ohne GPU (ADR 0013).
5. Die BWPM fehlt noch und ist zugleich das Verfahren, das am direktesten an `rtt-field` anschließt
   (Gitterfeld in einer Ebene, beide Richtungen). Sie wird als erstes „fusionsfähiges“ Modul von
   Anfang an mit den Konventionen und Ports von Raytatouille gebaut (Abschnitt 6.6).

## 3. Phasen im Überblick

| Phase | Inhalt | Ort | Voraussetzung | Abnahme |
| --- | --- | --- | --- | --- |
| **P0 Vorbereitung** | Inventar der Löser, ADR-Entwurf „Fusion“, Konformitätsprüfung, Prototyp-Repository, Pilotfälle | Fieldtracing / neues Repo | – (jetzt, parallel zu M8) | Entwürfe freigegeben |
| **P1 Port `rtt-field`** | I1 Vektor-Debye-PSF aus Strahlen, I2 Felder in Ebenen, I3 hybride Feldverfolgung | Raytatouille | v1.0, #285 entschieden | Sollwertdateien grün |
| **P2 Blackbox-Prototyp** | Schnittstelle, Pilotlöser FourierGrating und CliffordBem, Darstellungen, Umwandlungen, Pilotfälle K1–K3 in Python; parallel BWPM B0–B3 | Prototyp-Repo | P0; P1 nur für den Strahl-Feld-Teil | Pilotfälle gegen unabhängige Referenz |
| **P3 Kern der Fusion** | Darstellungen und Ports in C++, Gebiete im Modell, Tabellenkopplung im Tracer, Ausführung und Cache | Raytatouille | P1, P2, ADR Fusion | Pilotfälle aus P2 bitgleich bzw. in Toleranz |
| **P4 Pfadplanung** | Topologie, Darstellungswahl, Abtastung, Fehlerbudget, Kohärenz; Abstimmung mit M9 | Raytatouille | P3; M9-Entwurf | Entscheidungstabellen, End-to-End gegen Brute-Force |
| **P5 Ausbau** | BWPM (B4, K4), BioOptoMD (K5), SpeckleBem mit M9 (K6), GPU, Optimierung durch gekoppelte Systeme, GUI | Raytatouille und Adapter | P4 | je Modul eigene Abnahme |

P0 kann sofort beginnen. P1 und P2 laufen nach v1.0 parallel, weil P2 nur für den
Strahl-Feld-Übergang auf P1 wartet und sich bis dahin mit der Python-Referenz aus Fieldtracing
behilft. Die Abfolge folgt dem Muster, das sich bei Field Tracing bewährt hat: Prototyp in Python
mit unabhängigen Referenzen und Sollwertdateien, dann Port.

## 4. Phase P0: Vorbereitung (jetzt bis v1.0)

1. **Inventar der Löser** – erledigt (Abschnitt 2.1). Offen sind dort Laufzeit und Speicher je
   Pilotfall und die Determinismus-Angaben der GPU-Pfade.
2. **Pilotfälle** (Abschnitt 8): FourierGrating (K1, K2) und CliffordBem (K3) für P2, weil beide
   reif, validiert und ohne Änderung der Konvention (e^(−iωt)) koppelbar sind. Die BWPM wird
   parallel als neues Modul entwickelt (6.6) und ist der Pilot der Feldkopplung (K4).
3. **ADR-Entwurf „Fusion“** (Abschnitt 9): Scope, Darstellungen, Schnittstelle, Ausführung,
   Modell, Planer, Kohärenz, Determinismus.
4. **Konformitätsprüfung** für Löser spezifizieren (Abschnitt 6.5), bevor ein Adapter entsteht.
5. **Prototyp-Repository** anlegen (Vorschlag: `PhotonPower/Fusion`, Python, gleiche Regeln wie
   Fieldtracing: Test zuerst, unabhängige Sollwerte, Sollwertdateien im Format
   `raytatouille-result`). Fieldtracing bleibt die Referenz für Feldbausteine und wird dort als
   Abhängigkeit genutzt, nicht erweitert.
6. **M9 im Blick:** Der Entwurf von M9 (BVH, Ray Splitting, Streumodelle) fällt in dieselbe Zeit.
   Wie beim GRIN-ADR (X1, „ADR früh, damit … nicht fest einbauen“) sollte der Fusion-Entwurf vor
   dem M9-ADR vorliegen, damit Pfade und Strahlteilung Gebiete mit Ports nicht ausschließen.

## 5. Phase P1: Port von Field Tracing

Nach `docs/rtt-port/issues.md` und dem angenommenen ADR, in dieser Reihenfolge:

1. **Start-Checkliste** (E8 in `docs/roadmap.md`): `rtt-commit.txt` auf v1.0 setzen,
   `pytest -m rtt` auf dem Rechner des Maintainers, Abweichungen durch M8 (Flächentypen,
   Polygon-Aperturen, Obskurationen) als Issues; #285 entscheiden (Delaunay, KD-Baum,
   Interpolation, Minimax-LP, Bessel); Meilenstein „Field Tracing“ in Raytatouille.
2. **Grundlagen:** `fourier`, `propagation`, `geometric`, `mesh`, `pft`, `sft`, `edge`
   (einschließlich der gleichmäßigen 2D-Lösung, Sollwertdatei `edge_uniform_rectangle`),
   `lpia`, `field`, `selection` – je Header ein Issue mit Sollwertdatei.
3. **I1 Vektor-Debye-PSF aus Strahlen** (erster Integrationspunkt, schließt H7 aus ADR 0033).
4. **I2 Felder in beliebigen Ebenen** mit SFT/FFT nach den Auswahlregeln, Plan als Ergebnis.
5. **I3 hybride Feldverfolgung** durch das System (`chain`), Strahlen von `rtt-trace`.

Für die Fusion wichtig sind aus P1 vor allem: `rays.hpp` (Strahlen → hybrides Feld), die
umgekehrte Richtung (Feld → Strahlen über die lokalen Gradienten, in Fieldtracing als PFT und
`GridlessField` vorhanden), die Auswahlregeln als reine Funktion und der Plan als Datenobjekt.
Diese Teile sollten in P1 so geschnitten werden, dass P3 sie ohne Umbau wiederverwendet.

## 6. Die wellenoptische Blackbox

### 6.1 Gebiete und Ports

Ein **Gebiet** ist ein räumlicher Bereich, in dem genau ein Verfahren rechnet: ein Abschnitt
sequenzieller Flächen (Raytracing), freie Ausbreitung (Field Tracing), oder ein externer Löser
(BWPM, FourierGrating, hp-FEM, CliffordBem, SpeckleBem). Gebiete tauschen Licht nur über **Ports** aus. Ein Port hat:

- eine Geometrie: Ebene mit lokalem Koordinatensystem und Ausdehnung, oder geschlossene Fläche
  (Kugel, Hüllquader) bei Streukörpern;
- eine Richtung (ein, aus, beides) und die Medien auf beiden Seiten;
- eine **Darstellung** (6.2) mit Abtastung, Wellenlänge und Polarisationsbasis;
- eine Kohärenzangabe: kohärentes Feld oder inkohärente Leistungsgröße.

### 6.2 Darstellungen

Der Kern kennt eine kleine, abgeschlossene Menge von Darstellungen; jede Umwandlung zwischen
ihnen ist ein geprüfter Operator. Neue Darstellungen brauchen ein ADR.

| Darstellung | Inhalt | Kohärenz | Herkunft |
| --- | --- | --- | --- |
| **R Strahlen** | `RayBatch`: Ort, Richtung, OPL, P, Gewicht, Status | kohärent über OPL und P, oder inkohärent | Raytatouille |
| **H hybrides Feld** | `GridlessField`: U ∈ ℂ³, ψ, ∇ψ, Hesse-Matrix | kohärent | `rtt-field` |
| **G Gitterfeld** | `GriddedField`: Werte ∈ ℂ³ auf Gitter mit analytischer Trägerphase | kohärent | `rtt-field` |
| **S Ebene-Wellen-Spektrum** | diskrete Liste (k∥, Amplitude ∈ ℂ², Basis s/p) oder Gitterspektrum; für periodische Ports Floquet-Ordnungen | kohärent | neu |
| **M Streumatrix / Tabelle** | linearer Operator eines Gebiets: S-Matrix in Ordnungen, Jones- oder Mueller-Tabelle über (Einfallsrichtung, λ), T-Matrix in Kugelwellen | kohärent (Jones, S, T) bzw. inkohärent (Mueller) | neu, Tabellen in X2 vorgesehen |
| **P radiometrisch** | Strahlen oder Gitter mit Stokes-Vektor bzw. Strahldichte, ohne Phase | inkohärent | neu, erst mit M9 (Streumodelle, SpeckleBem-BSDF) |
| **Q Quellen** | Dipol-Emitter: Ort, elektrisches/magnetisches Moment, Wellenlänge; Ensembles | inkohärent über Emitter, kohärent je Emitter | neu (BioOptoMD → CliffordBem, hp-FEM) |
| **N Struktur** | Brechzahlprofil n(z) bzw. Schichtstapel, mit Streuung über Frames | – (Daten, kein Licht) | neu (BioOptoMD → Coatings, FourierGrating, CliffordBem) |

Q und N sind **Datenports**: Sie verbinden kein Licht, sondern liefern Quellen und Strukturen
einseitig an ein Gebiet. Der Planer behandelt sie als Abhängigkeiten im Ausführungs-DAG, nicht als
Pfadabschnitte. Kugelwellen (T-Matrix) braucht nach dem Inventar zunächst kein Modul, weil beide
BEM-Codes Winkelspektren als Einfall nehmen; sie bleiben ein späterer Inhalt von M.

**Vorgeschlagen (D13):** zwei weitere Datenports für Biosensoren über der Zeit, **F Feldkarte**
(komplexes E auf geschachtelten Volumengittern je mW, optional Absorptionsdichte und ΔT; Richtung
Licht → Biologie) und **B Belegung** (Anteile der Rezeptorzustände über Ort und Zeit, Ereignisliste
für Zellen; gewichtet die Q-Ensembles je Zustand), dazu die Antwortkarte eines Emitters als Tabelle
in M (Zerfalls- und Detektionstensoren, symmetrisch 3×3 je Ort). Entwurf: `PhotonPower/BindOptiX`, `docs/adr-entwurf-datenports.md`.

### 6.3 Umwandlungen

| von → nach | Operator | Gültigkeit (Kennzahl) | Quelle |
| --- | --- | --- | --- |
| R → H | Phase k₀·OPL, Amplitude nach Fieldtracing ADR 0004, Jacobi aus dem Netz | Strahlnetz dicht genug, keine Kaustik am Port | `rtt-field` I3 |
| H → R | Strahlen entlang ∇ψ, Gewicht aus \|U\|² und Jacobi | lokale Fresnel-Zahl ≫ 1 | Fieldtracing `field` |
| H ↔ G | Umtasten nur des Restfelds (ADR 0003 von Fieldtracing) | Abtastung nach den Auswahlregeln | `rtt-field` |
| G ↔ S | FFT, SFT, PFT | Tabelle 7.1 des Papers | `rtt-field` |
| S → R | eine Strahlfamilie je Ordnung bzw. Spektralbereich | Ordnungen getrennt, Bündel groß gegen λ | neu |
| S ↔ M | Anwenden der Matrix auf Ein-Ordnungen; Tabellen über k∥ interpolieren (nur glatte Größen, analog ADR 0003) | Tabelle fein genug, keine Resonanzen zwischen Stützstellen | neu |
| H/G → S für BEM | Winkelspektrum des Einfalls um den Streukörper (PFT bzw. FFT, dann Stützstellen auf der Richtungskugel) | Spektrum vollständig erfasst (Apertur, Evaneszenz vernachlässigbar) | `rtt-field` + Adapter (CliffordBem `BeamField`, SpeckleBem `AngularSpectrumBeam`) |
| Fernfeld → S → H | Fernfeld E∞(k̂) eines Streuers als Spektrum, dann inverse PFT bzw. Strahlen | Abstand ≫ Streuergröße²/λ oder Nahfeldkorrektur | neu |
| N → Struktur | n(z) als Schichtfolge (Diskretisierung, Mittel und Varianz) | Schichtdicke ≪ λ/n oder Profil-RK in FourierGrating | neu |
| R ↔ P | Gewicht und Stokes aus P; Rückweg ohne Phase | nur inkohärente Weiterrechnung | neu |

Jeder Operator bekommt eine Sollwertdatei mit unabhängiger Referenz und, bei Näherungen, einen
Konvergenztest – dieselben Regeln wie in Fieldtracing.

### 6.4 Schnittstelle eines Lösers

Ein Löser wird über einen **Adapter** eingebunden, der vier Aufrufe bereitstellt:

1. `describe()` – Metadaten: Name, Version, Verfahren, unterstützte Darstellungen je Port,
   Linearität, Reziprozität, Konventionen (Zeitvorzeichen, Fourier, Polarisationsbasis,
   Einheiten), Determinismus (bitgleich / in Toleranz), Ressourcen (GPU, Speicher).
2. `prepare(structure, wavelengths)` – Struktur einlesen und vorrechnen (etwa Eigenmoden der
   RCWA-Schichten, BEM-Matrix); liefert einen Griff und einen Inhalts-Hash.
3. `apply(handle, port_in, field) → {port_out: field}` – ein Einfallsfeld durchrechnen.
4. `operator(handle, port_in, port_out, basis)` – optional für lineare Löser: die Streumatrix bzw.
   Tabelle in einer gewählten Basis. Damit kann der Kern vorrechnen, zwischenspeichern und
   dieselbe Rechnung für viele Strahlen, Felder und Optimierungsschritte nutzen.
5. `estimate(handle, request) → (Laufzeit, Speicher, Fehler)` – Kostenmodell für den Planer;
   anfangs grob, aus Messungen kalibriert.

Die Daten an der Grenze sind **Werte, keine Objekte**: Darstellungen nach 6.2, kodiert wie
`raytatouille-result` (ADR 0023). Damit ist dieselbe Schnittstelle in drei Ausführungsarten
nutzbar:

| Ausführung | Wann | Eigenschaften |
| --- | --- | --- |
| **in-process Python** | Prototyp, vorhandene Python-Anbindungen | schnell zu bauen; GIL und Fremdabhängigkeiten im selben Prozess |
| **Unterprozess** | Standard für Produktion | Isolation (Abstürze, Lizenzen, GPU-Treiber), Austausch über Dateien oder Shared Memory, Cache auf Platte |
| **entfernter Dienst** | große FEM-/BEM-Läufe auf dem Rechner mit viel Speicher | gleiche Daten, Übertragung über Netz; erst bei Bedarf |

Der Kern von Raytatouille hängt von keinem Löser ab und linkt keinen. Adapter liegen in eigenen
Paketen (`raytatouille-fusion-rcwa` usw.) mit eigener Lizenzprüfung. Das hält den MIT-Kern sauber
und erfüllt die Scope-Aussage „rigorose Löser nicht im Kern“.

### 6.5 Konformitätsprüfung

Ein Adapter wird erst aufgenommen, wenn er einen gemeinsamen Testsatz besteht. Er prüft vor
allem die Konventionen, denn hier entstehen die teuersten Fehler (viele RCWA- und FEM-Codes
rechnen mit exp(−i(k·r − ωt)) oder anderen Polarisationsbasen):

- **Ebene Welle durch eine homogene Platte:** Fresnel nach Convention A und Airy-Formel; prüft
  Zeitvorzeichen, s/p-Basis, Normierung der Leistung (Raytatouille ADR 0021: leistungsnormierte
  P) und die Phase der Durchgänge.
- **Energieerhaltung** über alle Ports bei verlustfreier Struktur, **Reziprozität** bei
  reziproker Struktur.
- **Verfahrensspezifisch, unabhängig:** RCWA dünnes Gitter gegen skalare Beugungseffizienz und
  Gittergleichung (ADR 0025); BEM und FEM Kugel gegen Mie; BPM im homogenen Medium gegen das
  Winkelspektrum, im Parabolprofil gegen den analytischen Gaußmodus.
- **Determinismus:** gleiche Eingabe → gleiche Ausgabe; wo ein Löser auf der GPU nicht bitgleich
  ist, eine dokumentierte Toleranz statt Bitgleichheit.

Die Konventionsanpassung (Vorzeichen, Basiswechsel, Einheiten) geschieht nur im Adapter, nie im
Kern.

### 6.6 BWPM als erstes fusionsfähiges Modul

Grundlage: M. Wende, J. Drozella, A. M. Herkommer, „Fast bidirectional vector wave propagation
method …“, Opt. Express 31(18), 28874 (2023), doi:10.1364/OE.497244, und die dort zugrunde gelegte
FPWPM (ebd. Ref. 7). Das Verfahren propagiert die Querkomponenten des elektrischen Felds
schrittweise längs z mit FFT, zerlegt die Brechzahlverteilung je Schritt in M diskrete Werte
(Gewichtsfunktion I_r,s) und wendet in der Frequenzdarstellung Transmissions- und
Reflexionsmatrizen T, R an; vor- und rücklaufende Felder werden in Iterationen hin und zurück
überlagert. Der Aufwand je Schritt ist O(M² N log N); das Paper speichert die reflektierten Felder
im ganzen Volumen und rechnet deshalb in 2D.

Warum die BWPM in der Fusion besonders gut passt:

- Ihre Ports sind **Gitterfelder in Ebenen senkrecht zu z, in beiden Richtungen** – genau die
  Darstellung G von `rtt-field`; die Kopplung braucht nur Umtasten und das Abspalten der
  Trägerphase, keine neue Darstellung.
- Sie schließt die Lücke zwischen Field Tracing (glatte Felder, große Systeme) und den rigorosen
  Lösern (kleine Volumina): 3D-gedruckte Mikrooptik, Faserfacetten, Rückreflexe.
- Sie kann von Beginn an mit den Konventionen von Raytatouille gebaut werden (Zeitfaktor, Fresnel
  Convention A, s/p-Basis, leistungsnormierte Ports), statt sie nachträglich anzupassen.

Vorgehen in eigenem Repository (Vorschlag `PhotonPower/BWPM`, MIT), mit den Regeln von
Fieldtracing:

| Stufe | Inhalt | Referenz und Abnahme |
| --- | --- | --- |
| B0 | Spezifikation: Formeln der FPWPM und BWPM mit Herleitung im Volltext geprüft, Konventionsabgleich mit Raytatouille (insbesondere Vorzeichen der Reflexionsmatrix gegen Convention A) | Dokument, Freigabe |
| B1 | Python-Prototyp der FPWPM in 2D, aufbauend auf den Fieldtracing-Bausteinen (Winkelspektrum, Fresnel-Operator) | homogenes Medium = Winkelspektrum exakt; Grenzfläche = Fresnel; Konvergenz in dz |
| B2 | BWPM in 2D, Iterationen | Platte schräg beleuchtet gegen Airy-Formel, TE und TM, Energie R + T = 1; Konvergenz über die Iterationen wie im Paper (Abschnitt 2.1.3) |
| B3 | 3D mit Speicherung der Reflexe nur an Grenzflächen (im Paper als Ausblick genannt) | 3D-Platte und Linse; Speicher und Laufzeit gemessen, Ziel: Rechner des Maintainers (1 TB, RTX 6000 Ada) |
| B4 | C++/CUDA-Kern mit Python-Anbindung, Adapter nach 6.4, Konformitätstests nach 6.5 | Sollwerte aus B1–B3 bitgleich bzw. in Toleranz |
| B5 | Fusionspilot K4 (Abschnitt 8) | Werte aus dem Paper als Plausibilität, Platte und Fresnel als Referenz |

## 7. Pfadplanung

### 7.1 Fusionspfade

Raytatouille beschreibt einen Pfad als explizite Folge von Ereignissen an Flächen. Ein
**Fusionspfad** verallgemeinert das: eine Folge von Abschnitten, jeder in genau einem Gebiet, und
Übergängen an Ports. Ein Abschnitt in einem Löser-Gebiet trägt, analog zur Beugungsordnung am
Ereignis (ADR 0025), den gewählten **Kanal**: Beugungsordnung, Reflexion oder Transmission,
Streurichtung bzw. Ausgangsport. Damit bleiben Pfade Daten, Strahlteilung durch einen Löser
erzeugt mehrere Pfade, und der Tracer braucht keinen Sonderfall.

### 7.2 Planer in vier Stufen

Der Planer ist eine reine Funktion vom Modell und den Anforderungen (Ausgaben, Toleranz ε,
Ressourcen) auf einen **Plan als Datenobjekt**: welche Abschnitte, welche Darstellung und
Abtastung je Port, welcher Operator je Übergang, welches Fehlerbudget, mit Begründung. Der
Nutzer kann jeden Punkt überschreiben (wie `plan_propagation` in Fieldtracing).

1. **Topologie.** Alle Fusionspfade von den Quellen zu den Ausgaben aufzählen: Folge der Gebiete,
   Kanäle je Löser (Ordnungen, Reflexion/Transmission), Mehrfachdurchgänge bis zu einer Tiefe;
   Pfade mit geschätzter Leistung unter einer Schwelle werden mit Diagnose verworfen (wie das
   Ghost-Ranking, ADR 0027). Ergebnis: ein Baum bzw. DAG gemeinsamer Präfixe, damit gemeinsame
   Abschnitte nur einmal gerechnet werden.
2. **Darstellungswahl.** Je Abschnitt und Port die Darstellung, je Übergang den Operator. Formal
   ein kürzester Weg im Graphen der Zustände (Abschnitt, Darstellung): Kanten sind Operatoren mit
   Kosten aus `estimate` und Fehlerbeiträgen, unzulässige Kanten fallen über die
   Gültigkeitskennzahlen weg. Kennzahlen: lokale Fresnel-Zahlen und Kantenverhältnis (Kapitel 7
   des Papers), Strukturgröße zu λ (Strahl- vs. Wellenoptik), Größenparameter x = 2πa/λ
   (Partikel: Rayleigh, BEM, FEM, Mie), Winkelbreite des Spektrums gegen die Paraxialität des BPM,
   Zahl der Ordnungen, Kohärenzlänge zu Wegunterschieden.
3. **Abtastung und Budget.** Gitter, Strahlzahlen und Tabellendichte je Port aus dem
   Fehlerbudget; das Budget wird nach dem Prinzip des Papers verteilt (der größte Term dominiert,
   die übrigen deutlich darunter). Speichergrenze als harte Nebenbedingung (Beispiel 7.5:
   ≈ 1,3 GB je vektoriellem Feld).
4. **Ausführungsplan.** DAG aus Löser-Aufrufen und Umwandlungen, parallel über Wellenlängen,
   Feldpunkte und Polarisationszustände, mit Cache je Inhalts-Hash (Struktur, Löser-Version, λ,
   Einfall). Lineare Löser mit `operator` werden einmal tabelliert und dann angewandt.

### 7.3 Kohärenz und Überlagerung

Regeln, die der Plan festhält und die Ausgabe dokumentiert:

- Amplituden addieren sich innerhalb einer **Kohärenzgruppe**: gleiche Wellenlänge, gleicher
  Quellpunkt bzw. gleiche Mode, gleicher Polarisationszustand der Quelle, Wegunterschiede unter
  der Kohärenzlänge. Über Gruppen hinweg addieren sich Intensitäten bzw. Stokes-Vektoren.
- Zwei Pfade, die in denselben Ausgangsport münden, werden nur kohärent überlagert, wenn beide in
  einer kohärenten Darstellung ankommen; sonst warnt der Planer.
- Emitter-Ensembles aus BioOptoMD (Q) strahlen untereinander inkohärent: Je Emitter wird kohärent
  gerechnet, über Emitter und Orientierungen werden Intensitäten gemittelt. Streuung an rauen
  Flächen (SpeckleBem) liefert kohärente Speckle-Felder je Realisierung; erst die Mittelung über
  Realisierungen ergibt eine radiometrische Größe (P, BSDF) für M9.

### 7.4 Beziehung zu M9

M9 bringt BVH, Quellen, Detektoren, Strahlteilung und Streumodelle. Die Fusionspfade sollten
dieselbe Pfadbeschreibung nutzen: Ein Löser-Gebiet ist aus Sicht des nicht-sequenziellen Tracers
ein Körper mit Ports, an dem Strahlen enden (R → Einfall) und neue Strahlfamilien entstehen
(Ausgang → R). Ein Streumodell in M9 (Lambert, Gauß, ABg) und eine tabellierte Mueller- oder
BSDF-Antwort eines Lösers sind dann zwei Formen derselben Interaktion. Abstimmung über den
M9-ADR (Entscheidung D6).

## 8. Pilotfälle

Jeder Fall prüft eine Kopplungsart gegen eine Referenz, die nicht aus den gekoppelten Teilen
stammt.

| Nr. | Fall | Module | Kopplung | Referenz | Phase |
| --- | --- | --- | --- | --- | --- |
| K1 | Gitterspektrometer mit Effizienzen je Ordnung und Polarisation | Raytatouille + FourierGrating | Tabelle (M) als Interaktion im sequenziellen Tracer | Gittergleichung, dünnes Gitter; Energie über die Ordnungen | P2/P3 |
| K2 | Metafläche bzw. DOE als Linse, Tabelle über k∥ | `rtt-field` + FourierGrating | Tabelle mit Phase im Feldpfad (H) | Fokus gegen Richards-Wolf bei idealer Phase | P2/P3 |
| K3 | Partikel im Fokus eines Objektivs, Rückweg durch das Objektiv | Raytatouille + `rtt-field` + CliffordBem | R → H → S (`BeamField`) → Fernfeld → S → H/R | Mie (CliffordBem dagegen validiert); Fokusfeld gegen Richards-Wolf und `BeamField.focused` | P2/P3 |
| K4 | Bildleiter mit 3D-gedruckter Index-Matching-Kappe: Objektiv → Fokus → Facette mit Rückreflex → Objektiv → Lochblende | Raytatouille + `rtt-field` + BWPM | R → H → G → BWPM → G (rücklaufend) → H → R | Platte und Fresnel exakt; Werte aus Wende et al. 2023 (Fokusversatz, Rückreflex, Strehl) als Plausibilität | P5 (nach B4) |
| K5 | Biosensor: Rezeptorschicht aus BioOptoMD, Fluorophore als Emitter, Abbildung durch ein Objektiv | BioOptoMD + CliffordBem bzw. hp-FEM + `rtt-field` + Raytatouille | N → Schicht, Q → Emitter, Fernfeld → S → H → R | Emitterrate gegen Reihenlösung (in CliffordBem validiert); ohne Schicht = freier Dipol | P5 |
| K6 | Raue Fläche als Streuer in einem Abbildungssystem | Raytatouille (M9) + SpeckleBem | S → SpeckleBem → Fernfeld/BSDF → P bzw. R | Fresnel für glatte Fläche (in SpeckleBem validiert), Energiebilanz | mit M9 |
| K7 | Brute-Force-Vergleich eines kleinen Gesamtsystems (2D oder kleine Apertur) | gesamte Kette gegen hp-FEM | ganze Kette | eine rigorose hp-FEM-Rechnung über alles | P4 |
| K8 | Biosensor über der Zeit: Einfang von Bakterien, Rezeptorzustände, Fluoreszenz durch gedruckte SIL mit Spiegelrand und Relay-Optik auf die Kamera | BioOptoMD + hp-FEM bzw. CliffordBem + `rtt-field` + Raytatouille + Voxelspaetzle + BindOptiX | Q je Zustand, F, Antwortkarte (M) → BindOptiX → B → Signal | Teilprüfungen gegen Chance–Prock–Silbey, Lukosz–Kunz, Transfer-Matrix (DP1–DP8 im ADR-Entwurf); ohne Analyt = Untergrund allein | P5 (nach K5) |

K7 ist die Abnahme des Planers: Der Plan muss im vorhergesagten Budget liegen, wie es F2 und F3 in
Fieldtracing für die Feldbausteine gezeigt haben.

## 9. Entscheidungen, die vorher fallen müssen

| Nr. | Frage | Vorschlag |
| --- | --- | --- |
| D1 | Scope von Raytatouille: Kopplung rigoroser Löser aufnehmen? | Ja, als Kopplung über Ports; Löser bleiben extern (ADR „Fusion“, ändert den Satz „Nicht geplant“ in `architecture.md`) |
| D2 | Wo liegt der Prototyp? | Neues Repository (Python), Fieldtracing als Abhängigkeit; Regeln wie Fieldtracing |
| D3 | Ausführung der Löser | Unterprozess als Standard, in-process Python für den Prototyp, Dienst später |
| D4 | Wo läuft der Planer? | erst Python im Prototyp; im Kern als reine C++-Funktion in einer neuen Bibliothek der Schicht Workflows (`rtt-fusion`), die Ausführung in C++ mit Adapter-Aufrufen über die Unterprozess-Grenze |
| D5 | Modell und Format | Gebiete und Ports als neuer Abschnitt im Systemmodell mit Verweis auf Struktur-Dateien und Löser-Kennung (F, Schemasprung); Löser-Parameter als Daten, kein Code |
| D6 | Abstimmung mit M9 | Fusion-ADR vor dem M9-ADR; Löser-Gebiete als Körper mit Ports im nicht-sequenziellen Tracer |
| D7 | GPU | Kern weiter ohne GPU (ADR 0013), bis ein eigenes ADR folgt; Löser dürfen GPU nutzen. Erster GPU-Kandidat im Kern ist `rtt-field` (FFT/SFT), nach P1 und mit Messungen |
| D8 | Determinismus mit externen Lösern | Bitgleichheit nur für den Kern; je Löser die in `describe()` erklärte Toleranz, im Ergebnis vermerkt |
| D9 | Lizenzen | Adapter-Pakete getrennt vom MIT-Kern; Lizenz je Löser im Inventar |
| D10 | Reihenfolge der Löser | FourierGrating und CliffordBem zuerst (K1–K3); BWPM parallel neu (B0–B4), dann K4; BioOptoMD über Datenports (K5); hp-FEM als Gegenprobe (K7) und für Nahfelder; SpeckleBem mit M9 (K6) |
| D11 | Ort der BWPM | eigenes Repository mit den Regeln von Fieldtracing, Konventionen von Raytatouille von Anfang an (6.6) |
| D12 | Erweiterungen in den Modulen (Spektrum als Einfall in CliffordBem und SpeckleBem, absolute Leistung) | als Issues in den jeweiligen Repositories, nach Freigabe des Maintainers; dieses Projekt liest dort nur |
| D13 | Datenports für Biosensoren über der Zeit (F, B; Antwortkarte als M) | ADR-Entwurf aus BindOptiX (`PhotonPower/BindOptiX`, `docs/adr-entwurf-datenports.md`); F und B als Datenports wie N und Q, Antwortkarte ohne neue Darstellung |

## 10. Risiken

- **Konventionen:** Zeitvorzeichen, Polarisationsbasis und Normierung unterscheiden sich zwischen
  Codes (SpeckleBem e^(+jωt), FourierGrating mit Schichtung in y und lokalem Phasenbezug). Gegenmittel: Konformitätstest vor jedem Adapter, Anpassung nur im Adapter.
- **Abtastung an Ports:** Periodische Ports (RCWA), Gitter (BPM, Felder) und Strahlnetze passen
  nicht von selbst zusammen; schlechte Wahl verschiebt den Fehler an die Übergänge. Gegenmittel:
  Umwandlungen als geprüfte Operatoren mit Fehlermodell, Abtastung aus dem Budget.
- **Kosten in Schleifen:** Optimierung und Toleranzierung rufen das System hunderte Male auf.
  Gegenmittel: lineare Löser einmal tabellieren (`operator`), Cache nach Inhalts-Hash,
  Ersatzmodelle für Parameterabhängigkeit erst in P5.
- **Leistungsnormierung:** Raytatouilles P ist leistungsnormiert, Felder sind Amplituden.
  Gegenmittel: Umrechnung nach Fieldtracing ADR 0004 an jedem Port, als Test im Konformitätssatz.
- **Umfang:** Fünf Löser, sechs Darstellungen und ein Planer sind ein großes Vorhaben neben M9
  und M10. Gegenmittel: Phasen mit eigener Abnahme, Pilotfälle statt Vollständigkeit, jede
  Darstellung nur mit einem konkreten Nutzer.
- **Wartung vieler Adapter:** Gegenmittel: schmale Schnittstelle (6.4), Konformitätstests in der
  CI der Adapter, Löser-Version im Cache-Schlüssel.

## 11. Nächste Schritte

1. Entscheidungen D1, D2, D10–D12 treffen; danach den ADR-Entwurf „Fusion“ schreiben.
2. BWPM: Repository anlegen, Stufe B0 (Spezifikation mit Konventionsabgleich), dann B1.
3. Prototyp-Repository der Fusion anlegen und K1 als ersten Fall aufsetzen (FourierGrating-Tabelle
   im Tracer, über die Python-APIs von Raytatouille und `fouriergrating`).
4. Issues für D12 in CliffordBem und SpeckleBem vorschlagen.
5. Nach v1.0: P1 nach `docs/rtt-port/`; P2 parallel mit K2 und K3.

## Anhang A: Vorlage für das Inventar weiterer Löser

Ausgefüllt für die vorhandenen Module in Abschnitt 2.1; für neue Module dieselben Felder.

| Feld | Inhalt |
| --- | --- |
| Name, Repository, Version | |
| Verfahren, Gültigkeitsbereich | z. B. RCWA für 1D/2D-periodische Schichten, Zahl der Ordnungen |
| Sprache, Python-Anbindung | |
| Lizenz | |
| Eingaben | Struktur (Format), Einfall (ebene Welle, Feld, Moden) |
| Ausgaben | Ordnungen, Nahfeld, Fernfeld, S-/T-Matrix, Stokes/Mueller |
| Konventionen | Zeitvorzeichen, Fourier, Polarisationsbasis, Einheiten |
| Linear / reziprok / deterministisch | |
| Ressourcen | CPU/GPU, Speicher, typische Laufzeit für einen Pilotfall |
| Validierung | vorhandene Referenzfälle |

### A.1 BindOptiX (ergänzt 2026-11-09)

| Feld | Inhalt |
| --- | --- |
| Name, Repository, Version | BindOptiX, `PhotonPower/BindOptiX`, 0.0 (Spezifikation Z0 als PR #1) |
| Verfahren, Gültigkeitsbereich | kein Lichtlöser; Transport (Langevin-Partikel für Zellen, Konvektion–Diffusion für Moleküle), Oberflächenkinetik, Assay-Abläufe, Signal- und Kameramodell aus optischen Karten; quasistatisch: Optik einmal je Design, Biologie über s bis h |
| Sprache, Python-Anbindung | Python-Referenz (NumPy, SciPy), danach C++/CUDA-Kern mit Python-Anbindung |
| Lizenz | MIT |
| Eingaben | Geometrie (Nurbforge), Q je Rezeptorzustand (BioOptoMD), F (hp-FEM, CliffordBem), Antwortkarte als M (Kette von K5), Strömung, Assay-Protokoll |
| Ausgaben | B (Belegung über der Zeit), Photonenraten je Spot und Pixel, synthetische Kamerabilder, Kennzahlen (Nachweisgrenze, Zeit bis zum Ergebnis, Konfusionsmatrix) |
| Konventionen | keine Wellenoptik intern; F mit e^(−iωt), je mW; Sensorsystem z = 0, Lösung z > 0 (wie BioOptoMD); µm, s, K intern, mm im Austausch; Teilchenzahlen statt mol |
| Linear / reziprok / deterministisch | nichtlinear (Kinetik), stochastisch mit Seed; Ergebnisse mit Standardfehler |
| Ressourcen | Referenz auf der CPU; GPU (RTX 6000 Ada, RTX 3090) für > 10⁶ Partikel |
| Validierung | geplant: BindOptiX `docs/testfaelle.md` (T1–T26), Konformität DP1–DP8 im ADR-Entwurf |
