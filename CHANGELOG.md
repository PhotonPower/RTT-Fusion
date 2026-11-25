# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/), Versionen nach SemVer.

## [Unreleased]

## [0.5.2] - 2026-11-25

Die Kamera hinter der Tubuslinse im FrObSy-Belichtungskopf blieb in Raytatouille Fusion leer: ein
Detektor im Plugin-Prozess kannte den Plan seiner Verbindung nicht (#137).

### Behoben

- **Felder in einen Detektor im Plugin-Prozess** (#137): Raytatouille Fusion rechnet die `rtt.*`-
  Blöcke im Raytracer-Plugin, hinter dem JSON-Kanal eines `BlockProvider`. Der Plan der
  Verbindungen lag nur in der Engine, und so hielt der Detektor dort eine Verbindung auf H für R. Er
  ließ die Felder, die auf der Verbindung entstanden, unbeachtet und gab ein leeres Bild ohne
  Meldung (Kamera hinter der Tubuslinse im FrObSy-Belichtungskopf). Jede Anfrage nennt jetzt die
  geplante Darstellung der Verbindungen am Block (`"plan"`), und der Anbieter gibt sie dem Block als
  `Context.plan`. Ein Block vor einer H-Verbindung setzt damit auch dort den Brechungsindex hinter
  sich (#131). Kommen am Detektor Felder an, ohne dass ein Plan H sagt, rechnet er trotzdem mit
  ihnen und vermerkt das an der Verbindung.

## [0.5.1] - 2026-11-25

Rückmeldungen aus dem Modell der Zweigitter-Interferenzlithografie (#119–#125, #88) und dem
FrObSy-Belichtungskopf (#132, #133): Gitter im Durchlicht als Schritt eines Astes, Bündelindex und
tote Bündel an Gittern, dispersive Port-Medien, Startstrahlen einer gedrehten Blende, Feldbilder mit
einem Fokus neben der Achse, der Index hinter der letzten Fläche vor einem ersetzten Block; ein
knotenweiser Detektor auf „auto“ entscheidet am ersten Bündel, das ihn erreicht; `replace_block` auf
Modellen und alte Graphen mit Gittern im Durchlicht als Modell, sodass kein Skript mehr Version 1
schreiben muss.

### Hinzugefügt

- **`replace_block` auf einem Modell (Version 2)**: Verbindungen und Ast-Schritte ziehen auf die
  Ports des neuen Blocks – gleicher Name, sonst der einzige Port, sonst nach den Ästen: wo das Licht
  eintritt `front` (auch der eine Port eines Spiegels), wo es austritt der Port, in den der Typ von
  dort streut (`back`); nicht eindeutig: ValueError. Der neue Block behält die Raytatouille-Daten
  des alten (`rtt`, `runs` nach den neuen Ports). Dieselbe Regel nutzt der Editor von Raytatouille
  Fusion. Damit braucht kein Skript mehr den Socket-Graphen (Version 1).
- **Alte Graphen mit Gittern im Durchlicht als Modell** (#120, Folge): `model_from_graph` liest eine
  Kette `out` → `in` von Gittern auf `side = transmitted` (wie Skripte sie schreiben) über den
  hinteren Port `back_out`, statt mit „a port has one connection“ abzubrechen; die Äste laufen
  front → back. Ein Gitter auf „reflected“, das Licht an einen anderen Block weitergibt, bleibt
  abgelehnt.
- **Gitter im Durchlicht als Schritt eines Astes** (#120): `fg.grating` und `hpfem.grating` haben
  einen zweiten Port `back` mit dem Ausgang `back_out`, der immer die transmittierten Ordnungen
  liefert (Streutabelle front → front und front → back, wie `fg.stack`). Damit liegt ein
  Transmissionsgitter zwischen zwei Nachbarn im Modell (Version 2), etwa zwei Gitter
  hintereinander. `out` bleibt wie bisher (die Ordnungen der Seite `side`); alte Dokumente sind
  unverändert gültig. Eine Tabelle nur reflektierter Ordnungen am Rückport meldet das klar. Beim
  Ersetzen eines Blocks geht ein Socket, für den mehrere passen, an den der gleichen Seite, sonst an
  den vorderen Port.

### Behoben

- **Workflow-Dokument eines frischen Modells**: `workflow.model_document` setzte erst die Sockets
  eines Knotens und leitete daraus dann die Ports ab – Typen mit deklarierten Ports (die Gitter:
  `front`, `back`) bekamen die Namen `in`/`out`, und das Dokument ließ sich nicht wieder lesen
  („no port 'front'“). Betraf Modelle ohne gespeicherte Ports, z. B. aus Skripten oder aus der
  Umwandlung eines alten Graphen; von CAD geschriebene Dokumente tragen die Ports schon.
- **Knotenweiser Detektor auf „auto“ urteilt über Licht, das ihn erreicht**: er entschied zwischen
  Strahlen und Feldern nur am ersten ankommenden Bündel. Verfehlt das den Detektor (eine
  Beugungsordnung neben der Optik), blieb er bei Strahlen, auch wo die anderen Bündel auf ihn
  fokussieren – im FrObSy-Modell zeigte die Kamera hinter Objektiv und Tubuslinse nur einen
  geometrischen Strich statt der Streifen des Substrats (×100, 33,8 µm). Jetzt entscheidet das erste
  Bündel, dessen Licht ankommt.
- **Index hinter der letzten Fläche vor einem ersetzten Block** (#88, Befund 2): im geteilten Pfad
  (ein Fremdblock wie `fg.grating` an der Stelle der Probe) trugen die Strahlen den Index des
  Mediums, durch das sie an der letzten Fläche ankamen – hinter einem Objektiv dessen letztes Glas
  (1,886) statt Luft. Jeder folgende Weg (`Plane.intersect`, der Block, Teil B) zählte (n − 1)·t
  zu viel; bei einem geneigten Feldpunkt ist das ein Phasenkeil, und das Feld (H) fokussierte
  neben seinen Strahlen: an der Pupillenkamera des 100×-Modells verschmolzen die Paare der
  Beugungsordnungen zu je einem Fleck dazwischen, mit Kopien außerhalb der NA. Jetzt gilt das
  Medium hinter der Fläche (`medium_front`/`medium_back`); die Flecken liegen wie im Strahlmodell
  (0,199 · 0,316 · 0,590 · 0,707 · 1,098 · 1,215 · 1,488 mm), mit 8 und 16 Ringen gleich.
- **Feldbilder mit einem Fokus neben der Achse des Felds** (#88, Befund 1 und die Fälle A und B):
  Rayleigh–Sommerfeld und Debye legen ihr Ergebnis um die Achse der Feldkoordinaten, die im
  Schwerpunkt des Bündels liegt. Ein Bündel, das daneben fokussiert (außeraxialer Feldpunkt, eine
  von mehreren Beugungsordnungen, ein steiles Teilbündel), landete auf einer periodischen Kopie,
  am Rand oder gar nicht im Bild. `field.propagate_to` rechnet den Fleck jetzt dort, wo die
  Strahlen landen: Rayleigh–Sommerfeld um diesen Punkt, Debye um die eigene Achse mit einem
  Fenster, das Detektorfenster und Fleck enthält, bei einem Fleck weiter als einen halben
  Pupillenradius daneben um den Fokus. Zwei kohärente ±1-Ordnungen im Fokus einer Linse geben den
  Fleck in der Mitte mit Streifen der Periode λf/(2·Versatz) und Kontrast 1; ein Teilbündel unter
  37° ist nicht mehr leer. Felder, die ganz mehr als eine Fensterbreite neben dem Fenster landen,
  werden mit einer Protokollzeile übersprungen (bisher z. T. riesige Debye-Gitter). Offen: Debye
  um den Fokus mit einer Pupille neben der Achse liefert in Fieldtracing einige Prozent zu viel
  Spitze (Entwurf `docs/issue-entwurf-fieldtracing-debye-pupille.md`).
- **Gitter: Index der gebeugten Bündel** (#121): `fg.grating` und `hpfem.grating` setzen den Index
  einer Ordnung auf das Medium, in dem sie läuft (transmittiert: Substrat, reflektiert:
  Einfallsmedium), wie `fg.stack`. Bisher behielt das Bündel den Index des Einfallsmediums; der
  nächste Block rechnete den optischen Weg mit n = 1, ein `fg.stack` dahinter meldete den falschen
  Index. Kommen Strahlen in einem anderen Index an, als die Tabelle erwartet, nennt das Protokoll es.
- **Gitter und Stapel mit einem Bündel ohne lebende Strahlen** (#122), z. B. die gesperrte
  Polarisationsgruppe hinter einem idealen Polarisator: `fg.grating`, `hpfem.grating` und
  `fg.stack` überspringen es (bisher Absturz in der Tabelle bzw. eine unnötige Lösung).
- **Dispersives Superstrat oder Substrat im FourierGrating-Adapter** (#123): `{"n": {"sellmeier":
  …}}` (und die anderen Modelle von FourierGrating M22) am Port wertet FourierGrating bei der
  Wellenlänge aus (`material_eps`, FourierGrating 0.20); ein älteres FourierGrating meldet es klar
  statt eines `TypeError`.
- **Startstrahlen bei einer gedrehten ersten Blende** (#125): kann `make_rays` einen Pfad nicht
  zielen (gefaltet, gedreht), laufen die Ersatz-Startstrahlen jetzt längs der Feldrichtung im
  Rahmen der ersten Fläche – ein Laseraustritt längs +x vor einem Umlenkspiegel wird getroffen.
  Erreicht kein Startstrahl das Ende eines Pfades, nennt das Protokoll es (bisher ein Bild der
  Leistung 0 ohne Hinweis).
- **Zwei ersetzte Raytatouille-Blöcke ohne `rtt`-Daten** (#124): ihre Plätze heißen nach ihrem Ort
  (`(removed /root/children/1)`); bisher zweimal `(removed)`, und das System kompilierte nicht.
- **Feldausbreitung mit der Methode „mixed“** (#119): wählt Fieldtracings Planer eine Methode ohne
  Operatorkette (`mixed`: F_min ≤ 1 und F_max ≥ 1/ε, z. B. viele Ringe in einen Fokus), nimmt
  `propagate_to` die der nächsten Auswahlregel (nahe einer Kaustik Debye oder Rayleigh–Sommerfeld,
  sonst SFT oder FFT) und nennt das im Protokoll, statt mit einem `ValueError` abzubrechen;
  `near_focus` folgt derselben Wahl.

## [0.5.0] - 2026-11-24

hp-FEM als zweiter rigoroser Gitterlöser im Workflow (R35): Block `hpfem.grating`, Jones-Tabellen
aus jedem Adapter mit Port-Wellen.

### Hinzugefügt

- **hp-FEM-Gitter im Workflow** (Stufe R35): Block `hpfem.grating` – das rigorose Gitter von
  `fg.grating` mit hp-FEM statt FourierGrating, Sockets eines Spiegels, Beugungsordnungen als
  Kanäle. Das Gitter ist ein hp-FEM-Job (Schema v1, Aufgabe „grating“: Zelle, Netz,
  Materialien, Stapel, Löser), als Datei (`job`) oder im Block (`job_data`, „Check workflow“
  bettet ihn ein; eine Netzdatei bleibt außerhalb und wird gemeldet). Der Block baut um die
  Strahlen einer Wellenlänge eine Jones-Tabelle aus dem hp-FEM-Adapter (eine Richtung: genau
  dort; sonst Knoten alle `sigma_step`, höchstens `max_nodes`, je Knoten zwei FEM-Lösungen) und
  hält sie im Operator-Speicher. Felder und Spektren kommen als Strahlen an. Neu im Adapter:
  `structure_from_job` (Zelle aus dem Job mit hp-FEMs eigenem Leser, Ports auf Ober- und
  Unterkante von Streuern und Stapel) und `front_index`, `back_index`, `period_um` am Handle.
  Tests W99–W101 (Fresnel, hp-FEM direkt und Energie, Tabelle gegen die direkte Lösung). In
  Raytatouille Fusion bietet das hp-FEM-Plugin den Block an.

### Geändert

- `couplings.jones_table` baut Tabellen aus jedem Adapter mit Port-Wellen (FourierGrating,
  hp-FEM), nicht mehr nur aus FourierGrating; der Strahlengang von `fg.grating` steckt in
  `grating.rays_through` und gilt für beide Gitterblöcke (Ergebnisse unverändert).

## [0.4.0] - 2026-11-20

Biosensoren in der Fusion (RTT-Fusion#106: D13, Piloten K5 und K8): Datenports für BioOptoMD und
BindOptiX, Emitter an der Sensoroberfläche durch einen Schichtstapel, ideale Abbildung, Antwortkarte
für BindOptiX und die Zeitschleife.

### Hinzugefügt

- **Zeitschleife `bindoptix.time_loop`** (D13 Schritt 5, Pilot K8): ein Meta-Block über dem
  optischen Graphen. Die Biologie rechnet BindOptiX und liefert die Belegung B
  (`data.occupancy`); der Block setzt den Zustand des Emitters (`emitter.dipole_layer`,
  `state`) auf jeden Zustand von B, wertet den beobachteten Ausgang (Kamerabild) je Zustand
  einmal aus (gecacht, Optik einmal je Design) und setzt die Bilder über der Zeit exakt als
  Σ_s f_s(t)·I_s zusammen (Emitter sind inkohärent). Zustände ohne Fluorophor nennt
  `dark_states`. Ergebnis `TimeSeries` (`image(t, spot)`, `measures` über der Zeit, Kurven in
  `display.plot`). Tests W98: Bild bei f = (½, ½) = halbes inkohärentes Bild aller Emitter
  (1e−9), Optik genau einmal je Zustand, dunkle und fehlende Zustände, Kette mit B aus
  BindOptiX. Ein Szenario-Format, mit dem der Block die Biologie selbst starten kann, ist in
  BindOptiX vorgeschlagen (PhotonPower/BindOptiX#44).
- **Antwortkarte M aus dem Schichtstapel** (Pilot K5, Stufe K5d): Block `emitter.response_map`
  (Stapel wie `emitter.dipole_layer`, Adschicht aus N, Höhen `heights_nm`, Detektorbereiche
  `regions` als NA-Bänder im Substrat oder in der Lösung) gibt BindOptiX' `EmitterResponseTable`
  aus – A_tot, A_rad und A_det je Höhe als diag(a∥, a∥, a⊥) aus einem x- und einem z-Dipol
  (`adapters.bindoptix.response_tensors`, `response_table`). Damit läuft die Kette K5 → K8 ohne
  hp-FEM: `emitter.response_map` → `bindoptix.signal`. `couplings.dipole_layers.band_power`: die
  Leistung in einem NA-Band einer Seite. Überlappende Bereiche einer Seite sind ein Fehler.
  Tests W97: Tensoren = quadratische Form gekippter Dipole (1e−9, auch mit Gold), Spiegel =
  Spiegeldipol (DP4), BindOptiX nimmt die Tabelle an (DP6), Quantenausbeute und Sammelanteil
  = direkte Raten (DP5), Signal der Kette = ausgeschriebenes Modell, Diagnosen.
- **Block `optics.ideal_imaging`** (Pilot K5, Stufe K5c): ein ideales aplanatisches Mikroskop
  (Objektiv mit NA im Medium des Spektrums, Tubuslinse, Vergrößerung M, Bild umgekehrt) für
  Winkelspektren S nach Richards–Wolf: Spektrum auf die Objektebene, je Welle n′ sin θ′ = n sin θ/M,
  Feld in der Meridionalebene gedreht, Amplitude leistungserhaltend, Phase erhalten; Ausgang am
  Bildpunkt mit Punktquellen-Abtastung, so dass der Detektor die Summe der Wellen direkt bildet.
  Wellen jenseits der NA fallen weg (Anteil im Protokoll). Kette Emitter → Mikroskop → Detektor
  mit S. Tests W96: Leistung in der NA erhalten (1e−9), Bild eines Dipols bei NA 0,1 = Airy
  (1,8e−4), Verschiebung um −M·Δ (1e−9), z-Dipol als Ring mit Minimum auf der Achse.
- **Block `emitter.dipole_layer`** (Pilot K5, Stufe K5b): Fluorophore an der Sensoroberfläche als
  Lichtquelle. Eingänge Q (Dipol-Reihe, Zustand `state`, Stichprobe `samples`/`seed`) und N
  (Adschicht über der Oberfläche) sind optional; ohne Q ein Dipol aus `dipole` und `height_nm`.
  Der Stapel (Substrat, `layers` unter der Oberfläche wie ein Goldfilm, Adschicht, Lösung) wird je
  Emitter mit `couplings.dipole_layers` gelöst, Nahfeld und überkritisches Licht eingeschlossen.
  Jeder Emitter gibt ein kontinuierliches Spektrum S mit eigener Kohärenzgruppe, mit der
  Abtastung einer Punktquelle am Fußpunkt, verschoben auf die Port-Ebene (`distance`), so dass
  S → H → R dahinter gehen. Gesammelt wird im Substrat oder in der Lösung (`collect`); das
  Sensorsystem wird nach der Spiegelregel von `in_rtt_frame` (DP8) in die Pose gelegt. Leistung:
  `power_w` (freier Dipol) mal Anteil in der NA. `couplings.dipole_layers.angular_spectrum_grid`
  löst den Stapel einmal je θ für drei Einheitsdipole. Tests W95: Leistung = Anteil aus
  `rates` (2e−4, auch mit Gold), Wellen = gespiegeltes Winkelspektrum (1e−9), Emitter inkohärent
  und seitlich verschoben um e^(−ik∥·Δ), Profil gleich dem Lösungsmittel ändert nichts.
- **Dipol in einem planaren Schichtstapel** (Pilot K5, Stufe K5a, `couplings.dipole_layers`): je
  transversaler Wellenzahl der Stapel mit dem Dipol als Primärfeld (Weyl) gelöst, TE/TM über
  stetige Tangentialfelder, nur abklingende Exponentialfaktoren. `rates` gibt die Gesamtrate
  (Chance–Prock–Silbey), die Leistung nach oben und unten, den Anteil in einer NA und die
  Absorption, normiert auf den freien Dipol; `angular_spectrum` das Winkelspektrum je Raumwinkel
  im Substrat oder oben, überkritisches Licht (SAF, SPCE durch Gold) eingeschlossen. Der Dipol
  darf in einer Schicht sitzen (Adschicht aus n(z)). Tests W94: homogenes Medium = freier Dipol
  (Fernfeld 1e−12), Spiegel = Spiegeldipol (Fehler wie 1/|n|, < 1e−4 bei n = i·10⁵), Nahfeld-Rate
  = abgestrahlte Leistung ohne Verluste (1e−9), Gold absorbiert, Quadratur konvergiert (1e−6),
  Leistung des Winkelspektrums = Fluss ins Glas (1e−6).


- **BindOptiX als Datenknoten** (D13 Schritt 3, RTT-Fusion#106): Block `bindoptix.signal`
  (Eingänge `excitation` F, `response` M, `dipoles` Q, `occupancy` B; Ausgang `rates`) rechnet mit
  dem Signalmodell von BindOptiX (`bindoptix.signal.spot_rates_s`) die mittlere erkannte
  Photonenrate je Rezeptor, Spot, Zeit und Detektorbereich (`adapters.bindoptix.SignalResult`,
  in `display.plot` als Kurven über der Zeit). Q trägt alle Zustände, ein Eingang genügt; Q geht
  unverändert im Sensorsystem hinein (DP8). Detektorbereiche der Antwortkarte werden über den
  Parameter `regions` („label=bereich, …“) gebunden, unbekannte oder doppelte Namen sind
  `exec.data`; eine Karte mit einem Bereich braucht keine Bindung. Eine Feldkarte, die auf I₀
  normiert ist, wird mit `port_irradiance_w_m2_per_mw` und Status L/K/P je mW am Eingang
  umgerechnet. Dazu die Dateiblöcke `data.field_map`, `data.response_map`, `data.occupancy`
  (Umschlag `raytatouille-result`, gelesen von BindOptiX). BindOptiX wird erst bei Bedarf
  importiert (ohne C++-Kern nur NumPy und SciPy); fehlt es, ist das `exec.solver`. Tests W93
  (`tests/test_adapter_bindoptix.py`, ohne BindOptiX übersprungen): Rate gegen das ausgeschriebene
  Signalmodell (1e−12), Mischung nach der Belegung (DP7), Q im Sensorsystem bei einem Feld linear
  in z (DP8), Bindung der Bereiche, Diagnosen, Kurven.
- **Datenblock `biooptomd.layer` rechnet** (D13 Schritt 2, RTT-Fusion#106): liest ein
  Laufverzeichnis von BioOptoMD (Manifest `biooptomd-layer.json`: Profil-CSV mit Δz, n_Lösung und
  λ, je Rezeptorzustand eine `.npy` der Tracker-Reihe (n_replicas, n_frames, 7) in nm) und gibt
  `profile` (N) und `dipoles` (Q) aus. Q (`adapters.biooptomd.DipoleSeries`) bleibt im
  Sensorsystem von BioOptoMD; ins RTT-System rechnet erst ein Lichtknoten um
  (`DipoleSeries.in_rtt_frame`, Spiegelregel DP8). `read_layer`/`write_layer` im Adapter
  (`write_layer` ist die Vorlage für einen Export in BioOptoMD). Die Dateien gehen mit Größe und
  Änderungszeit in den Cache-Schlüssel ein. Neue Diagnose `exec.data` für fehlende oder falsche
  Datendateien. Tests W92: Rundreise bitgleich mit Replikaten und Frames, Spiegelregel
  ausgeschrieben, neu geschriebener Lauf wird neu gelesen, Diagnosen.

## [0.3.3] - 2026-11-19

### Geändert

- **Vektorielles Strahlbild (R32).** `imaging.coherent_image` (die Darstellung R an
  `optics.detector`, `rtt.detector` und `display.beam`) interpoliert je Strahl den komplexen
  Feldvektor a·prt·jones statt |a| mit einer mittleren Phase je Bündel; Felder einer
  Kohärenzgruppe addieren sich als Vektoren, das Bild ist Σ |E|². Licht orthogonaler
  Polarisation interferiert damit nicht mehr (vorher Streifen mit Sichtbarkeit 1), elliptisches
  und zirkulares Licht interferiert mit dem Überlapp seiner Polarisationen – wie in den Feld- (H)
  und Spektrumsbildern (S). Für Licht einer gemeinsamen linearen Polarisation bleibt das Bild
  gleich. `Image.ray_power` nimmt |a·prt·jones|². Tests W91: Zweistrahlinterferenz gegen
  I = cos θ (|e₁|² + |e₂|² + 2 Re(e₁·ē₂ e^(2ikx sin θ))) für parallele, gekreuzte, gleich und
  gegenläufig zirkulare und elliptische Polarisation (1e−9), prt im Bild.
- Die Beispielsysteme mit Raytatouille (Michelson, Linnik, Geisterplatten, K1, Kalkspat) geben
  dieselben Bilder (≤ 9e−16). W52 (Geist eines Singletts knotenweise gegen den ganzen
  Raytatouille-Pfad, der keine Polarisation trägt) weicht um 1,7e−5 ab statt < 1e−9: Die
  Polarisation dreht sich über die Pupille, und zwischen Nachbarstrahlen wird sie linear
  interpoliert (Toleranz 1e−4, `docs/plan-r31-polarisation.md`, Nachtrag R32). W43 (Michelson mit
  gekipptem Spiegel, R) weicht um 1,8e−6 ab: Die Felder der beiden Arme stehen 4 mrad schräg und
  interferieren nicht mehr vollständig (Toleranz 1e−5; mit H wie bisher 1e−9).

## [0.3.2] - 2026-11-19

### Hinzugefügt

- `optics.collimated_source.ellipticity_deg` χ (−45° … 45°, Vorgabe 0): elliptische und
  zirkulare Polarisation als gewöhnlicher Zahlparameter, damit auch in JSON-Dokumenten und im
  Knotendialog von Raytatouille Fusion. `polarisation` ist dann die große Halbachse p̂, das Feld
  cos χ·p̂ + i sin χ·(d × p̂) (tan χ = kleine/große Halbachse; χ > 0: S3 > 0, das Feld dreht von p̂
  nach d × p̂, Zeitkonvention e^(−iωt)). Ein komplexes `polarisation` zusammen mit χ ≠ 0 oder
  |χ| > 45° ist `exec.param`. Tests W90: Stokes-Parameter (1, cos 2χ, 0, sin 2χ) am schrägen
  Strahl, χ = 45° gleich (1, i, 0)/√2, zirkulares Licht durch den Fresnel-Rhombus wird linear.

## [0.3.1] - 2026-11-19

### Korrigiert

- Referenz für Agenten, Abschnitt 0.1: fieldtracing (`photonpower-fieldtracing`), fouriergrating
  und raytatouille gibt es auf PyPI; nur die übrigen Module kommen von GitHub.

### Hinzugefügt

- Befehlsreferenz für Sprachmodelle und Agenten im Paket: `src/rtt_fusion/REFERENCE.md`
  (englisch, kompakt: Orientierung, Denkmodell mit Darstellungen, Planer samt R30f, Umwandlungen,
  Kohärenz, Einheiten, Vorzeichen und Polarisation nach R31; lauffähige Rezepte; Blockkatalog;
  öffentliche API je Modul; Diagnosecodes mit Abhilfe; Stolperfallen und Stellschrauben der
  Laufzeit; Erkundung zur Laufzeit). Zugriff ohne Repository: `python -m rtt_fusion.reference`
  (`--list`, `--section NAME`, `--blocks [TYP …]`, `--blocks --live`) und
  `rtt_fusion.reference(abschnitt)`; `help(rtt_fusion)` nennt beides in der ersten Zeile, das
  README oben (damit auch die Beschreibung auf PyPI).
- Der Blockkatalog darin wird aus der Registry erzeugt (`tools/gen_llm_reference.py`,
  `rtt_fusion._reference.catalogue_markdown`): je Blocktyp Kategorie, benötigte optionale Module,
  nativ gerechnete Darstellungen, Sockets, Ports, Abhängigkeiten und Parameter mit Typ,
  Vorgabe, Einheit und Bereich. `tests/test_llm_reference.py` prüft, dass er aktuell ist, dass
  jeder Blocktyp seine Voraussetzungen nennt, dass jedes mit `python run` markierte Beispiel
  läuft (mit `needs=…` ohne das Modul übersprungen, mit `RTT_FUSION_REQUIRE_MODULES=1` Fehler), dass
  die Datei Paketdaten ist und die Kommandozeile sie ausgibt. Der Wheels-Workflow prüft sie im
  installierten Wheel.

### Geändert

- Planer, Stufe R30f (Entscheidung des Maintainers, `docs/plan-s-spektren.md`): „auto“ wählt S
  ohne Kostenvergleich auf Verbindungen zwischen zwei Blöcken, die S nativ nehmen, wenn das Licht
  am Quellblock ein Spektrum ist – hinter einer Verbindung mit S, hinter `fusion.convert` mit
  `to = S` und durch eine Kette nativer Blöcke (Gitter, Stapel, verdichtete Gruppe, Teilchen)
  dahinter. Der Plan sagt damit S, wo ein Spektrum fließt (bisher R mit der Anzeige „R→S“), und
  plant keine Umwandlung S → R → S. Strahlen werden nicht von selbst zu Spektren; bei gemischten
  Eingängen, vor Blöcken ohne S und wo die Gültigkeit S ablehnt, gilt die Kostenrangfolge wie
  bisher. Ergebnisse ändern sich nicht; Pläne mit Strahlen (Linnik, Beispiele, K3) bleiben gleich.
- Neu: `Requirements.produces` (was ein Block abgibt, was auch immer ankommt), aus dem Attribut
  `produces` der Auswerter (`fusion.convert`: sein `to`), `Executor.produces`,
  `workflow.produces_of`; `workflow.checked` gibt es dem Planer mit.
- Tests mit der beabsichtigten Planänderung angepasst: `fusion.convert(S)` → Detektor,
  Gitter → Stapel → Detektor hinter einem Spektrum und K3 mit S in das Teilchen planen S.
- Komplexe Polarisation der Strahlen, Stufe R31 (Plan `docs/plan-r31-polarisation.md`): die
  Polarisation eines Strahls ê = prt·jones ist ein komplexer Einheitsvektor, elliptische
  Polarisation übersteht jede Umwandlung. Wo das Feld eines Strahls in Skalar und Rest geteilt
  wird (`light.split_field`: |a| = |f|, Phase entlang der großen Halbachse), gibt es keine
  Projektion auf eine reelle Richtung mehr: H → R (`field.field_to_rays`, damit S → R und das
  Umtasten auf gekippte Detektoren), Strahlweg von `fusion.particle`, Gitter und Stapel
  (`grating.split_prt`), geteilte Raytatouille-Pfade. `optics.interface` und `fusion.condensed`
  tragen die Fresnel- bzw. Tabellenmatrix in prt (`light.physical_split`), gemischtes s/p-Licht
  und der Phasensprung der Totalreflexion sind exakt; der skalare Faktor bleibt für s- und
  p-Licht der bisherige. Der Debye-Weg des Teilchens nimmt die Phase der Strahlamplitude mit.
  Lineare Polarisation rechnet wie bisher (bis auf Rundung). `RayBundle` nimmt reelle `prt` und
  `jones` (Port-Daten von vor R31) als complex128; neu `RayBundle.field_vectors()`.
  Abnahme W86–W89 (`tests/test_complex_polarisation.py`, `tests/test_particle_spectra.py`): die
  Rückstreuung der Kugel auf einem um 1,5 mrad gekippten Detektor wie ungekippt (bisher bis
  13 % weniger an den Diagonalen), Fresnel-Rhombus gegen die Randbedingungen, Rundreise
  H → R → H eines elliptischen Gauß-Strahls mit gleichen Stokes-Parametern.

## [0.3.0] - 2026-11-19

### Hinzugefügt

- Spektren ebener Wellen (S) auf Verbindungen, Stufe R30a (`graph/spectrum.py`, Plan
  `docs/plan-s-spektren.md`): Lichttyp `SpectrumBundle` (Kohärenzgruppe, Port-Ebene mit
  Phasenbezug im Ursprung, `ports.PlaneWaves` mit Raumwinkeln dΩ, Abtastung für den Rückweg) und die
  Umwandlungen H → S (punktweise Fourier-Transformation von Fieldtracing, wo die lokale
  Fresnel-Zahl groß ist, sonst FFT auf einem Gitter), G → S (Gitterfeld über
  `ports.to_plane_waves`), S → H (inverse PFT bzw. FFT zurück auf die Stützstellen), R → S und
  S → R (über H). Die Ausführung wandelt an den Ports (`Executor._convert`): eine Verbindung mit S
  trägt Spektren, ein Block, der S nicht nimmt, bekommt ein Feld (S → H) oder Strahlen (S → R).
  Diagnosen als Warnung im Protokoll: `exec.spectrum_cut` (Anteil der Leistung, der evaneszent
  wegfällt), `exec.spectrum_rays` (Bündel nicht groß gegen λ). Abnahme W76–W78
  (`tests/test_spectrum_links.py`): Gauß-Strahl gegen das analytische Winkelspektrum, Rundreisen,
  Leistung, Fokus gegen Richards–Wolf.
- `execute.WITH_SPECTRA` = {R, H, S}; `workflow.available()` gibt es mit fieldtracing zurück.
- Planer: `Requirements.native` (was jeder Block nativ nimmt); „auto“ wählt S nur zwischen zwei
  Blöcken, die beide S nativ nehmen (Entscheidung 3). Bisherige Pläne bleiben unverändert.
  `workflow.native_of` und `workflow.checked(..., evaluators)` liefern es für das Dokument.
- `fusion.convert` rechnet `to = S` (G und P melden weiter `exec.convert_unsupported`).
- Detektor und freie Ausbreitung mit S, Stufe R30b: `optics.detector` und `rtt.detector` (ganze
  Pfade und knotenweise) nehmen Spektren nativ und bilden das Bild aus der Summe der ebenen Wellen
  auf den Pixeln (`spectrum.spectrum_image`) – exakt für diskrete Wellen (Streifen zweier ebener
  Wellen), für kontinuierliche Spektren die Quadratur des Winkelspektrums ohne Fresnel- oder
  Debye-Näherung, solange die Wellen für das Fenster dicht genug liegen (Prüfung des
  Phasenschritts; sonst S → H und Ausbreitung wie bei H, mit Vermerk im Protokoll).
  Gitterspektren auf parallelen Ebenen als zwei Matrixprodukte (`spectrum.grid_sum`). Felder
  und Spektren auf einem Eingang interferieren. Freie Ausbreitung im Spektrum
  `spectrum.propagate` (Plan S4): Phase e^(ik·Δr) je Welle, gekippte Ebenen mit Drehung in den
  neuen Rahmen, evaneszente und nicht vorwärts laufende Wellen fallen mit `exec.spectrum_cut`
  weg; die PFT-Abtastung wandert mit, sodass S → H auf der neuen Ebene gilt. Abnahme W79–W80
  (`tests/test_spectrum_detector.py`).
- Gitter, Schichtstapel und verdichtete Gruppen mit S, Stufe R30c: `fg.grating`, `fg.stack` und
  `fusion.condensed` nehmen Spektren nativ und geben Spektren ab. Jede Welle geht für sich durch
  FourierGrating (derselbe Weg wie der Feldoperator, `apply_many`; wahlweise mit verifizierter
  Interpolation auf dem Gitter eines Gitterspektrums) bzw. durch die Ebene-Wellen-Tabelle der
  Gruppe; Ordnung n liegt bei k∥ + n·K, jede Ordnung ein eigenes Spektrum, kontinuierliche
  Spektren bleiben kontinuierlich (Δ²k unverändert), diskrete bleiben diskret (Entscheidung 1).
  Die Abtastung geht mit (`spectrum.order_support`), sodass S → H hinter dem Block gilt; ein
  Gitterspektrum behält sie auch bei der Ausbreitung zwischen parallelen Ebenen. Mit H am Eingang
  rechnen die Blöcke wie bisher. Diagnose `exec.spectrum_rays` bei S → R für Ordnungen, deren
  Spektren sich überlappen. Abnahme W81–W83 (`tests/test_spectrum_blocks.py`): zwei Gitter als
  Kette gegen FourierGrating mit beiden Gittern in einer Struktur und gegen die H-Kette,
  Σ Ordnungen = R + T, Stapel und verdichtete Platte gegen die Airy-Formel je Welle.
- `SpectrumBundle` als Port-Daten (`portdata`, Blöcke anderer Prozesse): die Abtastung eines
  Gitterspektrums (`GridSupport`) ist jetzt reine Daten (Gitter, Bins, Stützstellen) statt einer
  Funktion.
- Streukörper mit S, Stufe R30d: `fusion.particle` nimmt Spektren nativ (Auswerter R und S) und
  gibt das rückgestreute Fernfeld als Spektrum ab. Das Einfallsspektrum geht ohne inneres R → S
  an CliffordBem: im Spektrum auf dessen Port verschoben (`particle.incident_waves`, exakte Phase
  je Welle, Achse = Normale des Spektrums). Das Fernfeld wird ein kontinuierliches Spektrum im
  gespiegelten Rahmen mit der Abtastung einer Punktquelle im Teilchen, 100 Wellenlängen davor
  (`particle.scattered_spectrum`), sodass S → H dort die Kugelwelle E∞ e^(ikr)/r ist und S → R
  Strahlen vom Teilchen weg liefert. Strahlen am Eingang rechnen wie bisher (Debye je Strahl,
  Strahlen vom Mittelpunkt). Abnahme W84–W85 (`tests/test_particle_spectra.py`): W29–W31 mit
  S-Verbindung, Fernfeld eines Dipols, K3 ohne Strahlen zwischen Fokus und Teilchen gegen
  CliffordBem allein.
- Diagnose `exec.spectrum_support` (`spectrum.SupportError`): S → H oder S → R eines
  kontinuierlichen Spektrums ohne Abtastung (hinter einer verdichteten Gruppe) meldet den Block
  und den Ausweg, statt als „diskrete Wellen“ zu scheitern (offener Punkt aus R30c).
- `spectrum.subset` (Teilmenge der Wellen mit ihrer Abtastung) und `spectrum.towards`.
- Genutzte Darstellung je Verbindung, Stufe R30e: `Executor._convert` trägt in `Context.used` ein,
  was eine Lichtverbindung im Lauf getragen hat (`representation`, z. B. S auf einer mit R
  geplanten Verbindung, wenn ein Spektrum aus `fusion.convert` einen Detektor erreicht), die
  Umwandlungen vom Ausgang der Quelle bis zu dem, womit der Verbraucher rechnet (`chain`, z. B.
  „H→S“ oder „S→H→R“), und einen Satz (`reason`) – überall, wo nicht nur Strahlen über eine
  R-Verbindung gehen, und auf jeder S-Verbindung. `workflow.with_used` gibt die Kette als
  `links[key]["used_chain"]` und je Verbindung eines Modells als `connections_used_reason` weiter;
  die CAD-GUI zeigt R→S an der Verbindung und die Kette im Tooltip.

### Geändert

- Standard von `fusion.convert` ist `to = H` statt S (Plan S6).
- Die Licht-Sockets der Blöcke `optics.*` nehmen R, H und S (vorher R und H); die Ausführung wandelt
  S an ihren Ports.
- `rtt.detector` mit S auf der Eingangsverbindung rechnet im Modus ganzer Pfade nicht mehr still mit
  Strahlen, sondern wandelt die Felder hinter der letzten Fläche in Spektren (offener Punkt aus
  R30a). Die Licht-Sockets importierter Raytatouille-Systeme (`importer`, `topology.unroll`,
  `parts`) nehmen R, H und S; vor den übrigen Raytatouille-Blöcken wandelt die Ausführung S in
  Strahlen. Pläne mit „auto“ bleiben unverändert.
- `field.field_image` aufgeteilt in `field_sums` (Feld je Kohärenzgruppe auf den Pixeln) und
  `image_of`; `field.pixel_centres`.
- Licht-Sockets von `fg.grating`, `fg.stack` und verdichteten Gruppen nehmen R, H und S (vorher R
  und H); ihre Auswerter nehmen S nativ (`fusion.condensed`: R und S, Felder als Strahlen). „auto“
  lässt S zwischen ihnen zu, die Kostenrangfolge behält aber R, wo R zulässig ist: bestehende
  Pläne (Linnik mit Gitterprobe, `auflicht_si_gitter`) bleiben unverändert.
- Ein Spektrum auf einer gegen die Blockebene gekippten Ebene mit Abtastung wird vor Gitter und
  Stapel über H auf die Blockebene gebracht (S → H → S), damit die Abtastung erhalten bleibt;
  parallele Ebenen und diskrete Wellen exakt im Spektrum.
- `couplings.grating_field`: die interpolierte Übertragung als `interpolated_transfer` (auch für
  Spektren), `transfer_waves` und `spectrum_through` für Wellen ohne Gitter.
- Licht-Sockets von `fusion.particle` nehmen R, H und S (vorher R und H). Bestehende Pläne mit
  „auto“ bleiben unverändert (R, keine Umwandlungen); nur die Begründung nennt S als zulässig
  zwischen Teilchen und nativen Blöcken („among RHS“).
- Detektor mit S: bevor ein PFT-Spektrum für das Bild über H geht (die Summe der Wellen würde
  falten), fallen Wellen weg, deren Strahlen mehr als `spectrum.REACH` = 10 halbe Diagonalen des
  Fensters an dessen Mitte vorbeigehen (mit Vermerk und Leistungsanteil). Sie erreichen das
  Fenster nicht; ohne das ließen die fast streifenden Richtungen eines Fernfelds den Planer von
  Fieldtracing Wellenausbreitung über die ganze Ebene wählen (Speicher).
- `field._onto` (Detektor gekippt gegen die Feldebene): Stützstellen, deren Strahlen die Ebene
  nicht vorwärts kreuzen, fallen mit Vermerk weg, statt das Bild mit „rays do not all cross the
  port plane forward“ abzubrechen.

## [0.2.1] - 2026-11-18

### Hinzugefügt

- Konvertierungsblock `fusion.convert` rechnet: R → H (Strahlen zu einem Feld auf der Portebene) und
  H → R (`graph/conversion.py`, mit fieldtracing). Licht in der Zieldarstellung geht unverändert
  durch; G, S und P melden `exec.convert_unsupported` mit dem Hinweis auf eine spätere Stufe.
- `exec.no_evaluator` erklärt, dass die Installation den Block nicht rechnen kann (geplanter Typ
  oder Plugin nicht geladen).

## [0.2.0] - 2026-11-18

### Hinzugefügt

- Gitterblock in C++: Die Schleife über die Strahlen (Tabellenbereich, Rayleigh-Anomalie,
  Interpolation der Jones-Tabelle, PRT-Matrix) läuft in der C++-Beugungstabelle aus Stufe R12, als
  optionale Erweiterung `rtt_fusion._diffraction` (pybind11, `setup.py`). Ohne Compiler wird das
  Paket ohne sie installiert, und der Block nimmt weiter die Python-Schleife. Die Ergebnisse sind
  gleich, auch welche Strahlen verworfen werden (`tests/test_diffraction_cpp.py`). Die Wheels sind
  jetzt plattformabhängig (cibuildwheel für Linux, Windows und macOS).
- Glaskataloge im Dokument: `embedded.embed_catalogs` speichert die AGF-Dateien, die der
  Raytracer in die Quelle schreibt (`rtt_catalogs`), als `rtt_catalog_data` im Workflow (Teil von
  `embed_files`). `embedded.load_catalogs` baut daraus und aus den übrigen Dateien die
  Materialbibliothek. Ein Katalog, der sich nicht lesen lässt und nicht im Dokument steht, ist eine
  Warnung (`graph.file_unreadable`) mit dem Hinweis, wie man ihn nachträgt, statt eines Abbruchs.
- Wheels und PyPI: `.github/workflows/wheels.yml` baut das Wheel (reines Python) und die sdist und
  testet das installierte Wheel unter Linux, Windows und macOS. Ein Tag `v<Version>` hängt sie an
  ein GitHub-Release; *Run workflow* auf dem Tag mit `pypi` lädt sie auf PyPI (Trusted
  Publishing, Projekt `rtt-fusion`, Umgebung `pypi`).

- Was ein Block zur Laufzeit auf einem Link „auto“ benutzt hat, steht jetzt in `Context.used`
  (Link-Schlüssel → Darstellung und Grund), z. B. „H“, wenn ein Detektor nahe einem Fokus mit
  Feldern rechnet, sonst „R“. Entfernte Blöcke schicken es mit ihrer Antwort zurück.
  `workflow.with_used(plan, used, graph, model, reg)` trägt es in den Plan ein (`links[k]["used"]`,
  `["used_reason"]`, `connections_used`), damit die GUI am Link „R→H“ zeigen kann.
- Automatische Wellenoptik nahe einem Fokus (`shared["auto_fields"]`, im CAD das Workflow-Panel):
  Ein Detektor, dessen Link auf „auto“ steht, prüft vorab mit einem Feldpunkt und einer
  Wellenlänge je Pfad, ob die Strahlen sich nahe der Detektorebene kreuzen oder fokussieren. Dafür
  nutzt er die Kaustik-Regel des Fieldtracing-Planers (1/F_target > ε). Trifft das zu, rechnet er
  ab der letzten Fläche mit Feldern, als stünde der Link auf H, und schreibt Grund und Kennzahlen
  ins Log. Sonst bleibt er bei Strahlen. Ein Link auf R oder H behält seine Einstellung. Das gilt
  für ganze Pfade, Split-Pfade (Gitter als Probe) und knotenweise Ausführung.
- Feldpropagation: Die Fehlergrenze ε des Planers von Fieldtracing ist jetzt eine Einstellung des
  Laufs (`shared["field_tolerance"]`, im CAD das Workflow-Panel). Ohne Angabe gilt weiter der
  Standard 0,02. Nahe einem Fokus lässt ε = 0,05 Debye statt Rayleigh–Sommerfeld zu.
- `fg.grating`: Die Auflösung der Jones-Tabelle für den Strahlenweg ist jetzt ein Parameter des
  Blocks (Gruppe „Ray path“): `sigma_step` (Knotenabstand in den Richtungskosinus σ = n·sin θ der
  einfallenden Strahlen, bisher fest 0,01) und `sigma_margin` (wie weit die Tabelle über die
  σ-Spanne hinausreicht; war schon lesbar, aber nicht deklariert).

### Behoben

- Ein System, das Raytatouille nicht übersetzen kann, bricht den Lauf nicht mehr mit einem
  `CompileError`-Traceback ab. Es meldet `exec.compile` mit dem, was fehlt. Fehlen Gläser, weil ein
  Glaskatalog nicht geladen ist, nennt die Meldung den Katalog, die betroffenen Gläser und wo man
  ihn nachträgt.
- Feldpropagation mit Rayleigh–Sommerfeld rechnet in Fieldtracing auf 64² Punkte im Abstand
  λ/(8 sin θ) um die Achse. Ein größeres Bildfenster wurde deshalb stillschweigend auf diesen
  Bereich beschnitten (Pupillenkamera: 218 µm von 11 mm). Deckt das RS-Gitter das Fenster nicht
  ab, nimmt `propagate_to` jetzt Debye, das das Fenster berücksichtigt, und schreibt den geschätzten
  Fehler neben ε ins Log.
- Abtasten eines Feldes auf die Pixel (`field._sample`): Der kubische Spline wird nur noch über
  den Teil des Gitters gelegt, den die Pixel brauchen (`RectBivariateSpline` statt
  `RegularGridInterpolator`). Bei einem Fokusgitter von 2016² dauerte das vorher 26–38 s pro Feld,
  jetzt unter 1 s. Die Genauigkeit liegt in derselben Größenordnung. Beim Auflichtmikroskop
  mit Si-Gitter rechnet die Pupillenkamera (Link H, ε = 0,05) damit in gut 3 Minuten statt gar nicht.
- Entfernte Blöcke (CAD ADR 0010): Ein Detektor mit einem Block eines anderen Plugins auf seinem
  Pfad (R19b, z. B. `fg.grating` als Probe) brach im Prozess des Raytracer-Plugins mit „no
  evaluator for fg.grating“ ab, weil dieser Prozess das Gitter nicht kennt. `register_remote`
  nimmt jetzt die eigenen Evaluatoren der Engine (`local=`); ein Block, dessen Evaluator mit
  `calls(ctx, node)` Blöcke eines anderen Anbieters nennt (`backend_rtt.foreign_blocks` für den
  Detektor), wird damit in der Engine gerechnet, die beide Plugins erreicht.
- `fg.grating` mit Konfiguration, Strahlenweg: Ein Bündel mit nur einer Einfallsrichtung (eine
  ebene Welle, z. B. ein Feldpunkt einer Köhler-Beleuchtung) bekommt seine Jones-Matrizen jetzt
  direkt an seinem σ (`jones_table.point_table`, 4 Lösungen) statt aus einer 4×4×4-Tabelle (128).
  Liegt eine Rayleigh-Anomalie näher als der Tabellenrand, wurde vorher das ganze Bündel verworfen
  („rays outside the table or near a Rayleigh anomaly“). Beim Auflichtmikroskop mit Si-Gitter
  (Λ = 1 µm, λ = 546 nm, σ = 0,447 neben der Anomalie bei 0,454) war das Kamerabild fast dunkel
  und brauchte 6 Minuten; jetzt 3 Sekunden.
- Gebaute Jones-Tabellen liegen im Operator-Speicher statt in `shared`. Im Prozess eines
  Plugins (entfernte Blöcke) ging der Cache nach jeder Anfrage verloren; Bild- und Pupillenast
  bauten dieselbe Tabelle zweimal.
- Die Meldung über verworfene Strahlen nennt die Abhilfe (σ-Rand und σ-Schritt im Dialog).

- `system_from_model` (die Zeichnung des Raytracers nach jeder Bearbeitung im Node-Editor) brach
  mit „branch …: no events for Probe front>front|“ ab, wenn ein Fremdblock (z. B. ein Gitter statt
  des Probenspiegels) für seine Ports keine Ereignisse hält. Ein solcher Schritt beendet jetzt den
  Pfad des Astes wie in `system_from_graph`; bei Raytatouille-Blöcken bleibt es ein Fehler.

- „Check workflow“ warnte bei jedem Block, der einen Raytatouille-Block ersetzt hat (z. B. ein
  `fg.grating` als Probe statt des Spiegels), mit `graph.param_unknown` vor `rtt` und `runs`.
  Diese Daten übernimmt `replace_block` absichtlich; das erzeugte System hält das Element weiter.
  Sie gelten jetzt bei jedem Block als bekannt, andere unbekannte Parameter werden weiter gemeldet.

- Ein Workflow-Dokument rechnet ohne weitere Dateien (W75): `fg.grating` und `fg.stack` nannten
  ihre FourierGrating-Konfiguration bzw. Jones-Tabelle nur über einen Pfad; auf einem anderen
  Rechner (oder ohne Zugriff auf den Ordner) brachen alle Anzeigen hinter dem Block mit „cannot
  read …“ ab. Die Blöcke halten den Inhalt jetzt selbst (`config`, `table_data`;
  `rtt_fusion.graph.embedded`) und lesen eine Datei nur noch, wenn sie keine Kopie haben.
  `workflow.embed_files` übernimmt lesbare Dateien einmalig in die Blöcke (auch in Gruppen) und
  leert die Pfade; unlesbare meldet es als `graph.file_unreadable`, STEP-Dateien von CAD-Schichten
  bleiben außen (`graph.file_external`). Fehlt beides, nennt der Fehler den Weg zur Abhilfe.

- Außeraxiale Bildleistung im Feldpfad lag bis zu 4,7 % über dem physikalischen Maximum, wenn
  eine Blende die Pupille beschneidet (Rest aus #69): Der Debye-Schritt von Fieldtracing füllte
  jede Pupille zur vollen Kreisscheibe auf. Behoben in Fieldtracing #79/#80; der Pin steht jetzt
  auf `ea6da7b`.
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

- Die Löser-Module mit festen Commits stehen jetzt in `requirements-modules.txt` statt im Extra
  `modules` (PyPI nimmt keine Git-Verweise in den Metadaten an): `pip install -r
  requirements-modules.txt` statt `pip install -e ".[modules]"`. Die Lizenzangabe ist jetzt ein
  SPDX-Ausdruck (setuptools ≥ 77).

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

