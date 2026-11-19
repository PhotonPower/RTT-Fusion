# Plan R31: Komplexe Polarisation der Strahlen

> Stufe R31 in [roadmap.md](roadmap.md). Anlass: offener Punkt „H → R und elliptische
> Polarisation“ aus dem Nachtrag R30d ([plan-s-spektren.md](plan-s-spektren.md)).

## Befund

`field.field_to_rays` (H → R) teilt das Feld u jeder Stützstelle in eine komplexe Amplitude a und
eine **reelle** Polarisationsrichtung e (Phase der größten Komponente heraus, Realteil normiert)
und setzt `prt` so, dass prt·jones = e. Von einem elliptisch polarisierten Feld bleibt nur der
lineare Anteil a·e = (u·e)e. Ebenso teilen `grating.split_prt` (Gitter, Stapel),
`backend_rtt._scalar_along` (geteilte Pfade) und `particle.spectrum_to_rays` (Strahlweg des
Teilchens) das Feld entlang einer reellen Richtung; die Grenzfläche `optics.interface` und die
verdichtete Gruppe (`condense.apply`) nehmen nur die Projektion auf die geometrisch mitgeführte
Polarisation („exakt für s- oder p-Licht“).

Folgen: Die Rückstreuung der Kugel (`fusion.particle`, CliffordBem) ist neben den Achsen
elliptisch; ein gegen die Spektrumsebene um 1,5 mrad gekippter Detektor (H → R → H in
`field._onto`) verliert an den Diagonalen bis 13 %. Der Strahlweg des Teilchens behält die
Leistung, macht die Polarisation aber linear. Gemischte s/p-Polarisation an einer Grenzfläche
(Totalreflexion mit Phasensprung) wird falsch.

## Entscheidung des Maintainers (2026-11-19)

Strahlen tragen eine komplexe Polarisation (Jones-/E-Feld-Vektor), damit elliptische Polarisation
die Umwandlungen Feld → Strahlen übersteht. Reelle Daten bleiben lesbar.

## Entwurf

1. **Datenmodell ohne neues Feld.** `RayBundle.prt` ist seit R18 komplex (N, 3, 3), `jones`
   komplex (3,); die Polarisation eines Strahls ist **ê = prt·jones**, ein komplexer Vektor. Neu
   ist die Regel: wer das Feld f eines Strahls in Skalar und Polarisation aufteilt, nimmt
   f = a·ê mit einem **komplexen Einheitsvektor** ê (ê†ê = 1) und |a| = |f|
   (`light.split_field`). Phase von a: die des Felds entlang der großen Halbachse der Ellipse,
   φ = ½ arg(f·f), Vorzeichen so, dass die größte Komponente von Re ê positiv ist; bei
   (nahezu) zirkularer Polarisation die Phase der größten Komponente. Für lineare Polarisation
   ist ê reell und a wie bisher (bis auf Rundung, < 10⁻¹⁵). Ein eigenes Feld (N, 3) wäre neben
   prt·jones redundant und müsste an jeder Wechselwirkung mit prt mitgeführt werden; prt bleibt
   nötig, weil im Raytatouille-Backend ein aufgezeichneter Pfad zwei Eingangspolarisationen
   (`jones`) bedient. Speicher und Rechenzeit: keine neuen Felder, prt war schon complex128.
2. **Erzeuger von Strahlen** nehmen ê komplex: `field_to_rays` (prt bildet die Achsen der
   Feldebene (u, v, n) auf (ê, k × ê*, k) ab, unitär), damit S → R und das gekippte
   Umtasten (`_onto`), `particle.spectrum_to_rays` (Strahlweg des Teilchens),
   `grating.split_prt` (Gitter-Strahlwege, `fg.stack`), `backend_rtt._scalar_along` (Pfade, die
   ein Block teilt).
3. **Wechselwirkungen:** `optics.interface` (Fresnel) und `fusion.condensed` tragen die
   physikalische Matrix c_s·ssᵀ + c_p·p_out p_inᵀ in prt statt nur der geometrischen; der
   skalare Faktor ist |P·ê|/|G·ê| mit der Phase der Projektion auf die geometrisch mitgeführte
   Polarisation G·ê – für s- oder p-Licht und senkrechten Einfall genau der bisherige Faktor.
   Gitter (Jones-Matrix je Ordnung) und Raytatouille-PRT stehen schon ganz in prt; Spiegel
   (Householder, skalarer Reflexionsfaktor) bleiben.
4. **Verbraucher:** `rays_to_field` (R → H, damit R → S) und `particle.rays_to_spectrum`
   (Debye) rechnen mit dem ganzen Feld a·prt·jones; das Teilchen nimmt dabei neu auch die Phase
   von a mit (bisher nur |a|; gleich für die bisherigen Fälle mit reellem a). Das skalare
   R-Bild (`imaging.coherent_image`, R15) nutzt |a|, also jetzt die volle Feldstärke einer
   elliptischen Polarisation statt ihres linearen Anteils; vektoriell wird es nicht (offen).
5. **Rückwärtskompatibel:** `RayBundle` wandelt reelle `prt`/`jones` (ältere Port-Daten,
   Aufrufer mit float-Feldern) in complex128 mit Imaginärteil 0; das Format der Port-Daten
   bleibt Version 1.
6. **Quelle:** `polarisation` der `optics.collimated_source` nimmt komplexe Werte (in Python;
   die Pose dreht Real- und Imaginärteil seit R26).

## Abnahme

| Nr. | Prüfung | Referenz |
| --- | --- | --- |
| W86 | R30d-Fall: Rückstreuung der Kugel (K3 mit S) auf einen um 1,5 mrad gekippten Detektor (H → R → H) | ungekippter Detektor bzw. \|E∞\|²/r² von CliffordBem allein, innerhalb des Abtastfehlers (bisher −13 % an den Diagonalen); dazu Dipol mit p = (1; 0,3i; 0) |
| W87 | Fresnel-Rhombus: zirkulare Polarisation aus 45°-linearer durch zwei Totalreflexionen (n = 1,51, 54,6°) | analytische Jones-Rechnung (Phasensprung der Totalreflexion, Born/Wolf 1.5.4) |
| W88 | Rundreise H → R → H eines elliptisch polarisierten Gauß-Strahls | Stokes-Parameter je Stützstelle und integriert unverändert |
| W89 | Lineare Polarisation unverändert; reelle prt aus älteren Port-Daten | Ergebnisse vor R31 (1e−12); komplex mit Imaginärteil 0 |

## Nachtrag R31 (Umsetzung, 2026-11-19)

**Gebaut.** Wie im Entwurf: `light.split_field` (f = a·ê, |a| = |f|, ê komplex mit ê†ê = 1;
Phase entlang der großen Halbachse, Vorzeichen an der größten Komponente von f wie bisher, bei
zirkularer Polarisation die Phase der größten Komponente), `light.split_along` und
`light.physical_split` (Grenzfläche, verdichtete Gruppe), `RayBundle.field_vectors()`
(a·prt·jones) und die Umwandlung reeller `prt`/`jones` in `RayBundle.__post_init__`. Benutzt in
`field.field_to_rays`, `particle.spectrum_to_rays`, `particle.rays_to_spectrum`/`ray_powers`,
`grating.split_prt` (auch `fg.stack`), `backend_rtt._scalar_along`, `optics._fresnel`,
`condense.apply`. Abnahme W86–W89 (`tests/test_complex_polarisation.py`, W86 (b) in
`tests/test_particle_spectra.py`):

| Prüfung | Referenz | Toleranz | vor R31 | gemessen |
|---|---|---|---|---|
| W86 (b) Kugel (K3 mit S), Detektor 1 mm davor um 1,5 mrad gekippt: Bild | \|E∞\|²/r² von CliffordBem allein | 1e−2 der Spitze | 7,0e−2 (Diagonale 0,867) | 4,4e−3 (Pixel ≥ 0,995) |
| W86 (b) dasselbe gegen den ungekippten Detektor | Bild ungekippt | 1e−2 | 7,0e−2 | 7,8e−3 |
| W86 (a) Dipol p = (1; 0,3i; 0), ungekippt / gekippt | Kugelwelle | 1e−2 | 3,9e−3 / 8,7e−2 | 3,9e−3 / 5,0e−3 (linear p: 4,9e−3) |
| W87 Fresnel-Rhombus, Feld je Strahl | Maxwell-Randbedingungen | 1e−12 | 0,78 (Leistung 2 %, S3 = 0) | 5e−16 (S3/S0 = 1) |
| W88 elliptischer Gauß-Strahl H → R → H, Stokes je Stützstelle | vorher | 1e−12 | S3 → 0, Leistung 0,876 | < 1e−12 |
| W89 lineare Polarisation: Aufteilung, Grenzflächenfaktor | Formeln von R18–R30 | 1e−15 / 1e−14 | – | erfüllt |
| W89 bestehende Tests (Ausgänge aller Auswertungen und Umwandlungen, 236 Tests der betroffenen Dateien) | Lauf auf `main` (6cc2f48) | 1e−12 | – | linear polarisierte Strahlen ≤ 2,7e−11 (Restelliptizität ~1e−7), Bilder rein linearer Fälle ≤ 1,3e−15; geändert nur mit elliptischem Licht (Rückstreuung der Kugel, S → R hinter Gittern und Fokus, Tabelle zwischen den Gitterpunkten: 1,5e−5) |

**Entscheidungen und Abweichungen.**

1. **Kein eigenes Feld für ê.** prt (N, 3, 3) war schon komplex; die Polarisation ist prt·jones.
   Ein Feld (N, 3) neben prt wäre redundant und an jeder Wechselwirkung doppelt zu führen.
2. **Phase des Skalars.** Das skalare R-Bild nimmt |a| und eine Phase je Bündel; für lineare
   Polarisation bleibt a bis auf Rundung gleich (auch das Vorzeichen bei zwei gleich großen
   Komponenten), für elliptische folgt die Phase der großen Halbachse (stetig über das Bündel,
   außer bei zirkularer Polarisation). An Grenzflächen und in verdichteten Gruppen bleibt der
   Faktor auf die geometrisch mitgeführte Polarisation bezogen, damit Teilsummen mehrerer
   Bündel (W50) ihre Vorzeichen behalten.
3. **Das R-Bild bleibt skalar** (R15). Es nutzt jetzt die volle Feldstärke einer elliptischen
   Polarisation, addiert die Bündel einer Kohärenzgruppe aber weiter ohne ihre Polarisation;
   vektoriell rechnen H und S.
4. **Debye-Weg des Teilchens mit der Phase der Amplitude.** Bisher ging nur |a| ein; für die
   Strahlen einer Linse (reelle Amplituden) gleich, nach einem Spiegel oder H → R jetzt richtig.
5. **Bestehende Tests unverändert** (nur ein Kommentar in `test_particle_spectra.py` zur
   gekippten Ebene und die neue Prüfung W86 (b) im Fixture `k3`). Geprüft zusätzlich mit einer
   Aufzeichnung aller Ausgänge von `Executor.evaluate`, `field_to_rays`, `rays_to_field`,
   `field_image`, `coherent_image` und der Umwandlungen des Teilchens in 236 Tests vor und nach
   R31 (Werte in der Tabelle). Abweichungen gibt es nur, wo Licht elliptisch ist; dort ist das
   neue Ergebnis das richtige (die Projektion fehlt).
6. **Laufzeit:** keine neuen Felder; `split_field` kostet so viel wie die alte Aufteilung
   (10⁵ Strahlen: 53 statt 55 ms), Grenzfläche und verdichtete Gruppe ein Produkt P·prt statt
   G·prt. Gesamtlaufzeiten der Tests sind auf dem belasteten Rechner nicht aussagekräftig
   (Teilchen und Spektren: 2904 s vorher, 2897 s nachher).

**Offen.**

- Vektorielles R-Bild (Bündel einer Kohärenzgruppe mit verschiedener Polarisation).
- Komplexe `polarisation` der Quelle im JSON-Format des Modells (bisher nur aus Python).
