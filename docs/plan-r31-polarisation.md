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
