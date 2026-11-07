# Konventionen

Verbindlich für alle Ports und Adapter; sie sind die von Raytatouille (`docs/architecture.md`,
„Konventionen“, Commit `6078e55`) und gleich denen von Fieldtracing und BWPM.

## Einheiten, Koordinaten, Felder

- Längen in mm, Vakuumwellenlänge in µm, Wellenzahlen in 1/mm, Winkel in rad.
- Rechtshändig, Hauptrichtung +z. Felder ∝ exp(i(k·r − ωt)) (Zeitfaktor e^(−iωt)).
- k_z = √(ñ²k₀² − |k⊥|²) mit Im k_z ≥ 0 (Re k_z ≥ 0, falls Im k_z = 0); Regel für streifenden
  Einfall wie BWPM: |ñ²k₀² − |k⊥|²| ≤ 10⁻¹²·|ñ²k₀²| gilt als k_z = 0.
- Fresnel Convention A (Byrnes, arXiv:1603.02720v5, Gl. (6)): s = k_ein × N/|…| mit N ins
  zweite Medium, p = k̂ × s je Welle; TE = s, TM = p; bei senkrechtem Einfall r_p = −r_s.

## Ports

Ein Port ist eine Ebene z = const mit

- Richtung +1 (Welle läuft in +z) oder −1,
- Medium ñ am Port (homogen in der Ebene),
- Vakuumwellenlänge,
- einer Darstellung des Felds:
  - **G Gitterfeld:** Querfeld (E_x, E_y) auf dem zentrierten Gitter
    x_m = (m − n_x//2)·dx, y_n = (n − n_y//2)·dy, Feld periodisch über das Fenster;
  - **S Ebene-Wellen-Spektrum:** Liste ebener Wellen (k_x, k_y, E_x, E_y), jede mit Phase 1 bei
    (x, y) = (0, 0) in der Port-Ebene; optional Ordnungen (n, m) bei periodischen Strukturen.

E_z folgt aus der Transversalität: E_z = −Richtung·(k_x E_x + k_y E_y)/k_z.

**Kontinuierliche Spektren** (Strahl, Fernfeld eines Streukörpers) werden auf Quadraturrichtungen
diskretisiert; jede Welle trägt dann A(k̂)·dΩ und die Port-Daten das Gewicht dΩ (`solid_angle`).
Ihre Leistung ist (2π)²/k Σ |A|² dΩ (`ports.spectrum_power`).

**Leistung:** Poynting-Fluss durch den Port in z, bis auf den Faktor 1/(2ωμ₀):
je ebener Welle Re(k_z)·|E|² (|E|² mit E_z), evaneszente Wellen tragen nichts. Für G summiert
über das FFT-Spektrum (Parseval) und multipliziert mit dx·dy/N; für S ist es die Leistungsdichte
je Welle (Fläche 1). Umrechnung auf leistungsnormierte Strahlen von Raytatouille an der
Port-Grenze nach Fieldtracing ADR 0004 (folgt mit dem Strahlen-Port).

## Umwandlung G ↔ S

- G → S: FFT auf dem zentrierten Gitter, auf (0, 0) bezogen: Ẽ(k) = (1/N) Σ E(r) e^(−i k·r);
  Frequenzen 2π·numpy.fft.fftfreq (Rundung wie BWPM, wichtig am streifenden Einfall).
- S → G: E(r) = Σ Ẽ e^(i k·r) für Wellen auf dem Gitter; Wellen zwischen den Gitterfrequenzen
  sind im Fenster nicht periodisch und werden abgelehnt.

## Konventionen der Löser (Umrechnung nur im Adapter)

| Löser | Unterschied | Umrechnung |
| --- | --- | --- |
| BWPM | keiner | – |
| hp-FEM (Gitter) | SI-Meter; Rahmen x = Periode, y = Normale (Licht von +y), z = Rillen; Jones-Matrizen leistungsnormiert in ŝ = k̂ × ŷ (ẑ bei senkrechtem Einfall), p̂ = ŝ × k̂; Bezug am linken Zellrand und an den Stapelgrenzen | wie FourierGrating gedreht; kartesische Felder mit hp-FEMs Basis; Leistungsnormierung entfernen; Phasen auf x = 0 und die Port-Ebenen verschieben |
| SpeckleBem | Zeitfaktor e^(+jωt) (passive Medien Im ε ≤ 0), SI-Meter, Phasenbezug der ebenen Welle am Netzursprung, Fernfeld E_s ≈ F e^(−jkr)/r | ε, Einfallsamplitude und Fernfeld komplex konjugieren; F[mm] = F[m]·10³; Amplitude am Körpermittelpunkt |
| CliffordBem | Längen in einer freien Einheit (hier 1 µm), ω = k₀; Einfall als Winkelspektrum mit Brennpunkt; Fernfeld E∞ mit E_s ≈ E∞ e^(ikr)/r | Einheiten; Brennpunkt = Port-Ursprung; Fernfeld → Port-Spektrum über Weyl, transversal projiziert |
| FourierGrating | Schichtung in y, Licht von oben (−y); x′ = x, y′ = z, z′ = −y (Drehung um x); Amplituden je Ordnung als (E·s, E·p) mit s = (γ, 0, −α)/|κ|, bei κ = 0 aber s = (0, 0, 1); Einfall A = cos δ·s_ein + sin δ·p_ein mit s_ein = (sin ψ, 0, −cos ψ) auch bei φ = 0; Phasenbezug lokal an der oberen bzw. unteren Grenzfläche; Längen in µm, Winkel in Grad | Felder kartesisch rekonstruieren, drehen; Einfall als Überlagerung zweier Läufe (δ = 0°, 90°) mit komplexen Gewichten |
