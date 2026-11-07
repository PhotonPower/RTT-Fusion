# ADR 00xx (Entwurf): Tabellierte Beugungsinteraktion (Jones-Matrix je Ordnung)

> **Entwurf aus RTT-Fusion (Stufe R6) für Raytatouille**, Präzisierung von Punkt 4 des
> freigegebenen Fusion-ADR-Entwurfs (`docs/adr-entwurf-fusion.md`, E1) und Erweiterung von
> ADR 0025. Nummer: direkt nach dem Fusion-ADR. Stand Raytatouille: Commit `708b6dd` (0.9.0).
> Prototyp: `rtt_fusion.couplings.jones_table`, Tests `tests/test_jones_table.py` (J1–J5 in
> `docs/testfaelle.md`).

**Status:** vorgeschlagen (für Raytatouille; in RTT-Fusion freigegeben am 2026-11-07, Entscheidung E2)

## Kontext

ADR 0025 gibt einer Gitterfläche eine Effizienz η je Ordnung: eine Zahl, unabhängig von
Wellenlänge, Einfallswinkel und Polarisation, die nur in `weight` eingeht und nach der eigenen
Interaktion der Fläche (Fresnel der Ordnung 0) wirkt. Pilot K1 (RTT-Fusion) zeigt an einem
Lamellengitter (Λ = 0,7 µm, n = 2,1/Luft, λ = 0,6 µm), dass das für rigoros gerechnete Gitter
nicht reicht: Effizienzen ändern sich schon bei 4,6° Einfall um bis zu 12 %, s und p unterscheiden
sich um bis zu 70 %, und die Phase je Ordnung fehlt ganz.

Ein rigoroser Löser (FourierGrating) liefert je Ordnung das vollständige Feld. Zu entscheiden
ist, in welcher Form Raytatouille das als Tabelle übernimmt: Basis, Normierung, Stützstellen,
Interpolation, Phasenbezug, Verhalten an Rayleigh-Anomalien und außerhalb der Tabelle.

## Entscheidung

### 1. Feld und Datei

- Neues optionales Flächenfeld `diffraction_table` mit dem Pfad einer Tabellendatei, nur an
  Flächen mit genau einer Phasenschicht `linear_grating` (der Gittervektor legt die lokale x-Achse
  der Tabelle fest). Nicht zusammen mit `diffraction_efficiency` (Diagnosecode
  `surface.table_conflict`).
- Datei im Format `raytatouille-result` (ADR 0023), Typ `DiffractionTable`: Periode Λ (µm), reale
  Brechzahlen des Einfalls- und Austrittsmediums, Achsen, Ordnungen je Seite (durchgelassen,
  reflektiert), die Jones-Matrizen und die Herkunft (Löser, Version, Inhalts-Hash, Konvergenz).
- Erzeugt wird die Datei außerhalb von Raytatouille durch einen Adapter (Fusion-ADR Punkt 5).

### 2. Basis: transversale kartesische Komponenten im lokalen Flächensystem

- Lokales System der Fläche: x = Gittervektor, y = Rillenrichtung, z = Normale in
  Lichtrichtung.
- Eine einfallende Welle mit tangentialen Richtungskosinus σ = (σ_x, σ_y) = k⊥/k₀ (also n·sin θ)
  und transversalem Feld (E_x, E_y) erzeugt in Ordnung m die Welle σ_m = σ + m λ/Λ x̂ mit dem
  transversalen Feld J_m (E_x, E_y); E_z folgt aus der Transversalität.
- **Begründung:** Die s/p-Basis ist bei senkrechtem Einfall nicht definiert; Raytatouille wählt
  s dort über eine Achsenregel (`rtt/polar/prt.hpp`). Für ein Gitter ist J dort anisotrop (TE ≠ TM,
  J1/J2: |J_xx| und |J_yy| unterscheiden sich deutlich), das Ergebnis hinge also von der
  Achsenregel ab. Die kartesische Basis hat keine solche Stelle.

### 3. Normierung und PRT-Matrix

- J_m ist **leistungsnormiert**: Feld-Jones-Matrix mal √(Re k_z,m / Re k_z,ein). Für ein
  Einheitsfeld (|E| = 1, mit E_z) ist |E_m|² der Leistungsanteil der Ordnung – die Bedeutung von P
  in ADR 0021.
- PRT-Matrix der Ordnung (k_ein, k_m Einheitsvektoren):

      P_m = L_m J_m Π (I − k_ein k_einᵀ) + k_m k_einᵀ

  mit Π der Projektion auf (x, y) und L_m der Ergänzung (E_x, E_y) → (E_x, E_y, ∓(σ_m·E)/(n cos θ_m))
  (oberes Vorzeichen durchgelassen, unteres reflektiert). Damit P_m k_ein = k_m und P_m E = E_m.
- Die Tabelle **ersetzt** die Interaktion der Fläche (Fresnel, Coating) und η: Grenzflächen
  zählen nicht doppelt. `weight` folgt wie in ADR 0021 aus P.
- Energie: Für ein verlustfreies Gitter ist Σ_m |P_m E|² über durchgelassene und reflektierte
  Ordnungen = 1 für jede Polarisation (J1: 10⁻⁹).

### 4. Phasenbezug

J_m bezieht sich auf den Ursprung der Phasenschicht (φ(0, 0) = 0, ADR 0025 Punkt 1). Ein um x₀
verschobenes Gitter ändert J_m um den Faktor e^(−i m 2π x₀/Λ) (J4: 10⁻¹⁰); an einem Auftreffpunkt
(x, y) kommt die Phase m·φ(x, y) der Phasenschicht wie bisher hinzu. Ein Test im Port prüft das
Vorzeichen gegen den OPL-Beitrag von ADR 0025 (Abnahme 5).

### 5. Stützstellen und Interpolation

- Gleichmäßige Achsen in λ (µm), σ_x, σ_y; Interpolation **kubisch** (Lagrange, 4 Punkte je Achse,
  Tensorprodukt) der komplexen Einträge.
- Gemessene Genauigkeit fern von Anomalien (J5, Lamellengitter, Bereich λ 0,58–0,62 µm,
  σ_x ±0,1): größter Fehler eines Eintrags

  | Schritt σ_x | linear | kubisch |
  | --- | --- | --- |
  | 0,05 | 9,0·10⁻³ | 2,5·10⁻³ |
  | 0,025 | 1,6·10⁻³ | 2,5·10⁻⁴ |
  | 0,0125 | 5,9·10⁻⁴ | 1,4·10⁻⁵ |

  (kubisch ≈ h⁴, linear ≈ h²). Die Datei gibt die Genauigkeit an, mit der sie erzeugt wurde
  (Kontrollpunkte zwischen den Stützstellen).

### 6. Rayleigh-Anomalien

- An |σ + m′λ/Λ x̂| = n (Einfalls- oder Austrittsmedium) wird eine Ordnung m′ streifend; dort hat
  **jede** Ordnung einen Wurzel-Verzweigungspunkt. Kreuzt eine solche Linie die Interpolationszelle,
  ist der Fehler groß und konvergiert nur wie √h (gemessen 0,11 → 0,087 → 0,031 bei Halbierung).
- Raytatouille berechnet diese Linien selbst aus Λ, n und λ (nichts in der Datei). Liegt eine in
  der Zelle eines Strahls, rechnet der Tracer weiter und zählt den Strahl unter dem Diagnosecode
  `trace.table_rayleigh_anomaly` (Warnung mit Anzahl und Fläche, ADR 0022/0023).
- Abhilfe für Nutzer: feinere Tabelle nahe der Anomalie oder Teiltabellen, deren Grenzen auf der
  Anomalie liegen (spätere Erweiterung).

### 7. Außerhalb der Tabelle und evaneszente Ordnungen

- Keine Extrapolation: Liegt (λ, σ) außerhalb der Achsen, endet der Strahl mit dem neuen Status
  `OUT_OF_TABLE` (Strahlverlust als Daten, ADR 0023).
- Eine Ordnung, die am Auftreffpunkt nicht propagiert, endet wie bisher mit `EVANESCENT`
  (ADR 0025 Punkt 7); Knoten, an denen sie nicht propagiert, enthalten NaN und machen die Zelle zu
  einer Anomaliezelle (Punkt 6).

## Verworfen

- **Effizienz je Polarisation (η_s, η_p):** keine Phase, keine Kopplung s ↔ p bei konischem
  Einfall, keine Winkelabhängigkeit, Fresnel weiter doppelt.
- **s/p-Basis nach Convention A:** bei senkrechtem Einfall von der Achsenregel abhängig (Punkt 2).
- **Interpolation von Betrag und Phase:** Phasensprünge um 2π und Nullstellen machen sie
  unstetig; die komplexen Einträge sind glatt außerhalb von Anomalien (J5).
- **Rechnung des Gitters zur Laufzeit im Tracer:** gehört in einen Löser außerhalb des Kerns
  (Fusion-ADR Punkt 1).

## Folgen

- `rtt-model`: Feld `diffraction_table`, Typ `DiffractionTable` im Ergebnisformat, Codes
  `surface.table_conflict`, `surface.table_invalid` (Achsen nicht gleichmäßig, weniger als 4 Punkte
  je Achse, Ordnungen doppelt, Periode ≠ Gitterkonstante der Phasenschicht).
- `rtt-compile`: Tabelle laden und an `CompiledSurface` hängen.
- `rtt-trace`/`rtt-polar`: Interaktion nach Punkt 3 in `apply_event`; Status `OUT_OF_TABLE`;
  Code `trace.table_rayleigh_anomaly`.
- **Abnahme:**
  1. Energie Σ|P_m E|² = 1 für zufällige Polarisationen (verlustfreies Gitter), P k = k_m;
  2. senkrechter Einfall: Ergebnis unabhängig von einer Drehung des globalen Systems um die
     Normale (keine Achsenregel);
  3. Pfadtransmission je Ordnung = Leistungsanteil aus FourierGrating (wie K1, jetzt mit
     Polarisation und Winkel);
  4. Interpolation an Kontrollpunkten innerhalb der angegebenen Genauigkeit; Anomaliezellen
     gemeldet;
  5. Phasenbezug: verschobenes Gitter (Punkt 4) gegen den OPL-Beitrag der Phasenschicht.
- Sollwertdateien aus RTT-Fusion (Stufe R7).
