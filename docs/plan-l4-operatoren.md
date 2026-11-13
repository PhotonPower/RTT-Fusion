# Plan L4: Wiederverwendbare Operatoren der Löser

Stand: Entwurf zur Freigabe (2026-11-13). Grundlage: Beschluss L4 in
[konzept-node-workflow.md](konzept-node-workflow.md) und eine Bestandsaufnahme der fünf Löser
(FourierGrating, hp-FEM, SpeckleBem, CliffordBem, BWPM) und der Stellen in RTT-Fusion, die heute
je Einfall rechnen. Dateiangaben beziehen sich auf den Stand vom 2026-11-13.

## Ziel

Jeder Adapter bietet

- **`apply_many(handle, incidents)`** (Pflicht): viele Einfälle auf einmal. Der Adapter darf intern
  schleifen, nutzt aber, was der Löser wiederverwenden kann (Faktorisierung, Systemmatrix,
  Recycling).
- **`operator(handle, basis_in, basis_out)`** (optional): eine Block-Streumatrix zwischen den Basen
  der Ports (ebene Wellen bzw. Beugungsordnungen × 2 Polarisationen) bei fester Wellenlänge, mit
  Inhalts-Hash für den Operator-Speicher.

**Nutzen:** Parameterstudien, Feldpunkte × Wellenlängen × Polarisationen und Optimierung rechnen
dieselbe Struktur oft hundertfach. Bei FourierGrating, hp-FEM und SpeckleBem dominiert der Aufbau
bzw. die Faktorisierung; er lässt sich über die Einfälle teilen.

## Was die Löser heute können

| Löser | Adapter heute | Wiederverwendung im Löser | Streumatrix | Änderung nötig |
|---|---|---|---|---|
| **FourierGrating** | je Welle zwei volle Lösungen (s und p getrennt) | Schichtmodelle bleiben, Eigenproblem und Sternprodukt je Lösung neu | intern vorhanden (S11, S12, S21, S22), nicht in Python | **Löser:** S-Matrix in Python freigeben |
| **hp-FEM** | `grating.jones` je (λ, θ, φ), s und p je eigene Faktorisierung | `ScatteringOperator` mit `solve_many`, aber nur für `Scattering`, nicht für die konische Streuung der Gitter | Jones je Ordnung, keine S-Matrix | **Löser:** konischer Operator |
| **SpeckleBem** | Aufbau einmal je λ, je Welle `set_excitation` + GMRES | Operator, Nahmatrizen und Vorkonditionierer bleiben | keine | **nur Adapter**; optional Löser (Mehrfach-RHS bei `direct`) |
| **CliffordBem** | eine Lösung je `apply` (Wellen summiert) | `use_recycling`, `use_gcrodr`, HODLR – im Adapter ungenutzt | keine (nur Mie) | **nur Adapter** |
| **BWPM** | ein Feld je Lauf, Löser-Kern je Aufruf neu gebaut | nur Cache der Ausbreitungsmatrizen | keine (Basis G zu groß) | **nur Adapter**; optional Löser (Stapel von Feldern) |

Kostenangaben der Löser selbst: FourierGrating – das Eigenproblem dominiert (etwa 30 ZGEMM);
hp-FEM – Faktorisierung dominiert (100 rechte Seiten gebündelt 0,78 s gegenüber je einer
Faktorisierung); SpeckleBem – Aufbau 179–792 s gegen 0,2–22 s für `set_excitation`; CliffordBem –
Aufbau und Lösung ähnlich teuer, Recycling bis 3,6×; BWPM – linear in Punkten × Durchläufen, kein
Aufbau.

## Stufen

### O1 Schnittstelle und Schleifen-`apply_many` in allen Adaptern (nur RTT-Fusion)

- `apply_many` in allen fünf Adaptern, zunächst als Schleife über `apply` mit dem, was ohne
  Löser-Änderung geht:
  - **FourierGrating:** Einfälle nach Richtung (λ, φ, ψ) gruppieren; beide Polarisationen einer
    Richtung aus einer gemeinsamen Lösung, sobald O3 da ist (bis dahin wie heute).
  - **SpeckleBem:** wie heute (`set_excitation` bleibt), Einfall von hinten zulassen (der Löser
    kann es, nur der Adapter lehnt ab).
  - **CliffordBem:** `use_recycling` bzw. `use_gcrodr` einschalten, je Einfall `solve_rhs`;
    Einfall von hinten zulassen.
  - **hp-FEM:** Cache je (λ, θ, φ) bleibt; Schleife.
  - **BWPM:** Löser-Kern im Handle halten statt je Aufruf neu bauen (falls ohne Löser-Änderung
    möglich, sonst Schleife).
- **Abnahme:** `apply_many` gleich `apply` je Einfall auf 10⁻¹² für 20 zufällige Einfälle (Seed);
  CliffordBem mit Recycling gegen ohne auf die Toleranz des Lösers.

### O2 Die Aufrufer nutzen `apply_many` (nur RTT-Fusion)

- Heute je Einfall: `graph/particle.py` (Schleife über Bündel), `couplings/jones_table.py`
  (λ × σ-Raster), `graph/condense.py` (Proben je Gitterpunkt), `graph/grating.py`.
- Der Ausführer ruft je Block einmal mit allen Eingängen; dort sammeln und gebündelt anwenden.
- **Abnahme:** Ergebnisse unverändert (bestehende Tests, bitgleich wo der Löser deterministisch
  ist); Zahl der Löser-Aufrufe im Protokoll sinkt (z. B. Jones-Tabelle: Aufbau einmal je λ).

### O3 Operatoren aus den Lösern (braucht Freigaben in den Löser-Repos)

Je ein Issue im Löser-Repo, Text als Entwurf:

- **FourierGrating:** „Die Stapel-S-Matrix für ein (λ, α₀, γ₀) in Python freigeben, z. B.
  `Simulation.smatrix(lambda_um, phi_deg, psi_deg)` mit S11/S12/S21/S22, den Ordnungsindizes (n, m)
  und der s/p-Basis. Damit liefert ein Eigenproblem beide Polarisationen und den Einfall von
  hinten.“
- **hp-FEM:** „Einen `ConicalScatteringOperator` (oder `grating.operator(...)`) ergänzen:
  einmal je (kx, β, ω) faktorisieren, dann `solve_many` für s/p von Deck- und Substratseite, mit den
  Jones- bzw. R/T-Matrizen der Ordnungen für alle vier Blöcke.“
- **SpeckleBem** (optional): „Die dichte LU (`DenseLu`) bei `solver='direct'` über
  `set_excitation` hinweg behalten und `solve_many(excitations)` ergänzen.“
- **BWPM** (optional): „Einen Stapel von Einfallsfeldern `(n, ny, nx, 2)` in `apply`/`core.run`
  annehmen, die FFTs bündeln und den Löser-Kern im Handle behalten.“
- Danach `operator(...)` im Adapter für FourierGrating und hp-FEM.
- **Abnahme (wie L4):** `operator` und `apply` liefern für 20 zufällige Einfälle dasselbe auf
  10⁻¹²; Reziprozität Sᵀ = S in symmetrischer Basis, wo der Löser reziprok ist; Einfall von hinten
  gegen gespiegelte Struktur (L5).

### O4 Verdichten mit allgemeinem Operator und Umschaltschwelle

- `graph/condense.py` nimmt heute nur seitlich invariante, isotrope Gruppen (Skalare c_s, c_p je
  Gitterpunkt). Mit O3 wird der Tabelleneintrag ein 2×2-Block (s/p ein → s/p aus) je Ausgangsrichtung
  bzw. Ordnung: anisotrope Stapel und Gitter lassen sich dann verdichten.
- Der Planer (`graph/planner.py`) hat noch keine Kosten aus `estimate` und keine Schwelle: er nimmt
  den Operator, sobald die Zahl der Einfälle die Schwelle aus `estimate` übersteigt.
- **Abnahme:** verdichtetes Gitter = aufgelöstes auf die Toleranz der Tabelle (wie W55); der Planer
  wählt bei wenigen Einfällen `apply`, bei vielen den Operator (Testfall mit festen `estimate`).

## Was L4 nicht löst

**Linsengruppen aus Raytatouille** (z. B. ein Objektiv als Gruppe) lassen sich mit L4 nicht
verdichten. Ihr Ergebnis hängt davon ab, wo das Licht eintritt; in einer Basis ebener Wellen wäre
der Operator dicht und groß, und die Strahlrechnung durch die Gruppe ist ohnehin der schnelle Weg.
Sie bleiben aufgelöst. (In der Übersicht vom 2026-11-13 stand irrtümlich, L4 mache das möglich.)

## Reihenfolge und Freigaben

1. **O1 und O2** – nur RTT-Fusion, kein Eingriff in Löser-Repos. Vorschlag: zuerst.
2. **O3** – je Löser ein Issue (Texte oben); Umsetzung im Löser-Repo nur mit deiner Freigabe je Repo.
   Wichtigste: FourierGrating und hp-FEM (dort liegt der Gewinn).
3. **O4** – nach O3.

Entscheidungen für den Maintainer:

- **E1:** Reihenfolge wie oben?
- **E2:** Issues für FourierGrating und hp-FEM anlegen? SpeckleBem und BWPM (optional) auch?
- **E3:** Einfall von hinten in den Adaptern von SpeckleBem und CliffordBem zulassen (O1), da die
  Löser es können?
