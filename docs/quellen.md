# Quellen

| Kürzel | Quelle | Status | Verwendet für |
| --- | --- | --- | --- |
| byrnes | S. J. Byrnes, „Multilayer optical calculations“, arXiv:1603.02720v5 | von Raytatouille übernommen (dort Volltext geprüft) | Fresnel Convention A |
| fg-conv | FourierGrating `docs/conventions.md`, `src/physics/incidence.cpp`, `src/physics/efficiency.cpp` (Commit `959c93b`; bis `174ef93` unverändert) | gelesen | Umrechnung im FourierGrating-Adapter |
| bwpm | BWPM `bwpm.adapter` (Commit auf `main`, B4e) | gelesen | Querprobe, BWPM-Adapter |
| rtt-adr0025 | Raytatouille ADR 0025 „Gitter als Phasenschicht“ und `tests/reference/m4/grating_transmission.rtt.json` (Commit `708b6dd`) | gelesen | Pilot K1: Ordnungen, Effizienz, Bench |
| bohren-huffman | C. F. Bohren, D. R. Huffman, Absorption and Scattering of Light by Small Particles, Wiley (1983), Kap. 4 | Standardformeln (a_n, b_n, S1, S2, Q_ext), unabhängig nachgerechnet gegen CliffordBem | Mie-Referenz `tests/mie_reference.py` |
| cbem | CliffordBem `docs/python.md`, `docs/results_beams.md` (Commit `41e64e6`) | gelesen | Adapter, Richards-Wolf-Strahl |
| sbem | SpeckleBem `docs/06_conventions.md`, `docs/10_python_api.md` (Commit `4f04efe`) | gelesen | SpeckleBem-Adapter |
| ft-adr0004 | Fieldtracing ADR 0004 „Feldamplitude aus Raytatouilles leistungsnormierter PRT-Matrix und Gewicht“, `rtt_bridge.py`, `mesh.py`, `chain.plan_propagation` (Commit `6da9983`) | gelesen | R18: Strahlen → Feld, Ausbreitung auf H-Bändern |
| richards-wolf | B. Richards, E. Wolf, „Electromagnetic diffraction in optical systems II“, Proc. R. Soc. A 253, 358 (1959); als Bessel-Integral in Fieldtracing `focus.richards_wolf` (Commit `6da9983`) | Formel aus Fieldtracing übernommen (dort gegen die Arbeit geprüft); Originalarbeit hier nicht gelesen | W13: Fokus der aplanatischen Linse |
| gauss-beam | Gauß-Strahl w(z) = w₀ sqrt(1 + (z/z_R)²), z_R = π w₀²/λ, I = I₀ (w₀/w)² exp(−2r²/w²) | Standardformel | W14, W15 |
| rtt-paraboloid | Raytatouille `tests/reference/m6/paraboloid_f50.rtt.json` (MIT, Commit `61d1e9c`), kopiert nach `tests/data/paraboloid_f50.rtt.json` | übernommen | W17: Fokus eines Paraboloids (Airy) |
| airy | Airy-Muster [2 J₁(v)/v]², v = k NA r, Spitze P π NA²/λ² | Standardformel | W17 |
| raman-nath | Dünnes Phasengitter: e^{iφ₀ sin Kx} = Σ J_m(φ₀) e^{imKx} (Jacobi–Anger), Beugungseffizienzen J_m(φ₀)² | Standardformel | W23 |
| rtt-ideal-lens | Raytatouille ADR 0031 (ideale Linse: kollineare Abbildung, OPL exakt für die Auslegungs-Konjugierte `object_distance`), `docs/architecture.md` | gelesen | Testsystem `linnik_mini`, Grenze bei großen Winkeln (nicht aplanatisch) |
| rtt-michelson | Raytatouille `tests/reference/m0/michelson.rtt.json` (MIT), kopiert nach `tests/data/michelson.rtt.json` | übernommen | Pilot W9 über das Raytatouille-Backend |
| hpfem | hp-FEM `docs/python.md` (Gitter, Konventionen), `python/hpfem/grating.py` (`jones`, Commit `dce07d5`) | gelesen | hp-FEM-Adapter |
