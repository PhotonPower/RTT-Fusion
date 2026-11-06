# Quellen

| Kürzel | Quelle | Status | Verwendet für |
| --- | --- | --- | --- |
| byrnes | S. J. Byrnes, „Multilayer optical calculations“, arXiv:1603.02720v5 | von Raytatouille übernommen (dort Volltext geprüft) | Fresnel Convention A |
| fg-conv | FourierGrating `docs/conventions.md`, `src/physics/incidence.cpp`, `src/physics/efficiency.cpp` (Commit `959c93b`) | gelesen | Umrechnung im FourierGrating-Adapter |
| bwpm | BWPM `bwpm.adapter` (Commit auf `main`, B4e) | gelesen | Querprobe, BWPM-Adapter |
| rtt-adr0025 | Raytatouille ADR 0025 „Gitter als Phasenschicht“ und `tests/reference/m4/grating_transmission.rtt.json` (Commit `708b6dd`) | gelesen | Pilot K1: Ordnungen, Effizienz, Bench |
| bohren-huffman | C. F. Bohren, D. R. Huffman, Absorption and Scattering of Light by Small Particles, Wiley (1983), Kap. 4 | Standardformeln (a_n, b_n, S1, S2, Q_ext), unabhängig nachgerechnet gegen CliffordBem | Mie-Referenz `tests/mie_reference.py` |
| cbem | CliffordBem `docs/python.md`, `docs/results_beams.md` (Commit `41e64e6`) | gelesen | Adapter, Richards-Wolf-Strahl |
