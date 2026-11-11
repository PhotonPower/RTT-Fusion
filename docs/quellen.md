# Quellen

| Kürzel | Quelle | Status | Verwendet für |
| --- | --- | --- | --- |
| byrnes | S. J. Byrnes, „Multilayer optical calculations“, arXiv:1603.02720v5 | von Raytatouille übernommen (dort Volltext geprüft) | Fresnel Convention A |
| fg-conv | FourierGrating `docs/conventions.md`, `src/physics/incidence.cpp`, `src/physics/efficiency.cpp` (Commit `959c93b`) | gelesen | Umrechnung im FourierGrating-Adapter |
| bwpm | BWPM `bwpm.adapter` (Commit auf `main`, B4e) | gelesen | Querprobe, BWPM-Adapter |
| rtt-adr0025 | Raytatouille ADR 0025 „Gitter als Phasenschicht“ und `tests/reference/m4/grating_transmission.rtt.json` (Commit `708b6dd`) | gelesen | Pilot K1: Ordnungen, Effizienz, Bench |
| bohren-huffman | C. F. Bohren, D. R. Huffman, Absorption and Scattering of Light by Small Particles, Wiley (1983), Kap. 4 | Standardformeln (a_n, b_n, S1, S2, Q_ext), unabhängig nachgerechnet gegen CliffordBem | Mie-Referenz `tests/mie_reference.py` |
| cbem | CliffordBem `docs/python.md`, `docs/results_beams.md` (Commit `41e64e6`) | gelesen | Adapter, Richards-Wolf-Strahl |
| sbem | SpeckleBem `docs/06_conventions.md`, `docs/10_python_api.md` (Commit `4f04efe`) | gelesen | SpeckleBem-Adapter |
| rtt-michelson | Raytatouille `tests/reference/m0/michelson.rtt.json` (MIT), kopiert nach `tests/data/michelson.rtt.json` | übernommen | Pilot W9 über das Raytatouille-Backend |
| hpfem | hp-FEM `docs/python.md` (Gitter, Konventionen), `python/hpfem/grating.py` (`jones`, Commit `dce07d5`) | gelesen | hp-FEM-Adapter |
