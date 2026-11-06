# Regeln für Programmieragenten

Diese Datei liest jeder Agent zu Beginn jeder Sitzung. Sie gilt zusammen mit `docs/roadmap.md`,
`docs/konventionen.md`, `docs/schnittstelle.md` und `docs/quellen.md`.

## Projekt in einem Satz

RTT-Fusion ist der Prototyp der „Raytatouille Fusion“: Er bindet wellenoptische Löser
(BWPM, FourierGrating, später CliffordBem, hp-FEM, SpeckleBem, BioOptoMD) über eine gemeinsame
Port-Schnittstelle an Raytatouille und Field Tracing an und prüft die Kopplung an Pilotfällen,
bevor sie in Raytatouille portiert wird (Plan: Fieldtracing `docs/rtt-fusion/plan.md`).

## Bauen und testen

```bash
pip install -r requirements-dev.txt && pip install -e ".[modules]"
pytest
ruff check . && ruff format --check . && mypy
```

Sollwertdateien für den C++-Port: `python tools/export_reference.py` schreibt
`reference/*.result.json` (Format `raytatouille-result`); `tests/test_reference_files.py`
prüft, dass sie aktuell sind.

Das Extra `modules` installiert die Löser aus ihren Repositories auf festen Commits
(`pyproject.toml`). Ohne ein Modul werden dessen Tests übersprungen, außer mit
`RTT_FUSION_REQUIRE_MODULES=1` (CI), dann schlagen sie fehl.

## Arbeitsregeln

1. **Ein Auftrag = ein Ziel einer Stufe** aus `docs/roadmap.md`.
2. **Erst Test, dann Code.** Sollwerte aus analytischen Lösungen oder unabhängigen Rechnungen
   (Fresnel, Airy, Transfer-Matrix, Gittergleichung, Energieerhaltung, ein zweiter Löser), nie aus
   dem geprüften Adapter. Tests und Sollwerte nur mit Freigabe des Maintainers ändern.
3. **Konventionen sind Gesetz** (`docs/konventionen.md`): die von Raytatouille. Umrechnungen der
   Konventionen eines Lösers geschehen **nur in dessen Adapter**, mit Konformitätstest.
4. **Löser-Repositories nur lesen** (BWPM, FourierGrating, Fieldtracing, Raytatouille,
   CliffordBem, hp-FEM, SpeckleBem, BioOptoMD). Wird dort etwas gebraucht, entsteht ein Vorschlag
   als Dokument oder Issue-Entwurf unter `docs/`.
5. **Formeln mit Quelle** (`docs/quellen.md`); Gleichungsnummern fremder Quellen nur aus dem
   Volltext.
6. **Determinismus:** Zufall nur mit explizitem Seed.
7. **Keine Geheimnisse** (Tokens, Passwörter, persönliche Pfade) in Dateien, Commits, Remote-URLs.
8. **Unklarheit = Frage**: Issue mit Label `frage`; Abweichung von einer Quelle: Label `paper`.

## Sprache und Stil

- Code, Kommentare, Docstrings, Commit-Messages englisch; Projektdokumentation deutsch.
- ruff (Zeilenlänge 100), mypy strict; Commits im Imperativ, ein Thema pro Commit.

## Definition of Done (PR-Checkliste)

- [ ] CI grün
- [ ] Referenztests mit unabhängigem Sollwert, Konformitätstest je Adapter
- [ ] `docs/testfaelle.md` aktualisiert (Status, Messwerte)
- [ ] neue Quellen in `docs/quellen.md`
- [ ] Docstrings mit Einheiten, Vorzeichen und Verweis
- [ ] `CHANGELOG.md` unter „Unreleased“
