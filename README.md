# RTT-Fusion

Prototyp der **Raytatouille Fusion**: wellenoptische Löser als Gebiete mit Ports, gekoppelt mit
Raytracing (Raytatouille) und Field Tracing, in den Konventionen von Raytatouille.

- Gemeinsame Darstellungen an Ports: Gitterfeld (G) und Ebene-Wellen-Spektrum (S), Umwandlung
  zwischen beiden, Leistung (Poynting-Fluss).
- Adapter je Löser mit `describe`, `prepare`, `apply`, `estimate` und Konformitätstests.
- Pilotfälle, die die Kopplung gegen unabhängige Referenzen prüfen.

Plan: Fieldtracing `docs/rtt-fusion/plan.md`; Stufen dieses Repositorys:
[docs/roadmap.md](docs/roadmap.md).

## Installieren und testen

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt
pip install -e ".[modules]"     # BWPM, FourierGrating, Fieldtracing aus ihren Repositories
pytest
ruff check . && ruff format --check . && mypy
```

FourierGrating braucht zum Bauen einen C++-Compiler und Eigen 3.4 (wird sonst geladen).

## Dokumentation

- [Stufen und Entscheidungen](docs/roadmap.md)
- [Konventionen](docs/konventionen.md)
- [Schnittstelle der Adapter](docs/schnittstelle.md)
- [Testfälle](docs/testfaelle.md)
- [ADR-Entwurf für Raytatouille](docs/adr-entwurf-fusion.md)
- [ADR-Entwurf tabellierte Beugungsinteraktion](docs/adr-entwurf-beugungstabelle.md)
- [Quellen](docs/quellen.md)
- [Regeln für Programmieragenten](AGENTS.md)

## Lizenz

MIT, siehe [LICENSE](LICENSE).
