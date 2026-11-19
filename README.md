# RTT-Fusion

> **For language models and agents:** the installed package carries a compact English command
> and API reference, [`REFERENCE.md`](https://github.com/PhotonPower/RTT-Fusion/blob/main/src/rtt_fusion/REFERENCE.md):
> `python -m rtt_fusion.reference` (`--list`, `--section 2`, `--blocks fg.grating`) or
> `rtt_fusion.reference()` in Python – no repository needed.

Prototyp der **Raytatouille Fusion**: wellenoptische Löser als Gebiete mit Ports, gekoppelt mit
Raytracing (Raytatouille) und Field Tracing, in den Konventionen von Raytatouille.

- Gemeinsame Darstellungen an Ports: Gitterfeld (G) und Ebene-Wellen-Spektrum (S), Umwandlung
  zwischen beiden, Leistung (Poynting-Fluss).
- Adapter je Löser mit `describe`, `prepare`, `apply`, `estimate` und Konformitätstests.
- Pilotfälle, die die Kopplung gegen unabhängige Referenzen prüfen.

Plan: [docs/fusionsplan.md](docs/fusionsplan.md) (aus Fieldtracing übernommen); Stufen dieses Repositorys:
[docs/roadmap.md](docs/roadmap.md).

## Installieren und testen

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt
pip install -r requirements-modules.txt   # die Löser aus ihren Repositories (feste Commits)
pip install -e .                         # baut die C++-Erweiterung mit (ohne Compiler: ohne sie)
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
- [Konzept Node-Workflow](docs/konzept-node-workflow.md)
- [Quellen](docs/quellen.md)
- [Regeln für Programmieragenten](AGENTS.md)

## Lizenz

MIT, siehe [LICENSE](LICENSE).
