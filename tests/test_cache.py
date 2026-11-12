"""Results kept across runs (R23): only what changed is computed again.

Definitions checked on the analytic Michelson of W4 with a display and a sweep:
- W37 a second run with the same cache computes nothing and every block is up to date; after a
  parameter of mirror B changes, exactly B and the blocks behind it are out of date and only they
  are computed again; a changed setting (rays) makes every block out of date.
- W38 the cache keeps to its limit, drops the least recently used results and says so.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from rtt_fusion.graph import Node, display, meta
from rtt_fusion.graph.execute import Executor, ResultCache

sys.path.insert(0, str(Path(__file__).parent))
from test_execute import WAVELENGTH_MM, michelson


def workflow(cache: ResultCache, rings: int = 30) -> Executor:
    base = michelson()
    g = base.graph
    g.nodes["src"].params["rings"] = rings
    ex = Executor(g, base.registry, base.evaluators, cache=cache)
    ex.evaluators.update(display.register(ex.registry))
    ex.evaluators.update(meta.register(ex.registry))
    g.add(Node("show", "display.image"))
    g.connect(("cam", "image"), ("show", "in"))
    g.add(
        Node("ps", "meta.sweep", params={"measure": "centre", "values": [0.0, WAVELENGTH_MM / 8]})
    )
    g.connect(("ps", "drive"), ("mA", "drive"))
    g.connect(("cam", "image"), ("ps", "observe"))
    g.add(Node("curve", "display.plot"))
    g.connect(("ps", "results"), ("curve", "in"))
    return ex


def run(ex: Executor) -> None:
    ex.evaluate("show", "view")
    ex.evaluate("curve", "view")


def test_w37_only_what_changed_is_computed() -> None:
    cache = ResultCache()
    first = workflow(cache)
    assert not any(first.status().values())  # nothing computed yet
    run(first)
    assert sum(first.calls.values()) > 0
    again = workflow(cache)
    assert all(again.status().values())  # the same workflow in a new executor: all up to date
    run(again)
    assert sum(again.calls.values()) == 0
    # a parameter of mirror B: B and everything behind it are out of date
    changed = workflow(cache)
    changed.graph.nodes["mB"].params["reflectance_amplitude"] = -0.5
    status = changed.status()
    assert not status["mB"] and not status["bs"] and not status["cam"]
    assert not status["show"] and not status["ps"] and not status["curve"]
    assert status["src"] and status["mA"]
    run(changed)
    assert ("src", "out") not in changed.calls and ("mA", "out") not in changed.calls
    assert changed.calls[("mB", "out")] >= 1
    # a different setting (rays) changes every key
    other = workflow(cache, rings=20)
    assert not any(other.status().values())


def test_w38_cache_limit() -> None:
    cache = ResultCache(limit_bytes=200_000)
    for k in range(5):
        cache[f"k{k}"] = np.zeros(10_000)  # 80 kB each
    assert cache.size_bytes <= 200_000
    assert "k0" not in cache and "k4" in cache
    assert cache.warnings and "dropped" in cache.warnings[-1]
    _ = cache["k3"]  # used: stays
    cache["k5"] = np.zeros(10_000)
    assert "k3" in cache


def test_w56_marks_count_what_the_displays_need() -> None:
    """A splitter feeds two cameras, only one is shown: after the run the splitter is up to date
    (its way to the other camera is not needed), the other camera gets no mark; without
    `wanted` every linked output counts (the splitter is out of date)."""
    from rtt_fusion.graph import Graph, builtin_registry, optics

    reg = builtin_registry()
    ev = dict(optics.register(reg))
    ev.update(display.register(reg))
    g = Graph()
    g.add(
        Node(
            "src",
            "optics.collimated_source",
            params={"center": [0, 0, -10], "radius": 1.0, "rings": 4, "polarisation": [0, 1, 0]},
        )
    )
    g.add(Node("split", "optics.splitter", params={"point": [0, 0, 0], "normal": [1, 0, -1]}))
    g.add(Node("camT", "optics.detector", params={"point": [0, 0, 10], "normal": [0, 0, -1]}))
    g.add(Node("camR", "optics.detector", params={"point": [10, 0, 0], "normal": [-1, 0, 0]}))
    g.add(Node("shown", "display.image"))
    g.connect(("src", "out"), ("split", "front_a_in"))
    g.connect(("split", "back_a_out"), ("camT", "in"))
    g.connect(("split", "front_b_out"), ("camR", "in"))
    g.connect(("camT", "image"), ("shown", "in"))
    ex = Executor(g, reg, ev)
    ex.evaluate("shown", "view")
    marks = ex.status(wanted=[("shown", "view")])
    assert marks == {"shown": True, "camT": True, "split": True, "src": True}
    assert ex.needed([("shown", "view")])["split"] == {"back_a_out"}
    assert ex.status()["split"] is False  # every linked output: the way to camR is missing
    ex.graph.nodes["camT"].params["half_size"] = [2.0, 2.0]
    assert ex.status(wanted=[("shown", "view")])["camT"] is False  # an edit: out of date
