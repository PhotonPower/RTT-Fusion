# Issue-Entwurf für Fieldtracing: Debye mit einer Pupille neben der Achse

Nicht eingereicht. Ein Vorschlag für das Fieldtracing-Repository, entstanden bei RTT-Fusion#88.

## Befund

`fieldtracing.chain.debye` rechnet für eine ideale Kugelwelle mit einer Pupillenscheibe, die nicht um
die Achse des Feldkoordinatensystems liegt, zu viel Leistung und eine zu hohe Spitze. Der Fokus
liegt dabei auf der Achse, wie es der Docstring verlangt („reference sphere centred on the axis in
the target plane“).

Pupillenradius 1,2 mm, Abstand 58 mm, λ = 546 nm, 24 Ringe, `n = 128`:

| Versatz der Pupille | sin_max | Spitze | Leistung im Bild | Leistung der Pupille |
|---|---|---|---|---|
| 0 | 0,02069 | 20 399 | 4,499 | 4,524 |
| 0,6 mm | 0,03102 | 21 677 (+6,3 %) | 4,620 (+2,1 %) | 4,524 |
| 1,2 mm | 0,04134 | 21 682 (+6,3 %) | 4,610 (+1,9 %) | 4,524 |

Mehr Stützstellen ändern daran nichts. Der Fleck ist zu schmal; das spricht für Licht unter zu
großen Winkeln, also eine zu große k-Domäne (`_support`, der Rand der Scheibe `|k⊥| ≤ k·sin_max`
oder die Füllung „nearest“ außerhalb der Samples).

## Nachstellen

```python
import math
import numpy as np
from fieldtracing import chain as ch
from fieldtracing.field import GridlessField

lam_um, d, R = 0.546, 58.0, 1.2
k = 2 * np.pi / (lam_um * 1e-3)


def field(s, rings=24):
    px, py = [0.0], [0.0]
    for r in range(1, rings + 1):
        for j in range(6 * r):
            a = 2 * np.pi * j / (6 * r)
            px.append(r / rings * np.cos(a)); py.append(r / rings * np.sin(a))
    p = np.stack([s + R * np.array(px), R * np.array(py)], -1)
    w = np.sqrt((p ** 2).sum(1) + d * d)
    hess = np.array([-k * (np.eye(2) / wi - np.outer(pi, pi) / wi ** 3) for pi, wi in zip(p, w)])
    amp = np.zeros((len(p), 3), complex); amp[:, 1] = 1
    return GridlessField(p, -k * (w - d), -k * p / w[:, None], hess, amp, lam_um)


for s in (0.0, 0.6, 1.2):
    f = field(s)
    a = float(np.hypot(*f.points.T).max())
    out, _ = ch.debye(f, d, a / math.hypot(a, d), n=128)
    v = (np.abs(out.values) ** 2).sum(-1)
    print(s, v.max(), v.sum() * out.step ** 2)
```

## Erwartet

Spitze und Leistung unabhängig vom Versatz (bis auf die Neigung des Flecks, hier < 0,1 %).

## Umgang in RTT-Fusion

`field.propagate_to` lässt die Pupille für Debye auf ihrer Achse, solange der Fleck höchstens einen
halben Pupillenradius daneben liegt (RTT-Fusion#88). Nur weiter außen, etwa bei einem steilen
Teilbündel, zentriert es auf den Fokus und nimmt diesen Fehler in Kauf. Dort lag die Spitze 6 bis
10 % über der direkten Summe ebener Wellen (RTT-Fusion#88, Fall B).

Stand: Fieldtracing 0.2.0 (Pin `ea6da7b`).
