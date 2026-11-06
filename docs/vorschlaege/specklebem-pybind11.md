# Issue-Entwurf für SpeckleBem: Bau mit pybind11 3.x

> Entwurf aus RTT-Fusion (R8); in SpeckleBem nur lesend. Stand SpeckleBem `4f04efe`.

**Titel:** Python-Bindings bauen nicht mit pybind11 3.x (`bind_excitation.cpp:151`)

**Beobachtung:** `pip install ./python` mit pybind11 3.1 (GCC 13, Ubuntu 24.04):

```
python/src/bind_excitation.cpp:151:38: error: inconsistent types
  'pybind11::typing::Tuple<Eigen::Matrix<std::complex<double>, 3, 1>, Eigen::Matrix<std::complex<double>, 3, 1>>'
  and 'pybind11::typing::Tuple<Eigen::Matrix<std::complex<double>, -1, 3, 1>, ...>'
  deduced for lambda return type
```

Eine Lambda gibt je nach Zweig das Feld an einem Punkt (3er-Vektoren) oder an N Punkten
(N×3-Matrizen) als `pybind11::typing::Tuple` zurück; mit pybind11 3.x sind das zwei
verschiedene Typen, und die Ableitung des Rückgabetyps scheitert. Mit pybind11 2.13 baut es.

**Vorschlag:** expliziter Rückgabetyp der Lambda (`pybind11::tuple` oder ein gemeinsamer
`pybind11::typing::Tuple<pybind11::array, pybind11::array>`), oder zwei Überladungen für einen
Punkt und für N Punkte; dazu `pybind11>=2.13,<4` in `python/pyproject.toml` und ein CI-Lauf mit
pybind11 3.x.

**Umgehung bis dahin** (RTT-Fusion CI): pybind11 < 3 installieren und SpeckleBem mit
`--no-build-isolation` bauen; eine Beschränkung über `PIP_CONSTRAINT` allein greift nicht, wenn
pybind11 3.x schon in der Umgebung liegt.
