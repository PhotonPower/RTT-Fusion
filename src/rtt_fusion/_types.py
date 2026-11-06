"""Shared array type aliases (inputs from np.linspace etc. are `floating[Any]`)."""

from __future__ import annotations

from typing import Any, TypeAlias

import numpy as np
from numpy.typing import NDArray

FloatArray: TypeAlias = NDArray[np.floating[Any]]
ComplexArray: TypeAlias = NDArray[np.complexfloating[Any, Any]]
IntArray: TypeAlias = NDArray[np.integer[Any]]
