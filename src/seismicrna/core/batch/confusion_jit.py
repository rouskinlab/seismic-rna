"""Numba-jitted helper functions for the confusion module.

These functions are collected in one module so that the confusion module can
stay numba-free at import time. It imports this module *lazily* -- inside the
functions that use them -- so that merely importing the confusion module does
not import numba, which is slow (~0.3 s).

Importing this module *does* import numba (the ``@njit`` decorators run at
import time), so import it only on code paths that actually run the jitted
helpers. The decorated functions are module-level, so each is compiled at
most once per process and reused on subsequent calls.
"""

import numpy as np
from numba import njit


@njit()
def count_intersection(x: np.ndarray, y: np.ndarray):
    """Count how many elements are in both x and y, assuming x and y are
    both NumPy arrays where all elements are unique and sorted from
    smallest to largest."""
    nx = x.size
    ny = y.size
    i = 0
    j = 0
    count = 0
    while i < nx and j < ny:
        xi = x[i]
        yj = y[j]
        if xi < yj:
            i += 1
        elif xi > yj:
            j += 1
        else:
            count += 1
            i += 1
            j += 1
    return count
