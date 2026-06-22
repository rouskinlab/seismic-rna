"""Numba-jitted array helpers.

These functions are kept separate from ``array`` so that importing the
(numba-free) ``array`` module does not import numba, which is slow.  Import
from this module only on code paths that actually need the jitted helpers
(e.g. ``calc_inverse`` imports them only when filling missing indexes).
"""

import numpy as np
from numba import njit

from .array import MISSING


@njit()
def fill_inverse_fwd(inverse: np.ndarray, default: int):
    """Fill missing indexes in `inverse` in forward order."""
    fill = default
    for i in range(inverse.size):
        inv = inverse[i]
        if inv == MISSING:
            inverse[i] = fill
        else:
            fill = inv


@njit()
def fill_inverse_rev(inverse: np.ndarray, default: int):
    """Fill missing indexes in `inverse` in reverse order."""
    fill = default
    for i in range(inverse.size - 1, -1, -1):
        inv = inverse[i]
        if inv == MISSING:
            inverse[i] = fill
        else:
            fill = inv
