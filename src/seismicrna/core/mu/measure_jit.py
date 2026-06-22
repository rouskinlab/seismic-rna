"""Numba-jitted helper for core.mu.measure.

Kept separate so that importing ``measure`` does not import numba (slow);
``calc_gini`` imports this lazily only when it actually computes.
"""

import numpy as np
from numba import njit


@njit()
def calc_sum_abs_diff(x: np.ndarray):
    """Sum the absolute difference along axis 0."""
    n = x.shape[0]
    sum_abs_diff = np.zeros(x.shape[1:])
    for i in range(n):
        for j in range(i + 1, n):
            sum_abs_diff += np.abs(x[i] - x[j])
    return sum_abs_diff
