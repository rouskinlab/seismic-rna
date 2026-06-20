"""Numba-jitted helper functions for the muts module.

These functions are collected in one module so that the muts module can
stay numba-free at import time. It imports this module *lazily* -- inside the
functions that use them -- so that merely importing the muts module does
not import numba, which is slow (~0.3 s).

Importing this module *does* import numba (the ``@jit`` decorators run at
import time), so import it only on code paths that actually run the jitted
helpers. The decorated functions are module-level, so each is compiled at
most once per process and reused on subsequent calls.
"""

import numpy as np
from numba import jit

from ..rel.pattern import MATCH


@jit()
def fill_matches(
    matrix: np.ndarray,
    index5s: np.ndarray,
    index3s: np.ndarray,
    unmasked_read_indexes: np.ndarray,
):
    """Fill all covered positions with matches."""
    for i, read_index in enumerate(unmasked_read_indexes):
        matrix[read_index, index5s[i] : index3s[i]] = MATCH
