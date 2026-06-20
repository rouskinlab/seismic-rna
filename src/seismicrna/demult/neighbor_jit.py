"""Numba-jitted neighbor generators for demultiplexing.

Kept separate so that importing ``neighbor`` (and thus ``demult``) does not
import numba, which is slow.  ``neighbor._get_neighbors_*`` import these
lazily, only when neighbors are actually generated.
"""

from numba import njit, types
from numba.typed import List


@njit(cache=False)
def rec_neighbors_2bit(
    orig: int,
    length: int,
    max_mismatches: int,
    pos: int,
    mismatches: int,
    current: int,
    out: List,
):
    """
    Recursively generate neighbor integers using 2-bit encoding.
    """
    if pos == length:
        out.append(current)
        return
    shift = 2 * (length - pos - 1)
    orig_base = (orig >> shift) & 0b11
    # Option 1: Keep the original base.
    candidate = current | (orig_base << shift)
    rec_neighbors_2bit(
        orig, length, max_mismatches, pos + 1, mismatches, candidate, out
    )
    # Option 2: Substitute with any alternative (A, C, G, T).
    if mismatches < max_mismatches:
        for base in range(4):
            if base != orig_base:
                candidate = current | (base << shift)
                rec_neighbors_2bit(
                    orig,
                    length,
                    max_mismatches,
                    pos + 1,
                    mismatches + 1,
                    candidate,
                    out,
                )


@njit(cache=False)
def generate_neighbors_2bit(orig: int, length: int, max_mismatches: int):
    """
    Generate all neighbor integers using 2-bit encoding.
    """
    out = List.empty_list(types.int64)
    rec_neighbors_2bit(orig, length, max_mismatches, 0, 0, 0, out)
    return out


@njit(cache=False)
def rec_neighbors_3bit(
    orig: int,
    length: int,
    max_mismatches: int,
    pos: int,
    mismatches: int,
    current: int,
    out: List,
):
    """
    Recursively generate neighbor integers using 3-bit encoding.
    This version allows substitutions to 'N' (encoded as 4).
    """
    if pos == length:
        out.append(current)
        return
    shift = 3 * (length - pos - 1)
    orig_base = (orig >> shift) & 0b111  # Will be one of 0-3 (barcode has no 'N').
    # Option 1: Keep the original base.
    candidate = current | (orig_base << shift)
    rec_neighbors_3bit(
        orig, length, max_mismatches, pos + 1, mismatches, candidate, out
    )
    # Option 2: Substitute with any alternative base (A, C, G, T, or N).
    if mismatches < max_mismatches:
        for base in range(5):
            if base != orig_base:
                candidate = current | (base << shift)
                rec_neighbors_3bit(
                    orig,
                    length,
                    max_mismatches,
                    pos + 1,
                    mismatches + 1,
                    candidate,
                    out,
                )


@njit(cache=False)
def generate_neighbors_3bit(orig: int, length: int, max_mismatches: int):
    """
    Generate all neighbor integers using 3-bit encoding.
    """
    out = List.empty_list(types.int64)
    rec_neighbors_3bit(orig, length, max_mismatches, 0, 0, 0, out)
    return out
