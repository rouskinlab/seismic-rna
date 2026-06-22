"""Numba-jitted FASTQ read/quality generators for simulation.

Kept separate so that importing ``fastq`` (and thus ``sim``) does not import
numba, which is slow.  ``fastq.generate_fastq_record`` imports these lazily.
"""

import numpy as np
from numba import njit

from ..core.rel.pattern import MATCH, SUB_A, SUB_C, SUB_G, SUB_T, DELET
from ..core.seq.region import BASEA, BASEC, BASEG, BASET, BASEN


@njit()
def _complement(base: str):
    """JIT-compiled function to get the complement of a base."""
    if base == BASEA:
        return BASET
    if base == BASEC:
        return BASEG
    if base == BASEG:
        return BASEC
    if base == BASET:
        return BASEA
    return BASEN


@njit()
def generate_fastq_read_qual(
    rels: np.ndarray,
    refseq: str,
    adapter: str,
    read_length: int,
    revcomp: bool,
    hi_qual: str,
    lo_qual: str,
):
    """Generate a FASTQ line for a read."""
    (ref_length,) = rels.shape
    # Map each type of relationship to a type of base in the read.
    if revcomp:
        ref_pos = ref_length - 1
        ref_inc = -1
        sub_a = BASET
        sub_c = BASEG
        sub_g = BASEC
        sub_t = BASEA
    else:
        ref_pos = 0
        ref_inc = 1
        sub_a = BASEA
        sub_c = BASEC
        sub_g = BASEG
        sub_t = BASET
    # Fill the read with high-quality Gs (the default base in Illumina).
    read = np.full(read_length, BASEG)
    qual = np.full(read_length, hi_qual)
    # Add bases to the read.
    read_pos = 0
    while read_pos < read_length and 0 <= ref_pos < ref_length:
        rel = rels[ref_pos]
        if rel == MATCH:
            ref_base = refseq[ref_pos]
            read[read_pos] = _complement(ref_base) if revcomp else ref_base
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel == SUB_T:
            read[read_pos] = sub_t
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel == SUB_G:
            read[read_pos] = sub_g
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel == SUB_C:
            read[read_pos] = sub_c
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel == SUB_A:
            read[read_pos] = sub_a
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel == DELET:
            # Deletion: do not add any base to the read.
            pass
        else:
            # Assume a low-quality base call: add an N to the read.
            read[read_pos] = BASEN
            qual[read_pos] = lo_qual
            read_pos += 1
        ref_pos += ref_inc
    # Add the adapter to the end of the read.
    adapter_pos = 0
    adapter_length = len(adapter)
    while read_pos < read_length and adapter_pos < adapter_length:
        read[read_pos] = adapter[adapter_pos]
        read_pos += 1
        adapter_pos += 1
    return "".join(read), "".join(qual)
