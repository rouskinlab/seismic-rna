"""

Simulate SAM Files Module
========================================================================

"""

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from ..core.rand import rng
from ..core.rel import (MIN_QUAL, MAX_QUAL, NOCOV)
from ..core.relaux import relvec_to_read, validate_relvec
from ..core.sect import index_to_seq
from ..core.seq import DNA
from ..core.xam import (as_sam, sam_header, SAM_NOREF, FLAG_PAIRED, FLAG_PROPER,
                        FLAG_FIRST, FLAG_SECOND, FLAG_REVERSE, FLAG_MREVERSE)


def _relvec_to_sam_line(read: str,
                        relvec: np.ndarray,
                        ref: str,
                        refseq: DNA, *,
                        mapq: int,
                        flag: int = 0,
                        hi_qual: str = MAX_QUAL,
                        lo_qual: str = MIN_QUAL,
                        ins_len: int | Sequence[int] = 1):
    seq, qual, cig, end5, end3 = relvec_to_read(refseq,
                                                relvec,
                                                hi_qual,
                                                lo_qual,
                                                ins_len)
    return as_sam(read, flag, ref, end5, mapq, cig, SAM_NOREF, 0, 0, seq, qual)


def _find_blank_range(side3: bool, length: int, end5: int, end3: int):
    if not length:
        length = end3 - end5 + 1
    if length < 1:
        raise ValueError(f"Length of read must be ≥ 1, but got {length}")
    return (
        (end5 - 1, max(end5 - 1, end3 - length))
        if side3 else
        (min(end3, end5 - 1 + length), end3)
    )


def _mask_relvec(relvec: np.ndarray, mask5: int, mask3: int):
    masked = relvec.copy()
    masked[mask5: mask3] = NOCOV
    return masked


def _relvec_to_sam_pair(read: str,
                        relvec: np.ndarray,
                        ref: str,
                        refseq: DNA, *,
                        len1: int = 0,
                        len2: int = 0,
                        **kwargs):
    # Find the 5' and 3' coordinates of the relation vector.
    end5, end3 = validate_relvec(relvec)
    # Determine whether the first read is reversed.
    rev1 = bool(rng.integers(2))
    # Set the flags of the first and second reads accordingly.
    flag1 = FLAG_PAIRED + FLAG_PROPER + FLAG_FIRST
    flag2 = FLAG_PAIRED + FLAG_PROPER + FLAG_SECOND
    if rev1:
        flag1 += FLAG_REVERSE
        flag2 += FLAG_MREVERSE
    else:
        flag1 += FLAG_MREVERSE
        flag2 += FLAG_REVERSE
    # Split the relation vector into parts for read 1 and read 2.
    rv1 = _mask_relvec(relvec, *_find_blank_range(rev1, len1, end5, end3))
    rv2 = _mask_relvec(relvec, *_find_blank_range(not rev1, len2, end5, end3))
    # Convert read 1 and read 2 into one SAM line each.
    line1 = _relvec_to_sam_line(read, rv1, ref, refseq, flag=flag1, **kwargs)
    line2 = _relvec_to_sam_line(read, rv2, ref, refseq, flag=flag2, **kwargs)
    return line1, line2


def relvecs_to_sam_lines(relvecs: pd.DataFrame,
                         ref: str,
                         paired: bool,
                         **kwargs):
    refseq = index_to_seq(relvecs.columns)
    yield sam_header(ref, refseq)
    reads = zip(relvecs.index, relvecs.values, strict=True)
    if paired:
        for read, relvec in reads:
            yield from _relvec_to_sam_pair(read, relvec, ref, refseq, **kwargs)
    else:
        for read, relvec in reads:
            yield _relvec_to_sam_line(read, relvec, ref, refseq, **kwargs)


def relvecs_to_sam_file(file: Path, *args, overwrite: bool = False, **kwargs):
    with open(file, 'w' if overwrite else 'x') as f:
        for line in relvecs_to_sam_lines(*args, **kwargs):
            f.write(line)
