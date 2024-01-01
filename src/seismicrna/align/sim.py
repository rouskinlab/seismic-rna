"""

Simulate SAM Files Module
========================================================================

"""

from os import linesep
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from ..core.ngs import (HI_QUAL,
                        LO_QUAL,
                        FLAG_PAIRED,
                        FLAG_PROPER,
                        FLAG_FIRST,
                        FLAG_SECOND,
                        FLAG_REVERSE,
                        FLAG_MREVERSE,
                        MAX_FLAG,
                        SAM_DELIM,
                        SAM_NOREF,
                        SAM_SEQLEN,
                        SAM_SEQLINE,
                        SAM_SEQNAME)
from ..core.rel import NOCOV
from ..core.seq import DNA, index_to_seq
from ..relate.aux.infer import infer_read


def sam_header(ref: str, length: int | DNA):
    if isinstance(length, DNA):
        length = len(length)
    if not isinstance(length, int):
        raise TypeError(f"Length must be int, but got {type(length).__name__}")
    return SAM_DELIM.join((SAM_SEQLINE,
                           f"{SAM_SEQNAME}{ref}",
                           f"{SAM_SEQLEN}{length}{linesep}"))


def as_sam(name: str, flag: int, ref: str, end5: int, mapq: int, cigar: str,
           rnext: str, pnext: int, tlen: int, read: DNA, qual: str):
    """
    Return a line in SAM format from the given fields.

    Parameters
    ----------
    name: str
        Name of the read.
    flag: int
        SAM flag. Must be in [0, MAX_FLAG].
    ref: str
        Name of the reference.
    end5: int
        Most 5' position to which the read mapped (1-indexed).
    mapq: int
        Mapping quality score.
    cigar: str
        CIGAR string. Not checked for compatibility with the read.
    rnext: str
        Name of the mate's reference (if paired-end).
    pnext: int
        Most 5' position of the mate (if paired-end).
    tlen: int
        Length of the template.
    read: DNA
        Base calls in the read. Must be equal in length to `read`.
    qual: str
        Phred quality score string of the base calls. Must be equal in
        length to `read`.

    Returns
    -------
    str
        A line in SAM format containing the given fields.
    """
    if not name:
        raise ValueError("Read name is empty")
    if not 0 <= flag <= MAX_FLAG:
        raise ValueError(f"Invalid SAM flag: {flag}")
    if not ref:
        raise ValueError("Reference name is empty")
    if not end5 >= 1:
        raise ValueError(f"Invalid 5' mapping position: {end5}")
    if not cigar:
        raise ValueError("CIGAR string is empty")
    if not rnext:
        raise ValueError("Next reference name is empty")
    if not pnext >= 0:
        raise ValueError(f"Invalid next 5' mapping position: {pnext}")
    if not len(read) == len(qual):
        raise ValueError(
            f"Lengths of read ({len(read)}) and qual ({len(qual)}) disagree")
    return SAM_DELIM.join(map(str, (name, flag, ref, end5, mapq, cigar, rnext,
                                    pnext, tlen, read, f"{qual}{linesep}")))


def _relvec_to_sam_line(read: str,
                        relvec: np.ndarray,
                        ref: str,
                        refseq: DNA, *,
                        mapq: int,
                        flag: int = 0,
                        hi_qual: str = HI_QUAL,
                        lo_qual: str = LO_QUAL,
                        ins_len: int | Sequence[int] = 1):
    seq, qual, cig, end5, end3 = infer_read(refseq,
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
    end5, end3 = find_relvec_ends(relvec)
    # Determine whether the first read is reversed.
    rev1 = bool(np.random.default_rng().integers(2))
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

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
