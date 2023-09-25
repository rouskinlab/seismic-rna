from typing import Sequence

import numpy as np

from .cigar import CIG_ALIGN, CIG_MATCH, CIG_SUBST, CIG_DELET, CIG_INSRT
from .cigarcount import CigarOp
from .encode import BASE_DECODINGS
from ..core.qual import LO_QUAL, HI_QUAL
from ..core.relvect import (MATCH, DELET, INS_5, INS_3, INS_8,
                            SUB_A, SUB_C, SUB_G, SUB_T, ANY_N,
                            IRREC, NOCOV, REL_TYPE)
from ..core.seq import BASEN, DNA


def find_relvec_ends(relvec: np.ndarray):
    """ Find the 5' and 3' ends of a relation vector. """
    if not isinstance(relvec, np.ndarray):
        raise TypeError(f"Expected {np.ndarray}, but got {type(relvec)}")
    if relvec.dtype.type is not REL_TYPE:
        raise TypeError(f"Expected an array of type {REL_TYPE}, but got "
                        f"type {relvec.dtype.type}")
    if relvec.ndim != 1:
        raise ValueError("Expected an array with 1 dimension, but got "
                         f"{relvec.ndim} dimensions")
    called = np.flatnonzero(relvec != NOCOV)
    if called.size == 0:
        raise ValueError("Relation vector is blank")
    end5 = int(called[0]) + 1
    end3 = int(called[-1]) + 1
    if (nexp := end3 - end5 + 1) != called.size:
        raise ValueError(f"Expected {nexp} base calls "
                         f"({np.arange(end5, end3 + 1)}), but got "
                         f"{called.size} ({called})")
    return end5, end3


def inverse_relate(refseq: DNA, relvec: np.ndarray,
                   hi_qual: str = HI_QUAL,
                   lo_qual: str = LO_QUAL,
                   ins_len: int | Sequence[int] = 1):
    """
    Infer a read sequence and quality string from a reference sequence
    and relation vector.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    relvec: ndarray
        Relation vector.
    hi_qual: str = MAX_QUAL
        Character to put in the quality string at every position that is
        high-quality according to the relation vector.
    lo_qual: str = MIN_QUAL
        Character to put in the quality string at every position that is
        low-quality according to the relation vector.
    ins_len: int | Sequence[int] = 1
        Number of bases to insert into the read and quality strings upon
        finding an insertion in the relation vector. If an integer, then
        insert that number of bases for every insertion. If a sequence
        of integers, then the ith insertion gets a number of bases equal
        to the ith element of ins_len.

    Returns
    -------

    """
    # Validate the relation vector and reference sequence.
    end5, end3 = find_relvec_ends(relvec)
    if len(refseq) != len(relvec):
        raise ValueError(f"Reference sequence has {len(refseq)} nt, "
                         f"but relation vector has {len(relvec)} nt")
    if hi_qual < lo_qual:
        raise ValueError(f"The high-quality score ({hi_qual}) is less than "
                         f"than the low-quality score ({lo_qual})")
    # Validate the quality codes.
    if len(hi_qual) != 1:
        raise ValueError(f"hi_qual must be 1 character, but got {len(hi_qual)}")
    if len(lo_qual) != 1:
        raise ValueError(f"lo_qual must be 1 character, but got {len(hi_qual)}")
    # Build the read sequence, quality scores, and CIGAR string one base
    # at a time.
    read: list[str] = list()
    qual: list[str] = list()
    cigars: list[CigarOp] = list()
    ins3_next = False
    ins_count = 0

    def add_to_cigar(op: str):
        """ Add one base of the relation vector to the CIGAR string. """
        if cigars and cigars[-1].op == op:
            # The current operation matches that of the last CigarOp:
            # just increment its length.
            cigars[-1].lengthen()
        else:
            # Either no CigarOp objects exist or the current operation
            # does not match that of the last one: start a new CigarOp.
            cigars.append(CigarOp(op))

    for pos in range(end5, end3 + 1):
        ref_base = refseq[pos - 1]
        rel = relvec[pos - 1]
        if rel == NOCOV:
            raise ValueError(f"Position {pos} in {end5}-{end3} is not covered")
        if rel & INS_8:
            # Specially handle insertions because they may overlap any
            # other relation except for a deletion.
            if rel & DELET:
                # Being both a deletion and an insertion is forbidden.
                raise ValueError(f"Relation {rel} in {relvec} is del and ins")
            if rel & INS_3:
                if not pos > end5:
                    # Insertions cannot occur before the beginning of
                    # the read.
                    raise ValueError(f"Position {pos} in {end5}-{end3} cannot "
                                     f"be 3' of an insertion")
                if not ins3_next:
                    raise ValueError(f"Unexpected 3' ins at {pos} in {relvec}")
                # The current position is 5' of an insertion, so the
                # inserted base must be added before the base at the
                # current position is added. Insert a number of bases
                # given by ins_len.
                if isinstance(ins_len, int):
                    # Add the same number of bases for every insertion.
                    n_ins = ins_len
                else:
                    # Add a number of bases for this specific insertion.
                    n_ins = ins_len[ins_count]
                read.append(BASEN * n_ins)
                qual.append(hi_qual * n_ins)
                for _ in range(n_ins):
                    add_to_cigar(CIG_INSRT)
                ins_count += 1
                # Being 3' of an insertion is not allowed until the next
                # position 5' of an insertion is reached.
                ins3_next = False
            elif ins3_next:
                # Check if this position should be 3' of an insertion.
                raise ValueError(f"Missing 3' ins at {pos} in {relvec}")
            if rel & INS_5:
                if not pos < end3:
                    raise ValueError(f"Position {pos} in {end5}-{end3} cannot "
                                     f"be 5' of an insertion")
                # The current position is 5' of an insertion, so the
                # inserted base must be added after the base at the
                # current position is added. Defer the responsibility of
                # adding the inserted base to the base that lies 3' of
                # the insertion and indicate that the next base must be
                # 3' of an insertion by setting ins3_next to True.
                ins3_next = True
            # Switch off the insertion flag(s) in the relation so that
            # the base at this position can be added below.
            rel = rel & ~INS_8
        elif ins3_next:
            # Check if this position should be 3' of an insertion.
            raise ValueError(f"Missing 3' ins at {pos} in {relvec}")
        if rel == MATCH:
            # Match: Add the reference base to the read.
            read.append(ref_base)
            qual.append(hi_qual)
            add_to_cigar(CIG_MATCH)
        elif rel == DELET:
            # Deletion from the read.
            if not end5 < pos < end3:
                raise ValueError(
                    f"Deletion cannot be at position {pos} in {end5}-{end3}")
            add_to_cigar(CIG_DELET)
        elif rel ^ ANY_N in (SUB_A, SUB_C, SUB_G, SUB_T, IRREC):
            # Ambiguous substitution: Add any nucleotide as low quality.
            read.append(BASEN)
            qual.append(lo_qual)
            add_to_cigar(CIG_ALIGN)
        else:
            # Unambiguous substitution: Add the substituted nucleotide
            # as high quality.
            try:
                read_base = BASE_DECODINGS[rel]
            except KeyError:
                raise ValueError(f"Invalid relation {rel} in {relvec}")
            if ref_base == read_base:
                raise ValueError(
                    f"Cannot substitute {ref_base} to itself in {relvec}")
            read.append(read_base)
            qual.append(hi_qual)
            add_to_cigar(CIG_SUBST)
    # Check that the 5' end was found.
    if end5 == 0:
        raise ValueError(f"Relation vector had no 5' end: {relvec}")
    if len(read) != len(qual):
        raise ValueError(
            f"Lengths of read ({len(read)}) and qual ({len(qual)}) differed")
    if not read:
        raise ValueError("Read contained no bases")
    # Assemble and return the read, quality, and CIGAR strings.
    read_seq = DNA("".join(read))
    qual_string = "".join(qual)
    cigar_string = "".join(map(str, cigars))
    return read_seq, qual_string, cigar_string, end5, end3

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
