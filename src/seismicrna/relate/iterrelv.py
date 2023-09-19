"""

Relate Auxiliary Module
========================================================================

------------------------------------------------------------------------

"""

from itertools import (product, combinations,
                       combinations_with_replacement as cwr)
from typing import Sequence

import numpy as np

from .encode import encode_relate
from ..core.rel import DELET, INS_5, INS_3, NP_TYPE, NOCOV
from ..core.sect import Section
from ..core.seq import BASEA, BASEC, BASEG, BASET, DNA

MIN_INDEL_GAP = 1


def iter_relvecs_q53(refseq: DNA, low_qual: Sequence[int] = (),
                     end5: int | None = None, end3: int | None = None,
                     max_ins: int = 2):
    """
    For a given reference sequence, yield every possible unambiguous
    relation vector between positions end5 and end3 that follows the
    given low-quality positions and has at most two insertions.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    low_qual: Sequence[int]
        List of positions in the read that are low-quality.
    end5: int | None = None
        5' end of the read; 1-indexed with respect to `refseq`.
    end3: int | None = None
        3' end of the read; 1-indexed with respect to `refseq`.
    max_ins: int = 2
        Maximum number of insertions in the read. Must be in [0, 2].
    """
    low_qual = set(low_qual)
    # Determine the section of the reference sequence that is occupied
    # by the read.
    sect = Section("", refseq, end5=end5, end3=end3)
    if low_qual - set(sect.range_int):
        raise ValueError(f"Invalid positions in low_qual: "
                         f"{sorted(low_qual - set(sect.range))}")
    if max_ins not in range(3):
        raise ValueError(f"Invalid value for max_ins: {max_ins}")
    # Find the possible relationships at each position in the section,
    # not including insertions.
    rel_opts: list[tuple[int, ...]] = [(NOCOV,)] * len(refseq)
    for pos in sect.range_int:
        # Find the base in the reference sequence (pos is 1-indexed).
        ref_base = refseq[pos - 1]
        if low_qual:
            # The only option is low-quality.
            opts = encode_relate(ref_base, ref_base, '0', '1'),
        else:
            # The options are a match and three substitutions.
            opts = tuple(encode_relate(ref_base, read_base, '1', '1')
                         for read_base in (BASEA, BASEC, BASEG, BASET))
        # A deletion is an option at every position except the ends of
        # the covered region.
        if sect.end5 < pos < sect.end3:
            opts = opts + (DELET,)
        rel_opts[pos - 1] = opts
    # Iterate through all possible relationships at each position.
    margin = MIN_INDEL_GAP + 1
    for rel in product(*rel_opts):
        # Generate a relation vector from the relationships.
        relvec = np.array(rel, dtype=NP_TYPE)
        yield relvec
        if max_ins == 0:
            # Do not introduce any insertions into the read.
            continue
        # Allow up to two insertions in one read, at every position
        # except the ends of the covered region and near deletions.
        for ins1_pos5, ins2_pos5 in cwr(range(sect.end5, sect.end3), 2):
            if max_ins == 1 and ins1_pos5 != ins2_pos5:
                # If at most one insertion can be in the read, then skip
                # cases where ins1_pos5 and ins2_pos5 are not equal.
                continue
            # Skip insertions nearby any deletion.
            mrg15 = max(ins1_pos5 - margin, 0)
            mrg13 = min(ins1_pos5 + margin, relvec.size)
            if np.any(relvec[mrg15: mrg13] == DELET):
                continue
            mrg25 = max(ins2_pos5 - margin, 0)
            mrg23 = min(ins2_pos5 + margin, relvec.size)
            if np.any(relvec[mrg25: mrg23] == DELET):
                continue
            # Yield the relation vector with these insertions.
            relvec_ins = relvec.copy()
            if ins1_pos5 >= 1:
                relvec_ins[ins1_pos5 - 1] |= INS_5
            if ins1_pos5 < relvec.size:
                relvec_ins[ins1_pos5] |= INS_3
            if ins2_pos5 >= 1:
                relvec_ins[ins2_pos5 - 1] |= INS_5
            if ins2_pos5 < relvec.size:
                relvec_ins[ins2_pos5] |= INS_3
            yield relvec_ins


def iter_relvecs_all(refseq: DNA, max_ins: int = 2):
    """
    For a given reference sequence, yield every possible unambiguous
    relation vector that has at most two insertions.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    max_ins: int = 2
        Maximum number of insertions in a read.
    """
    # Use every possible pair of 5' and 3' end positions.
    for end5, end3 in cwr(range(1, len(refseq) + 1), 2):
        # Use every possible number of low-quality base calls.
        for n_low_qual in range((end3 - end5 + 1) + 1):
            # Use every possible set of low-quality positions.
            for low_qual in combinations(range(end5, end3 + 1), n_low_qual):
                # Yield every possible relation vector.
                yield from iter_relvecs_q53(refseq, low_qual, end5, end3,
                                            max_ins)
