"""

Read Iteration Module
========================================================================

------------------------------------------------------------------------

"""

from collections import defaultdict
from itertools import product

import numpy as np

from .cigar import CIG_INSRT
from .cigarcount import count_cigar_muts, find_cigar_op_pos
from .invert import inverse_relate
from .iterrelv import iter_relvecs_all
from ..core.relvect import INS_5, INS_3, REL_TYPE, NOCOV
from ..core.seq import DNA, expand_degenerate_seq


def ref_to_alignments(refseq: DNA,
                      max_ins: int = 2,
                      max_ins_len: int = 1,
                      max_ins_bases: int | None = None):
    """
    For a given reference sequence, return maps from every possible read
    to the CIGAR string(s) and (possibly ambiguous) relation vector for
    the read.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    max_ins: int = 2
        Maximum number of insertions in the read. Must be ≥ 0.
    max_ins_len: int = 1
        Maximum length of (i.e. number of bases in) one insertion.
        Must be ≥ 1.
    max_ins_bases: int | None = None
        Maximum total number of bases inserted. Must be ≥ `max_ins`.
        If `None`, there is no limit.
    """
    # Initialize maps of reads to CIGAR strings and relation vectors.
    quals = dict()
    cigars = defaultdict(lambda: defaultdict(list))
    relvecs = defaultdict(lambda: defaultdict(lambda: np.zeros(len(refseq),
                                                               dtype=REL_TYPE)))
    if max_ins < 0:
        raise ValueError(f"max_ins must be ≥ 0, but got {max_ins}")
    if max_ins > 0:
        if max_ins_len < 1:
            raise ValueError(f"max_ins_len must be ≥ 1, but got {max_ins_len}")
        if max_ins_bases is not None and max_ins_bases < max_ins:
            raise ValueError(f"max_ins_bases ({max_ins_bases}) "
                             f"must be ≥ max_ins ({max_ins})")
    # Iterate through all possible relation vectors.
    for relvec in iter_relvecs_all(refseq, max_ins):
        # Check if there are insertions in the relation vector.
        n_ins = int(np.count_nonzero(np.logical_and(relvec & INS_5,
                                                    relvec != NOCOV)))
        n_ins3 = int(np.count_nonzero(np.logical_and(relvec & INS_3,
                                                     relvec != NOCOV)))
        if n_ins != n_ins3:
            raise ValueError(f"Got {n_ins} 5' and {n_ins3} 3' insertions")
        if n_ins > max_ins:
            # This should not happen. Checking just in case.
            raise ValueError(f"Got {n_ins} insertions, but max_ins = {max_ins}")
        # Iterate over all possible combinations of insertion lengths.
        for ins_len in product(range(1, max_ins_len + 1), repeat=n_ins):
            if max_ins_bases is not None and sum(ins_len) > max_ins_bases:
                # Skip insertion lengths whose sum exceeds the limit.
                continue
            # Determine the read(s) corresponding to this relation vector.
            degen, qual, cigar, end5, end3 = inverse_relate(refseq, relvec,
                                                            ins_len=ins_len)
            if n_ins > 0:
                # If there are insertions, find their positions.
                ins_pos = list(find_cigar_op_pos(cigar, CIG_INSRT))
                # Remove quality codes of inserted bases because -- for
                # the purpose of aggregating reads based on sequence,
                # quality score, and position -- the "fake" quality
                # scores of inserted bases should not be considered.
                qual_no_ins = "".join(q for i, q in enumerate(qual, start=1)
                                      if i not in ins_pos)
            else:
                qual_no_ins = qual
            # Count the mutations in the CIGAR string.
            nmuts = count_cigar_muts(cigar)
            for read in expand_degenerate_seq(degen):
                key = read, qual_no_ins, end5, end3
                # Record the original quality score.
                quals[key] = qual
                # Gather every CIGAR string for the read.
                cigars[key][nmuts].append(cigar)
                # Accumulate the bitwise OR of all relation vectors.
                relvecs[key][nmuts] |= relvec
    # For every read-quality-end5-end3 key, keep only the CIGAR strings
    # and relation vector with the fewest mutations.
    cigars_best: dict[tuple[DNA, str, int, int], list[str]] = dict()
    relvec_best: dict[tuple[DNA, str, int, int], np.ndarray] = dict()
    for key, cig_key in cigars.items():
        # Find the minimum number of mutations for this key.
        min_nmuts = min(cig_key)
        # Export the results with the fewest mutations.
        cigars_best[key] = cigars[key][min_nmuts]
        relvec_best[key] = relvecs[key][min_nmuts]
    return quals, cigars_best, relvec_best


def iter_alignments(*args, **kwargs):
    """ For a given reference sequence, find every read that could come
    from the reference (with up to 2 bases inserted). For each read,
    yield the (possibly ambiguous) relation vector and every possible
    CIGAR string. """
    quals, cigars, relvecs = ref_to_alignments(*args, **kwargs)
    for key, relvec in relvecs.items():
        read, _, end5, end3 = key
        qual = quals[key]
        for cigar in cigars[key]:
            yield read, qual, cigar, end5, end3, relvec

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
