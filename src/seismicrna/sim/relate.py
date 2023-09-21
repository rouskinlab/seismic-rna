"""

Simulate Relation Vectors Module
========================================================================

"""

import numpy as np
import pandas as pd

from ..core.qual import HI_QUAL, LO_QUAL
from ..core.rand import rng
from ..core.rel import MATCH, DELET, ANY_N, SUB_A, SUB_C, SUB_G, SUB_T
from ..core.sect import BASE_NAME, POS_NAME
from ..core.seq import BASEA, BASEC, BASEG, BASET, BASEN, DNA
from ..relate.blank import blank_relvec
from ..relate.encode import encode_match


def sim_relvecs(refseq: DNA,
                reads: pd.Index | np.ndarray | list | int,
                ploq: pd.Series,
                pmut: pd.Series,
                ptri: float = 0.2,
                ptrc: float = 0.4,
                ptrn: float = 0.3):
    """
    Simulate relation vectors.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    reads: pd.Index | np.ndarray | list | int
        Names of the reads or number of reads (if int).
    ploq: Series
        Probability that each position is low-quality.
    pmut: Series
        Mutation rate of each position.
    ptri: float = 0.2
        Probability that a mutation is a transition.
    ptrc: float = 0.4
        Probability that a mutation is a complementary transversion.
    ptrn: float = 0.3
        Probability that a mutation is a non-complementary transversion.

    Returns
    -------
    """
    # Obtain the reference sequence as an array.
    seqarr = refseq.to_array()
    # Validate ploq and pmut.
    if not ploq.index.equals(pmut.index):
        raise ValueError(f"Indexes differ between ploq ({ploq.index}) "
                         f"and pmut ({pmut.index})")
    positions = ploq.index.get_level_values(POS_NAME)
    if positions.min() < 1 or positions.max() > len(refseq):
        raise ValueError(f"Positions must all be in [1, {len(refseq)}], "
                         f"but got {positions}")
    if not np.array_equal(seqarr[positions - 1],
                          ploq.index.get_level_values(BASE_NAME)):
        raise ValueError(f"Reference sequence {refseq} disagrees with index "
                         f"of probabilities {ploq.index}")
    if np.any(ploq < 0.) or np.any(ploq >= 1.):
        raise ValueError(f"All ploq must be in [0, 1), but got {ploq}")
    if np.any(pmut < 0.) or np.any(pmut >= 1.):
        raise ValueError(f"All pmut must be in [0, 1), but got {pmut}")
    # Validate ptri, ptrc, and ptrn.
    if not 0. <= ptri <= 1.:
        raise ValueError(f"ptri must be in [0, 1], but got {ptri}")
    if not 0. <= ptrc <= 1.:
        raise ValueError(f"ptrc must be in [0, 1], but got {ptrc}")
    if not 0. <= ptrn <= 1.:
        raise ValueError(f"ptrn must be in [0, 1], but got {ptrn}")
    # The sum of the probabilities cannot exceed 1.
    if (qdel := ptri + ptrc + ptrn) > 1.:
        raise ValueError(f"(ptri + ptrc + ptrn) must be ≤ 1, but got {qdel}")
    # Compute the probability that a mutation is a deletion.
    pdel = 1. - qdel
    # Initialize a DataFrame of blank relation vectors.
    relvecs = blank_relvec(refseq, reads)
    n_reads, n_pos = relvecs.shape
    # Initially, set all relationships within the section to matches.
    relvecs.loc[:, ploq.index] = MATCH
    # Simulate whether each position is low-quality. Force every N base
    # in the reference sequence to be low-quality.
    is_loq = np.logical_or(np.less(rng.random((n_reads, ploq.size)),
                                   ploq.values[np.newaxis, :]),
                           np.equal(seqarr, BASEN))
    # Find the indexes of the low-quality positions.
    loq_rows, loq_cols = np.nonzero(is_loq)
    n_loqs = loq_rows.size
    if n_loqs != loq_cols.size:
        raise ValueError(f"Counts differ for low-quality rows ({n_loqs}) "
                         f"and columns ({loq_cols.size})")
    # Mark every low-quality base.
    relvecs.values[loq_rows, loq_cols] = encode_match(BASEN, LO_QUAL, HI_QUAL)
    # Simulate whether each high-quality position is mutated.
    is_mut = np.logical_and(np.less(rng.random((n_reads, pmut.size)),
                                    pmut.values[np.newaxis, :]),
                            np.logical_not(is_loq))
    # Find the indexes of the mutated positions.
    mut_rows, mut_cols = np.nonzero(is_mut)
    n_muts = mut_rows.size
    if n_muts != mut_cols.size:
        raise ValueError(f"Counts differ for mutated rows ({n_muts}) "
                         f"and columns ({mut_cols.size})")
    # Simulate the type of each mutation.
    pmut_types = pdel, ptrc, ptri, ptrn
    is_mut_types = rng.choice(np.arange(len(pmut_types)), n_muts, p=pmut_types)
    # Mark every base with each type of mutation.
    mut_maps = [
        {BASEA: DELET, BASEC: DELET, BASEG: DELET, BASET: DELET, BASEN: DELET},
        {BASEA: SUB_T, BASEC: SUB_G, BASEG: SUB_C, BASET: SUB_A, BASEN: ANY_N},
        {BASEA: SUB_G, BASEC: SUB_T, BASEG: SUB_A, BASET: SUB_C, BASEN: ANY_N},
        {BASEA: SUB_C, BASEC: SUB_A, BASEG: SUB_T, BASET: SUB_G, BASEN: ANY_N},
    ]
    for mut_type, mut_map in enumerate(mut_maps):
        # Select only the mutations of the current type.
        is_mut_type = np.equal(is_mut_types, mut_type)
        mtype_rows = mut_rows[is_mut_type]
        mtype_cols = mut_cols[is_mut_type]
        if mut_type == 0:
            # The mutation type is a deletion, which is forbidden at the
            # first and last positions.
            is_not_end = np.logical_and(0 < mtype_cols, mtype_cols < n_pos - 1)
            mtype_rows = mtype_rows[is_not_end]
            mtype_cols = mtype_cols[is_not_end]
        # Mark each type of base.
        for base, mut in mut_map.items():
            # Select only the mutations where the reference base is the
            # current type of base.
            is_base = np.equal(seqarr[mtype_cols], base)
            # Mark those mutations.
            relvecs.values[mtype_rows[is_base], mtype_cols[is_base]] = mut
    return relvecs

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                              #
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
