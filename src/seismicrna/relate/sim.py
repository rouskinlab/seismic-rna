"""

Simulate Relation Vectors Module
========================================================================

"""

import numpy as np
import pandas as pd

from .c.encode import encode_match
from ..core.qual import HI_QUAL, LO_QUAL
from ..core.rand import rng
from ..core.rel import MATCH, DELET, ANY_N, SUB_A, SUB_C, SUB_G, SUB_T
from ..core.seq import BASEA, BASEC, BASEG, BASET, BASEN, POS_NAME, index_to_seq


def sim_relvecs(reads: pd.Index | np.ndarray | list | int,
                pmut: pd.DataFrame, *,
                pend: pd.DataFrame | None = None,
                ploq: pd.Series | None = None,
                ptri: float = 0.2,
                ptrv: float = 0.4,
                ptrn: float = 0.3):
    """
    Simulate relation vectors.

    Parameters
    ----------
    reads: pd.Index | np.ndarray | list | int
        Names of the reads or number of reads (if int).
    pmut: DataFrame
        Mutation rate of each position in each cluster.
    pend: DataFrame | None = None
        Probability that a read starts at the position in the index and
        ends at the position in the column. Only the upper triangle of
        the matrix is used.
    ploq: Series | None = None
        Probability that each position is low-quality.
    ptri: float = 0.2
        Probability that a mutation is a transition.
    ptrv: float = 0.4
        Probability that a mutation is a complementary transversion.
    ptrn: float = 0.3
        Probability that a mutation is a non-complementary transversion.

    Returns
    -------
    """
    npos, nclust = pmut.shape
    if ploq is not None:
        if not pmut.index.equals(ploq.index):
            raise ValueError(f"Indexes differ between pmut ({pmut.index}) "
                             f"and ploq ({ploq.index})")
    # Validate the positions and assemble the reference sequence.
    positions = pmut.index.get_level_values(POS_NAME)
    if not np.array_equal(positions, np.arange(1, npos + 1)):
        raise ValueError(f"Expected positions 1 - {npos}, but got "
                         f"{positions}")
    refseq = index_to_seq(pmut.index)
    # Validate the probabilities.
    if np.any(pmut < 0.) or np.any(pmut > 1.):
        raise ValueError(f"All pmut must be in [0, 1], but got {pmut}")

    if np.any(ploq < 0.) or np.any(ploq > 1.):
        raise ValueError(f"All ploq must be in [0, 1], but got {ploq}")
    # Validate ptri, ptrv, and ptrn.
    if not 0. <= ptri <= 1.:
        raise ValueError(f"ptri must be in [0, 1], but got {ptri}")
    if not 0. <= ptrv <= 1.:
        raise ValueError(f"ptrv must be in [0, 1], but got {ptrv}")
    if not 0. <= ptrn <= 1.:
        raise ValueError(f"ptrn must be in [0, 1], but got {ptrn}")
    # The sum of the probabilities cannot exceed 1.
    if (pivn := ptri + ptrv + ptrn) > 1.:
        raise ValueError(f"(ptri + ptrc + ptrn) must be ≤ 1, but got {pivn}")
    # Compute the probability that a mutation is a deletion.
    pdel = 1. - pivn
    # Initialize a DataFrame of blank relation vectors.
    relvecs = blank_relvec(refseq, reads)
    nreads, npos_ = relvecs.shape
    if npos_ != npos:
        raise ValueError(f"Inconsistent number of positions: {npos} ≠ {npos_}")
    if not relvecs.columns.equals(pmut.index):
        raise ValueError(
            f"Inconsistent positions: {pmut.index} ≠ {relvecs.columns}")
    if pend is None:
        # All positions in the relation vectors are covered.
        relvecs.loc[:, :] = MATCH
    else:
        # Validate the 5' and 3' end probabilities.
        if not pend.index.equals(pend.columns):
            raise ValueError(f"Got mismatched indexes ({pend.index}) and "
                             f"columns ({pend.columns}) of pend")
        if not pend.index.equals(pmut.index):
            raise ValueError(f"Indexes differ between pmut ({pmut.index}) "
                             f"and pend ({pend.index})")
        if np.any(pend < 0.) or np.any(pend > 1.):
            raise ValueError(f"All pend must be in [0, 1], but got {pend}")
        if np.sum(pend) > 1.:
            raise ValueError(
                f"The sum of pend must be in [0, 1], but got {np.sum(pend)}")
        # Choose 5' and 3' ends for the reads.


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
    pmut_types = pdel, ptrv, ptri, ptrn
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
