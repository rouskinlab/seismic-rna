from functools import cached_property
from typing import Callable

import numpy as np
import pandas as pd

from ..core.array import calc_inverse, check_naturals, get_length
from ..core.batch import (AllReadBatch,
                          SectionMutsBatch,
                          simulate_muts,
                          simulate_segment_ends)
from ..core.seq import Section, index_to_pos, index_to_seq
from ..core.rel import RelPattern


def format_read_name(batch: int, read: int):
    """ Format a read name. """
    return f"batch-{batch}_read-{read}"


class QnamesBatch(AllReadBatch):

    @classmethod
    def simulate(cls,
                 batch: int,
                 num_reads: int,
                 formatter: Callable[[int, int], str] = format_read_name,
                 **kwargs):
        """ Simulate a batch.

        Parameters
        ----------
        batch: int
            Batch number.
        num_reads: int
            Number of reads in the batch.
        formatter: Callable[[int, int], str]
            Function to generate the name of each read: must accept the
            batch number and the read number and return a string.
        """
        return cls(batch=batch,
                   names=[formatter(batch, read) for read in range(num_reads)],
                   **kwargs)

    def __init__(self, *, names: list[str] | np.ndarray, **kwargs):
        super().__init__(**kwargs)
        self.names = np.asarray(names, dtype=str)

    @cached_property
    def num_reads(self):
        return get_length(self.names, "read names")


class RelateBatch(SectionMutsBatch, AllReadBatch):

    @classmethod
    def simulate(cls,
                 batch: int,
                 ref: str,
                 pmut: pd.DataFrame,
                 uniq_end5s: np.ndarray,
                 uniq_end3s: np.ndarray,
                 pends: np.ndarray,
                 paired: bool,
                 read_length: int,
                 p_rev: float,
                 min_mut_gap: int,
                 num_reads: int,
                 **kwargs):
        """ Simulate a batch.

        Parameters
        ----------
        batch: int
            Batch number.
        ref: str
            Name of the reference.
        pmut: pd.DataFrame
            Rate of each type of mutation at each position.
        uniq_end5s: np.ndarray
            Unique read 5' end coordinates.
        uniq_end3s: np.ndarray
            Unique read 3' end coordinates.
        pends: np.ndarray
            Probability of each set of unique end coordinates.
        paired: bool
            Whether to simulate paired-end or single-end reads.
        read_length: int
            Length of each read segment (paired-end reads only).
        p_rev: float
            Probability that mate 1 is reversed (paired-end reads only).
        min_mut_gap: int
            Minimum number of positions between two mutations.
        num_reads: int
            Number of reads in the batch.
        """
        check_naturals(index_to_pos(pmut.index), "positions")
        section = Section(ref, index_to_seq(pmut.index))
        # Simulate a batch, ignoring min_mut_gap.
        seg_end5s, seg_end3s = simulate_segment_ends(uniq_end5s,
                                                     uniq_end3s,
                                                     pends,
                                                     num_reads,
                                                     (read_length if paired
                                                      else 0),
                                                     p_rev)
        # Drop any reads with zero coverage.
        has_coverage = np.less_equal(seg_end5s, seg_end3s).any(axis=1)
        seg_end5s = seg_end5s[has_coverage]
        seg_end3s = seg_end3s[has_coverage]
        simulated_all = cls(batch=batch,
                            section=section,
                            seg_end5s=seg_end5s,
                            seg_end3s=seg_end3s,
                            muts=simulate_muts(pmut, seg_end5s, seg_end3s),
                            **kwargs)
        if min_mut_gap == 0:
            return simulated_all
        # Remove reads with two mutations too close.
        reads_noclose = simulated_all.reads_noclose_muts(RelPattern.muts(),
                                                         min_mut_gap)
        reads_exclude = np.setdiff1d(simulated_all.read_nums,
                                     reads_noclose,
                                     assume_unique=True)
        renumber = calc_inverse(reads_noclose, what="reads_noclose")
        return cls(batch=batch,
                   section=section,
                   seg_end5s=seg_end5s[reads_noclose],
                   seg_end3s=seg_end3s[reads_noclose],
                   muts={pos: {rel: renumber[np.setdiff1d(reads,
                                                          reads_exclude,
                                                          assume_unique=True)]
                               for rel, reads in rels.items()}
                         for pos, rels in simulated_all.muts.items()},
                   **kwargs)

    @property
    def read_weights(self):
        return None

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
