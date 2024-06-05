from functools import cached_property
from typing import Callable

import numpy as np
import pandas as pd

from ..core.array import check_naturals, get_length
from ..core.batch import (AllReadBatch,
                          SectionMutsBatch,
                          simulate_muts,
                          simulate_segment_ends)
from ..core.seq import Section, index_to_pos, index_to_seq


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
        num_reads: int
            Number of reads in the batch.
        """
        check_naturals(index_to_pos(pmut.index), "positions")
        seg_end5s, seg_end3s = simulate_segment_ends(uniq_end5s,
                                                     uniq_end3s,
                                                     pends,
                                                     num_reads,
                                                     (read_length if paired
                                                      else 0),
                                                     p_rev)
        return cls(batch=batch,
                   section=Section(ref, index_to_seq(pmut.index)),
                   seg_end5s=seg_end5s,
                   seg_end3s=seg_end3s,
                   muts=simulate_muts(pmut, seg_end5s, seg_end3s),
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
