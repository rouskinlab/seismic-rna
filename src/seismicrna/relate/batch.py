from abc import ABC, abstractmethod
from functools import cached_property
from typing import Callable, Self

import numpy as np
import pandas as pd

from ..core.array import calc_inverse, check_naturals, get_length
from ..core.batch import (ReadBatch,
                          MutsBatch,
                          RegionMutsBatch,
                          simulate_muts,
                          simulate_segment_ends)
from ..core.rel import RelPattern
from ..core.seq import Region, index_to_pos, index_to_seq


def format_read_name(batch_num: int, read_num: int):
    """ Format a read name. """
    return f"batch-{batch_num}_read-{read_num}"


class FullReadBatch(ReadBatch, ABC):

    @classmethod
    @abstractmethod
    def simulate(cls, *args, **kwargs) -> Self:
        """ Simulate a batch. """

    @cached_property
    def read_nums(self):
        return np.arange(self.num_reads, dtype=self.read_dtype)

    @cached_property
    def max_read(self):
        return self.num_reads - 1

    @property
    def read_indexes(self):
        return self.read_nums


class ReadNamesBatch(FullReadBatch):

    @classmethod
    def simulate(cls,
                 branches: dict[str, str],
                 batch: int,
                 num_reads: int,
                 formatter: Callable[[int, int], str] = format_read_name,
                 **kwargs):
        """ Simulate a batch.

        Parameters
        ----------
        branches: dict[str, str]
            Branches of the workflow.
        batch: int
            Batch number.
        num_reads: int
            Number of reads in the batch.
        formatter: Callable[[int, int], str]
            Function to generate the name of each read: must accept the
            batch number and the read number and return a string.
        """
        return cls(branches=branches,
                   batch=batch,
                   names=[formatter(batch, read) for read in range(num_reads)],
                   **kwargs)

    def __init__(self, *, names: list[str] | np.ndarray, **kwargs):
        super().__init__(**kwargs)
        self.names = np.asarray(names, dtype=str)

    @cached_property
    def num_reads(self):
        return get_length(self.names, "read names")


class RelateMutsBatch(FullReadBatch, MutsBatch, ABC):

    @property
    def read_weights(self):
        read_weights = None
        if self.masked_reads_bool.any():
            read_weights = np.ones(self.num_reads)
            read_weights[self.masked_reads_bool] = 0
            read_weights = pd.DataFrame(read_weights)
        return read_weights


class RelateRegionMutsBatch(RelateMutsBatch, RegionMutsBatch):

    @classmethod
    def simulate(cls,
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
        region = Region(ref, index_to_seq(pmut.index))
        # Simulate a batch, ignoring min_mut_gap.
        seg_end5s, seg_end3s = simulate_segment_ends(uniq_end5s,
                                                     uniq_end3s,
                                                     pends,
                                                     num_reads,
                                                     (read_length if paired
                                                      else 0),
                                                     p_rev)
        simulated_all = cls(region=region,
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
        return cls(region=region,
                   seg_end5s=seg_end5s[reads_noclose],
                   seg_end3s=seg_end3s[reads_noclose],
                   muts={pos: {rel: renumber[np.setdiff1d(reads,
                                                          reads_exclude,
                                                          assume_unique=True)]
                               for rel, reads in rels.items()}
                         for pos, rels in simulated_all.muts.items()},
                   **kwargs)
