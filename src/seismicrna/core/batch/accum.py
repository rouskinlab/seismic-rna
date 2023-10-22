from typing import Iterable

import numpy as np
import pandas as pd

from .index import REL_NAME, add_to_rel, get_rel_index
from .muts import MutsBatch
from ..rel import RelPattern
from ..seq import POS_INDEX, DNA, seq_pos_to_index


def accumulate(batches: Iterable[MutsBatch],
               refseq: DNA,
               patterns: dict[str, RelPattern],
               clusters: Iterable[tuple[int, int]] | None = None, *,
               pos_nums: np.ndarray | None = None,
               per_read: bool = True,
               get_info: bool = True):
    dtype = float if clusters is not None else int
    poc_index = get_rel_index(patterns, clusters)
    # Initialize the total read counts.
    if clusters is not None:
        num_reads_index = poc_index.droplevel(REL_NAME).drop_duplicates()
        num_reads = pd.Series(0, index=num_reads_index)
    else:
        num_reads = 0
    # Initialize the counts per position.
    if pos_nums is not None:
        index_per_pos = seq_pos_to_index(refseq, pos_nums, POS_INDEX)
        fits_per_pos = pd.DataFrame(dtype(0), index_per_pos, poc_index)
        info_per_pos = (pd.DataFrame(dtype(0), index_per_pos, poc_index)
                        if get_info else None)
    else:
        fits_per_pos = None
        info_per_pos = None
    # Initialize the counts per read.
    if per_read:
        fits_per_read_per_batch = list()
        info_per_read_per_batch = list() if get_info else None
    else:
        fits_per_read_per_batch = None
        info_per_read_per_batch = None
    # Accumulate the counts from the batches.
    for batch in batches:
        if pos_nums is not None and not np.array_equal(batch.pos_nums,
                                                       pos_nums):
            raise ValueError(f"Positions of {batch} ({batch.pos_nums}) "
                             f"differ from pos_nums ({pos_nums})")
        # Count the total number of reads in the batch.
        num_reads += batch.num_reads
        # Count the positions and/or reads matching each pattern.
        if fits_per_read_per_batch is not None:
            fits_per_read_per_batch.append(pd.DataFrame(dtype(0),
                                                        batch.batch_read_index,
                                                        poc_index))
        if info_per_read_per_batch is not None:
            info_per_read_per_batch.append(pd.DataFrame(dtype(0),
                                                        batch.batch_read_index,
                                                        poc_index))
        for column, pattern in patterns.items():
            if fits_per_pos is not None:
                # Count the matching reads per position.
                ipp, fpp = batch.count_per_pos(refseq, pattern)
                add_to_rel(fits_per_pos, fpp, column)
                if info_per_pos is not None:
                    add_to_rel(info_per_pos, ipp, column)
            if fits_per_read_per_batch is not None:
                # Count the matching positions per read.
                ipr, fpr = batch.count_per_read(refseq, pattern)
                fits_per_read_per_batch[-1].loc[:, column] = fpr.values
                if info_per_read_per_batch is not None:
                    info_per_read_per_batch[-1].loc[:, column] = ipr.values

    def get_data_per_read(data_per_read_per_batch: pd.DataFrame | None):
        if data_per_read_per_batch is not None:
            if data_per_read_per_batch:
                # Concatenate the per-read counts for the batches.
                return pd.concat(data_per_read_per_batch, axis=0)
            return pd.DataFrame(columns=poc_index, dtype=dtype)
        return None

    return (num_reads,
            (fits_per_pos,
             info_per_pos),
            (get_data_per_read(fits_per_read_per_batch),
             get_data_per_read(info_per_read_per_batch)))


def accum_per_pos(batches: Iterable[MutsBatch],
                  refseq: DNA,
                  pos_nums: np.ndarray,
                  patterns: dict[str, RelPattern],
                  clusters: Iterable[tuple[int, int]] | None = None,
                  ) -> tuple[int | pd.Series, pd.DataFrame, pd.DataFrame]:
    """ Count reads with each relationship at each position in a section
    over multiple batches. """
    num_reads, (fpp, ipp), _ = accumulate(batches,
                                          refseq,
                                          patterns,
                                          clusters,
                                          pos_nums=pos_nums,
                                          per_read=False)
    return num_reads, fpp, ipp


def accum_fits(batches: Iterable[MutsBatch],
               refseq: DNA,
               pos_nums: np.ndarray,
               patterns: dict[str, RelPattern],
               clusters: Iterable[tuple[int, int]] | None = None,
               ) -> tuple[int | pd.Series, pd.DataFrame, pd.DataFrame]:
    """ Count positions and reads fitting each relationship. """
    num_reads, (fpp, _), (fpr, __) = accumulate(batches,
                                                refseq,
                                                patterns,
                                                clusters,
                                                pos_nums=pos_nums,
                                                get_info=False)
    return num_reads, fpp, fpr
