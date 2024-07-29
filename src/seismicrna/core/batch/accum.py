from typing import Iterable

import numpy as np
import pandas as pd

from .ends import END_COORDS
from .muts import SectionMutsBatch
from ..header import make_header
from ..rel import RelPattern
from ..seq import DNA, seq_pos_to_index


def _add_to_rel(added: pd.Series | pd.DataFrame, frame: pd.DataFrame, rel: str):
    """ Add the values in `added` to the column `rel` of `frame`. """
    frame_rel = frame[rel]
    if not frame_rel.index.equals(added.index):
        raise ValueError(f"Got different indexes for frame {frame_rel.index} "
                         f"and added values {added.index}")
    if (isinstance(added, pd.DataFrame)
            and not frame_rel.columns.equals(added.columns)):
        raise ValueError(f"Got different columns for frame {frame_rel.columns} "
                         f"and added values {added.columns}")
    frame[rel] = (frame_rel + added).values


def accumulate(batches: Iterable[SectionMutsBatch],
               refseq: DNA,
               patterns: dict[str, RelPattern], *,
               ks: Iterable[int] | None = None,
               pos_nums: np.ndarray | None = None,
               per_read: bool = True,
               get_info: bool = True,
               count_ends: bool = True):
    header = make_header(rels=list(patterns), ks=ks)
    end_counts_index = pd.MultiIndex.from_arrays([np.array([], dtype=int)
                                                  for _ in END_COORDS],
                                                 names=END_COORDS)
    # Initialize the total read counts and end coordinate counts.
    if header.clustered():
        dtype = float
        rel_header = header.get_rel_header()
        clust_index = header.get_clust_header().index
        num_reads = pd.Series(0, index=clust_index)
        end_counts = (pd.DataFrame(index=end_counts_index,
                                   columns=clust_index,
                                   dtype=dtype)
                      if count_ends else None)
    else:
        dtype = int
        rel_header = header
        num_reads = 0
        end_counts = (pd.Series(index=end_counts_index,
                                dtype=dtype)
                      if count_ends else None)
    zero = dtype(0)
    # Initialize the counts per position.
    if pos_nums is not None:
        index_per_pos = seq_pos_to_index(refseq, pos_nums, 1)
        fits_per_pos = pd.DataFrame(zero, index_per_pos, header.index)
        if get_info:
            info_per_pos = pd.DataFrame(zero, index_per_pos, header.index)
        else:
            info_per_pos = None
    else:
        fits_per_pos = None
        info_per_pos = None
    # Initialize the counts per read.
    if per_read:
        fits_per_read_per_batch = list()
        if get_info:
            info_per_read_per_batch = list()
        else:
            info_per_read_per_batch = None
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
        # Count the end coordinates.
        if end_counts is not None:
            end_counts = end_counts.add(batch.read_end_counts,
                                        fill_value=zero).astype(dtype,
                                                                copy=False)
        # Count the positions and/or reads matching each pattern.
        if fits_per_read_per_batch is not None:
            fits_per_read_per_batch.append(pd.DataFrame(0,
                                                        batch.batch_read_index,
                                                        rel_header.index))
        if info_per_read_per_batch is not None:
            info_per_read_per_batch.append(pd.DataFrame(0,
                                                        batch.batch_read_index,
                                                        rel_header.index))
        for column, pattern in patterns.items():
            if fits_per_pos is not None:
                # Count the matching reads per position.
                ipp, fpp = batch.count_per_pos(pattern)
                _add_to_rel(fpp, fits_per_pos, column)
                if info_per_pos is not None:
                    _add_to_rel(ipp, info_per_pos, column)
            if fits_per_read_per_batch is not None:
                # Count the matching positions per read.
                ipr, fpr = batch.count_per_read(pattern)
                fits_per_read_per_batch[-1].loc[:, column] = fpr.values
                if info_per_read_per_batch is not None:
                    info_per_read_per_batch[-1].loc[:, column] = ipr.values

    def get_data_per_read(data_per_read_per_batch: pd.DataFrame | None):
        if data_per_read_per_batch is not None:
            if data_per_read_per_batch:
                # Concatenate the per-read counts for the batches.
                return pd.concat(data_per_read_per_batch, axis=0)
            return pd.DataFrame(columns=rel_header.index, dtype=dtype)
        return None

    return (num_reads,
            (fits_per_pos,
             info_per_pos),
            (get_data_per_read(fits_per_read_per_batch),
             get_data_per_read(info_per_read_per_batch)),
            end_counts)


def accum_per_pos(batches: Iterable[SectionMutsBatch],
                  refseq: DNA,
                  pos_nums: np.ndarray,
                  patterns: dict[str, RelPattern],
                  ks: Iterable[int] | None = None):
    """ Count reads with each relationship at each position in a section
    over multiple batches. """
    num_reads, (fpp, ipp), (_, __), ___ = accumulate(batches,
                                                     refseq,
                                                     patterns,
                                                     ks=ks,
                                                     pos_nums=pos_nums,
                                                     per_read=False,
                                                     count_ends=False)
    return num_reads, fpp, ipp


def accum_fits(batches: Iterable[SectionMutsBatch],
               refseq: DNA,
               pos_nums: np.ndarray,
               patterns: dict[str, RelPattern],
               ks: Iterable[int] | None = None):
    """ Count positions and reads fitting each relationship. """
    num_reads, (fpp, _), (fpr, __), ends = accumulate(batches,
                                                      refseq,
                                                      patterns,
                                                      ks=ks,
                                                      pos_nums=pos_nums,
                                                      get_info=False)
    return num_reads, fpp, fpr, ends

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
