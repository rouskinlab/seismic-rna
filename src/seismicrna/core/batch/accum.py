from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd

from .ends import END_COORDS
from ..header import make_header
from ..logs import logger
from ..rel import RelPattern
from ..seq import DNA, seq_pos_to_index
from ..task import as_list_of_tuples, dispatch


def accumulate_counts(batch_counts: Iterable[tuple[Any, Any, Any, Any]],
                      refseq: DNA,
                      pos_nums: np.ndarray,
                      patterns: dict[str, RelPattern],
                      ks: Iterable[int] | None = None, *,
                      count_ends: bool = True,
                      count_pos: bool = True,
                      count_read: bool = True,
                      validate: bool = True):
    """ """
    logger.routine(f"Began accumulating patterns {patterns} in {batch_counts}")
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
    count_per_pos = (pd.DataFrame(zero,
                                  seq_pos_to_index(refseq, pos_nums, 1),
                                  header.index)
                     if count_pos else None)
    # Initialize the counts per read.
    count_per_batch_read = list() if count_read else None
    logger.detail(
        "\n".join(["Initialized accumulated",
                   f"num_reads = {num_reads}",
                   f"end_counts = {end_counts}",
                   f"count_per_pos = {count_per_pos}",
                   f"count_per_batch_read = {count_per_batch_read}"])
    )
    # Accumulate the counts from the batches.
    for (i, (num_reads_i,
             end_counts_i,
             count_per_pos_i,
             count_per_read_i)) in enumerate(batch_counts):
        if not isinstance(num_reads_i, type(num_reads)):
            raise TypeError(
                f"num_reads_i must be {type(num_reads).__name__}, "
                f"but got {type(num_reads_i).__name__}"
            )
        if (validate
                and isinstance(num_reads, pd.Series)
                and not num_reads_i.index.equals(num_reads.index)):
            raise ValueError("Got different indexes for "
                             f"num_reads_i ({num_reads_i.index}) "
                             f"and num_reads ({num_reads.index})")
        num_reads += num_reads_i
        if end_counts is not None:
            if not isinstance(end_counts_i, type(end_counts)):
                raise TypeError(
                    f"end_counts_i must be {type(end_counts).__name__}, "
                    f"but got {type(end_counts_i).__name__}"
                )
            if (validate
                    and isinstance(end_counts, pd.DataFrame)
                    and not end_counts_i.columns.equals(end_counts.columns)):
                raise ValueError("Got different columns for "
                                 f"end_counts_i ({end_counts_i.columns}) "
                                 f"and end_counts ({end_counts.columns})")
            end_counts = end_counts.add(end_counts_i,
                                        fill_value=zero).astype(dtype, copy=False)
        if count_per_pos is not None:
            if not isinstance(count_per_pos_i, pd.DataFrame):
                raise TypeError(f"count_per_pos_i must be DataFrame, "
                                f"but got {type(count_per_pos_i).__name__}")
            if (validate
                    and not count_per_pos_i.index.equals(count_per_pos.index)):
                raise ValueError(
                    "Got different indexes for "
                    f"count_per_pos_i ({count_per_pos_i.index}) "
                    f"and count_per_pos ({count_per_pos.index})"
                )
            if (validate and
                    not count_per_pos_i.columns.equals(count_per_pos.columns)):
                raise ValueError(
                    "Got different columns for "
                    f"count_per_pos_i ({count_per_pos_i.columns}) "
                    f"and count_per_pos ({count_per_pos.columns})"
                )
            count_per_pos += count_per_pos_i
        if count_per_batch_read is not None:
            if not isinstance(count_per_read_i, pd.DataFrame):
                raise TypeError(f"count_per_read_i must be DataFrame, "
                                f"but got {type(count_per_read_i).__name__}")
            if (validate
                    and not count_per_read_i.columns.equals(rel_header.index)):
                raise ValueError(
                    "Got different columns for "
                    f"count_per_read_i ({count_per_read_i.columns}) "
                    f"and header ({rel_header.index})"
                )
            count_per_batch_read.append(count_per_read_i)
        logger.detail(
            "\n".join([f"After batch {i}, accumulated",
                       f"num_reads = {num_reads}",
                       f"end_counts = {end_counts}",
                       f"count_per_pos = {count_per_pos}",
                       f"count_per_batch_read = {count_per_batch_read}"])
        )
    # Concatenate the per-read counts for the batches.
    if count_per_batch_read:
        count_per_read = pd.concat(count_per_batch_read, axis=0)
    elif count_per_batch_read is not None:
        count_per_read = pd.DataFrame(columns=rel_header.index, dtype=dtype)
    else:
        count_per_read = None
    logger.routine(f"Ended accumulating patterns {patterns} in {batch_counts}")
    return num_reads, end_counts, count_per_pos, count_per_read


def accumulate_batches(
        get_batch_count_all: Callable[[int], tuple[Any, Any, Any, Any]],
        num_batches: int,
        refseq: DNA,
        pos_nums: np.ndarray,
        patterns: dict[str, RelPattern],
        ks: Iterable[int] | None = None, *,
        count_ends: bool = True,
        count_pos: bool = True,
        count_read: bool = True,
        validate: bool = True,
        max_procs: int = 1
):
    # Generate the counts for the batches in parallel.
    batch_counts = dispatch(get_batch_count_all,
                            max_procs=max_procs,
                            pass_n_procs=False,
                            raise_on_error=True,
                            args=as_list_of_tuples(range(num_batches)),
                            kwargs=dict(patterns=patterns,
                                        ks=ks,
                                        count_ends=count_ends,
                                        count_pos=count_pos,
                                        count_read=count_read))
    # Accumulate the counts for all batches.
    return accumulate_counts(batch_counts,
                             refseq,
                             pos_nums,
                             patterns,
                             ks,
                             count_ends=count_ends,
                             count_pos=count_pos,
                             count_read=count_read,
                             validate=validate)

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
