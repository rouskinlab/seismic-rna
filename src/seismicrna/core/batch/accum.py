from typing import Any, Iterable

import numpy as np
import pandas as pd

from .ends import END_COORDS
from .muts import SectionMutsBatch
from ..header import make_header
from ..logs import logger
from ..rel import RelPattern
from ..seq import DNA, seq_pos_to_index


def accumulate_counts(batch_counts: Iterable[tuple[Any, Any, Any, Any]],
                      refseq: DNA,
                      pos_nums: np.ndarray,
                      patterns: dict[str, RelPattern],
                      ks: Iterable[int] | None = None,
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
        end_counts = pd.DataFrame(index=end_counts_index,
                                  columns=clust_index,
                                  dtype=dtype)
    else:
        dtype = int
        rel_header = header
        num_reads = 0
        end_counts = pd.Series(index=end_counts_index,
                               dtype=dtype)
    zero = dtype(0)
    # Initialize the counts per position.
    count_per_pos = pd.DataFrame(zero,
                                 seq_pos_to_index(refseq, pos_nums, 1),
                                 header.index)
    # Initialize the counts per read.
    count_per_batch_read = list()
    logger.detail(
        "\n".join(["Initialized accumulated",
                   f"num_reads = {num_reads}",
                   f"end_counts = {end_counts}",
                   f"count_per_pos = {count_per_pos}",
                   f"count_per_batch_read = {count_per_batch_read}"])
    )
    # Accumulate the counts from the batches.
    for (i, (batch_num_reads,
             batch_end_counts,
             batch_count_per_pos,
             batch_count_per_read)) in enumerate(batch_counts):
        if not isinstance(batch_num_reads, type(num_reads)):
            raise TypeError(
                f"batch_num_reads must be {type(num_reads).__name__}, "
                f"but got {type(batch_num_reads.__name__)}"
            )
        if (validate
                and isinstance(num_reads, pd.Series)
                and not batch_num_reads.index.equals(num_reads.index)):
            raise ValueError("Got different indexes for "
                             f"batch_num_reads ({batch_num_reads.index}) "
                             f"and num_reads ({num_reads.index})")
        num_reads += batch_num_reads
        if not isinstance(batch_end_counts, type(end_counts)):
            raise TypeError(
                f"batch_end_counts must be {type(end_counts).__name__}, "
                f"but got {type(batch_end_counts.__name__)}"
            )
        if (validate
                and isinstance(end_counts, pd.DataFrame)
                and not batch_end_counts.columns.equals(end_counts.columns)):
            raise ValueError("Got different columns for "
                             f"batch_end_counts ({batch_end_counts.columns}) "
                             f"and end_counts ({end_counts.columns})"
                             )
        end_counts = end_counts.add(batch_end_counts,
                                    fill_value=zero).astype(dtype, copy=False)
        if not isinstance(batch_count_per_pos, pd.DataFrame):
            raise TypeError(f"batch_count_per_pos must be DataFrame, "
                            f"but got {type(batch_count_per_pos.__name__)}")
        if (validate
                and not batch_count_per_pos.index.equals(count_per_pos.index)):
            raise ValueError(
                "Got different indexes for "
                f"batch_count_per_pos ({batch_count_per_pos.index}) "
                f"and count_per_pos ({count_per_pos.index})"
            )
        if (validate and
                not batch_count_per_pos.columns.equals(count_per_pos.columns)):
            raise ValueError(
                "Got different columns for "
                f"batch_count_per_pos ({batch_count_per_pos.columns}) "
                f"and count_per_pos ({count_per_pos.columns})"
            )
        count_per_pos += batch_count_per_pos
        if not isinstance(batch_count_per_read, pd.DataFrame):
            raise TypeError(f"batch_count_per_read must be DataFrame, "
                            f"but got {type(batch_count_per_read.__name__)}")
        if (validate
                and not batch_count_per_read.columns.equals(rel_header.index)):
            raise ValueError(
                "Got different columns for "
                f"batch_count_per_read ({batch_count_per_read.columns}) "
                f"and header ({rel_header.index})"
            )
        count_per_batch_read.append(batch_count_per_read)
        logger.detail(
            "\n".join([f"After batch {i}, accumulated"
                       f"num_reads = {num_reads}"
                       f"end_counts = {end_counts}"
                       f"count_per_pos = {count_per_pos}",
                       f"count_per_batch_read = {count_per_batch_read}"])
        )
    # Concatenate the per-read counts for the batches.
    if count_per_batch_read:
        count_per_read = pd.concat(count_per_batch_read, axis=0)
    else:
        count_per_read = pd.DataFrame(columns=rel_header.index, dtype=dtype)
    logger.routine(f"Ended accumulating patterns {patterns} in {batch_counts}")
    return num_reads, count_per_pos, count_per_read, end_counts


def accumulate_batches(batches: Iterable[SectionMutsBatch],
                       refseq: DNA,
                       pos_nums: np.ndarray,
                       patterns: dict[str, RelPattern],
                       ks: Iterable[int] | None = None,
                       validate: bool = True):
    return accumulate_counts((batch.calc_all_counts(patterns, ks)
                              for batch in batches),
                             refseq,
                             pos_nums,
                             patterns,
                             ks,
                             validate)

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
