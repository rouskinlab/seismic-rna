from abc import ABC
from logging import getLogger
from typing import Iterable

import numpy as np

from ..core.batch import (RefseqMutsBatch,
                          PartialMutsBatch,
                          PartialReadBatch,
                          has_mids,
                          sanitize_pos)

logger = getLogger(__name__)


class MaskReadBatch(PartialReadBatch):

    def __init__(self, *, read_nums: np.ndarray, **kwargs):
        self._read_nums = read_nums
        super().__init__(**kwargs)

    @property
    def read_nums(self):
        return self._read_nums


class MaskMutsBatch(MaskReadBatch, RefseqMutsBatch, PartialMutsBatch, ABC):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.num_discontiguous_reads:
            raise ValueError(
                f"{self} has {self.num_discontiguous_reads} discontiguous "
                "paired-end reads, which are not (yet) supported for masking"
            )


def apply_mask(batch: RefseqMutsBatch,
               reads: Iterable[int] | None = None,
               positions: Iterable[int] | None = None,
               clip5: int | None = None,
               clip3: int | None = None):
    # Determine masked read numbers.
    masked_reads = (np.setdiff1d(batch.read_nums, reads)
                    if reads is not None
                    else None)
    # Clean and validate the selection.
    if positions is not None:
        positions = sanitize_pos(positions, batch.max_pos)
    # Select mutations at each position.
    muts = dict()
    for pos in positions if positions is not None else batch.pos_nums:
        if clip5 is not None and pos < clip5:
            logger.warning(f"Skipped clipped position {pos} (< {clip5})")
            continue
        if clip3 is not None and pos > clip3:
            logger.warning(f"Skipped clipped position {pos} (> {clip3})")
            continue
        muts[pos] = dict()
        # Remove masked reads with each type of mutation at this position.
        for mut, pos_mut_reads in batch.muts.get(pos, dict()).items():
            muts[pos][mut] = (np.setdiff1d(pos_mut_reads,
                                           masked_reads,
                                           assume_unique=True)
                              if reads is not None
                              else pos_mut_reads)
    if reads is not None:
        # Select specific read indexes.
        read_nums = np.asarray(reads, dtype=batch.read_dtype)
        read_indexes = batch.read_indexes[read_nums]
        end5s = batch.end5s[read_indexes]
        end3s = batch.end3s[read_indexes]
        if has_mids(batch.mid5s, batch.mid3s):
            mid5s = batch.mid5s[read_indexes]
            mid3s = batch.mid3s[read_indexes]
        else:
            mid5s = None
            mid3s = None
    else:
        # Use all reads.
        read_nums = batch.read_nums
        end5s = batch.end5s
        end3s = batch.end3s
        mid5s = batch.mid5s
        mid3s = batch.mid3s
    # Clip the 5' and 3' end and middle coordinates.
    if clip5 is not None or clip3 is not None:
        if clip5 is not None and clip3 is not None and clip5 > clip3:
            raise ValueError("Must have clip5 ≤ clip3, "
                             f"but got {clip5} > {clip3}")
        if clip5 is not None:
            end5s = np.maximum(end5s, clip5)
        if clip3 is not None:
            end3s = np.minimum(end3s, clip3)
        if clip5 is not None and mid5s is not None:
            mid5s = np.minimum(end3s, np.maximum(mid5s, clip5))
        if clip3 is not None and mid3s is not None:
            mid3s = np.maximum(end5s, np.minimum(mid3s, clip3))
        sanitize = True
    else:
        sanitize = False
    return MaskMutsBatch(batch=batch.batch,
                         refseq=batch.refseq,
                         muts=muts,
                         read_nums=read_nums,
                         end5s=end5s,
                         mid5s=mid5s,
                         mid3s=mid3s,
                         end3s=end3s,
                         sanitize=sanitize)

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
