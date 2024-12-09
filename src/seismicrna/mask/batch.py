from abc import ABC
from functools import cached_property

import numpy as np
import pandas as pd

from ..core.array import calc_inverse, get_length
from ..core.batch import ReadBatch, RegionMutsBatch
from ..core.seq import Region


class PartialReadBatch(ReadBatch, ABC):

    @cached_property
    def max_read(self):
        return self.read_nums.max(initial=0)

    @cached_property
    def read_indexes(self):
        return calc_inverse(self.read_nums, what="read_nums", verify=False)


class PartialRegionMutsBatch(PartialReadBatch, RegionMutsBatch, ABC):
    pass


class MaskReadBatch(PartialReadBatch):

    def __init__(self, *, read_nums: np.ndarray, **kwargs):
        if not isinstance(read_nums, np.ndarray):
            raise TypeError(
                f"read_nums must be ndarray, but got {type(read_nums).__name__}"
            )
        self._read_nums = read_nums
        super().__init__(**kwargs)

    @property
    def read_nums(self):
        return self._read_nums

    @cached_property
    def num_reads(self):
        return get_length(self.read_nums, "read_nums")


class MaskMutsBatch(MaskReadBatch, PartialRegionMutsBatch):

    @property
    def read_weights(self):
        read_weights = None
        if self.masked_reads_bool.any():
            read_weights = np.ones(self.num_reads)
            read_weights[self.masked_reads_bool] = 0
            read_weights = pd.DataFrame(read_weights)
        return read_weights


def apply_mask(batch: RegionMutsBatch,
               read_nums: np.ndarray | None = None,
               region: Region | None = None,
               sanitize: bool = False):
    # Determine which reads to use.
    if read_nums is not None:
        # Select specific read indexes.
        read_indexes = batch.read_indexes[read_nums]
        seg_end5s = batch.seg_end5s[read_indexes]
        seg_end3s = batch.seg_end3s[read_indexes]
        # Determine which reads were masked.
        masked_reads = np.setdiff1d(batch.read_nums, read_nums)
    else:
        # Use all reads.
        read_nums = batch.read_nums
        seg_end5s = batch.seg_end5s
        seg_end3s = batch.seg_end3s
        masked_reads = None
    # Determine the region of the new batch.
    if region is not None:
        # Clip the read coordinates to the region bounds.
        seg_end5s = seg_end5s.clip(region.end5, region.end3 + 1)
        seg_end3s = seg_end3s.clip(region.end5 - 1, region.end3)
    else:
        # Use the same region as the given batch.
        region = batch.region
    # Select only the given positions and reads.
    muts = {pos: ({mut: np.setdiff1d(pos_mut_reads,
                                     masked_reads,
                                     assume_unique=True)
                   for mut, pos_mut_reads in batch.muts[pos].items()}
                  if masked_reads is not None
                  else batch.muts[pos])
            for pos in region.unmasked_int}
    return MaskMutsBatch(batch=batch.batch,
                         read_nums=read_nums,
                         region=region,
                         seg_end5s=seg_end5s,
                         seg_end3s=seg_end3s,
                         muts=muts,
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
