from functools import cached_property
from logging import getLogger

import numpy as np

from ..core.array import get_length
from ..core.batch import (SectionMutsBatch,
                          PartialMutsBatch,
                          PartialReadBatch)
from ..core.seq import Section

logger = getLogger(__name__)


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


class MaskMutsBatch(MaskReadBatch, SectionMutsBatch, PartialMutsBatch):

    @property
    def read_weights(self):
        return None


def apply_mask(batch: SectionMutsBatch,
               read_nums: np.ndarray | None = None,
               section: Section | None = None,
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
    # Determine the section of the new batch.
    if section is not None:
        # Clip the read coordinates to the section bounds.
        seg_end5s = seg_end5s.clip(section.end5, section.end3 + 1)
        seg_end3s = seg_end3s.clip(section.end5 - 1, section.end3)
    else:
        # Use the same section as the given batch.
        section = batch.section
    # Select only the given positions and reads.
    muts = {pos: ({mut: np.setdiff1d(pos_mut_reads,
                                     masked_reads,
                                     assume_unique=True)
                   for mut, pos_mut_reads in batch.muts[pos].items()}
                  if masked_reads is not None
                  else batch.muts[pos])
            for pos in section.unmasked_int}
    return MaskMutsBatch(batch=batch.batch,
                         read_nums=read_nums,
                         section=section,
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
