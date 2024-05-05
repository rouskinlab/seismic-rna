from abc import ABC
from logging import getLogger

import numpy as np

from ..core.batch import (SectionMutsBatch,
                          PartialMutsBatch,
                          PartialReadBatch)
from ..core.seq import Section

logger = getLogger(__name__)


class MaskReadBatch(PartialReadBatch):

    def __init__(self, *, read_nums: np.ndarray, **kwargs):
        self._read_nums = read_nums
        super().__init__(**kwargs)

    @property
    def read_nums(self):
        return self._read_nums


class MaskMutsBatch(MaskReadBatch, SectionMutsBatch, PartialMutsBatch, ABC):
    pass


def apply_mask(batch: SectionMutsBatch,
               read_nums: np.ndarray | None = None,
               section: Section | None = None,
               sanitize: bool = False):
    # Determine which reads to use.
    if read_nums is not None:
        # Select specific read indexes.
        read_indexes = batch.read_indexes[read_nums]
        ends = [(e5[read_indexes], e3[read_indexes]) for e5, e3 in batch.ends]
        # Determine which reads were masked.
        masked_reads = np.setdiff1d(batch.read_nums, read_nums)
    else:
        # Use all reads.
        read_nums = batch.read_nums
        ends = batch.ends
        masked_reads = None
    # Determine the section of the new batch.
    if section is not None:
        # Clip the read coordinates to the section bounds.
        ends = [(np.clip(e5, section.end5, section.end3 + 1),
                 np.clip(e3, section.end5 - 1, section.end3))
                for e5, e3 in ends]
    else:
        # Use the same section as the given batch.
        section = batch.section
    # Select only the given positions and reads.
    muts = dict()
    for pos in section.unmasked_int:
        muts[pos] = dict()
        # Remove masked reads.
        for mut, pos_mut_reads in batch.muts[pos].items():
            muts[pos][mut] = (np.setdiff1d(pos_mut_reads,
                                           masked_reads,
                                           assume_unique=True)
                              if masked_reads is not None
                              else pos_mut_reads)
    return MaskMutsBatch(batch=batch.batch,
                         read_nums=read_nums,
                         section=section,
                         ends=ends,
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
