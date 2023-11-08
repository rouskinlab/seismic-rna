from abc import ABC
from typing import Iterable

import numpy as np

from ..core.batch import (RefseqMutsBatch,
                          PartialMutsBatch,
                          PartialReadBatch,
                          sanitize_pos)


class MaskReadBatch(PartialReadBatch):

    def __init__(self, *, read_nums: np.ndarray, **kwargs):
        self._read_nums = read_nums
        super().__init__(**kwargs)

    @property
    def read_nums(self):
        return self._read_nums


class MaskMutsBatch(MaskReadBatch, RefseqMutsBatch, PartialMutsBatch, ABC):
    pass


def apply_mask(batch: RefseqMutsBatch,
               reads: Iterable[int] | None = None,
               positions: Iterable[int] | None = None):
    # Clean and validate the selection.
    if positions is not None:
        positions = sanitize_pos(positions, batch.max_pos)
    # Select mutations at each position.
    muts = dict()
    for pos in positions if positions is not None else batch.pos_nums:
        muts[pos] = dict()
        # Select reads with each type of mutation at this position.
        for mut, pos_mut_reads in batch.muts.get(pos, dict()).items():
            muts[pos][mut] = (np.intersect1d(pos_mut_reads,
                                             reads,
                                             assume_unique=True)
                              if reads is not None
                              else pos_mut_reads)
    if reads is not None:
        read_nums = np.asarray(reads, dtype=batch.read_dtype)
        read_indexes = batch.read_indexes[read_nums]
        end5s = batch.end5s[read_indexes]
        mid5s = batch.mid5s[read_indexes]
        mid3s = batch.mid3s[read_indexes]
        end3s = batch.end3s[read_indexes]
    else:
        read_nums = batch.read_nums
        end5s = batch.end5s
        mid5s = batch.mid5s
        mid3s = batch.mid3s
        end3s = batch.end3s
    return MaskMutsBatch(batch=batch.batch,
                         refseq=batch.refseq,
                         muts=muts,
                         read_nums=read_nums,
                         end5s=end5s,
                         mid5s=mid5s,
                         mid3s=mid3s,
                         end3s=end3s,
                         sanitize=False)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
