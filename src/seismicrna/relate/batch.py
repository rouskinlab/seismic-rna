from abc import ABC
from functools import cached_property

import numpy as np

from ..core.batch import (AllReadBatch,
                          MutsBatch,
                          ReflenMutsBatch,
                          RefseqMutsBatch,
                          get_length)
from ..core.seq import POS_INDEX


class QnamesBatch(AllReadBatch):

    def __init__(self, *, names: list[str] | np.ndarray, **kwargs):
        super().__init__(**kwargs)
        self.names = np.asarray(names, dtype=str)

    @cached_property
    def num_reads(self):
        return get_length(self.names, "read names")


class RelateBatch(MutsBatch, AllReadBatch, ABC):

    @cached_property
    def num_reads(self):
        nums_reads = list(set(map(get_length, (self.end5s,
                                               self.mid5s,
                                               self.mid3s,
                                               self.end3s))))
        if len(nums_reads) != 1:
            raise ValueError(f"Got inconsistent numbers of reads: {nums_reads}")
        return nums_reads[0]

    @cached_property
    def pos_nums(self):
        return np.arange(POS_INDEX,
                         POS_INDEX + self.max_pos,
                         dtype=self.pos_dtype)


class RelateReflenBatch(RelateBatch, ReflenMutsBatch):
    pass


class RelateRefseqBatch(RelateBatch, RefseqMutsBatch):
    pass

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
