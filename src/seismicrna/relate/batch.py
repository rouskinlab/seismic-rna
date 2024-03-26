from abc import ABC
from functools import cached_property

import numpy as np

from ..core.batch import (AllReadBatch,
                          MutsBatch,
                          ReflenMutsBatch,
                          RefseqMutsBatch,
                          ensure_same_length,
                          get_length)


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
        return ensure_same_length(self.end5s, self.end3s, "end5s", "end3s")

    @cached_property
    def pos_nums(self):
        return np.arange(1, self.max_pos + 1, dtype=self.pos_dtype)


class RelateReflenBatch(RelateBatch, ReflenMutsBatch):
    pass


class RelateRefseqBatch(RelateBatch, RefseqMutsBatch):
    pass

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
