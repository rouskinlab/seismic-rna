from functools import cached_property

import pandas as pd

from ..core.batch import PartialMutsBatch, PartialReadBatch, RefseqMutsBatch


class ClustReadBatch(PartialReadBatch):

    def __init__(self, *, resps: pd.DataFrame, **kwargs):
        self.resps = resps
        super().__init__(**kwargs)

    @cached_property
    def num_reads(self):
        return self.resps.sum(axis=0)

    @cached_property
    def read_nums(self):
        return self.resps.index.values


class ClustMutsBatch(ClustReadBatch, PartialMutsBatch, RefseqMutsBatch):

    @property
    def read_weights(self):
        return self.resps

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
