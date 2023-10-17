from functools import cached_property

from .io import MaskReadBatchIO
from .report import MaskReport
from ..core.batch import MaskMutsBatch
from ..core.io import (CountMutsF,
                       CountRefsF,
                       MinMutGapF,
                       PosKeptF,
                       BatchedLoadedDataset,
                       BatchedMergedDataset,
                       MergedMutsDataset)
from ..core.rel import RelPattern
from ..relate.data import RelateLoader
from ..relate.io import RelateBatchIO


class MaskLoader(BatchedLoadedDataset):
    """ Load batches of masked relation vectors. """

    @classmethod
    def get_report_type(cls):
        return MaskReport

    @classmethod
    def get_data_type(cls):
        return MaskReadBatchIO

    @property
    def min_mut_gap(self):
        return self.report.get_field(MinMutGapF)

    @property
    def pos_kept(self):
        return self.report.get_field(PosKeptF)

    @cached_property
    def pattern(self):
        return RelPattern(self.report.get_field(CountMutsF),
                          self.report.get_field(CountRefsF))


class MaskMerger(BatchedMergedDataset, MergedMutsDataset):
    """ Merge mutation data with masked reads. """

    MASK_NAME = "mask"

    @classmethod
    def get_data_type(cls):
        return MaskMutsBatch

    @classmethod
    def get_data1_type(cls):
        return RelateLoader

    @classmethod
    def get_data2_type(cls):
        return MaskLoader

    @property
    def min_mut_gap(self):
        return self.data2.min_mut_gap

    @property
    def pattern(self):
        return self.data2.pattern

    @cached_property
    def section(self):
        section = super().section
        section.add_mask(self.MASK_NAME, self.data2.pos_kept, invert=True)
        return section

    def _merge(self, batch1: RelateBatchIO, batch2: MaskReadBatchIO):
        return batch1.mask(positions=self.data2.pos_kept,
                           reads=batch2.read_nums)

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
