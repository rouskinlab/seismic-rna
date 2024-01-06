from functools import cached_property

from .batch import apply_mask
from .io import MaskBatchIO
from .report import MaskReport
from ..core.data import ChainedMutsDataset, LoadedDataset, LoadFunction
from ..core.rel import RelPattern
from ..core.report import CountMutsF, CountRefsF, MinMutGapF, PosKeptF
from ..pool.data import load_relate_pool_dataset
from ..relate.batch import RelateRefseqBatch


class MaskReadDataset(LoadedDataset):
    """ Load batches of masked relation vectors. """

    @classmethod
    def get_report_type(cls):
        return MaskReport

    @classmethod
    def get_batch_type(cls):
        return MaskBatchIO

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


class MaskMutsDataset(ChainedMutsDataset):
    """ Chain mutation data with masked reads. """

    MASK_NAME = "mask"

    @classmethod
    def get_dataset1_load_func(cls):
        return load_relate_pool_dataset

    @classmethod
    def get_dataset2_type(cls):
        return MaskReadDataset

    @property
    def min_mut_gap(self):
        return getattr(self.data2, "min_mut_gap")

    @property
    def pattern(self):
        return self.data2.pattern

    @cached_property
    def section(self):
        section = super().section
        section.add_mask(self.MASK_NAME,
                         getattr(self.data2, "pos_kept"),
                         invert=True)
        return section

    def _chain(self, batch1: RelateRefseqBatch, batch2: MaskBatchIO):
        return apply_mask(batch1,
                          batch2.read_nums,
                          getattr(self.data2, "pos_kept"))


load_mask_dataset = LoadFunction(MaskMutsDataset)

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
