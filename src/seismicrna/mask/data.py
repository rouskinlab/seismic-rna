from functools import cached_property

from .batch import MaskMutsBatch, apply_mask
from .io import MaskBatchIO
from .report import MaskReport
from ..core.data import (ArrowDataset,
                         LoadedDataset,
                         LoadFunction,
                         MergedUnbiasDataset,
                         UnbiasDataset)
from ..core.rel import RelPattern
from ..core.report import (CountMutsF,
                           CountRefsF,
                           MinMutGapF,
                           PosKeptF,
                           QuickUnbiasF,
                           QuickUnbiasThreshF)
from ..joinbase.data import (BATCH_NUM,
                             READ_NUMS,
                             SEG_END5S,
                             SEG_END3S,
                             MUTS,
                             JoinMutsDataset)
from ..joinbase.report import JoinMaskReport
from ..pool.data import load_relate_dataset
from ..relate.batch import RelateBatch


class MaskReadDataset(LoadedDataset, UnbiasDataset):
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
    def quick_unbias(self):
        return self.report.get_field(QuickUnbiasF)

    @property
    def quick_unbias_thresh(self):
        return self.report.get_field(QuickUnbiasThreshF)

    @property
    def pos_kept(self):
        """ Positions kept after masking. """
        return self.report.get_field(PosKeptF)

    @cached_property
    def pattern(self):
        return RelPattern(self.report.get_field(CountMutsF),
                          self.report.get_field(CountRefsF))


class MaskMutsDataset(ArrowDataset, UnbiasDataset):
    """ Chain mutation data with masked reads. """

    MASK_NAME = "mask"

    @classmethod
    def get_dataset1_load_func(cls):
        return load_relate_dataset

    @classmethod
    def get_dataset2_type(cls):
        return MaskReadDataset

    @property
    def pattern(self):
        return self.data2.pattern

    @property
    def min_mut_gap(self):
        return getattr(self.data2, "min_mut_gap")

    @property
    def quick_unbias(self):
        return getattr(self.data2, "quick_unbias")

    @property
    def quick_unbias_thresh(self):
        return getattr(self.data2, "quick_unbias_thresh")

    @cached_property
    def section(self):
        # Mask the positions that were not kept.
        section = super().section
        section.add_mask(self.MASK_NAME,
                         getattr(self.data2, "pos_kept"),
                         complement=True)
        return section

    def _integrate(self, batch1: RelateBatch, batch2: MaskBatchIO):
        return apply_mask(batch1,
                          batch2.read_nums,
                          self.section,
                          sanitize=False)


class JoinMaskMutsDataset(JoinMutsDataset, MergedUnbiasDataset):

    @classmethod
    def get_report_type(cls):
        return JoinMaskReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_mask_dataset

    @classmethod
    def get_batch_type(cls):
        return MaskMutsBatch

    @classmethod
    def name_batch_attrs(cls):
        return [BATCH_NUM, READ_NUMS, SEG_END5S, SEG_END3S, MUTS]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self._clusts is not None:
            raise TypeError(f"{self} has no clusters, but got {self._clusts}")


load_mask_dataset = LoadFunction(MaskMutsDataset, JoinMaskMutsDataset)

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
