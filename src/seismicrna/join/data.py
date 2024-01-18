from functools import cached_property
from typing import Iterable

import numpy as np

from .report import JoinMaskReport
from ..core.data import LoadFunction, JoinedDataset
from ..mask.batch import MaskReadBatch
from ..mask.data import MaskReadDataset, MaskMutsDataset


class JoinMaskReadDataset(JoinedDataset):

    @classmethod
    def get_report_type(cls):
        return JoinMaskReport

    @classmethod
    def get_dataset_load_func(cls):
        return LoadFunction(MaskReadDataset)

    @cached_property
    def min_mut_gap(self):
        return self._get_common_attr("min_mut_gap")

    @cached_property
    def pos_kept(self):
        # Take the union of all positions kept among all datasets.
        pos_kept = None
        for data_pos in self._list_dataset_attr("pos_kept"):
            if pos_kept is not None:
                pos_kept = np.union1d(pos_kept, data_pos)
            else:
                pos_kept = data_pos
        if pos_kept is None:
            raise ValueError("Got no datasets to determine positions kept")
        return pos_kept

    def _join(self, batches: Iterable[MaskReadBatch]):
        batch_num = None
        read_nums = None
        for batch in batches:
            if batch_num is None:
                batch_num = batch.batch
            elif batch.batch != batch_num:
                raise ValueError(
                    f"Inconsistent batch number: {batch_num} ≠ {batch.batch}"
                )
            if read_nums is None:
                read_nums = batch.read_nums
            else:
                read_nums = np.union1d(read_nums, batch.read_nums)
        return MaskReadBatch(batch=batch_num, read_nums=read_nums)


class JoinMaskMutsDataset(MaskMutsDataset):

    @classmethod
    def get_dataset2_type(cls):
        return JoinMaskReadDataset

    @property
    def sects(self):
        return getattr(self.data2, "sects")


load_mask_join_data = LoadFunction(MaskMutsDataset, JoinMaskMutsDataset)

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
