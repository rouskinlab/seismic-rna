from abc import ABC
from typing import Iterable

import numpy as np
import pandas as pd

from .report import JoinReport
from ..cluster.batch import ClusterReadBatch
from ..cluster.data import ClusterReadDataset
from ..core.data import LoadFunction, JoinedMutsDataset
from ..mask.batch import MaskReadBatch
from ..mask.data import MaskReadDataset


class JoinDataset(JoinedMutsDataset, ABC):
    """ Mask or Cluster batches joined together. """

    @classmethod
    def get_report_type(cls):
        return JoinReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_mask_cluster_dataset


class JoinMaskDataset(JoinDataset):

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


class JoinClusterDataset(JoinDataset):

    def _join(self, batches: Iterable[ClusterReadBatch]):
        batch_num = None
        resps = None
        for batch in batches:
            if batch_num is None:
                batch_num = batch.batch
            elif batch.batch != batch_num:
                raise ValueError(
                    f"Inconsistent batch number: {batch_num} ≠ {batch.batch}"
                )
            if resps is None:
                resps = batch.resps
            else:
                # FIXME: use the cluster linkage table to join clusters.
                resps = pd.concat([resps, batch.resps])
        return ClusterReadBatch(batch=batch_num, resps=resps)


load_mask_cluster_dataset = LoadFunction(MaskReadDataset, ClusterReadDataset)

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
