from abc import ABC
from functools import cached_property
from typing import Iterable

import numpy as np
import pandas as pd

from .report import JoinMaskReport, JoinClusterReport
from ..cluster.batch import ClusterReadBatch
from ..cluster.data import ClusterReadDataset, ClusterMutsDataset
from ..core.data import LoadFunction, JoinedDataset, ChainedMutsDataset
from ..core.header import (ClustHeader,
                           index_orders_clusts,
                           list_clusts,
                           list_orders)
from ..mask.batch import MaskReadBatch
from ..mask.data import MaskReadDataset, MaskMutsDataset


class JoinMutsDataset(ChainedMutsDataset, ABC):

    @property
    def sects(self):
        return getattr(self.data2, "sects")


class JoinMaskReadDataset(JoinedDataset):

    @classmethod
    def get_report_type(cls):
        return JoinMaskReport

    @classmethod
    def get_dataset_load_func(cls):
        return LoadFunction(MaskReadDataset)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self._clusts is not None:
            raise TypeError(f"{self} has no clusters, but got {self._clusts}")

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

    def _join(self, batches: Iterable[tuple[str, MaskReadBatch]]):
        batch_num = None
        read_nums = None
        for _, batch in batches:
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
        if batch_num is None:
            raise ValueError("No batches were given to join")
        return MaskReadBatch(batch=batch_num, read_nums=read_nums)


class JoinMaskMutsDataset(JoinMutsDataset, MaskMutsDataset):

    @classmethod
    def get_dataset2_type(cls):
        return JoinMaskReadDataset


class JoinClusterReadDataset(JoinedDataset):

    @classmethod
    def get_report_type(cls):
        return JoinClusterReport

    @classmethod
    def get_dataset_load_func(cls):
        return LoadFunction(ClusterReadDataset)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not self._clusts:
            self._clusts = {sect: {order: {clust: clust
                                           for clust in list_clusts(order)}
                                   for order in list_orders(self.max_order)}
                            for sect in self.sects}
        if sorted(self._clusts) != sorted(self.sects):
            raise ValueError(f"{self} expected clusters for {self.sects}, "
                             f"but got {self._clusts}")

    @cached_property
    def max_order(self):
        return self._get_common_attr("max_order")

    @cached_property
    def clusts(self):
        """ Index of order and cluster numbers. """
        return index_orders_clusts(self.max_order)

    def _sect_cols(self, sect: str):
        """ Get the columns for a section's responsibilities. """
        clusts = self._clusts[sect]
        return pd.MultiIndex.from_tuples(
            [(order, clusts[order][clust])
             for order, clust in self.clusts],
            names=ClustHeader.level_names()
        )

    def _sect_resps(self, sect: str, resps: pd.DataFrame):
        """ Get the cluster responsibilities for a section. """
        # Reorder the columns.
        reordered = resps.loc[:, self._sect_cols(sect)]
        # Rename the columns by increasing order and cluster.
        reordered.columns = self.clusts
        return reordered

    def _join(self, batches: Iterable[tuple[str, ClusterReadBatch]]):
        batch_num = None
        resps = None
        for sect, batch in batches:
            if batch_num is None:
                batch_num = batch.batch
            elif batch.batch != batch_num:
                raise ValueError(
                    f"Inconsistent batch number: {batch_num} ≠ {batch.batch}"
                )
            if resps is None:
                resps = self._sect_resps(sect, batch.resps)
            else:
                resps = resps.add(self._sect_resps(sect, batch.resps),
                                  fill_value=0.)
        if batch_num is None:
            raise ValueError("No batches were given to join")
        return ClusterReadBatch(batch=batch_num,
                                resps=resps.fillna(0.).sort_index())


class JoinClusterMutsDataset(JoinMutsDataset, ClusterMutsDataset):

    @classmethod
    def get_dataset1_load_func(cls):
        return LoadFunction(JoinMaskMutsDataset)

    @classmethod
    def get_dataset2_type(cls):
        return JoinClusterReadDataset

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
