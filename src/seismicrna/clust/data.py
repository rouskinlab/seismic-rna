from functools import cached_property
from logging import getLogger

from .batch import ClustMutsBatch
from .io import ClustBatchIO
from .report import ClustReport
from ..core.data import (BatchedLoadedDataset,
                         BatchedMergedDataset,
                         MergedMutsDataset)
from ..core.header import index_orders_clusts
from ..core.report import NumClustsF
from ..mask.batch import MaskMutsBatch
from ..mask.data import MaskMerger

logger = getLogger(__name__)


class ClustLoader(BatchedLoadedDataset):
    """ Load clustering results. """

    @classmethod
    def get_report_type(cls):
        return ClustReport

    @classmethod
    def get_data_type(cls):
        return ClustBatchIO

    @cached_property
    def max_order(self):
        """ Number of clusters. """
        return self.report.get_field(NumClustsF)

    @property
    def pattern(self):
        return None


class ClustMerger(BatchedMergedDataset, MergedMutsDataset):
    """ Merge cluster responsibilities with mutation data. """

    @classmethod
    def get_data_type(cls):
        return ClustMutsBatch

    @classmethod
    def get_dataset1_type(cls):
        return MaskMerger

    @classmethod
    def get_dataset2_type(cls):
        return ClustLoader

    @property
    def min_mut_gap(self):
        return self.data1.min_mut_gap

    @property
    def pattern(self):
        return self.data1.pattern

    @cached_property
    def section(self):
        return self.data1.section

    @cached_property
    def max_order(self):
        return self.data2.max_order

    @cached_property
    def clusters(self):
        return index_orders_clusts(self.max_order)

    def _merge(self, batch1: MaskMutsBatch, batch2: ClustBatchIO):
        return self.get_data_type()(batch=batch1.batch,
                                    refseq=batch1.refseq,
                                    muts=batch1.muts,
                                    end5s=batch1.end5s,
                                    mid5s=batch1.mid5s,
                                    mid3s=batch1.mid3s,
                                    end3s=batch1.end3s,
                                    resps=batch2.resps,
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
