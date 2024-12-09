from abc import ABC
from functools import cached_property

import pandas as pd

from .data import ClusterMutsDataset
from ..core import path
from ..core.header import ClustHeader, RelClustHeader, make_header, parse_header
from ..core.table import (AbundanceTable,
                          RelTypeTable,
                          PositionTableWriter,
                          AbundanceTableWriter)
from ..mask.table import (PartialTable,
                          PartialPositionTable,
                          PartialTabulator,
                          PartialDatasetTabulator)
from ..relate.table import TableLoader, PositionTableLoader


class ClusterTable(RelTypeTable, ABC):

    @classmethod
    def kind(cls):
        return path.CMD_CLUST_DIR

    @classmethod
    def header_type(cls):
        return RelClustHeader


class ClusterPositionTable(ClusterTable, PartialPositionTable, ABC):
    pass


class ClusterAbundanceTable(AbundanceTable, PartialTable, ABC):

    @classmethod
    def kind(cls):
        return path.CMD_CLUST_DIR

    @classmethod
    def header_type(cls):
        return ClustHeader

    @classmethod
    def index_depth(cls):
        return cls.header_depth()

    def _get_header(self):
        return parse_header(self.data.index)


class ClusterPositionTableWriter(PositionTableWriter, ClusterPositionTable):
    pass


class ClusterAbundanceTableWriter(AbundanceTableWriter, ClusterAbundanceTable):
    pass


class ClusterPositionTableLoader(PositionTableLoader, ClusterPositionTable):
    """ Load cluster data indexed by position. """


class ClusterAbundanceTableLoader(TableLoader, ClusterAbundanceTable):
    """ Load cluster data indexed by cluster. """

    @cached_property
    def data(self) -> pd.Series:
        data = pd.read_csv(self.path,
                           index_col=self.index_cols()).squeeze(axis=1)
        if not isinstance(data, pd.Series):
            raise ValueError(f"{self} must have one column, but got\n{data}")
        # Any numeric data in the header will be read as strings and
        # must be cast to integers using parse_header.
        header = parse_header(data.index)
        # The index must be replaced with the header index for the
        # type casting to take effect.
        data.index = header.index
        return data


class ClusterTabulator(PartialTabulator, ABC):

    @classmethod
    def table_types(cls):
        return [ClusterPositionTableWriter, ClusterAbundanceTableWriter]

    def __init__(self, *, ks: list[int], **kwargs):
        super().__init__(**kwargs)
        if ks is None:
            raise ValueError(
                f"{type(self).__name__} requires clusters, but got ks={ks}"
            )
        self.ks = ks

    @cached_property
    def clust_header(self):
        """ Header of the per-cluster data. """
        return make_header(ks=self.ks)

    @cached_property
    def data_per_clust(self):
        """ Number of reads in each cluster. """
        n_rels, n_clust = self._adjusted
        n_clust.name = "Number of Reads"
        return n_clust


class ClusterDatasetTabulator(ClusterTabulator, PartialDatasetTabulator):

    def __init__(self, *, dataset: ClusterMutsDataset, **kwargs):
        super().__init__(dataset=dataset, ks=dataset.ks, **kwargs)

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
