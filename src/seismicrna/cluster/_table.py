from abc import ABC
from functools import cached_property

import pandas as pd

from ..core import path
from ..core.header import (NUM_CLUSTS_NAME,
                           ClustHeader,
                           RelClustHeader,
                           make_header,
                           parse_header)
from ..core.table import (BatchTabulator,
                          CountTabulator,
                          AbundanceTable,
                          RelTypeTable,
                          PositionTableWriter,
                          AbundanceTableWriter)
from ..mask.table import (PartialTable,
                          PartialPositionTable,
                          PartialTabulator)
from ..relate.table import TableLoader, PositionTableLoader


class ClusterTable(RelTypeTable, ABC):

    @classmethod
    def kind(cls):
        return path.CLUSTER_STEP

    @classmethod
    def header_type(cls):
        return RelClustHeader


class ClusterPositionTable(ClusterTable, PartialPositionTable, ABC):
    pass


class ClusterAbundanceTable(AbundanceTable, PartialTable, ABC):

    @classmethod
    def kind(cls):
        return path.CLUSTER_STEP

    @classmethod
    def header_type(cls):
        return ClustHeader

    @classmethod
    def index_depth(cls):
        return cls.header_depth()

    def _get_header(self):
        return parse_header(self.data.index)

    @cached_property
    def proportions(self):
        return self.data / self.data.groupby(level=NUM_CLUSTS_NAME).sum()


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


class ClusterBatchTabulator(BatchTabulator, ClusterTabulator):
    pass


class ClusterCountTabulator(CountTabulator, ClusterTabulator):
    pass
