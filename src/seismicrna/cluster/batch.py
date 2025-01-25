from functools import cached_property

import pandas as pd

from ..mask.batch import PartialRegionMutsBatch, PartialReadBatch


class ClusterReadBatch(PartialReadBatch):

    def __init__(self, *, resps: pd.DataFrame, **kwargs):
        self.resps = resps
        super().__init__(**kwargs)

    @cached_property
    def num_reads(self) -> pd.Series:
        return self.resps.sum(axis=0)

    @cached_property
    def read_nums(self):
        return self.resps.index.values


class ClusterMutsBatch(ClusterReadBatch, PartialRegionMutsBatch):

    @property
    def read_weights(self):
        self.resps.loc[self.masked_reads_bool] = 0
        return self.resps
