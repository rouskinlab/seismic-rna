from abc import ABC
from functools import cached_property

import pandas as pd

from ..core.batch import PartialMutsBatch, PartialReadBatch


class ClustReadBatch(PartialReadBatch):

    def __init__(self, *, resps: pd.DataFrame, **kwargs):
        self.resps = resps
        super().__init__(**kwargs)

    @cached_property
    def num_reads(self):
        return self.resps.sum(axis=0)

    @cached_property
    def read_nums(self):
        return self.resps.index.values

    @property
    def clusters(self):
        return self.resps.columns


class ClustMutsBatch(ClustReadBatch, PartialMutsBatch, ABC):

    @property
    def read_weights(self):
        return self.resps
