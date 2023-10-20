from functools import cached_property

import pandas as pd

from .names import ORD_CLS_NAME, READ_NAME
from ..core.batch import PartialReadBatch


class ClusterReadBatch(PartialReadBatch):

    def __init__(self, *, resps: pd.DataFrame, **kwargs):
        # Rename the index.
        resps.index.rename(READ_NAME, inplace=True)
        # Rename the column levels.
        resps.columns.rename(ORD_CLS_NAME, inplace=True)
        self.resps = resps
        super().__init__(**kwargs)

    @cached_property
    def read_nums(self):
        return self.resps.index.values

    @property
    def clusters(self):
        """ Order and number of each cluster. """
        return self.resps.columns
