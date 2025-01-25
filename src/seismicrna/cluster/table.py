from ._table import *
from .dataset import load_cluster_dataset
from ..mask.table import PartialDatasetTabulator


class ClusterDatasetTabulator(PartialDatasetTabulator, ClusterTabulator):

    @classmethod
    def load_function(cls):
        return load_cluster_dataset

    @classmethod
    def init_kws(cls):
        return super().init_kws() + ["ks"]
