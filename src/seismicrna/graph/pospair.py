from abc import ABC, abstractmethod
from functools import cached_property

import pandas as pd

from .cgroup import make_tracks
from .dataset import DatasetGraph, DatasetWriter, DatasetRunner
from .trace import get_pairwise_position_trace
from ..core.batch import POSITION_A, POSITION_B, accumulate_confusion_matrices
from ..core.dataset import UnbiasDataset
from ..core.header import NO_CLUST, NO_CLUSTS, NUM_CLUSTS_NAME, CLUST_NAME


class PositionPairGraph(DatasetGraph, ABC):
    """ Function of pairs of positions. """

    @classmethod
    @abstractmethod
    def get_pair_func(cls):
        """ Function to compare each pair of positions. """

    @property
    def x_title(self):
        return POSITION_A

    @property
    def y_title(self):
        return POSITION_B

    @cached_property
    def row_tracks(self):
        return make_tracks(self.dataset, self.k, self.clust)

    @cached_property
    def data(self):
        if self.row_tracks is not None and self.row_tracks != NO_CLUSTS:
            clusters = pd.MultiIndex.from_tuples(
                self.row_tracks,
                names=[NUM_CLUSTS_NAME, CLUST_NAME]
            )
        else:
            clusters = None
        pair_func = self.get_pair_func()
        data = pair_func(*accumulate_confusion_matrices(
            self.dataset.get_batch,
            self.dataset.num_batches,
            self.pattern,
            self.dataset.region.unmasked,
            clusters,
            min_gap=(self.dataset.min_mut_gap
                     if isinstance(self.dataset, UnbiasDataset)
                     else 0),
            num_cpus=self.num_cpus
        ))
        if clusters is None:
            assert isinstance(data, pd.Series)
            data = data.to_frame(name=NO_CLUST)
        return data

    def get_traces(self):
        for row, (_, values) in enumerate(self.data.items(), start=1):
            trace = get_pairwise_position_trace(values,
                                                self.dataset.region.end5,
                                                self.dataset.region.end3)
            yield (row, 1), trace


class PositionPairWriter(DatasetWriter):

    @classmethod
    @abstractmethod
    def graph_type(cls):
        """ Type of graph. """
        return type[PositionPairGraph]

    def get_graph(self, rel, **kwargs):
        graph_type = self.graph_type()
        return graph_type(dataset=self.dataset, rel=rel, **kwargs)


class PositionPairRunner(DatasetRunner, ABC):

    @classmethod
    @abstractmethod
    def get_writer_type(cls) -> type[PositionPairWriter]:
        pass
