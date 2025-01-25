from abc import ABC, abstractmethod
from functools import cached_property
from itertools import combinations

import numpy as np
import pandas as pd

from .cgroup import make_tracks
from .dataset import DatasetGraph, DatasetWriter, DatasetRunner
from .trace import get_pairwise_position_trace
from ..core.header import NO_CLUSTS, NUM_CLUSTS_NAME, CLUST_NAME

POSITION_A = "Position A"
POSITION_B = "Position B"


def calc_phi(n: int | float | np.ndarray | pd.Series,
             a: int | float | np.ndarray | pd.Series | pd.DataFrame,
             b: int | float | np.ndarray | pd.Series | pd.DataFrame,
             a_and_b: int | float | np.ndarray | pd.Series | pd.DataFrame):
    """ Calculate the phi correlation coefficient for a 2x2 matrix.

    +----+----+
    | AB | AO | A.
    +----+----+
    | OB | OO | O.
    +----+----+
      .B   .O   ..

    where
      A. = AB + AO
      .B = AB + OB
      .. = A. + O. = .B + .O

    Parameters
    ----------
    n: int | float | np.ndarray | pd.Series
        Observations in total (..)
    a: int | float | np.ndarray | pd.Series | pd.DataFrame
        Observations for which A is true, regardless of B (A.)
    b: int | float | np.ndarray | pd.Series | pd.DataFrame
        Observations for which B is true, regardless of A (.B)
    a_and_b: int | float | np.ndarray | pd.Series | pd.DataFrame
        Observations for which A and B are both true (AB)

    Returns
    -------
    float | np.ndarray | pd.Series | pd.DataFrame
        Phi correlation coefficient
    """
    a_x_b = a * b
    with np.errstate(divide="ignore", invalid="ignore"):
        return (n * a_and_b - a_x_b) / np.sqrt(a_x_b * (n - a) * (n - b))


class PositionPairGraph(DatasetGraph, ABC):
    """ Function of pairs of positions. """

    @classmethod
    @abstractmethod
    def _pair_func(cls):
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
        positions = pd.MultiIndex.from_tuples(
            combinations(self.dataset.region.unmasked_int, 2),
            names=[POSITION_A, POSITION_B]
        )
        clusters = pd.MultiIndex.from_tuples(self.row_tracks,
                                             names=[NUM_CLUSTS_NAME,
                                                    CLUST_NAME])
        # Initialize the confusion matrix.
        if self.row_tracks == NO_CLUSTS:
            # The dataset has no clusters.
            zero = 0
            n = zero
        else:
            # The dataset has clusters.
            zero = 0.
            n = pd.Series(zero, clusters)
        a_accum = pd.DataFrame(zero, positions, clusters)
        b_accum = pd.DataFrame(zero, positions, clusters)
        ab_accum = pd.DataFrame(zero, positions, clusters)
        # Fill the confusion matrix, accumulating over the batches.
        for batch in self.dataset.iter_batches():
            reads_per_pos = {
                pos: batch.read_indexes[reads] for pos, reads
                in batch.reads_per_pos(self.pattern).items()
            }
            # Count the reads in the batch.
            if batch.read_weights is None:
                n += batch.num_reads
            elif isinstance(batch.read_weights, pd.DataFrame):
                for clust in clusters:
                    n += batch.read_weights[clust].sum()
            else:
                raise TypeError(batch.read_weights)
            # For each pair of positions a and b, count the reads that
            # fit the relationship pattern for positions a and b (ab),
            # only position a (ao), only position b (ob), and neither
            # position (oo).
            for pos_ab in positions:
                pos_a, pos_b = pos_ab
                reads_a = reads_per_pos[pos_a]
                reads_b = reads_per_pos[pos_b]
                reads_ab = np.intersect1d(reads_a, reads_b, assume_unique=True)
                if batch.read_weights is None:
                    # There is only one cluster, and every read has the
                    # same weight.
                    a = reads_a.size
                    b = reads_b.size
                    ab = reads_ab.size
                    a_accum.loc[pos_ab] += a
                    b_accum.loc[pos_ab] += b
                    ab_accum.loc[pos_ab] += ab
                elif isinstance(batch.read_weights, pd.DataFrame):
                    # There are multiple clusters, where each read gets
                    # a different weight in each cluster.
                    for clust in clusters:
                        weights = batch.read_weights[clust].values
                        a = weights[reads_a].sum()
                        b = weights[reads_b].sum()
                        ab = weights[reads_ab].sum()
                        a_accum.loc[pos_ab, clust] += a
                        b_accum.loc[pos_ab, clust] += b
                        ab_accum.loc[pos_ab, clust] += ab
                else:
                    raise TypeError(batch.read_weights)
        pair_func = self._pair_func()
        return pair_func(n, a_accum, b_accum, ab_accum)

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
