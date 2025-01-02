import os
from functools import cached_property
from itertools import combinations

import numpy as np
import pandas as pd
from click import command

from .base import get_action_name, make_tracks
from .onedataset import OneDatasetGraph, OneDatasetWriter, OneDatasetRunner
from .rel import OneRelGraph
from .trace import get_pairwise_position_trace
from ..core.header import NO_CLUSTS, NUM_CLUSTS_NAME, CLUST_NAME

COMMAND = __name__.split(os.path.extsep)[-1]

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


class PositionCorrelationGraph(OneDatasetGraph, OneRelGraph):
    """ Phi correlations between pairs of positions. """

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Phi correlation"

    @cached_property
    def action(self):
        return get_action_name(self.dataset)

    @property
    def x_title(self):
        return "Position A"

    @property
    def y_title(self):
        return "Position B"

    @cached_property
    def row_tracks(self):
        return make_tracks(self.dataset, self.k, self.clust)

    @cached_property
    def data(self):
        positions = pd.MultiIndex.from_tuples(
            combinations(self.dataset.region.unmasked_int, 2),
            names=[POSITION_A, POSITION_B]
        )
        # Initialize the confusion matrix.
        if self.row_tracks == NO_CLUSTS:
            # The dataset has no clusters.
            clusters = None
            n = 0
            a_accum = pd.Series(0, positions)
            b_accum = pd.Series(0, positions)
            ab_accum = pd.Series(0, positions)
        else:
            # The dataset has clusters.
            clusters = pd.MultiIndex.from_tuples(self.row_tracks,
                                                 names=[NUM_CLUSTS_NAME,
                                                        CLUST_NAME])
            n = pd.Series(0., clusters)
            a_accum = pd.DataFrame(0., positions, clusters)
            b_accum = pd.DataFrame(0., positions, clusters)
            ab_accum = pd.DataFrame(0., positions, clusters)
        # Fill the confusion matrix, accumulating over the batches.
        for batch in self.dataset.iter_batches():
            reads_per_pos = {
                pos: batch.read_indexes[reads] for pos, reads
                in batch.reads_per_pos(self.dataset.pattern).items()
            }
            # Count the reads in the batch.
            if clusters is not None:
                if not isinstance(batch.read_weights, pd.DataFrame):
                    raise TypeError(
                        "batch.read_weights must be a DataFrame, "
                        f"but got {type(batch.read_weights).__name__}"
                    )
                for clust in clusters:
                    n += batch.read_weights[clust].sum()
            else:
                if batch.read_weights is not None:
                    raise TypeError(
                        "batch.read_weights must be None, "
                        f"but got {type(batch.read_weights).__name__}"
                    )
                n += batch.num_reads
            # For each pair of positions a and b, count the reads that
            # fit the relationship pattern for positions a and b (ab),
            # only position a (ao), only position b (ob), and neither
            # position (oo).
            for pos_ab in positions:
                pos_a, pos_b = pos_ab
                reads_a = reads_per_pos[pos_a]
                reads_b = reads_per_pos[pos_b]
                reads_ab = np.intersect1d(reads_a, reads_b, assume_unique=True)
                if clusters is not None:
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
                    # There is only one cluster, and every read has the
                    # same weight.
                    a = reads_a.size
                    b = reads_b.size
                    ab = reads_ab.size
                    a_accum.loc[pos_ab] += a
                    b_accum.loc[pos_ab] += b
                    ab_accum.loc[pos_ab] += ab
        return calc_phi(n, a_accum, b_accum, ab_accum)

    def get_traces(self):
        if isinstance(self.data, pd.Series):
            trace = get_pairwise_position_trace(self.data,
                                                self.dataset.region.end5,
                                                self.dataset.region.end3)
            yield (1, 1), trace
        else:
            assert isinstance(self.data, pd.DataFrame)
            for row, (_, values) in enumerate(self.data.items(), start=1):
                trace = get_pairwise_position_trace(values,
                                                    self.dataset.region.end5,
                                                    self.dataset.region.end3)
                yield (row, 1), trace


class PositionCorrelationWriter(OneDatasetWriter):

    def get_graph(self, rel, **kwargs):
        return PositionCorrelationGraph(dataset=self.dataset, rel=rel, **kwargs)


class PositionCorrelationRunner(OneDatasetRunner):

    @classmethod
    def get_writer_type(cls):
        return PositionCorrelationWriter


@command(COMMAND, params=PositionCorrelationRunner.params())
def cli(*args, **kwargs):
    """ Phi correlations between pairs of positions. """
    return PositionCorrelationRunner.run(*args, **kwargs)

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
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
