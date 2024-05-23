import os
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd
from click import command
from plotly import graph_objects as go

from .base import PosGraphWriter, PosGraphRunner
from .onetable import OneTableGraph, OneTableRunner, OneTableWriter
from .rel import OneRelGraph
from .roll import RollingGraph, RollingRunner
from .trace import iter_rolling_gini_traces
from ..core.header import format_clust_name
from ..core.mu import calc_signal_noise
from ..core.seq import BASEA, BASEC, BASE_NAME, iter_windows

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class RollingSNRGraph(OneTableGraph, OneRelGraph, RollingGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Rolling signal-to-noise ratio"

    @property
    def y_title(self):
        return "Signal-to-noise ratio"

    @cached_property
    def data(self):
        data = self._fetch_data(self.table,
                                order=self.order,
                                clust=self.clust)
        snr = pd.DataFrame(index=data.index, dtype=float)
        for cluster, cluster_data in data.items():
            cluster_snr = pd.Series(index=snr.index, dtype=float)
            for center, (window,) in iter_windows(cluster_data,
                                                  size=self._size,
                                                  min_count=self._min_count):
                seq = window.index.get_level_values(BASE_NAME)
                signal = np.logical_or(seq == BASEA, seq == BASEC)
                cluster_snr.loc[center] = calc_signal_noise(window, signal)
            if isinstance(cluster, tuple):
                _, order, clust = cluster
                label = format_clust_name(order, clust)
            else:
                label = cluster
            snr[label] = cluster_snr
        return snr

    def get_traces(self):
        for row, trace in enumerate(iter_rolling_gini_traces(self.data),
                                    start=1):
            yield (row, 1), trace

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class RollingSNRWriter(OneTableWriter, PosGraphWriter):

    def get_graph(self, rels_group: str, **kwargs):
        return RollingSNRGraph(table=self.table, rel=rels_group, **kwargs)


class RollingSNRRunner(RollingRunner, OneTableRunner, PosGraphRunner):

    @classmethod
    def get_writer_type(cls):
        return RollingSNRWriter


@command(COMMAND, params=RollingSNRRunner.params())
def cli(*args, **kwargs):
    """ Rolling signal-to-noise ratio. """
    return RollingSNRRunner.run(*args, **kwargs)

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
