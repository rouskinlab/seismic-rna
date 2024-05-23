import os
from functools import cached_property
from logging import getLogger

import pandas as pd
from click import command
from plotly import graph_objects as go

from .base import PosGraphWriter, PosGraphRunner
from .onetable import OneTableGraph, OneTableRunner, OneTableWriter
from .rel import OneRelGraph
from .roll import RollingGraph, RollingRunner
from .trace import iter_rolling_gini_traces
from ..core.header import format_clust_name
from ..core.mu import calc_gini
from ..core.seq import iter_windows

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class RollingGiniGraph(OneTableGraph, OneRelGraph, RollingGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Rolling Gini coefficient"

    @property
    def y_title(self):
        return "Gini"

    @cached_property
    def data(self):
        data = self._fetch_data(self.table,
                                order=self.order,
                                clust=self.clust)
        gini = pd.DataFrame(index=data.index, dtype=float)
        for cluster, cluster_data in data.items():
            cluster_gini = pd.Series(index=gini.index, dtype=float)
            for center, (window,) in iter_windows(cluster_data,
                                                  size=self._size,
                                                  min_count=self._min_count):
                cluster_gini.loc[center] = calc_gini(window)
            if isinstance(cluster, tuple):
                _, order, clust = cluster
                label = format_clust_name(order, clust)
            else:
                label = cluster
            gini[label] = cluster_gini
        return gini

    def get_traces(self):
        for row, trace in enumerate(iter_rolling_gini_traces(self.data),
                                    start=1):
            yield (row, 1), trace

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class RollingGiniWriter(OneTableWriter, PosGraphWriter):

    def get_graph(self, rels_group: str, **kwargs):
        return RollingGiniGraph(table=self.table, rel=rels_group, **kwargs)


class RollingGiniRunner(RollingRunner, OneTableRunner, PosGraphRunner):

    @classmethod
    def get_writer_type(cls):
        return RollingGiniWriter


@command(COMMAND, params=RollingGiniRunner.params())
def cli(*args, **kwargs):
    """ Rolling Gini coefficient. """
    return RollingGiniRunner.run(*args, **kwargs)

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
