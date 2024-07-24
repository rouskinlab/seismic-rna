from abc import ABC, abstractmethod
from functools import cached_property
from typing import Callable

import pandas as pd
from plotly import graph_objects as go

from .base import PosGraphWriter, PosGraphRunner
from .onetable import OneTableGraph, OneTableRunner, OneTableWriter
from .rel import OneRelGraph
from .roll import RollingGraph, RollingRunner
from .trace import iter_line_traces
from ..core.header import format_clust_name
from ..core.seq import iter_windows


class RollingStatGraph(OneTableGraph, OneRelGraph, RollingGraph, ABC):

    @classmethod
    @abstractmethod
    def stat_func(cls) -> Callable[[pd.Series], pd.Series]:
        """ Function to compute a statistic on the data. """

    @cached_property
    def data(self):
        stat_func = self.stat_func()
        data = self._fetch_data(self.table,
                                k=self.k,
                                clust=self.clust)
        stat = pd.DataFrame(index=data.index, dtype=float)
        for cluster, cluster_data in data.items():
            cluster_stat = pd.Series(index=stat.index, dtype=float)
            for center, (window,) in iter_windows(cluster_data,
                                                  size=self._size,
                                                  min_count=self._min_count):
                cluster_stat.loc[center] = stat_func(window)
            if isinstance(cluster, tuple):
                _, k, clust = cluster
                label = format_clust_name(k, clust)
            else:
                label = cluster
            stat[label] = cluster_stat
        return stat

    def get_traces(self):
        for row, trace in enumerate(iter_line_traces(self.data),
                                    start=1):
            yield (row, 1), trace

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class RollingStatWriter(OneTableWriter, PosGraphWriter, ABC):

    @classmethod
    @abstractmethod
    def get_graph_type(cls) -> type[RollingGraph]:
        """ Type of graph. """

    def get_graph(self, rels_group: str, **kwargs):
        graph_type = self.get_graph_type()
        return graph_type(table=self.table, rel=rels_group, **kwargs)


class RollingStatRunner(RollingRunner, OneTableRunner, PosGraphRunner, ABC):
    pass

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
