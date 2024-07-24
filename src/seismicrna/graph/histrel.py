from abc import ABC, abstractmethod
from functools import cached_property

import numpy as np
import pandas as pd

from .color import ColorMapGraph, RelColorMap
from .hist import HistogramGraph, HistogramRunner, get_edges_index
from .onetable import OneTableGraph, OneTableRunner, OneTableWriter
from .rel import MultiRelsGraph
from .trace import iter_hist_traces
from ..core.header import parse_header


class RelHistogramGraph(OneTableGraph,
                        HistogramGraph,
                        MultiRelsGraph,
                        ColorMapGraph,
                        ABC):
    """ Histogram of relationship(s) in one table. """

    @classmethod
    def get_cmap_type(cls):
        return RelColorMap

    @property
    def x_title(self):
        return self.data_kind.capitalize()

    @cached_property
    def data(self):
        # Fetch the raw data from the table.
        data = self._fetch_data(self.table,
                                k=self.k,
                                clust=self.clust)
        # Determine the edges of the bins.
        edges = self.get_edges(data)
        index = get_edges_index(edges, self.use_ratio)
        # Construct a histogram for each column.
        hist = dict()
        for col_name in data.columns:
            col_hist, _ = np.histogram(data[col_name], bins=edges)
            hist[col_name] = pd.Series(col_hist, index=index)
        return pd.DataFrame.from_dict(hist).reindex(columns=data.columns)

    @cached_property
    def data_header(self):
        """ Header of the selected data (not of the entire table). """
        return parse_header(self.data.columns)

    def get_traces(self):
        for row, index in zip(range(1, self.nrows + 1),
                              self.data_header.iter_clust_indexes(),
                              strict=True):
            for trace in iter_hist_traces(self.data.loc[:, index], self.cmap):
                yield (row, 1), trace


class RelHistogramWriter(OneTableWriter, ABC):

    @classmethod
    @abstractmethod
    def get_graph_type(cls) -> type[RelHistogramGraph]:
        """ Type of graph. """

    def get_graph(self, rels_group: str, **kwargs):
        graph_type = self.get_graph_type()
        return graph_type(table=self.table, rels=rels_group, **kwargs)


class RelHistogramRunner(OneTableRunner, HistogramRunner, ABC):
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
