import os
from functools import cached_property
from logging import getLogger

import pandas as pd
from click import command
from plotly import graph_objects as go

from .base import PosGraphWriter, PosGraphRunner
from .onestruct import (StructOneTableGraph,
                        StructOneTableRunner,
                        StructOneTableWriter)
from .roc import PROFILE_NAME, rename_columns
from .roll import RollingGraph, RollingRunner
from .trace import iter_rolling_auc_traces

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class RollingAUCGraph(StructOneTableGraph, RollingGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Rolling AUC-ROC"

    @property
    def y_title(self):
        return "AUC-ROC"

    @cached_property
    def data(self):
        # Collect the rolling AUC-ROC from every RNA state.
        data = dict()
        for state in self.iter_states():
            key = state.data_name, state.title
            if key in data:
                raise ValueError(f"Duplicate RNA state: {key}")
            data[key] = state.rolling_auc(self._size, self._min_count)
        if not data:
            raise ValueError(f"Got no data for {self}")
        # Covert the data into a DataFrame and rename the column levels.
        return rename_columns(pd.DataFrame.from_dict(data))

    @cached_property
    def profile_names(self):
        """ Names of the profiles as they appear in the data. """
        return self.data.columns.unique(PROFILE_NAME)

    def get_traces(self):
        for row, profile in enumerate(self.profile_names, start=1):
            for trace in iter_rolling_auc_traces(self.data.loc[:, profile],
                                                 profile):
                yield (row, 1), trace

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class RollingAUCWriter(StructOneTableWriter, PosGraphWriter):

    def get_graph(self, rels_group: str, **kwargs):
        return RollingAUCGraph(table=self.table, rel=rels_group, **kwargs)


class RollingAUCRunner(RollingRunner, StructOneTableRunner, PosGraphRunner):

    @classmethod
    def get_writer_type(cls):
        return RollingAUCWriter


@command(COMMAND, params=RollingAUCRunner.params())
def cli(*args, **kwargs):
    """ Rolling AUC-ROC comparing a profile to a structure. """
    return RollingAUCRunner.run(*args, **kwargs)

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
