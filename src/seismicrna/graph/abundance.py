import os
from collections import defaultdict
from functools import cached_property
from plotly import graph_objects as go

import pandas as pd
from click import command

from .onetable import OneTableGraph, OneTableWriter, OneTableRunner
from .table import AbundanceTableRunner
from .trace import iter_stack_bar_traces
from ..core.header import NUM_CLUSTS_NAME, CLUST_NAME
from ..core.run import log_command

COMMAND = __name__.split(os.path.extsep)[-1]


class ClusterAbundanceGraph(OneTableGraph):
    """ Abundance of each cluster. """

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Cluster abundance"

    @property
    def details(self):
        return list()

    @property
    def predicate(self):
        return list()

    @property
    def path_subject(self):
        return self.action

    @cached_property
    def _title_main(self):
        return [f"{self.what()} {self.data_kind}s "
                f"in {self.title_action_sample} "
                f"over reference {repr(self.ref)} "
                f"region {repr(self.reg)}"]

    @property
    def x_title(self):
        return "Number of clusters (K)"

    @property
    def y_title(self):
        return "Abundance of each cluster"

    @property
    def row_tracks(self):
        return None

    @property
    def col_tracks(self):
        return None

    @cached_property
    def data(self):
        if self.use_ratio:
            table_data = self.table.proportions
        else:
            table_data = self.table.data
        data_dict = defaultdict(dict)
        for k, clust in table_data.index:
            data_dict[k][clust] = table_data.at[(k, clust)]
        data_frame = pd.DataFrame.from_dict(data_dict, orient="index")
        data_frame.index.name = NUM_CLUSTS_NAME
        data_frame.columns.name = CLUST_NAME
        return data_frame

    def get_traces(self):
        for trace in iter_stack_bar_traces(self.data):
            yield (1, 1), trace

    def _figure_layout(self, figure: go.Figure):
        super()._figure_layout(figure)
        # Stack the bars at each position.
        figure.update_layout(barmode="stack")


class ClusterAbundanceWriter(OneTableWriter):

    def get_graph(self, **kwargs):
        return ClusterAbundanceGraph(table=self.table, **kwargs)

    def iter_graphs(self, **kwargs):
        yield self.get_graph(**kwargs)


class ClusterAbundanceRunner(OneTableRunner, AbundanceTableRunner):

    @classmethod
    def get_writer_type(cls):
        return ClusterAbundanceWriter
    
    @classmethod
    @log_command(COMMAND)
    def run(cls, *args, **kwargs):
        return super().run(*args, **kwargs)


@command(COMMAND, params=ClusterAbundanceRunner.params())
def cli(*args, **kwargs):
    """ Abundance of each cluster. """
    return ClusterAbundanceRunner.run(*args, **kwargs)
