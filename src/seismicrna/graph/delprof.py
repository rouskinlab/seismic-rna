import os

import numpy as np
from click import command
from plotly import graph_objects as go

from .color import ColorMapGraph, SeqColorMap
from .rel import OneRelGraph
from .table import PositionTableRunner
from .trace import iter_seq_base_bar_traces
from .twotable import (TwoTableMergedClusterGroupGraph,
                       TwoTableRelClusterGroupWriter,
                       TwoTableRelClusterGroupRunner)
from ..core.run import log_command
from ..core.seq import POS_NAME

COMMAND = __name__.split(os.path.extsep)[-1]


class DeltaProfileGraph(TwoTableMergedClusterGroupGraph,
                        OneRelGraph,
                        ColorMapGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Delta profile graph"

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    @classmethod
    def _trace_function(cls):
        return iter_seq_base_bar_traces

    @property
    def x_title(self) -> str:
        return POS_NAME

    @property
    def _trace_kwargs(self):
        return super()._trace_kwargs | dict(cmap=self.cmap)

    @property
    def y_title(self):
        return f"Difference between {self.data_kind}s 1 and 2"

    @property
    def _merge_data(self):
        """ Compute the difference between the profiles. """
        return np.subtract

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class DeltaProfileWriter(TwoTableRelClusterGroupWriter):

    @classmethod
    def get_graph_type(cls):
        return DeltaProfileGraph


class DeltaProfileRunner(TwoTableRelClusterGroupRunner, PositionTableRunner):

    @classmethod
    def get_writer_type(cls):
        return DeltaProfileWriter

    @classmethod
    @log_command(COMMAND)
    def run(cls, *args, **kwargs):
        return super().run(*args, **kwargs)


@command(COMMAND, params=DeltaProfileRunner.params())
def cli(*args, **kwargs):
    """ Bar graph of differences between two profiles per position. """
    return DeltaProfileRunner.run(*args, **kwargs)
