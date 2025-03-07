import os
from functools import cached_property, partial

from click import command
from plotly import graph_objects as go

from .rel import OneRelGraph
from .roll import RollingGraph, RollingRunner
from .table import PositionTableRunner
from .trace import iter_seq_line_traces
from .twotable import (TwoTableMergedClusterGroupGraph,
                       TwoTableRelClusterGroupWriter,
                       TwoTableRelClusterGroupRunner)
from ..core.arg import opt_metric
from ..core.mu import compare_windows, get_comp_name
from ..core.run import log_command

COMMAND = __name__.split(os.path.extsep)[-1]


class RollingCorrelationGraph(TwoTableMergedClusterGroupGraph,
                              OneRelGraph,
                              RollingGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Rolling correlation"

    @classmethod
    def _trace_function(cls):
        return iter_seq_line_traces

    def __init__(self, *, metric: str, **kwargs):
        super().__init__(**kwargs)
        self._metric = metric

    @cached_property
    def predicate(self):
        return super().predicate + [self._metric]

    @cached_property
    def details(self):
        return super().details + [f"metric = {self._metric.upper()}"]

    @cached_property
    def y_title(self):
        return f"{get_comp_name(self._metric)} of {self.data_kind}s"

    @cached_property
    def _merge_data(self):
        return partial(compare_windows,
                       method=self._metric,
                       size=self._size,
                       min_count=self._min_count)

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class RollingCorrelationWriter(TwoTableRelClusterGroupWriter):

    @classmethod
    def get_graph_type(cls):
        return RollingCorrelationGraph


class RollingCorrelationRunner(TwoTableRelClusterGroupRunner,
                               RollingRunner,
                               PositionTableRunner):

    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_metric]

    @classmethod
    def get_writer_type(cls):
        return RollingCorrelationWriter

    @classmethod
    @log_command(COMMAND)
    def run(cls, *args, **kwargs):
        return super().run(*args, **kwargs)


@command(COMMAND, params=RollingCorrelationRunner.params())
def cli(*args, **kwargs):
    """ Rolling correlation/comparison of two profiles. """
    return RollingCorrelationRunner.run(*args, **kwargs)
