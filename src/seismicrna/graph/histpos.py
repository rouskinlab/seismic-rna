import os

from click import command

from .table import TableWriter, PositionTableRunner
from .hist import HistogramGraph, HistogramWriter, HistogramRunner
from ..core.run import log_command

COMMAND = __name__.split(os.path.extsep)[-1]


class PositionHistogramGraph(HistogramGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Histogram per position"

    @property
    def y_title(self):
        return "Number of positions"


class PositionHistogramWriter(HistogramWriter, TableWriter):

    @classmethod
    def get_graph_type(cls):
        return PositionHistogramGraph


class PositionHistogramRunner(HistogramRunner, PositionTableRunner):

    @classmethod
    def get_writer_type(cls):
        return PositionHistogramWriter

    @classmethod
    @log_command(COMMAND)
    def run(cls, *args, **kwargs):
        return super().run(*args, **kwargs)


@command(COMMAND, params=PositionHistogramRunner.params())
def cli(*args, **kwargs):
    """ Histogram of relationship(s) per position. """
    return PositionHistogramRunner.run(*args, **kwargs)
