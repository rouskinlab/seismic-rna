import os

from click import command

from .table import TableWriter, ReadTableRunner
from .hist import HistogramGraph, HistogramWriter, HistogramRunner
from ..core.run import log_command

COMMAND = __name__.split(os.path.extsep)[-1]


class ReadHistogramGraph(HistogramGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Histogram per read"

    @property
    def y_title(self):
        return "Number of reads"


class ReadHistogramWriter(HistogramWriter, TableWriter):

    @classmethod
    def get_graph_type(cls):
        return ReadHistogramGraph


class ReadHistogramRunner(HistogramRunner, ReadTableRunner):

    @classmethod
    def get_writer_type(cls):
        return ReadHistogramWriter

    @classmethod
    @log_command(COMMAND)
    def run(cls, *args, **kwargs):
        return super().run(*args, **kwargs)


@command(COMMAND, params=ReadHistogramRunner.params())
def cli(*args, **kwargs):
    """ Histogram of relationship(s) per read. """
    return ReadHistogramRunner.run(*args, **kwargs)
