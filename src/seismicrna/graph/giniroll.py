import os

from click import command

from .statroll import RollingStatGraph, RollingStatRunner, RollingStatWriter
from ..core.mu import calc_gini
from ..core.run import log_command

COMMAND = __name__.split(os.path.extsep)[-1]


class RollingGiniGraph(RollingStatGraph):

    @classmethod
    def stat_func(cls):
        return calc_gini

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Rolling Gini coefficient"

    @property
    def y_title(self):
        return "Gini coefficient"


class RollingGiniWriter(RollingStatWriter):

    @classmethod
    def get_graph_type(cls):
        return RollingGiniGraph


class RollingGiniRunner(RollingStatRunner):

    @classmethod
    def get_writer_type(cls):
        return RollingGiniWriter

    @classmethod
    @log_command(COMMAND)
    def run(cls, *args, **kwargs):
        return super().run(*args, **kwargs)


@command(COMMAND, params=RollingGiniRunner.params())
def cli(*args, **kwargs):
    """ Rolling Gini coefficient. """
    return RollingGiniRunner.run(*args, **kwargs)
