import os
from abc import ABC

from click import command

from .pospair import PositionPairGraph, PositionPairWriter, PositionPairRunner
from ..core.batch import calc_confusion_phi
from ..core.run import log_command

COMMAND = __name__.split(os.path.extsep)[-1]


class PositionCorrelationGraph(PositionPairGraph, ABC):
    """ Phi correlations between pairs of positions. """

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Phi correlation between"

    @classmethod
    def get_pair_func(cls):
        """ Function to compare each pair of positions. """
        return calc_confusion_phi


class PositionCorrelationWriter(PositionPairWriter):

    @classmethod
    def graph_type(cls):
        return PositionCorrelationGraph


class PositionCorrelationRunner(PositionPairRunner):

    @classmethod
    def get_writer_type(cls):
        return PositionCorrelationWriter

    @classmethod
    @log_command(COMMAND)
    def run(cls, *args, **kwargs):
        return super().run(*args, **kwargs)


@command(COMMAND, params=PositionCorrelationRunner.params())
def cli(*args, **kwargs):
    """ Phi correlations between pairs of positions. """
    return PositionCorrelationRunner.run(*args, **kwargs)
