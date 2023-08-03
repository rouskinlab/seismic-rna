import os
from logging import getLogger

import pandas as pd
from click import command
from plotly import graph_objects as go

from .seqpair import SeqPairGraphRunner, SeqPairGraphWriter, SeqPairOneAxisGraph
from .traces import iter_seq_base_bar_traces

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class SeqDiffGraphRunner(SeqPairGraphRunner):

    @classmethod
    def writer_type(cls):
        return SeqDiffGraphWriter


@command(COMMAND, params=SeqPairGraphRunner.params)
def cli(*args, **kwargs):
    """ Create bar graphs of differences between pairs of samples at
    each position in a sequence. """
    return SeqDiffGraphRunner.run(*args, **kwargs)


class SeqDiffGraph(SeqPairOneAxisGraph):

    @classmethod
    def graph_type(cls):
        return COMMAND

    @property
    def y_title(self):
        return f"{self.quantity}-2 minus {self.quantity}-1"

    @classmethod
    def _trace_function(cls):
        return iter_seq_base_bar_traces

    @property
    def _merge_data(self):

        def diff(vals1: pd.Series, vals2: pd.Series):
            """ Compute the difference between the Series. """
            return vals2 - vals1

        return diff

    def _figure_layout(self, fig: go.Figure):
        super()._figure_layout(fig)
        fig.update_yaxes(gridcolor="#d0d0d0")


class SeqDiffGraphWriter(SeqPairGraphWriter):

    @property
    def graph_type(self):
        return SeqDiffGraph
