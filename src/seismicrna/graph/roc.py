import os
from abc import ABC
from logging import getLogger

from click import command

from .onetable import OneTableGraph, OneTableRunner, OneTableWriter
from .traces import iter_seq_base_bar_traces

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class ROCGraph(OneTableGraph, ABC):
    """ Graph of a receiver operating characteristic (ROC) curve. """

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "ROC curve"

    @property
    def x_title(self):
        return "False positive rate"

    @property
    def y_title(self):
        return "True positive rate"


class EndoROCGraph(ROCGraph):

    def get_traces(self):
        for row, (_, values) in enumerate(self.data.items(), start=1):
            for trace in iter_seq_base_bar_traces(values):
                yield (row, 1), trace


class ExoROCGraph(ROCGraph):

    def get_traces(self):
        for row, (_, values) in enumerate(self.data.items(), start=1):
            for trace in iter_seq_base_bar_traces(values):
                yield (row, 1), trace


class ROCWriter(OneTableWriter):

    @classmethod
    def get_graph_type(cls, ):
        return ROCGraph


class ROCRunner(OneTableRunner):

    @classmethod
    def writer_type(cls):
        return ROCWriter


@command(ROCGraph.graph_kind(),
         params=ROCRunner.params())
def cli(*args, **kwargs):
    """ Create bar graphs of positions in a sequence. """
    return ROCRunner.run(*args, **kwargs)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
