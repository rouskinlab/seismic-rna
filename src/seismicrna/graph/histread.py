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

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
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
