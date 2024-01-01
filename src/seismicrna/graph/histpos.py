import os
from logging import getLogger

from click import command

from .base import PosGraphWriter, PosGraphRunner
from .histrel import RelHistogramGraph, RelHistogramWriter, RelHistogramRunner

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


class PosHistogramGraph(RelHistogramGraph):

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Histogram per position"

    @property
    def y_title(self):
        return "Number of positions"


class PosHistogramWriter(RelHistogramWriter, PosGraphWriter):

    @classmethod
    def get_graph_type(cls):
        return PosHistogramGraph


class PosHistogramRunner(RelHistogramRunner, PosGraphRunner):

    @classmethod
    def get_writer_type(cls):
        return PosHistogramWriter


@command(COMMAND, params=PosHistogramRunner.params())
def cli(*args, **kwargs):
    """ Histogram of relationship(s) per position. """
    return PosHistogramRunner.run(*args, **kwargs)

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
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
