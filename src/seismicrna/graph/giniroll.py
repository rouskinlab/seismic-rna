import os

from click import command

from .statroll import RollingStatGraph, RollingStatRunner, RollingStatWriter
from ..core.mu import calc_gini

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


@command(COMMAND, params=RollingGiniRunner.params())
def cli(*args, **kwargs):
    """ Rolling Gini coefficient. """
    return RollingGiniRunner.run(*args, **kwargs)

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
