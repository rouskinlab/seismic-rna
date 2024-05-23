import os

import numpy as np
import pandas as pd
from click import command

from .statroll import RollingStatGraph, RollingStatRunner, RollingStatWriter
from ..core.mu import calc_signal_noise
from ..core.seq import BASEA, BASEC, BASE_NAME

COMMAND = __name__.split(os.path.extsep)[-1]


class RollingSNRGraph(RollingStatGraph):

    @classmethod
    def stat_func(cls):

        def calc_snr(data: pd.Series):
            bases = data.index.get_level_values(BASE_NAME)
            signal = np.logical_or(bases == BASEA, bases == BASEC)
            return calc_signal_noise(data, signal)

        return calc_snr

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "Rolling signal-to-noise ratio"

    @property
    def y_title(self):
        return "Signal-to-noise ratio"


class RollingSNRWriter(RollingStatWriter):

    @classmethod
    def get_graph_type(cls):
        return RollingSNRGraph


class RollingSNRRunner(RollingStatRunner):

    @classmethod
    def get_writer_type(cls):
        return RollingSNRWriter


@command(COMMAND, params=RollingSNRRunner.params())
def cli(*args, **kwargs):
    """ Rolling signal-to-noise ratio. """
    return RollingSNRRunner.run(*args, **kwargs)

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
