import os
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd
from click import command

from .onetable import OneTableStructureGraph, OneTableRunner, OneTableWriter
from .traces import iter_roc_traces

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]

# Index level names.
AXIS_NAME = "Axis"
PROFILE_NAME = "Profile"
STRUCT_NAME = "Structure"


class ROCGraph(OneTableStructureGraph):
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

    @cached_property
    def _roc(self):
        """ ROC curve as two DataFrames of false positive rates (FPR)
        and true positive rates (TPR). """
        # Gather the FPR and TPR data from each RNA state.
        fpr = dict()
        tpr = dict()
        for state in self.iter_states():
            key = state.data_name, state.title
            if key in fpr or key in tpr:
                raise ValueError(f"Duplicate RNA state: {key}")
            fpr[key], tpr[key] = state.roc
        # Consolidate the FPR and TPR data into two DataFrames.
        fpr = pd.DataFrame.from_dict(fpr)
        fpr.columns.names = PROFILE_NAME, STRUCT_NAME
        tpr = pd.DataFrame.from_dict(tpr)
        tpr.columns.names = PROFILE_NAME, STRUCT_NAME
        return fpr, tpr

    @property
    def fpr(self):
        """ False positive rate (FPR) of each RNA state. """
        fpr, _ = self._roc
        return fpr

    @property
    def tpr(self):
        """ True positive rate (TPR) of each RNA state. """
        _, tpr = self._roc
        return tpr

    @cached_property
    def data(self):
        # Join the FPR and TPR data horizontally.
        data = pd.concat([self.fpr, self.tpr], axis=1, join="inner")
        # Add the axis name as the first level of the columns.
        axes = np.hstack([np.repeat([self.x_title], self.fpr.columns.size),
                          np.repeat([self.y_title], self.tpr.columns.size)])
        names = [AXIS_NAME] + list(data.columns.names)
        data.columns = pd.MultiIndex.from_arrays(
            [(axes if name == AXIS_NAME
              else data.columns.get_level_values(name).values)
             for name in names],
            names=names
        )
        return data

    @cached_property
    def profile_names(self):
        """ Names of the profiles as they appear in the ROC data. """
        profile_names = self.fpr.columns.unique(PROFILE_NAME)
        if not profile_names.equals(self.tpr.columns.unique(PROFILE_NAME)):
            raise ValueError(f"Profile names differ: {profile_names} "
                             f"≠ {self.tpr.columns.unique(PROFILE_NAME)}")
        return profile_names

    def get_traces(self):
        for row, profile in enumerate(self.profile_names, start=1):
            for trace in iter_roc_traces(self.fpr.loc[:, profile],
                                         self.tpr.loc[:, profile]):
                yield (row, 1), trace


class ROCWriter(OneTableWriter):

    @classmethod
    def get_graph_type(cls, dummy):
        # FIXME: remove the dummy argument.
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
