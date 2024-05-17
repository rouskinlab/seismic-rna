import os
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd
from click import command

from .base import PosGraphRunner, PosGraphWriter
from .onestruct import (StructOneTableGraph,
                        StructOneTableRunner,
                        StructOneTableWriter)
from .trace import iter_roc_traces

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]

# Index level names.
AXIS_NAME = "Axis"
PROFILE_NAME = "Profile"
STRUCT_NAME = "Structure"
COL_NAMES = PROFILE_NAME, STRUCT_NAME

# Axis names.
TPR = "True positive rate"
FPR = "False positive rate"


def rename_columns(df: pd.DataFrame):
    """ Rename the levels of the columns. """
    # The DataFrame's columns must be a MultiIndex with two levels named
    # "Profile" and "Structure".
    if df.size > 0:
        # If the DataFrame has at least one column, then it will have a
        # two-level MultiIndex already: just rename the column levels.
        df.columns.names = COL_NAMES
    else:
        # If it is empty, then its columns will default to a RangeIndex
        # (which has one level), so they must be replaced.
        df.columns = pd.MultiIndex.from_arrays([[], []], names=COL_NAMES)
    return df


def _consolidate_pr(pr: dict):
    """ Consolidate a true or false positive rate (PR) forming half the
    ROC from a dict into a DataFrame. """
    return rename_columns(pd.DataFrame.from_dict(pr))


class ROCGraph(StructOneTableGraph):
    """ Graph of a receiver operating characteristic (ROC) curve. """

    @classmethod
    def graph_kind(cls):
        return COMMAND

    @classmethod
    def what(cls):
        return "ROC curve"

    @property
    def x_title(self):
        return FPR

    @property
    def y_title(self):
        return TPR

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
        if not fpr or not tpr:
            raise ValueError(f"Got no data for {self}")
        # Consolidate the FPR and TPR data into two DataFrames.
        return _consolidate_pr(fpr), _consolidate_pr(tpr)

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
        axes = np.hstack([np.repeat([FPR], self.fpr.columns.size),
                          np.repeat([TPR], self.tpr.columns.size)])
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
        """ Names of the profiles as they appear in the data. """
        profile_names = self.fpr.columns.unique(PROFILE_NAME)
        if not profile_names.equals(self.tpr.columns.unique(PROFILE_NAME)):
            raise ValueError(f"Profile names differ: {profile_names} "
                             f"≠ {self.tpr.columns.unique(PROFILE_NAME)}")
        return profile_names

    def get_traces(self):
        for row, profile in enumerate(self.profile_names, start=1):
            for trace in iter_roc_traces(self.fpr.loc[:, profile],
                                         self.tpr.loc[:, profile],
                                         profile):
                yield (row, 1), trace


class ROCWriter(StructOneTableWriter, PosGraphWriter):

    def get_graph(self, rels_group: str, **kwargs):
        return ROCGraph(table=self.table, rel=rels_group, **kwargs)


class ROCRunner(StructOneTableRunner, PosGraphRunner):

    @classmethod
    def get_writer_type(cls):
        return ROCWriter


@command(COMMAND, params=ROCRunner.params())
def cli(*args, **kwargs):
    """ ROC curve comparing a profile to a structure. """
    return ROCRunner.run(*args, **kwargs)

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
