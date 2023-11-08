import numpy as np
import pandas as pd

from ..seq import POS_NAME


def compute_roc_curve(paired: pd.Series, reacts: pd.Series):
    """ Compute the receiver operating characteristic (ROC) curve to
    compute how well chemical reactivities agree with a structure. """
    # Use only positions with non-missing reactivities.
    reacts_use = reacts.loc[np.logical_not(np.isnan(reacts.values))]
    # Count the number of positions with non-missing reactivities.
    n_use = reacts_use.size
    # Sort the positions in ascending order by reactivity.
    order = reacts_use.sort_values().index.get_level_values(POS_NAME)
    paired_use = paired.loc[order].astype(bool, copy=False)
    # Compute the cumulative number of paired bases at each threshold.
    paired_cumsum = np.hstack([[0], np.cumsum(paired_use)])
    # Count the paired positions.
    n_paired = paired_cumsum[-1]
    # Get the false positive rate: (paired and reactive) / paired.
    fpr = 1. - paired_cumsum / n_paired
    # Get the true positive rate: (unpaired and reactive) / unpaired.
    tpr = 1. - (np.arange(n_use + 1) - paired_cumsum) / (n_use - n_paired)
    # Traditionally, false positive rate is plotted on the x-axis and
    # true positive rate on the y-axis of an ROC curve.
    return fpr, tpr


def compute_auc_roc(fpr: np.ndarray, tpr: np.ndarray):
    """ Compute the area under the curve (AUC) of the receiver operating
    characteristic (ROC). """
    return -np.vdot(np.diff(fpr), tpr[1:])

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
