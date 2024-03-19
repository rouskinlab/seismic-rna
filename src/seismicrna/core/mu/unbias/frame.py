import numpy as np
import pandas as pd

from .algo import (calc_p_noclose_given_ends_numpy,
                   calc_params_numpy)
from ...seq import Section


def _data_to_matrix(data: pd.Series | pd.DataFrame, section: Section):
    return data.reindex(index=section.range,
                        fill_value=0.).values.reshape((section.length, -1))


def calc_p_noclose_given_ends_frame(section: Section,
                                    p_mut_given_span: pd.DataFrame | pd.Series,
                                    min_gap: int):
    """ Calculate the observed fraction of reads in each cluster given
    their mutation rates adjusted for observer bias.

    Parameters
    ----------
    p_mut_given_span: pd.DataFrame
        Adjusted fraction of mutated bits at each non-excluded position
        (index) in each cluster (column). All values must be in [0, 1).
    section: Section
        The section over which to compute mutation rates, including all
        positions that were excluded.
    min_gap: int
        Minimum number of non-mutated bases between two mutations.
        Must be ≥ 0.

    Returns
    -------
    np.ndarray
        For each cluster, the fraction of all bit vectors coming from
        that cluster that would be observed.
    """
    return calc_p_noclose_given_ends_numpy(
        _data_to_matrix(p_mut_given_span, section), min_gap
    )


def calc_params_frame(section: Section,
                      p_mut_given_span_noclose: pd.DataFrame | pd.Series,
                      p_ends_given_noclose: np.ndarray,
                      p_clust_given_noclose: pd.Series | float,
                      min_gap: int):
    """ Correct the mutation rates of a DataFrame for observer bias.

    Parameters
    ----------
    mu_obs: pd.DataFrame
        Fraction of mutated bits at each non-excluded position (index)
        in each cluster (column). All values must be in [0, 1).
    section: Section
        The section over which to compute mutation rates, including all
        positions that were excluded.
    min_gap: int
        Minimum number of non-mutated bases between two mutations.
        Must be ≥ 0.

    Returns
    -------
    pd.DataFrame | pd.Series
        Data frame of the adjusted mutation rates with the same index
        and columns as `mu_obs`.
    """
    p_mut_given_span, p_ends, p_clust = calc_params_numpy(
        _data_to_matrix(p_mut_given_span_noclose, section),
        p_ends_given_noclose,
        np.atleast_1d(p_clust_given_noclose),
        min_gap
    )
    if isinstance(p_mut_given_span_noclose, pd.DataFrame):
        p_mut_given_span_frame = pd.DataFrame(
            p_mut_given_span,
            section.range,
            p_mut_given_span_noclose.columns
        )
    elif isinstance(p_mut_given_span_noclose, pd.Series):
        p_mut_given_span_frame = pd.Series(
            p_mut_given_span.reshape(section.length),
            section.range
        )
    else:
        raise TypeError(
            "p_mut_given_span_noclose must be a Series or DataFrame, "
            f"but got {type(p_mut_given_span_noclose).__name__}"
        )
    p_mut_given_span_frame = (
        p_mut_given_span_frame.loc[p_mut_given_span_noclose.index]
    )
    return p_mut_given_span_frame, p_ends, p_clust

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
