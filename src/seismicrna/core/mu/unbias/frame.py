import pandas as pd

from .algo import calc_f_obs_numpy, calc_mu_adj_numpy, calc_prop_adj_numpy
from ...seq import Section


def _mus_to_matrix(mus: pd.Series | pd.DataFrame, section: Section):
    return mus.reindex(index=section.range,
                       fill_value=0.).values.reshape((section.length, -1))


def calc_f_obs_frame(mu_adj: pd.DataFrame | pd.Series,
                     section: Section,
                     min_gap: int):
    """ Calculate the observed fraction of reads in each cluster given
    their mutation rates adjusted for observer bias.

    Parameters
    ----------
    mu_adj: pd.DataFrame
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
    pd.Series | float
        For each cluster, the fraction of all bit vectors coming from
        that cluster that would be observed.
    """
    f_obs = calc_f_obs_numpy(_mus_to_matrix(mu_adj, section), min_gap)
    if isinstance(mu_adj, pd.DataFrame):
        return pd.Series(f_obs, index=mu_adj.columns)
    elif isinstance(mu_adj, pd.Series):
        return float(f_obs.item())
    raise TypeError(f"Expected mu_adj to be a Series or DataFrame, but got "
                    f"{type(mu_adj).__name__}")


def calc_mu_adj_frame(mu_obs: pd.DataFrame | pd.Series,
                      section: Section,
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
    mu_adj = calc_mu_adj_numpy(_mus_to_matrix(mu_obs, section), min_gap)
    if isinstance(mu_obs, pd.DataFrame):
        mu_adj_frame = pd.DataFrame(mu_adj, section.range, mu_obs.columns)
    elif isinstance(mu_obs, pd.Series):
        mu_adj_frame = pd.Series(mu_adj.reshape(section.length), section.range)
    else:
        raise TypeError(f"Expected mu_obs to be a Series or DataFrame, but got "
                        f"{type(mu_obs).__name__}")
    return mu_adj_frame.loc[mu_obs.index]


def calc_prop_adj_frame(prop_obs: pd.Series, f_obs: pd.Series):
    """ Calculate the adjusted proportion of the clusters given their
    observed proportions and the observer bias. """
    if not (clusters := prop_obs.index).equals(f_obs.index):
        raise ValueError(f"Got different clusters for observed proportions of "
                         f"clusters {clusters} and fraction of reads observed "
                         f"in each cluster {f_obs.index}")
    return pd.Series(calc_prop_adj_numpy(prop_obs.values, f_obs.values),
                     index=clusters)

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
