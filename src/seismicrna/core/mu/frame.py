"""

Mutation Rate Core Module

========================================================================

The functions in this module serve two main purposes:
 1. Adjust mutation rates to correct for observer bias.
 2. Normalize and winsorize mutation rates

------------------------------------------------------------------------

Adjust mutation rates to correct for observer bias

Our lab has found that pairs of mutations in DMS-MaPseq data rarely have
fewer than three non-mutated bases separating them. We suspect that the
reverse transcriptase is prone to falling off or stalling at locations
in the RNA where DMS has methylated bases that are too close.

Regardless of the reason, mutations on nearby bases are not independent,
which violates a core assumption of the Bernoulli mixture model that we
use in the expectation-maximization clustering algorithm. Namely, that
mutations occur independently of each other, such that the likelihood of
observing a bit vector is the product of the mutation rate of each base
that is mutated and one minus the mutation rate ("non-mutation rate") of
each base that is not mutated.

In order to use the Bernoulli mixture model for expectation-maximization
clustering, we modify it such that bases separated by zero, one, or two
other bases are no longer assumed to mutate independently. Specifically,
since pairs of close mutations are rare, we add the assumption that no
mutations separated by fewer than three non-mutated bases may occur, and
exclude the few bit vectors that have such close mutations.

When these bit vectors are assumed to exist in the original population
of RNA molecules but not appear in the observed data, the mutation rates
that are observed will differ from the real, underlying mutation rates,
which would include the unobserved bit vectors. The real mutation rates
must therefore be estimated from the observed mutation rates.

It is relatively easy to estimate the observed mutation rates given the
real mutation rates, but there is no closed-form solution that we know
of to estimate the real mutation rates from the observed mutation rates.
Thus, we use an iterative approach, based on Newton's method for finding
the roots of functions. We initially guess the real mutation rates. Each
iteration, we estimate the mutation rates that would have been observed
given our current guess of the real mutation rates, and then subtract
the mutation rates that were actually observed. This difference would be
zero if we had accurately guessed the real mutation rates. Thus, we use
Newton's method to solve for the real mutation rates that minimize this
difference. The details are described in the comments of this module and
in Tomezsko et al. (2020) (https://doi.org/10.1038/s41586-020-2253-5).

------------------------------------------------------------------------

Normalize and winsorize mutation rates

The mutation rates of an RNA may vary among different samples because of
variations in the chemical probing and mutational profiling procedure.
Thus, it is often helpful to compute "normalized" mutation rates that
can be compared directly across different samples and used to predict
secondary structures.

This module currently provides a simple method for normalizing mutation
rates. First, a specific quantile of the dataset is chosen, such as 0.95
(i.e. the 95th percentile). The mutation rate with this quantile is set
to 1.0, and all other mutation rates are scaled linearly.

If the chosen quantile is less than 1.0, then any mutation rates above
the quantile will be scaled to values greater than 1.0. Since these high
mutation rates may be exceptionally reactive bases, it is often helpful
to cap the normalized mutation rates to a maximum of 1.0. The winsorize
function in this module performs normalization and then sets any value
greater than 1.0 to 1.0.

"""

from logging import getLogger

import pandas as pd

from .unbias import calc_f_obs_numpy, calc_mu_adj_numpy, calc_prop_adj_numpy
from ..seq import Section

logger = getLogger(__name__)


def _mus_to_matrix(mus: pd.Series | pd.DataFrame, section: Section):
    return mus.reindex(index=section.range,
                       fill_value=0.).values.reshape((section.length, -1))


def calc_f_obs_frame(mu_adj: pd.DataFrame | pd.Series,
                     section: Section,
                     min_gap: int):
    """
    Calculate the observed fraction of reads in each cluster given their
    mutation rates adjusted for observer bias.

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
    """
    Correct the mutation rates of a DataFrame for observer bias.

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
