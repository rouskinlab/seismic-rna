import numpy as np

from ..core.array import find_dims
from ..core.unbias import (READS,
                           POSITIONS,
                           CLUSTERS,
                           calc_p_noclose_given_clust,
                           calc_p_noclose_given_ends_auto,
                           calc_p_ends_given_clust_noclose,
                           calc_p_clust_given_noclose,
                           triu_log)


def _zero_masked(p_mut: np.ndarray, unmasked: np.ndarray):
    """ Set mutation rates of masked positions to zero. """
    p_mut_unmasked = np.zeros_like(p_mut)
    p_mut_unmasked[unmasked] = p_mut[unmasked]
    return p_mut_unmasked


def _calc_logp_joint(p_mut: np.ndarray,
                     p_ends: np.ndarray,
                     p_clust: np.ndarray,
                     end5s: np.ndarray,
                     end3s: np.ndarray,
                     unmasked: np.ndarray,
                     muts_per_pos: list[np.ndarray],
                     min_mut_gap: int):
    # Validate the dimensions.
    find_dims([(POSITIONS, CLUSTERS),
               (POSITIONS, POSITIONS),
               (CLUSTERS,),
               (READS,),
               (READS,)],
              [p_mut, p_ends, p_clust, end5s, end3s],
              ["p_mut", "p_ends", "p_clust", "end5s", "end3s"],
              nonzero=True)
    # Ensure the mutation rates of masked positions are 0.
    p_mut = _zero_masked(p_mut, unmasked)
    # Calculate the end probabilities.
    p_noclose_given_ends = calc_p_noclose_given_ends_auto(p_mut, min_mut_gap)
    p_ends_given_clust_noclose = calc_p_ends_given_clust_noclose(
        p_ends,
        p_noclose_given_ends
    )
    # Calculate the cluster probabilities.
    p_noclose = calc_p_noclose_given_clust(p_ends, p_noclose_given_ends)
    p_clust_given_noclose = calc_p_clust_given_noclose(p_clust, p_noclose)
    # Compute the probability that a read would have no two mutations
    # too close given its end coordinates.
    logp_noclose_given_ends = triu_log(calc_p_noclose_given_ends_auto(p_mut, min_mut_gap))
    # Compute the logs of the parameters.
    with np.errstate(divide="ignore"):
        # Suppress warnings about taking the log of zero, which is a
        # valid mutation rate.
        logp_mut = np.log(p_mut)
    logp_not = np.log(1. - p_mut)
    logp_ends_given_clust_noclose = triu_log(p_ends_given_clust_noclose)
    logp_clust_given_noclose = np.log(p_clust_given_noclose)
    # For each cluster, calculate the probability that a read up to and
    # including each position would have no mutations.
    logp_nomut_incl = np.cumsum(logp_not, axis=0)
    # For each cluster, calculate the probability that a read up to but
    # not including each position would have no mutations.
    logp_nomut_excl = np.vstack([np.zeros_like(p_clust), logp_nomut_incl[:-1]])
    # For each unique read, calculate the probability that a random read
    # with the same end coordinates would have no mutations.
    logp_nomut_given_clust = logp_nomut_incl[end3s] - logp_nomut_excl[end5s]
    # For each unique read, compute the joint probability that a random
    # read would have the same end coordinates, come from each cluster,
    # and have no mutations; normalize by the fraction of all reads that
    # have no two mutations too close.
    # 2D (unique reads x clusters)
    logp_joint = (logp_clust_given_noclose[np.newaxis, :]
                  + logp_ends_given_clust_noclose[end5s, end3s]
                  - logp_noclose_given_ends[end5s, end3s]
                  + logp_nomut_given_clust)
    # For each unique read, compute the likelihood of observing it
    # (including its mutations) by adjusting the above likelihood
    # of observing the end coordinates with no mutations.
    for j, mut_reads in zip(unmasked, muts_per_pos, strict=True):
        logp_joint[mut_reads] += logp_mut[j] - logp_not[j]
    return logp_joint


def _calc_logp_marginal(logp_joint: np.ndarray):
    # For each unique observed read, the marginal probability that a
    # random read would have the end coordinates and mutations, no
    # matter which cluster it came from, is the sum of the joint
    # probability over all clusters (axis 1).
    # 1D (unique reads)
    return np.logaddexp.reduce(logp_joint, axis=1)


def calc_marginal(*args, **kwargs):
    return _calc_logp_marginal(_calc_logp_joint(*args, **kwargs))


def calc_marginal_resps(*args, **kwargs):
    logp_joint = _calc_logp_joint(*args, **kwargs)
    logp_marginal = _calc_logp_marginal(logp_joint)
    # Calculate the posterior probability that each read came from
    # each cluster by dividing the joint probability (observing the
    # read and coming from the cluster) by the marginal probability
    # (observing the read in any cluster).
    # 2D (unique reads x clusters)
    resps = np.exp(logp_joint - logp_marginal[:, np.newaxis])
    return logp_marginal, resps

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
