from functools import cached_property
from itertools import combinations, filterfalse
from logging import getLogger
from typing import Callable

import numpy as np
import pandas as pd

from .names import CLUST_PROP_NAME
from .uniq import UniqReads
from ..core.array import find_dims, get_length
from ..core.header import ClustHeader
from ..core.mu import (READS,
                       POSITIONS,
                       CLUSTERS,
                       calc_p_noclose,
                       calc_p_noclose_given_ends,
                       calc_params,
                       calc_params_observed,
                       calc_p_ends_given_noclose,
                       calc_p_clust_given_noclose,
                       calc_pearson,
                       calc_nrmsd,
                       triu_log)

logger = getLogger(__name__)

LOG_LIKE_PRECISION = 3  # number of digits to round the log likelihood


def _calc_bic(n_params: int,
              n_data: int,
              log_like: float,
              min_data_param_ratio: float = 10.):
    """ Compute the Bayesian Information Criterion (BIC) of a model.
    Typically, the model with the smallest BIC is preferred.

    Parameters
    ----------
    n_params: int
        Number of parameters that the model estimates
    n_data: int
        Number of data points in the sample from which the parameters
        were estimated
    log_like: float
        Natural logarithm of the likelihood of observing the data given
        the parameters
    min_data_param_ratio: float = 10.0
        In order for the BIC approximation to be valid, the sample size
        must be much larger than the number of estimated parameters.
        Issue a warning if the sample size is less than this ratio times
        the number of parameters, but still compute and return the BIC.

    Returns
    -------
    float
        Bayesian Information Criterion (BIC)
    """
    if log_like > 0.:
        raise ValueError(f"log_like must be ≤ 0, but got {log_like}")
    if n_params < 0:
        raise ValueError(f"n_params must be ≥ 0, but got {n_params}")
    if n_data < 0:
        raise ValueError(f"n_data must be ≥ 0, but got {n_data}")
    if n_data == 0:
        logger.warning("The Bayesian Information Criterion (BIC) is undefined "
                       "for 0 data points")
    if n_data < min_data_param_ratio * n_params:
        logger.warning(f"The Bayesian Information Criterion (BIC) uses an "
                       f"approximation that is valid only when the size of the "
                       f"sample (n = {n_data}) is much larger than the number "
                       f"of parameters being estimated (p = {n_params}). "
                       f"This model does not meet this criterion, so the BIC "
                       f"may not indicate the model's complexity accurately.")
    with np.errstate(divide="ignore", invalid="ignore"):
        return n_params * np.log(n_data) - 2. * log_like


def _zero_masked(p_mut: np.ndarray, unmasked: np.ndarray):
    """ Set mutation rates of masked positions to zero. """
    p_mut_unmasked = np.zeros_like(p_mut)
    p_mut_unmasked[unmasked] = p_mut[unmasked]
    return p_mut_unmasked


def _expectation(p_mut: np.ndarray,
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
    p_noclose_given_ends = calc_p_noclose_given_ends(p_mut, min_mut_gap)
    p_ends_given_noclose = calc_p_ends_given_noclose(p_ends,
                                                     p_noclose_given_ends)
    # Calculate the cluster probabilities.
    p_noclose = calc_p_noclose(p_ends, p_noclose_given_ends)
    p_clust_given_noclose = calc_p_clust_given_noclose(p_clust, p_noclose)
    # Compute the probability that a read would have no two mutations
    # too close given its end coordinates.
    logp_noclose = triu_log(calc_p_noclose_given_ends(p_mut, min_mut_gap))
    # Compute the logs of the parameters.
    with np.errstate(divide="ignore"):
        # Suppress warnings about taking the log of zero, which is a
        # valid mutation rate.
        logp_mut = np.log(p_mut)
    logp_not = np.log(1. - p_mut)
    logp_ends_given_noclose = triu_log(p_ends_given_noclose)
    logp_clust_given_noclose = np.log(p_clust_given_noclose)
    # For each cluster, calculate the probability that a read up to and
    # including each position would have no mutations.
    logp_nomut_incl = np.cumsum(logp_not, axis=0)
    # For each cluster, calculate the probability that a read up to but
    # not including each position would have no mutations.
    logp_nomut_excl = np.vstack([np.zeros_like(p_clust), logp_nomut_incl[:-1]])
    # For each unique read, calculate the probability that a random read
    # with the same end coordinates would have no mutations.
    logp_nomut = logp_nomut_incl[end3s] - logp_nomut_excl[end5s]
    # For each unique read, compute the joint probability that a random
    # read would have the same end coordinates, come from each cluster,
    # and have no mutations; normalize by the fraction of all reads that
    # have no two mutations too close, i.e. logp_noclose[end5s, end3s].
    # 2D (unique reads x clusters)
    logp_joint = (logp_clust_given_noclose[np.newaxis, :]
                  + logp_ends_given_noclose[end5s, end3s]
                  - logp_noclose[end5s, end3s]
                  + logp_nomut)
    # For each unique read, compute the likelihood of observing it
    # (including its mutations) by adjusting the above likelihood
    # of observing the end coordinates with no mutations.
    for j, mut_reads in zip(unmasked, muts_per_pos, strict=True):
        logp_joint[mut_reads] += logp_mut[j] - logp_not[j]
    # For each unique observed read, the marginal probability that a
    # random read would have the end coordinates and mutations, no
    # matter which cluster it came from, is the sum of the joint
    # probability over all clusters (axis 1).
    # 1D (unique reads)
    logp_marginal = np.logaddexp.reduce(logp_joint, axis=1)
    # Calculate the posterior probability that each read came from
    # each cluster by dividing the joint probability (observing the
    # read and coming from the cluster) by the marginal probability
    # (observing the read in any cluster).
    # 2D (unique reads x clusters)
    membership = np.exp(logp_joint - logp_marginal[:, np.newaxis])
    return logp_marginal, membership


def _calc_log_like(logp_marginal: np.ndarray, counts_per_uniq: np.ndarray):
    # Calculate the log likelihood of observing all the reads
    # by summing the log probability over all reads, weighted
    # by the number of times each read occurs. Cast to a float
    # explicitly to verify that the product is a scalar.
    return float(np.vdot(logp_marginal, counts_per_uniq))


class EMRun(object):
    """ Run expectation-maximization to cluster the given reads into the
    specified number of clusters. """

    def __init__(self,
                 uniq_reads: UniqReads,
                 k: int, *,
                 min_iter: int,
                 max_iter: int,
                 em_thresh: float):
        """
        Parameters
        ----------
        uniq_reads: UniqReads
            Container of unique reads
        k: int
            Number of clusters; must be a positive integer.
        min_iter: int
            Minimum number of iterations for clustering. Must be a
            positive integer no greater than `max_iter`.
        max_iter: int
            Maximum number of iterations for clustering. Must be a
            positive integer no less than `min_iter`.
        em_thresh: float
            Stop the algorithm when the difference in log likelihood
            between two successive iterations becomes smaller than the
            convergence threshold (and at least min_iter iterations have
            run). Must be a positive real number.
        """
        # Unique reads of mutations
        self.uniq_reads = uniq_reads
        # Number of clusters
        if not k >= 1:
            raise ValueError(f"k must be ≥ 1, but got {k}")
        self.k = k
        # Minimum number of iterations of EM
        if not min_iter >= 1:
            raise ValueError(f"min_iter must be ≥ 1, but got {min_iter}")
        self.min_iter = min_iter
        # Maximum number of iterations of EM
        if not max_iter >= min_iter:
            raise ValueError(f"max_iter must be ≥ min_iter ({min_iter}), "
                             f"but got {max_iter}")
        self.max_iter = max_iter
        # Cutoff for convergence of EM
        if not em_thresh >= 0.:
            raise ValueError(f"em_thresh must be ≥ 0, but got {em_thresh}")
        self.em_thresh = em_thresh
        # Mutation rates adjusted for observer bias.
        # 2D (all positions x clusters)
        self.p_mut = np.zeros((self.n_pos_total, self.k))
        # Read end coordinate proportions adjusted for observer bias.
        # 2D (all positions x all positions)
        self.p_ends = np.zeros((self.n_pos_total, self.n_pos_total))
        # Cluster proportions adjusted for observer bias.
        # 1D (clusters)
        self.p_clust = np.zeros(self.k)
        # Likelihood of each unique read coming from each cluster.
        # 2D (unique reads x clusters)
        self.resps = np.zeros((self.uniq_reads.num_uniq, self.k))
        # Marginal likelihoods of observing each unique read.
        # 1D (unique reads)
        self.logp_marginal = np.zeros(self.uniq_reads.num_uniq)
        # Trajectory of log likelihood values.
        self.log_likes: list[float] = list()
        # Number of iterations.
        self.iter = 0
        # Whether the algorithm has converged.
        self.converged = False

    @property
    def section_end5(self):
        """ 5' end of the section. """
        return self.uniq_reads.section.end5

    @cached_property
    def unmasked(self):
        """ Unmasked positions (0-indexed). """
        return self.uniq_reads.section.unmasked_zero

    @cached_property
    def masked(self):
        """ Masked positions (0-indexed). """
        return self.uniq_reads.section.masked_zero

    @property
    def n_pos_total(self):
        """ Number of positions, including those masked. """
        return self.uniq_reads.section.length

    @cached_property
    def n_pos_unmasked(self):
        """ Number of unmasked positions. """
        return get_length(self.unmasked, "unmasked positions")

    @cached_property
    def end5s(self):
        """ 5' end coordinates (0-indexed). """
        return self.uniq_reads.read_end5s_zero

    @cached_property
    def end3s(self):
        """ 3' end coordinates (0-indexed). """
        return self.uniq_reads.read_end3s_zero

    @cached_property
    def clusters(self):
        """ MultiIndex of k and cluster numbers. """
        return ClustHeader(ks=[self.k]).index

    @property
    def log_like(self):
        """ Return the current log likelihood, which is the last item in
        the trajectory of log likelihood values. """
        try:
            return self.log_likes[-1]
        except IndexError:
            # No log likelihood values have been computed.
            return np.nan

    @property
    def log_like_prev(self):
        """ Return the previous log likelihood, which is the penultimate
        item in the trajectory of log likelihood values. """
        try:
            return self.log_likes[-2]
        except IndexError:
            # Fewer than two log likelihood values have been computed.
            return np.nan

    @property
    def delta_log_like(self):
        """ Compute the change in log likelihood from the previous to
        the current iteration. """
        return self.log_like - self.log_like_prev

    @property
    def n_params(self):
        """ Number of parameters in the model. """
        # The parameters estimated by the model are
        # - the mutation rates for each position in each cluster
        # - the proportion of each cluster
        # - the proportion of every pair of read end coordinates
        # The cluster memberships are latent variables, not parameters.
        # The degrees of freedom of the mutation rates equals the total
        # number of unmasked positions times clusters.
        df_mut = self.n_pos_unmasked * self.k
        # The degrees of freedom of the coordinate proportions equals
        # the number of non-zero proportions (which are the only ones
        # that can be estimated) minus one because of the constraint
        # that the proportions of all end coordinates must sum to 1.
        df_ends = max(0, np.count_nonzero(
            self.p_ends[np.triu_indices_from(self.p_ends)]
        ) - 1)
        # The degrees of freedom of the cluster proportions equals the
        # number of clusters minus one because of the constraint that
        # the proportions of all clusters must sum to 1.
        df_clust = self.k - 1
        # The number of parameters is the sum of the degrees of freedom.
        return df_mut + df_ends + df_clust

    @property
    def n_data(self):
        """ Number of data points in the model. """
        # The number of data points is the total number of reads.
        return self.uniq_reads.num_nonuniq

    @property
    def bic(self):
        """ Bayesian Information Criterion of the model. """
        return _calc_bic(self.n_params, self.n_data, self.log_like)

    def _max_step(self):
        """ Run the Maximization step of the EM algorithm. """
        # Estimate the parameters based on observed data.
        (p_mut_observed,
         p_ends_observed,
         p_clust_observed) = calc_params_observed(
            self.n_pos_total,
            self.k,
            self.unmasked,
            self.uniq_reads.muts_per_pos,
            self.end5s,
            self.end3s,
            self.uniq_reads.counts_per_uniq,
            self.resps
        )
        if self.iter > 1:
            # If this iteration is not the first, then use the previous
            # values of the parameters as the initial guesses.
            guess_p_mut = self.p_mut
            guess_p_ends = self.p_ends
            guess_p_clust = self.p_clust
        else:
            # Otherwise, do not guess the initial parameters.
            guess_p_mut = None
            guess_p_ends = None
            guess_p_clust = None
        # Update the parameters.
        self.p_mut, self.p_ends, self.p_clust = calc_params(
            p_mut_observed,
            p_ends_observed,
            p_clust_observed,
            self.uniq_reads.min_mut_gap,
            guess_p_mut,
            guess_p_ends,
            guess_p_clust,
            prenormalize=False,
            quick_unbias=self.uniq_reads.quick_unbias,
            quick_unbias_thresh=self.uniq_reads.quick_unbias_thresh,
        )
        # Ensure all masked positions have a mutation rate of 0.
        if n_nonzero := np.count_nonzero(self.p_mut[self.masked]):
            p_mut_masked = self.p_mut[self.masked]
            logger.warning(
                f"{n_nonzero} masked position(s) have a mutation rate ≠ 0: "
                f"{p_mut_masked[p_mut_masked != 0.]}"
            )
            self.p_mut[self.masked] = 0.

    def _exp_step(self):
        """ Run the Expectation step of the EM algorithm. """
        # Update the marginal probabilities and cluster memberships.
        self.logp_marginal, self.resps = _expectation(
            self.p_mut,
            self.p_ends,
            self.p_clust,
            self.end5s,
            self.end3s,
            self.unmasked,
            self.uniq_reads.muts_per_pos,
            self.uniq_reads.min_mut_gap
        )
        # Calculate the log likelihood and append it to the trajectory.
        log_like = _calc_log_like(self.logp_marginal,
                                  self.uniq_reads.counts_per_uniq)
        self.log_likes.append(round(log_like, LOG_LIKE_PRECISION))

    def run(self, seed: int | None = None):
        """ Run the EM clustering algorithm.

        Parameters
        ----------
        seed: int | None = None
            Random number generator seed.

        Returns
        -------
        EMRun
            This instance, in order to permit statements such as
            ``return [em.run() for em in em_clusterings]``
        """
        logger.info(f"{self} began with {self.min_iter}-{self.max_iter} "
                    f"iterations")
        rng = np.random.default_rng(seed)
        # Erase the trajectory of log likelihood values (if any).
        self.log_likes.clear()
        # Choose the concentration parameters using a standard uniform
        # distribution so that the reads assigned to each cluster can
        # vary widely among the clusters (helping to make the clusters
        # distinct from the beginning) and among different runs of the
        # algorithm (helping to start the runs in very different states
        # so that they can explore much of the parameter space).
        # Use the half-open interval (0, 1] because the concentration
        # parameter of a Dirichlet distribution can be 1 but not 0.
        conc_params = 1. - rng.random(self.k)
        # Initialize cluster membership with a Dirichlet distribution.
        self.resps = rng.dirichlet(alpha=conc_params,
                                   size=self.uniq_reads.num_uniq)
        if self.uniq_reads.num_uniq == 0 or self.n_pos_unmasked == 0:
            logger.warning(f"{self} got 0 reads or positions: stopping")
            self.log_likes.append(0.)
            return self
        # Run EM until the log likelihood converges or the number of
        # iterations reaches max_iter, whichever happens first.
        self.converged = False
        for self.iter in range(1, self.max_iter + 1):
            # Update the parameters.
            self._max_step()
            # Update the cluster membership and log likelihood.
            self._exp_step()
            if not np.isfinite(self.log_like):
                raise ValueError(f"{self}, iteration {self.iter} returned a "
                                 f"non-finite log likelihood: {self.log_like}")
            logger.debug(f"{self}, iteration {self.iter}: "
                         f"log likelihood = {self.log_like}")
            # Check for convergence.
            if self.delta_log_like < 0.:
                # The log likelihood should not decrease.
                logger.warning(f"{self}, iteration {self.iter} returned a "
                               f"smaller log likelihood ({self.log_like}) than "
                               f"the previous iteration ({self.log_like_prev})")
            if (self.delta_log_like < self.em_thresh
                    and self.iter >= self.min_iter):
                # Converge if the increase in log likelihood is
                # smaller than the convergence cutoff and at least
                # the minimum number of iterations have been run.
                self.converged = True
                logger.info(f"{self} converged on iteration {self.iter}: "
                            f"last log likelihood = {self.log_like}")
                break
        else:
            # The log likelihood did not converge within the maximum
            # number of iterations.
            logger.warning(f"{self} failed to converge within {self.max_iter} "
                           f"iterations: last log likelihood = {self.log_like}")
        # Return this instance so that any code that runs multiple
        # EmClustering objects can create and run them in one line.
        return self

    @property
    def logn_exp(self):
        """ Log number of expected observations of each read. """
        return (np.log(self.uniq_reads.num_nonuniq) + self.logp_marginal
                if self.uniq_reads.num_nonuniq > 0
                else np.zeros(0))

    def get_props(self):
        """ Real and observed log proportion of each cluster. """
        return pd.DataFrame(self.p_clust[:, np.newaxis],
                            index=self.clusters,
                            columns=[CLUST_PROP_NAME])

    def get_mus(self):
        """ Log mutation rate at each position for each cluster. """
        return pd.DataFrame(self.p_mut[self.unmasked],
                            index=self.uniq_reads.section.unmasked_int,
                            columns=self.clusters)

    def get_resps(self, batch_num: int):
        """ Cluster memberships of the reads in the batch. """
        batch_uniq_nums = self.uniq_reads.batch_to_uniq[batch_num]
        return pd.DataFrame(self.resps[batch_uniq_nums],
                            index=batch_uniq_nums.index,
                            columns=self.clusters)

    def _calc_p_mut_pairs(self,
                          stat: Callable[[np.ndarray, np.ndarray], float],
                          summ: Callable[[list[float]], float]):
        """ Calculate a statistic for each pair of mutation rates and
        summarize them into one value. """
        stats = list(filterfalse(
            np.isnan,
            [stat(*p) for p in combinations(self.p_mut[self.unmasked].T, 2)]
        ))
        return summ(stats) if stats else np.nan

    def calc_max_pearson(self):
        """ Maximum Pearson correlation among any pair of clusters. """
        return self._calc_p_mut_pairs(calc_pearson, max)

    def calc_min_nrmsd(self):
        """ Minimum NRMSD among any pair of clusters. """
        return self._calc_p_mut_pairs(calc_nrmsd, min)

    def __str__(self):
        return f"{type(self).__name__} {self.uniq_reads}, {self.k} cluster(s)"

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
