from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .names import CLUST_PROP_NAME
from .uniq import UniqReads
from ..core.batch import get_length
from ..core.header import index_order_clusts
from ..core.mu import (calc_p_noclose_given_ends_numpy,
                       calc_p_ends_given_noclose,
                       calc_spanning_sum,
                       calc_params_numpy,
                       triu_log)

logger = getLogger(__name__)

LOG_LIKE_PRECISION = 3  # number of digits to round the log likelihood
FIRST_ITER = 1  # number of the first iteration


def calc_bic(n_params: int,
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
    if n_data < min_data_param_ratio * n_params:
        logger.warning(f"The Bayesian Information Criterion (BIC) uses an "
                       f"approximation that is valid only when the size of the "
                       f"sample (n = {n_data}) is much larger than the number "
                       f"of parameters being estimated (p = {n_params}). "
                       f"This model does not meet this criterion, so the BIC "
                       f"may not indicate the model's complexity accurately.")
    return n_params * np.log(n_data) - 2. * log_like


class EmClustering(object):
    """ Run expectation-maximization to cluster the given reads into the
    specified number of clusters. """

    def __init__(self,
                 uniq_reads: UniqReads,
                 order: int, *,
                 min_iter: int,
                 max_iter: int,
                 conv_thresh: float):
        """
        Parameters
        ----------
        uniq_reads: UniqReads
            Container of unique reads
        order: int
            Number of clusters into which to cluster the reads;
            must be a positive integer
        min_iter: int
            Minimum number of iterations for clustering. Must be a
            positive integer no greater than `max_iter`.
        max_iter: int
            Maximum number of iterations for clustering. Must be a
            positive integer no less than `min_iter`.
        conv_thresh: float
            Stop the algorithm when the difference in log likelihood
            between two successive iterations becomes smaller than the
            convergence threshold (and at least min_iter iterations have
            run). Must be a positive real number.
        """
        # Unique reads of mutations
        self.uniq_reads = uniq_reads
        # Number of clusters (i.e. order)
        if not order >= 1:
            raise ValueError(f"order must be ≥ 1, but got {order}")
        self.order = order
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
        if not conv_thresh >= 0.:
            raise ValueError(f"conv_thresh must be ≥ 0, but got {conv_thresh}")
        self.conv_thresh = conv_thresh
        # Number of reads with no two mutations too close.
        # self.n_reads_noclose = np.empty(self.order)
        # Probability that a read with each pair of end coordinates in
        # each cluster would have no two mutations too close.
        # 3D (all positions x all positions x clusters)
        self.logp_noclose = np.empty((self.n_pos_total,
                                      self.n_pos_total,
                                      self.order))
        # Mutation rates adjusted for observer bias.
        # 2D (all positions x clusters)
        self.p_mut = np.empty((self.n_pos_total, self.order))
        # Read end coordinate proportions adjusted for observer bias.
        # 2D (all positions x all positions)
        self.p_ends = np.empty((self.n_pos_total, self.n_pos_total))
        # Cluster proportions adjusted for observer bias.
        # 1D (clusters)
        self.p_clust = np.empty(self.order)
        # Likelihood of each unique read coming from each cluster.
        # 2D (unique reads x clusters)
        self.membership = np.empty((self.uniq_reads.num_uniq, self.order))
        # Marginal likelihoods of observing each unique read.
        # 1D (unique reads)
        self.log_marginals = np.empty(self.uniq_reads.num_uniq)
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
        return self.uniq_reads.end5s_zero

    @cached_property
    def end3s(self):
        """ 3' end coordinates (0-indexed). """
        return self.uniq_reads.end3s_zero

    @cached_property
    def clusters(self):
        """ MultiIndex of the order and cluster numbers. """
        return index_order_clusts(self.order)

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
    def bic(self):
        """ Bayesian Information Criterion of the model. """
        # The parameters estimated by the model are
        # - the mutation rates for each position in each cluster
        # - the proportion of each cluster
        # - the proportion of every pair of read end coordinates
        # The cluster memberships are latent variables, not parameters.
        # The degrees of freedom of the mutation rates equals the total
        # number of unmasked positions times clusters.
        df_mut = self.n_pos_unmasked * self.order
        # The degrees of freedom of the coordinate proportions equals
        # the number of non-zero proportions (which are the only ones
        # that can be estimated) minus one because of the constraint
        # that the proportions of all end coordinates must sum to 1.
        df_ends = np.count_nonzero(
            self.p_ends[np.triu_indices_from(self.p_ends)]
        ) - 1
        # The degrees of freedom of the cluster proportions equals the
        # number of clusters minus one because of the constraint that
        # the proportions of all clusters must sum to 1.
        df_clust = self.order - 1
        # The number of parameters is the sum of the degrees of freedom.
        n_params = df_mut + df_ends + df_clust
        # The number of data points is the total number of reads.
        n_data = self.uniq_reads.num_nonuniq
        return calc_bic(n_params, n_data, self.log_like)

    def _max_step(self):
        """ Run the Maximization step of the EM algorithm. """
        # Count each unique read in each cluster.
        # 2D (unique reads x clusters)
        n_each_read_each_clust = (self.uniq_reads.counts_per_uniq[:, np.newaxis]
                                  * self.membership)
        # Count the total number of reads in each cluster.
        # 1D (clusters)
        n_reads_per_clust = np.sum(n_each_read_each_clust, axis=0)
        # Calculate the observed proportion of reads in each cluster.
        # 1D (clusters)
        p_clust_given_noclose = n_reads_per_clust / n_reads_per_clust.sum()
        # Calculate the proportion of each unique read in each cluster.
        # 2D (unique reads x clusters)
        p_each_read_each_clust = n_each_read_each_clust / n_reads_per_clust
        # Calculate the proportion of observed reads with each pair of
        # 5' and 3' end coordinates.
        # 3D (all positions x all positions x clusters)
        p_ends_given_noclose = calc_p_ends_given_noclose(self.n_pos_total,
                                                         self.end5s,
                                                         self.end3s,
                                                         p_each_read_each_clust,
                                                         check_values=False)
        # Count the observed reads that cover each position.
        # 2D (all positions x clusters)
        n_reads_per_pos = calc_spanning_sum(p_ends_given_noclose) * n_reads_per_clust
        # Count the observed mutations at each position.
        # 2D (all positions x clusters)
        n_muts_per_pos = np.zeros((self.n_pos_total, self.order))
        for j, mut_reads in zip(self.unmasked,
                                self.uniq_reads.muts_per_pos,
                                strict=True):
            # Calculate the number of mutations at each position in each
            # cluster by summing the count-weighted likelihood that each
            # read with a mutation at (j) came from the cluster.
            n_muts_per_pos[j] = (
                    self.uniq_reads.counts_per_uniq[mut_reads]
                    @ self.membership[mut_reads]
            )
        # Calculate the observed mutation rate at each position.
        # 2D (all positions x clusters)
        p_mut_noclose = n_muts_per_pos / n_reads_per_pos
        if self.iter > FIRST_ITER:
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
        self.p_mut, self.p_ends, self.p_clust = calc_params_numpy(
            p_mut_noclose,
            p_ends_given_noclose,
            p_clust_given_noclose,
            self.uniq_reads.min_mut_gap,
            guess_p_mut,
            guess_p_ends,
            guess_p_clust,
            prenormalize=False,
        )
        # Update the log probability that a read with each pair of 5'/3'
        # end coordinates would have no two mutations too close.
        self.logp_noclose = triu_log(calc_p_noclose_given_ends_numpy(
            self.p_mut, self.uniq_reads.min_mut_gap
        ))

    def _exp_step(self):
        """ Run the Expectation step of the EM algorithm. """
        # Compute the logs of the parameters.
        with np.errstate(divide="ignore"):
            # Suppress warnings about taking the log of zero, which is a
            # valid mutation rate.
            logp_mut = np.log(self.p_mut)
        logp_not = np.log(1. - self.p_mut)
        logp_ends = triu_log(self.p_ends)
        logp_clust = np.log(self.p_clust)
        # For each cluster, calculate the probability that a read up to
        # and including each position would have no mutations.
        logp_nomut_incl = np.cumsum(logp_not, axis=0)
        # For each cluster, calculate the probability that a read up to
        # but not including each position would have no mutations.
        logp_nomut_excl = logp_nomut_incl - logp_nomut_incl[0]
        # For each unique read, compute the joint probability that a
        # random read would come from each cluster and have identical
        # end coordinates but no mutations.
        # To save memory, overwrite self.membership, which has the same
        # dimensions and is no longer needed.
        # 2D (unique reads x clusters)
        self.membership = (logp_clust[np.newaxis, :]
                           + (logp_ends[self.end5s, self.end3s, np.newaxis]
                              - self.logp_noclose[self.end5s, self.end3s])
                           + (logp_nomut_incl[self.end3s]
                              - logp_nomut_excl[self.end5s]))
        # For each unique read, compute the likelihood of observing it
        # (including its mutations) by adjusting the above likelihood
        # of observing the end coordinates with no mutations.
        for j, mut_reads in zip(self.unmasked,
                                self.uniq_reads.muts_per_pos,
                                strict=True):
            self.membership[mut_reads] += logp_mut[j] - logp_not[j]
        # For each unique observed read, the marginal probability that a
        # random read would have the end coordinates and mutations, no
        # matter which cluster it came from, is the sum of the joint
        # probability over all clusters (axis 1).
        # 1D (unique reads)
        self.log_marginals = np.logaddexp.reduce(self.membership, axis=1)
        # Calculate the posterior probability that each read came from
        # each cluster by dividing the joint probability (observing the
        # read and coming from the cluster) by the marginal probability
        # (observing the read in any cluster).
        # 2D (unique reads x clusters)
        self.membership = np.exp(self.membership
                                 - self.log_marginals[:, np.newaxis])
        # Calculate the log likelihood of observing all the reads
        # by summing the log probability over all reads, weighted
        # by the number of times each read occurs. Cast to a float
        # explicitly to verify that the product is a scalar.
        log_like = float(np.vdot(self.log_marginals,
                                 self.uniq_reads.counts_per_uniq))
        self.log_likes.append(round(log_like, LOG_LIKE_PRECISION))

    def run(self, props_seed: int | None = None, resps_seed: int | None = None):
        """ Run the EM clustering algorithm.

        Parameters
        ----------
        props_seed: int | None = None
            Random number generator seed for cluster proportions
        resps_seed: int | None = None
            Random number generator seed for read responsibilities

        Returns
        -------
        EmClustering
            This instance, in order to permit statements such as
            ``return [em.run() for em in em_clusterings]``
        """
        # Import scipy here instead of at the top of this module because
        # its import is slow enough to impact global startup time.
        from scipy.stats import dirichlet
        logger.info(f"{self} began with {self.min_iter} - {self.max_iter} "
                    f"iterations")
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
        conc_params = 1. - np.random.default_rng(props_seed).random(self.order)
        # Initialize cluster membership with a Dirichlet distribution.
        self.membership = dirichlet.rvs(alpha=conc_params,
                                        size=self.uniq_reads.num_uniq,
                                        random_state=resps_seed)
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
            elif (self.delta_log_like < self.conv_thresh
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
        return np.log(self.uniq_reads.num_nonuniq) + self.log_marginals

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

    def get_members(self, batch_num: int):
        """ Cluster memberships of the reads in the batch. """
        batch_uniq_nums = self.uniq_reads.batch_to_uniq[batch_num]
        return pd.DataFrame(self.membership[batch_uniq_nums],
                            index=batch_uniq_nums.index,
                            columns=self.clusters)

    def __str__(self):
        return f"{type(self).__name__} {self.uniq_reads} to order {self.order}"

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
