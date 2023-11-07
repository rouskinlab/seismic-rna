from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .uniq import UniqReads
from ..core.header import index_order_clusts
from ..core.mu import calc_f_obs_numpy, calc_mu_adj_numpy, calc_prop_adj_numpy
from ..core.rand import rng

logger = getLogger(__name__)

LOG_LIKE_PRECISION = 3  # number of digits to round the log likelihood
FIRST_ITER = 1  # number of the first iteration

# Names
OBSERVED = "Observed"
ADJUSTED = "Adjusted"
BITVECTOR = "Bit Vector"
LOGOBS = "Log Observed"
LOGEXP = "Log Expected"


def calc_bic(n_params: int,
             n_data: int,
             log_like: float,
             min_data_param_ratio: float = 10.):
    """
    Compute the Bayesian Information Criterion (BIC) of a model.
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
        The Bayesian Information Criterion (BIC)
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
    """ Run expectation-maximization to cluster the given bit vectors
    into the specified number of clusters."""

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
            Number of clusters into which to cluster the bit vectors;
            must be a positive integer
        min_iter: int
            Minimum number of iterations for clustering. Must be a
            positive integer no greater than max_iter.
        max_iter: int
            Maximum number of iterations for clustering. Must be a
            positive integer no less than min_iter.
        conv_thresh: float
            Stop the algorithm when the difference in log likelihood
            between two successive iterations becomes smaller than the
            convergence threshold (and at least min_iter iterations have
            run). Must be a positive real number (ideally close to 0).
        """
        # Unique bit vectors of mutations
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
        # Number of reads observed in each cluster, not adjusted
        self.nreads = np.empty(self.order, dtype=float)
        # Log of the fraction observed for each cluster
        self.log_f_obs = np.empty(self.order, dtype=float)
        # Mutation rates of used positions (row) in each cluster (col),
        # adjusted for observer bias
        self.mus = np.empty((self.uniq_reads.num_pos, self.order), dtype=float)
        # Positions of the section that will be used for clustering
        # (0-indexed from the beginning of the section)
        self.sparse_pos = self.uniq_reads.section.unmasked_zero
        # Likelihood of each vector (col) coming from each cluster (row)
        self.resps = np.empty((self.order, self.uniq_reads.num_uniq),
                              dtype=float)
        # Marginal probabilities of observing each bit vector.
        self.log_marginals = np.empty(self.uniq_reads.num_uniq, dtype=float)
        # Trajectory of log likelihood values.
        self.log_likes: list[float] = list()
        # Number of iterations.
        self.iter = 0
        # Whether the algorithm has converged.
        self.converged = False

    @cached_property
    def nreads_total(self):
        """ Total number of reads, including redundant ones. """
        return int(round(self.uniq_reads.counts_per_uniq.sum()))

    @property
    def sparse_mus(self):
        """ Mutation rates of all positions (row) in each cluster (col),
        including positions that are not used for clustering. The rate
        for every unused position always remains zero. """
        # Initialize an array of zeros with one row for each position in
        # the section and one column for each cluster.
        sparse_mus = np.zeros((self.uniq_reads.section.length, self.order),
                              dtype=float)
        # Copy the rows of self.mus that correspond to used positions.
        # The rows for unused positions remain zero (hence "sparse").
        sparse_mus[self.sparse_pos] = self.mus
        return sparse_mus

    @cached_property
    def clusters(self):
        """ MultiIndex of the order and cluster numbers. """
        return index_order_clusts(self.order)

    @property
    def prop_obs(self):
        """ Observed proportion of each cluster, without adjusting for
        observer bias. """
        return self.nreads / np.sum(self.nreads)

    @property
    def prop_adj(self):
        """ Proportion of each cluster, adjusted for observer bias. """
        return calc_prop_adj_numpy(self.prop_obs, np.exp(self.log_f_obs))

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
        """ Compute this model's Bayesian Information Criterion. """
        # The parameters estimated by the model are the mutation rates
        # for each position in each cluster and the proportion of each
        # cluster in the population. By contrast, the cluster membership
        # values are latent variables because each describes one item in
        # the sample (a bit vector), not a parameter of the population.
        # The number of data points is the number of unique bit vectors.
        return calc_bic(self.mus.size + self.prop_adj.size,
                        self.uniq_reads.num_uniq,
                        self.log_like)

    def _max_step(self):
        """ Run the Maximization step of the EM algorithm. """
        # Calculate the number of reads in each cluster by summing the
        # count-weighted likelihood that each bit vector came from the
        # cluster.
        self.nreads = self.resps @ self.uniq_reads.counts_per_uniq
        # If this iteration is not the first, then use mutation rates
        # from the previous iteration as the initial guess for this one.
        mus_guess = self.sparse_mus if self.iter > FIRST_ITER else None
        # Compute the observed mutation rate at each position (j).
        # To save memory and enable self.sparse_mus to use the observed
        # rates of mutation, overwrite self.mus rather than allocate a
        # new array.
        for j, muts_j in enumerate(self.uniq_reads.muts_per_pos):
            # Calculate the number of mutations at each position in each
            # cluster by summing the count-weighted likelihood that each
            # bit vector with a mutation at (j) came from the cluster,
            # then divide by the count-weighted sum of the number of
            # reads in the cluster to find the observed mutation rate.
            self.mus[j] = ((self.resps[:, muts_j]
                            @ self.uniq_reads.counts_per_uniq[muts_j])
                           / self.nreads)
        # Solve for the real mutation rates that are expected to yield
        # the observed mutation rates after considering read drop-out.
        self.mus = calc_mu_adj_numpy(self.sparse_mus,
                                     self.uniq_reads.min_mut_gap,
                                     mus_guess)[self.sparse_pos]

    def _exp_step(self):
        """ Run the Expectation step of the EM algorithm. """
        # Import scipy here instead of at the top of this module because
        # its import is slow enough to impact global startup time.
        from scipy.special import logsumexp
        # Update the log fraction observed of each cluster.
        self.log_f_obs = np.log(calc_f_obs_numpy(self.sparse_mus,
                                                 self.uniq_reads.min_mut_gap))
        # Compute the logs of the mutation and non-mutation rates.
        with np.errstate(divide="ignore"):
            # Suppress warnings about taking the log of zero, which is a
            # valid mutation rate. Thus, here we CAN compute log of 0,
            # just like Chuck Norris.
            log_mus = np.log(self.mus)
        log_nos = np.log(1. - self.mus)
        # Compute for each cluster and observed bit vector the joint
        # probability that a random read vector would both come from the
        # same cluster and have exactly the same set of bits. To save
        # memory, overwrite self.resps rather than allocate a new array.
        for k in range(self.order):
            # Compute the probability that a bit vector has no mutations
            # given that it comes from cluster (k), which is the sum of
            # all not-mutated log probabilities, sum(log_nos[:, k]),
            # minus the cluster's log fraction observed, log_f_obs[k].
            log_prob_no_muts_given_k = np.sum(log_nos[:, k]) - self.log_f_obs[k]
            # Initialize the log probability for all bit vectors in
            # cluster (k) to the log probability that an observed bit
            # vector comes from cluster (k) (log(prob_obs[k])) and has
            # no mutations given that it came from cluster (k).
            log_prob_no_muts_and_k = (np.log(self.prop_obs[k])
                                      + log_prob_no_muts_given_k)
            self.resps[k].fill(log_prob_no_muts_and_k)
            # Loop over each position (j); each iteration adds the log
            # PMF for one additional position in each bit vector to the
            # accumulating total log probability of each bit vector.
            # Using accumulation with this loop also uses less memory
            # than would holding the PMF for every position in an array
            # and then summing over the position axis.
            for j in range(self.uniq_reads.num_pos):
                # Compute the probability that each bit vector would
                # have the bit observed at position (j). The probability
                # is modeled using a Bernoulli distribution, with PMF:
                # log_mus[k, j] if muts[j, i] else log_nos[k, j].
                # This PMF could be computed explicitly, such as with
                # scipy.stats.bernoulli.logpmf(muts[j], mus[k, j]).
                # But few of the bits are mutated, so a more efficient
                # (and much, MUCH faster) way is to assume initially
                # that no bit is mutated -- that is, to initialize the
                # probability of every bit vector with the probability
                # that the bit vector has no mutations, as was done.
                # Then, for only the few bits that are mutated, the
                # probability is adjusted by adding the log of the
                # probability that the base is mutated minus the log of
                # the probability that the base is not mutated.
                log_adjust_jk = log_mus[j, k] - log_nos[j, k]
                self.resps[k][self.uniq_reads.muts_per_pos[j]] += log_adjust_jk
        # For each observed bit vector, the marginal probability that a
        # random read would have the same series of bits (regardless of
        # which cluster it came from) is the sum over all clusters of
        # the joint probability of coming from the cluster and having
        # the specific series of bits.
        self.log_marginals = logsumexp(self.resps, axis=0)
        # Calculate the posterior probability that each bit vector came
        # from each cluster by dividing the joint probability (observing
        # the bit vector and coming from the cluster) by the marginal
        # probability (observing the bit vector in any cluster).
        self.resps = np.exp(self.resps - self.log_marginals)
        # Calculate the log likelihood of observing all the bit vectors
        # by summing the log probability over all bit vectors, weighted
        # by the number of times each bit vector occurs. Cast to a float
        # explicitly to verify that the product is a scalar.
        log_like = float(np.vdot(self.log_marginals,
                                 self.uniq_reads.counts_per_uniq))
        self.log_likes.append(round(log_like, LOG_LIKE_PRECISION))

    def run(self):
        """ Run the EM clustering algorithm. """
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
        conc_params = 1. - rng.random(self.order)
        # Initialize cluster membership with a Dirichlet distribution.
        self.resps = dirichlet.rvs(conc_params, self.uniq_reads.num_uniq).T
        # Run EM until the log likelihood converges or the number of
        # iterations reaches max_iter, whichever happens first.
        self.converged = False
        for self.iter in range(1, self.max_iter + 1):
            # Update the mutation rates and cluster proportions.
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
        """ Log number of expected observations of each bit vector. """
        return np.log(self.nreads_total) + self.log_marginals

    def get_props(self):
        """ Real and observed log proportion of each cluster. """
        return pd.DataFrame(np.vstack([self.prop_obs, self.prop_adj]).T,
                            index=self.clusters,
                            columns=[OBSERVED, ADJUSTED])

    def get_mus(self):
        """ Log mutation rate at each position for each cluster. """
        return pd.DataFrame(self.mus,
                            index=self.uniq_reads.section.unmasked,
                            columns=self.clusters)

    def get_resps(self, batch_num: int):
        """ Responsibilities of the reads in the batch. """
        batch_uniq_nums = self.uniq_reads.batch_to_uniq[batch_num]
        return pd.DataFrame(self.resps.T[batch_uniq_nums],
                            index=batch_uniq_nums.index,
                            columns=self.clusters)

    '''
    
    @property
    def logp_contig(self):
        """ Log probability of every bit vector from the most likely to
        the rarest observed bit vector, including any bit vector with an
        intermediate probability that had an observed count of 0. """
        # Find the log probability of the rarest observed bit vector.
        min_logp = np.min(self.log_marginals)
        # Initialize a bit vector iterator for each cluster.
        bvecs = {clust: iter_all_bit_vectors(cmus,
                                             self.data.section,
                                             self.data.min_mut_gap)
                 for clust, cmus in self.get_mus().items()}
        # For each cluster, iterate through the bit vectors up to and
        # including the rarest observed bit vector; record the log
        # probability of each bit vector encountered.
        logps = dict()
        cump = 0.
        for clust, clust_bvecs in bvecs.items():
            # Initialize the map of bit vectors to log probabilities
            # for this cluster.
            logps[clust] = dict()
            while True:
                # Find the next bit vector; record its log probability.
                bvec, logp = next(clust_bvecs)
                logps[clust][bvec.tobytes().decode()] = logp
                cump += np.exp(logp)
                # Iterate until the log probability falls below that of
                # the rarest bit vector.
                if logp <= min_logp:
                    break
        # Find all "non-rare" bit vectors that had a higher probability
        # than the rarest bit vector in any cluster.
        bvecs_non_rare = set.union(*map(set, logps.values()))
        # For each cluster, compute the log probability of all non-rare
        # bit vectors that were not already encountered in the cluster.
        for clust, clust_bvecs in bvecs.items():
            # Find the non-rare bit vectors that were not encountered
            # already in the cluster.
            unseen_non_rare = bvecs_non_rare - set(logps[clust])
            # Iterate until all such bit vectors are encountered.
            while unseen_non_rare:
                # Find the next bit vector; record its log probability.
                bvec, logp = next(clust_bvecs)
                bvec_str = bvec.tobytes().decode()
                logps[clust][bvec_str] = logp
                if bvec_str in unseen_non_rare:
                    # This bit vector is no longer unseen.
                    unseen_non_rare.remove(bvec_str)
        # Now that all necessary bit vectors have been encountered, join
        # the maps into a single DataFrame (taking the intersection of
        # the indexes) of bit vectors (indexes) and clusters (columns).
        logps_df = pd.concat([pd.Series(clust_logps, name=clust)
                              for clust, clust_logps in logps.items()],
                             join="inner")
        # Compute the joint probability that a random bit vector would
        # match the bit vector corresponding to a row in logps_df and
        # come from the cluster corresponding to a column in logps_df.
        logps_joint = logps_df + np.log(self.get_props()[self.ADJUSTED])
        # Compute the marginal probability of observing each bit vector,
        # regardless of the cluster from which it came.
        logps_marginal = logsumexp(logps_joint, axis=0)
        # Sort the bit vectors by decreasing log probability.
        return logps_marginal.sort_values(ascending=False)

    def output_cdf(self):
        """ Return the cumulative fraction of reads in descending order
        of likelihood. """

        # Map each bit vector to its log probability and observed count.
        ns_obs: dict[str, int] = dict()
        logps: dict[int, dict[str, float]] = {c: dict() for c in bvecs}
        # Iterate through each unique bit vector that was observed.
        for bvec_arr, n_obs in zip(self.data.get_matrix(),
                                   self.data.counts_per_uniq,
                                   strict=True):
            # Convert the bit vector from an array to a str so that it
            # becomes hashable and internable (unlike bytes).
            bvec_str = bvec_arr.tobytes().decode()
            # Map the bit vector to its observed count.
            ns_obs[bvec_str] = n_obs
            # For each cluster, check if the bit vector has appeared.
            for clust, clogps in logps.items():
                # If the bit vector has not yet appeared, then iterate
                # through the bit vectors until finding it.
                while bvec_str not in clogps:
                    # Get the next bit vector and its log likelihood.
                    bvec, logp = next(bvecs[clust])
                    # Map the bit vector to its log probability.
                    clogps[bvec.tobytes().decode()] = logp
        # Find the bit vectors that are
        # Compute the adjusted proportion of each cluster.
        props = self.get_props()[self.ADJUSTED]
    
    '''

    def __str__(self):
        return f"Clustering {self.uniq_reads} to order {self.order}"

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
