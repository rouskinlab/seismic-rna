import time
from typing import Iterable

import numpy as np
from numba import jit

from .marginal import calc_marginal
from ..core.arg import MUT_COLLISIONS_DROP, MUT_COLLISIONS_MERGE
from ..core.array import find_dims, get_length
from ..core.logs import logger
from ..core.random import get_random_integer_generator, stochastic_round
from ..core.unbias import (CLUSTERS,
                           POSITIONS,
                           READS,
                           UNIQUE_READS,
                           calc_p_noclose_given_clust,
                           calc_p_noclose_given_ends_auto,
                           calc_p_clust_given_noclose,
                           calc_p_ends_given_clust_noclose,
                           calc_p_clust_given_ends_noclose)
from ..core.validate import require_atleast

SUM_EXP_PRECISION = 3


def linearize_ends_matrix(p_ends: np.ndarray):
    """ Linearize an N x N matrix of end coordinate probabilities/counts
    into parallel arrays of 5'/3' coordinates and their values.

    Only entries in the upper triangle (end5 ≤ end3) with non-zero
    values are included.

    Parameters
    ----------
    p_ends: np.ndarray
        2D (positions x positions) array of end coordinate
        probabilities or counts.

    Returns
    -------
    end5s: np.ndarray
        1D array of 5' end coordinates of non-zero entries.
    end3s: np.ndarray
        1D array of 3' end coordinates of non-zero entries.
    values: np.ndarray
        1D array of values at (end5s, end3s) in p_ends.
    """
    end5s, end3s = np.nonzero(np.triu(p_ends))
    return end5s, end3s, p_ends[end5s, end3s]


@jit()
def _assign_clusters_dropped(clusters: np.ndarray,
                     p_clust_given_read: np.ndarray,
                     n_reads_per_clust: np.ndarray,
                     rng: np.random.Generator):
    """ Assign one cluster to each read.

    Parameters
    ----------
    clusters: np.ndarray
        1D (reads) integer array into which the cluster assigned to each
        read will be written.
    p_clust_given_read: np.ndarray
        2D (reads x clusters) float array of the probability that each
        read belongs to each cluster.
    n_reads_per_clust: np.ndarray
        1D (clusters) integer array of the number of reads to assign to
        each cluster.
    rng: np.random.Generator
        Random number generator used to permute read order and draw
        random values for cluster selection.

    Additionally, this function makes the following assumptions; if any
    assumption is violated, this function can produce incorrect results
    or attempt to access invalid memory addresses, potentially causing
    a segmentation fault.
    - There must be at least 1 cluster.
    - Every element of `p_clust_given_read` is ≥ 0 and ≤ 1.
    - Every row of `p_clust_given_read` sums to 1.
    - All elements of `n_reads_per_clust` are non-negative integers.
    - `n_reads_per_clust` sums to the number of reads (i.e. the lengths
      of `clusters` and `p_clust_given_read`).
    """
    n_reads, n_clust = p_clust_given_read.shape
    # Iterate through the reads in random order to eliminate bias that
    # can arise from the order in which reads are assigned.
    read_order = rng.permutation(n_reads)
    for i in read_order:
        rand_i = rng.random()
        # Calculate the probability of selecting each cluster for this
        # read.
        p_clust_i = p_clust_given_read[i] * n_reads_per_clust
        p_clust_i /= p_clust_i.sum()
        # Select the cluster based on the probabilities.
        p_clust_i_sum = 0.
        for k in range(n_clust):
            p_clust_i_sum += p_clust_i[k]
            if rand_i < p_clust_i_sum:
                # Assign the read to the cluster.
                clusters[i] = k
                # Decrement the number of reads remaining in that cluster.
                n_reads_per_clust[k] -= 1
                break
    return clusters


def _sim_clusters_dropped(p_clust_given_read: np.ndarray, seed: int | None):
    """ Simulate a cluster assignment for each read (drop mode).

    Reads are assigned so that the number assigned to each cluster
    matches a stochastic rounding of the total probability mass for
    that cluster, and the assignment of each individual read respects
    its per-cluster probabilities.

    Parameters
    ----------
    p_clust_given_read: np.ndarray
        2D (reads x clusters) array of the probability that each read
        belongs to each cluster.  Every row must be positive and sum
        to 1.
    seed: int | None
        Random number generator seed.

    Returns
    -------
    np.ndarray
        1D (reads) integer array of the cluster index assigned to each
        read.
    """
    n_reads, n_clust = p_clust_given_read.shape
    if p_clust_given_read.size > 0 and p_clust_given_read.min() <= 0.:
        raise ValueError("All p_clust_given_read must be > 0, but got "
                         f"{np.count_nonzero(p_clust_given_read <= 0.)} "
                         "probabilities < 0")
    if not np.allclose(p_clust_given_read.sum(axis=1), 1.):
        p_clust_given_read_sum = p_clust_given_read.sum(axis=1)
        not_close = np.logical_not(np.isclose(p_clust_given_read_sum, 1.))
        raise ValueError("All rows of p_clust_given_read must sum to 1, "
                         f"but got {p_clust_given_read_sum[not_close]}")
    # Initially assign all reads to cluster 0.
    clusters = np.zeros(n_reads, dtype=int)
    if n_clust == 1:
        # There is only one cluster, so all reads must be in it.
        return clusters
    # Choose the number of reads for each cluster, ensuring that the sum
    # equals the total number of reads.
    seeds = get_random_integer_generator(seed)
    p_clust_given_read_sum = p_clust_given_read.sum(axis=0)
    while np.sum(n_reads_per_clust := stochastic_round(
            p_clust_given_read_sum,
            preserve_sum=True,
            seed=next(seeds),
    )) != n_reads:
        pass
    # Choose the cluster for each read.
    rng = np.random.default_rng(seed)
    _assign_clusters_dropped(clusters,
                     p_clust_given_read,
                     n_reads_per_clust,
                     rng)
    return clusters


def _sim_clusters_merged(p_clust: np.ndarray, n_reads: int, seed: int | None):
    """ Simulate a cluster assignment for each read (merge mode).

    Each cluster receives a stochastically rounded number of reads
    proportional to `p_clust`, and reads are then permuted randomly
    so that each read is independently assigned.

    Parameters
    ----------
    p_clust: np.ndarray
        1D (clusters) array of the probability that a read belongs to
        each cluster.  Must be positive and sum to 1.
    n_reads: int
        Total number of reads to assign.
    seed: int | None
        Random number generator seed.

    Returns
    -------
    np.ndarray
        1D (reads) integer array of the cluster index assigned to each
        read.
    """
    n_clust, = p_clust.shape
    if n_clust <= 0:
        raise ValueError("At least one cluster is required")
    if p_clust.min() <= 0.:
        raise ValueError("All p_clust must be > 0, but got "
                         f"{np.count_nonzero(p_clust <= 0.)} "
                         "probabilities ≤ 0")
    if not np.isclose(p_clust.sum(), 1.):
        raise ValueError("All p_clust must sum to 1, but got "
                         f"{p_clust} (sum = {p_clust.sum()})")
    if n_clust == 1:
        # There is only one cluster, so all reads must be in it.
        return np.zeros(n_reads, dtype=int)
    # Choose the number of reads for each cluster, ensuring that the sum
    # equals the total number of reads.
    seeds = get_random_integer_generator(seed)
    while np.sum(n_reads_per_clust := stochastic_round(
            p_clust * n_reads,
            preserve_sum=True,
            seed=next(seeds),
    )) != n_reads:
        pass
    # Choose the cluster for each read.
    rng = np.random.default_rng(seed)
    return rng.permutation(np.repeat(np.arange(n_clust), n_reads_per_clust))


@jit()
def _sim_muts_dropped_jit(
    muts: np.ndarray,
    end5s: np.ndarray,
    end3s: np.ndarray,
    clusts: np.ndarray,
    p_mut: np.ndarray,
    min_mut_gap: int,
    max_attempts: int,
    rng: np.random.Generator
):
    """ Simulate the array of mutations. Drop reads with mutations too
    close.

    This function is meant to be called only after the arguments have
    been validated and thus makes the following assumptions:

    - muts is a NumPy boolean array with shape (reads x positions),
      initially all False
    - end5s, end3s, and clusts are 1D NumPy integer arrays of the same 
      length.
    - p_mut is a 2D NumPy array of floats.
    - p_mut has at least one column.
    - max_attempts ≥ 1
    """
    n_reads, = clusts.shape
    # Initialize a flag for if any read has an error.
    error = False
    # Assign mutations to each read.
    for i in range(n_reads):
        k = clusts[i]
        end5 = end5s[i]
        end3 = end3s[i]
        # Attempt to assign mutations to the read.
        read_valid = False
        remaining_attempts = max_attempts
        while not read_valid and remaining_attempts > 0:
            read_valid = True
            remaining_attempts -= 1
            # Assign mutations one position j at a time to avoid needing
            # a temporary array.
            for j in range(end5, end3 + 1):
                if rng.random() < p_mut[j, k]:
                    # Attempt to mutate this position. First, check for
                    # existing mutations too close before this position.
                    if (min_mut_gap > 0 
                        and np.any(muts[i, max(j - min_mut_gap, 0): j])):
                        # An existing mutation is too close.
                        read_valid = False
                        # Erase the existing mutations so that the next
                        # attempt will start fresh.
                        muts[i, :j] = False
                        # Break the for loop to abort this attempt.
                        break
                    else:
                        # Mutate this position.
                        muts[i, j] = True
        if not read_valid:
            # Flag that a read had an error and abort the whole simulation.
            error = True
            break
    return error


@jit()
def _sim_muts_merged_jit(
    muts: np.ndarray,
    end5s: np.ndarray,
    end3s: np.ndarray,
    clusts: np.ndarray,
    p_mut: np.ndarray,
    min_mut_gap: int,
    rng: np.random.Generator
):
    """ Simulate the array of mutations. Merge mutations too close.

    This function is meant to be called only after the arguments have
    been validated and thus makes the following assumptions:

    - muts is a NumPy boolean array with shape (reads x positions),
      initially all False
    - end5s, end3s, and clusts are 1D NumPy integer arrays of the same 
      length.
    - p_mut is a 2D NumPy array of floats.
    - p_mut has at least one column.
    - min_mut_gap ≥ 0
    """
    n_reads, = clusts.shape
    # Assign mutations to each read.
    for i in range(n_reads):
        k = clusts[i]
        end5 = end5s[i]
        end3 = end3s[i]
        # Assign mutations one position j at a time to avoid needing
        # a temporary array.
        j = end3
        while j >= end5:
            if rng.random() < p_mut[j, k]:
                # Mutate this position.
                muts[i, j] = True
                # Skip over min_mut_gap positions, none of which can
                # have a mutation.
                j -= (min_mut_gap + 1)
            else:
                # Leave the position without a mutation and go to the
                # next position.
                j -= 1


def _sim_muts_dropped(
    end5s: np.ndarray,
    end3s: np.ndarray,
    clusts: np.ndarray,
    p_mut: np.ndarray,
    min_mut_gap: int,
    seed: int | None,
    max_attempts: int = 256
):
    """ Simulate mutations and write them into reads (drop mode).

    Mutations are simulated position by position for each read.  If a
    read ends up with two mutations closer than `min_mut_gap`, the read
    is discarded and resimulated from scratch.

    Parameters
    ----------
    end5s: np.ndarray
        1D (reads) array of the 5' end coordinate (0-indexed) of each
        read.
    end3s: np.ndarray
        1D (reads) array of the 3' end coordinate (0-indexed) of each
        read.
    clusts: np.ndarray
        1D (reads) integer array of the cluster to which each read
        belongs.
    p_mut: np.ndarray
        2D (positions x clusters) array of the probability that a base
        at each position is mutated given the read covers the position.
    min_mut_gap: int
        Minimum number of non-mutated positions required between two
        mutations.
    seed: int | None
        Random number generator seed.
    max_attempts: int
        Maximum number of simulation attempts per read before raising
        a RuntimeError.

    Returns
    -------
    np.ndarray
        2D (reads x positions) boolean array of whether each base in
        each read is mutated.
    """
    require_atleast("min_mut_gap", min_mut_gap, 0, classes=int)
    # Validate the dimensions.
    dims = find_dims([(READS,),
                      (READS,),
                      (READS,),
                      (POSITIONS, CLUSTERS)],
                     [end5s,
                      end3s,
                      clusts,
                      p_mut],
                     ["end5s",
                      "end3s",
                      "clusts",
                      "p_mut"],
                     nonzero=[CLUSTERS])
    n_reads = dims[READS]
    n_pos = dims[POSITIONS]
    rng = np.random.default_rng(seed)
    # Initialize an array for the mutated positions; all positions start
    # out as not mutated (False).
    muts = np.zeros((n_reads, n_pos), dtype=bool)
    error = _sim_muts_dropped_jit(
        muts, end5s, end3s, clusts, p_mut, min_mut_gap, max_attempts, rng
    )
    if error:
        raise RuntimeError(
            f"Failed to simulate mutations in {max_attempts} attempts per read"
        )
    return muts


def _sim_muts_merged(
    end5s: np.ndarray,
    end3s: np.ndarray,
    clusts: np.ndarray,
    p_mut: np.ndarray,
    min_mut_gap: int,
    seed: int | None,
):
    """ Simulate mutations and write them into reads (merge mode).

    Mutations are simulated from the 3' end to the 5' end.  When a
    position is mutated, the `min_mut_gap` positions immediately 5' of
    it are skipped, enforcing the minimum gap constraint.

    Parameters
    ----------
    end5s: np.ndarray
        1D (reads) array of the 5' end coordinate (0-indexed) of each
        read.
    end3s: np.ndarray
        1D (reads) array of the 3' end coordinate (0-indexed) of each
        read.
    clusts: np.ndarray
        1D (reads) integer array of the cluster to which each read
        belongs.
    p_mut: np.ndarray
        2D (positions x clusters) array of the probability that a base
        at each position is mutated given the read covers the position.
    min_mut_gap: int
        Minimum number of non-mutated positions required between two
        mutations.
    seed: int | None
        Random number generator seed.

    Returns
    -------
    np.ndarray
        2D (reads x positions) boolean array of whether each base in
        each read is mutated.
    """
    require_atleast("min_mut_gap", min_mut_gap, 0, classes=int)
    # Validate the dimensions.
    dims = find_dims([(READS,),
                      (READS,),
                      (READS,),
                      (POSITIONS, CLUSTERS)],
                     [end5s,
                      end3s,
                      clusts,
                      p_mut],
                     ["end5s",
                      "end3s",
                      "clusts",
                      "p_mut"],
                     nonzero=[CLUSTERS])
    n_reads = dims[READS]
    n_pos = dims[POSITIONS]
    rng = np.random.default_rng(seed)
    # Initialize an array for the mutated positions; all positions start
    # out as not mutated (False).
    muts = np.zeros((n_reads, n_pos), dtype=bool)
    _sim_muts_merged_jit(
        muts, end5s, end3s, clusts, p_mut, min_mut_gap, rng
    )
    return muts


def _sim_reads_dropped(
    end5s: np.ndarray,
    end3s: np.ndarray,
    p_clust_given_ends_noclose: np.ndarray,
    p_mut: np.ndarray,
    min_mut_gap: int,
    seed: int | None
):
    """ Simulate reads with cluster assignments and mutations (drop mode).

    Parameters
    ----------
    end5s: np.ndarray
        1D (reads) array of the 5' end coordinate (0-indexed) of each
        read.
    end3s: np.ndarray
        1D (reads) array of the 3' end coordinate (0-indexed) of each
        read.
    p_clust_given_ends_noclose: np.ndarray
        3D (positions x positions x clusters) array of the probability
        that a read with given end coordinates and no mutations too
        close belongs to each cluster.
    p_mut: np.ndarray
        2D (positions x clusters) array of the mutation rate at each
        position in each cluster.
    min_mut_gap: int
        Minimum number of non-mutated positions between two mutations.
    seed: int | None
        Random number generator seed.

    Returns
    -------
    reads: np.ndarray
        2D (reads x (positions + 2)) integer array encoding mutations
        (boolean columns) and end coordinates (last two columns).
    clusts: np.ndarray
        1D (reads) integer array of the cluster assigned to each read.
    """
    find_dims([(READS,),
               (READS,),
               (POSITIONS, POSITIONS, CLUSTERS,),
               (POSITIONS, CLUSTERS)],
              [end5s,
               end3s,
               p_clust_given_ends_noclose,
               p_mut],
              ["end5s",
               "end3s",
               "p_clust_given_ends_noclose",
               "p_mut_given_span_noclose"],
              nonzero=[CLUSTERS])
    # Simulate the clusters and the mutations in each read.
    clusts = _sim_clusters_dropped(p_clust_given_ends_noclose[end5s, end3s], seed)
    muts = _sim_muts_dropped(end5s,
                     end3s,
                     clusts,
                     p_mut,
                     min_mut_gap,
                     seed)
    # Merge the mutation data and 5'/3' ends into one array of reads.
    reads = np.hstack([muts,
                       end5s[:, np.newaxis],
                       end3s[:, np.newaxis]],
                      dtype=int)
    return reads, clusts


def _sim_reads_merged(
    end5s: np.ndarray,
    end3s: np.ndarray,
    p_clust: np.ndarray,
    p_mut: np.ndarray,
    min_mut_gap: int,
    seed: int | None
):
    """ Simulate reads with cluster assignments and mutations (merge mode).

    Parameters
    ----------
    end5s: np.ndarray
        1D (reads) array of the 5' end coordinate (0-indexed) of each
        read.
    end3s: np.ndarray
        1D (reads) array of the 3' end coordinate (0-indexed) of each
        read.
    p_clust: np.ndarray
        1D (clusters) array of the probability that a read belongs to
        each cluster.
    p_mut: np.ndarray
        2D (positions x clusters) array of the mutation rate at each
        position in each cluster.
    min_mut_gap: int
        Minimum number of non-mutated positions between two mutations.
    seed: int | None
        Random number generator seed.

    Returns
    -------
    reads: np.ndarray
        2D (reads x (positions + 2)) integer array encoding mutations
        (boolean columns) and end coordinates (last two columns).
    clusts: np.ndarray
        1D (reads) integer array of the cluster assigned to each read.
    """
    dims = find_dims([(READS,),
               (READS,),
               (CLUSTERS,),
               (POSITIONS, CLUSTERS)],
              [end5s,
               end3s,
               p_clust,
               p_mut],
              ["end5s",
               "end3s",
               "p_clust_given_ends_noclose",
               "p_mut_given_span_noclose"],
              nonzero=[CLUSTERS])
    # Simulate the clusters and the mutations in each read.
    n_reads = dims[READS]
    clusts = _sim_clusters_merged(p_clust, n_reads, seed)
    muts = _sim_muts_merged(end5s,
                     end3s,
                     clusts,
                     p_mut,
                     min_mut_gap,
                     seed)
    # Merge the mutation data and 5'/3' ends into one array of reads.
    reads = np.hstack([muts,
                       end5s[:, np.newaxis],
                       end3s[:, np.newaxis]],
                      dtype=int)
    return reads, clusts


def _calc_obs_exp(reads: np.ndarray,
                  clusts: np.ndarray,
                  p_mut: np.ndarray,
                  p_ends: np.ndarray,
                  p_clust: np.ndarray,
                  min_mut_gap: int,
                  mut_collisions: str,
                  unmasked: np.ndarray):
    """ Calculate observed and expected counts for simulated reads.

    Parameters
    ----------
    reads: np.ndarray
        2D (reads x (positions + 2)) integer array of simulated reads,
        as returned by `_sim_reads_dropped` or `_sim_reads_merged`.
    clusts: np.ndarray
        1D (reads) integer array of the cluster assigned to each read.
    p_mut: np.ndarray
        2D (positions x clusters) array of the mutation rate at each
        position in each cluster.
    p_ends: np.ndarray
        2D (positions x positions) array of the probability distribution
        of 5'/3' end coordinates.
    p_clust: np.ndarray
        1D (clusters) array of the probability that a read belongs to
        each cluster.
    min_mut_gap: int
        Minimum number of non-mutated positions between two mutations.
    mut_collisions: str
        How to handle mutations that are less than `min_mut_gap`
        positions apart.
    unmasked: np.ndarray
        1D array of unmasked position indices (0-indexed).

    Returns
    -------
    num_obs: np.ndarray
        1D (unique reads) integer array of the observed count of each
        unique read.
    log_exp: np.ndarray
        1D (unique reads) float array of the log expected count of each
        unique read.
    """
    dims = find_dims([(READS, POSITIONS), (READS,)],
                     [reads, clusts],
                     ["reads", "clusts"])
    n_reads = dims[READS]
    # Count each unique read.
    uniq_reads, num_obs = np.unique(reads,
                                    axis=0,
                                    return_counts=True)
    uniq_end5s = uniq_reads[:, -2]
    uniq_end3s = uniq_reads[:, -1]
    # Calculate the properties of the simulated reads.
    muts_per_pos = [np.flatnonzero(uniq_reads[:, j]) for j in unmasked]
    log_exp = np.log(n_reads) + calc_marginal(p_mut,
                                              p_ends,
                                              p_clust,
                                              uniq_end5s,
                                              uniq_end3s,
                                              unmasked,
                                              muts_per_pos,
                                              min_mut_gap,
                                              mut_collisions)
    return num_obs, log_exp


def sim_obs_exp(end5s: np.ndarray,
                end3s: np.ndarray,
                p_mut: np.ndarray,
                p_ends: np.ndarray,
                p_clust: np.ndarray,
                min_mut_gap: int,
                mut_collisions: str,
                unmasked: np.ndarray,
                seed: int | None):
    """ Simulate observed and expected counts (infinite generator).

    Yields one (num_obs, log_exp) pair per simulated bootstrap replicate
    indefinitely.  The caller is responsible for stopping iteration.

    Parameters
    ----------
    end5s: np.ndarray
        1D (reads) array of the 5' end coordinate (0-indexed) of each
        read.
    end3s: np.ndarray
        1D (reads) array of the 3' end coordinate (0-indexed) of each
        read.
    p_mut: np.ndarray
        2D (positions x clusters) array of the mutation rate at each
        position in each cluster.
    p_ends: np.ndarray
        2D (positions x positions) array of the probability distribution
        of 5'/3' end coordinates.
    p_clust: np.ndarray
        1D (clusters) array of the probability that a read belongs to
        each cluster.
    min_mut_gap: int
        Minimum number of non-mutated positions between two mutations.
    mut_collisions: str
        How to handle mutations that are less than `min_mut_gap`
        positions apart.
    unmasked: np.ndarray
        1D array of unmasked position indices (0-indexed).
    seed: int | None
        Seed for the random number generator stream.

    Yields
    ------
    num_obs: np.ndarray
        1D (unique reads) integer array of observed counts.
    log_exp: np.ndarray
        1D (unique reads) float array of log expected counts.
    """
    # Validate the dimensions.
    find_dims([(READS,),
               (READS,),
               (POSITIONS, CLUSTERS),
               (POSITIONS, POSITIONS),
               (CLUSTERS,)],
              [end5s, end3s, p_mut, p_ends, p_clust],
              ["end5s", "end3s", "p_mut", "p_ends", "p_clust"])
    # Generate a unique random seed for each simulation.
    seeds = get_random_integer_generator(seed)
    if mut_collisions == MUT_COLLISIONS_DROP:
        # Calculate the parameters of reads with no mutations too close.
        p_noclose_given_ends = calc_p_noclose_given_ends_auto(
            p_mut, min_mut_gap
        )
        p_noclose_given_clust = calc_p_noclose_given_clust(
            p_ends, p_noclose_given_ends
        )
        p_ends_given_clust_noclose = calc_p_ends_given_clust_noclose(
            p_ends, p_noclose_given_ends
        )
        p_clust_given_noclose = calc_p_clust_given_noclose(
            p_clust, p_noclose_given_clust
        )
        p_clust_given_ends_noclose = calc_p_clust_given_ends_noclose(
            p_ends_given_clust_noclose, p_clust_given_noclose
        )
        for s in seeds:
            # Simulate the reads and the clusters to which they belong.
            reads, clusts = _sim_reads_dropped(
                end5s,
                end3s,
                p_clust_given_ends_noclose,
                p_mut,
                min_mut_gap,
                seed=s
            )
            yield _calc_obs_exp(reads,
                                clusts,
                                p_mut,
                                p_ends,
                                p_clust,
                                min_mut_gap,
                                mut_collisions,
                                unmasked)
    elif mut_collisions == MUT_COLLISIONS_MERGE:
        for s in seeds:
            # Simulate the reads and the clusters to which they belong.
            reads, clusts = _sim_reads_merged(
                end5s,
                end3s,
                p_clust,
                p_mut,
                min_mut_gap,
                seed=s
            )
            yield _calc_obs_exp(reads,
                                clusts,
                                p_mut,
                                p_ends,
                                p_clust,
                                min_mut_gap,
                                mut_collisions,
                                unmasked)
    else:
        raise ValueError(
            f"Invalid value for mut_collisions: {repr(mut_collisions)}"
        )


def calc_semi_g_anomaly(num_obs: int | np.ndarray,
                        log_exp: float | np.ndarray):
    """ Calculate each read's semi-G-anomaly.

    The semi-G-anomaly is half of each read's contribution to the
    G-test statistic: O * log(O / E) = O * (log(O) - log(E)).

    Parameters
    ----------
    num_obs: int | np.ndarray
        Observed count(s) of each unique read.
    log_exp: float | np.ndarray
        Log expected count(s) of each unique read.

    Returns
    -------
    int | np.ndarray
        Semi-G-anomaly: num_obs * (log(num_obs) - log_exp).
    """
    with np.errstate(divide="ignore"):
        result = num_obs * (np.log(num_obs) - log_exp)
    # By convention, 0 * log(0/E) = 0 (the limit as O→0⁺).
    if not np.isscalar(num_obs):
        result = np.where(num_obs == 0, 0., result)
    elif num_obs == 0:
        result = 0.
    return result


def calc_jackpot_score(semi_g_anomalies: np.ndarray, n_reads: int):
    """ Calculate the jackpotting score.

    The jackpotting score is defined as the average log of the ratio of
    observed to expected counts for a read, weighted by the count of
    each read:

    JS = sum{O_i * log(O_i / E_i)} / sum{O_i}
    where i is each unique read.

    This formula is based on that of the G-test statistic, which is
    identical except that the latter is multiplied by 2 and not divided
    by the number of observations.

    Parameters
    ----------
    semi_g_anomalies: np.ndarray
        1D array of semi-G-anomalies, one per unique read.
    n_reads: int
        Total number of reads (sum of all observed counts).

    Returns
    -------
    float
        Jackpotting score.
    """
    if n_reads == 0:
        if semi_g_anomalies.size > 0:
            raise ValueError("If n_reads is 0, then semi_g_anomalies must have "
                             f"length 0, but got {semi_g_anomalies.size}")
        return 0.
    return semi_g_anomalies.sum() / n_reads


def calc_jackpot_score_ci(jackpot_scores: Iterable[float],
                          confidence_level: float):
    """ Calculate a confidence interval for the mean of jackpotting scores.

    Uses a Student's t-distribution assuming the scores are normally
    distributed.  Returns (nan, nan) if fewer than two scores are given.

    Parameters
    ----------
    jackpot_scores: Iterable[float]
        Collection of jackpotting scores (e.g. from bootstrap replicates).
    confidence_level: float
        Confidence level for the interval, e.g. 0.95.

    Returns
    -------
    ci_lo: float
        Lower bound of the confidence interval.
    ci_up: float
        Upper bound of the confidence interval.
    """
    # Ensure g_stats is a NumPy array of floats with no NaN/inf values.
    if not isinstance(jackpot_scores, (np.ndarray, list)):
        jackpot_scores = list(jackpot_scores)
    jackpot_scores = np.asarray_chkfinite(jackpot_scores, dtype=float)
    # Calculate the confidence interval of the mean, assuming that the
    # jackpotting scores are distributed normally.
    n = get_length(jackpot_scores, "jackpot_scores")
    if n > 1:
        # Import SciPy here instead of at the top of the module because
        # the latter would take so long as to slow down the start-up.
        from scipy.stats import t
        t_lo, t_up = t.interval(confidence_level, n - 1)
        mean = jackpot_scores.mean()
        std_err = np.sqrt(jackpot_scores.var() / n)
        ci_lo = mean + std_err * t_lo
        ci_up = mean + std_err * t_up
    else:
        ci_lo = np.nan
        ci_up = np.nan
    return ci_lo, ci_up


def calc_jackpot_quotient(real_jackpot_score: float,
                          null_jackpot_score: float):
    """ Calculate the jackpotting quotient.

    The jackpotting quotient indicates how much more overrepresented the
    average read is in the real dataset compared to the null dataset.

    Since the jackpotting score is the expected log-ratio of a read's
    observed to expected count, raising e to the power of it yields the
    observed-to-expected ratio, which measures jackpotting.

    JS = average{log(O / E)}
    exp(JS) = exp(average{log(O / E)})

    Then, the jackpotting quotient is the jackpotting score of the real
    dataset divided by that of the null dataset:

    JQ = exp(JS_real) / exp(JS_null)
       = exp(JS_real - JS_null)

    Parameters
    ----------
    real_jackpot_score: float
        Jackpotting score of the real dataset.
    null_jackpot_score: float
        Jackpotting score of the null model.

    Returns
    -------
    float
        Jackpotting quotient
    """
    return np.exp(real_jackpot_score - null_jackpot_score)


def bootstrap_jackpot_scores(uniq_end5s: np.ndarray,
                             uniq_end3s: np.ndarray,
                             counts_per_uniq: np.ndarray,
                             p_mut: np.ndarray,
                             p_ends: np.ndarray,
                             p_clust: np.ndarray,
                             min_mut_gap: int,
                             mut_collisions: str,
                             unmasked: np.ndarray,
                             real_jackpot_score: float,
                             confidence_level: float,
                             max_jackpot_quotient: float,
                             max_jackpot_sims: int,
                             seed: int | None):
    """ Bootstrap jackpotting scores from the null model.

    Parameters
    ----------
    uniq_end5s: np.ndarray
        1D (unique reads) array of the 5' ends (0-indexed) of the unique
        reads.
    uniq_end3s: np.ndarray
        1D (unique reads) array of the 3' ends (0-indexed) of the unique
        reads.
    counts_per_uniq: np.ndarray
        1D (unique reads) array of the number of times each unique read
        was observed.
    p_mut: np.ndarray
        2D (positions x clusters) array of the mutation rate at each
        position in each cluster.
    p_ends: np.ndarray
        2D (positions x positions) array of the probability distribution
        of 5'/3' end coordinates.
    p_clust: np.ndarray
        1D (clusters) probability that a read comes from each cluster.
    min_mut_gap: int
        Minimum number of non-mutated positions between two mutations.
    mut_collisions: str
        How to handle mutations that are less than min_mut_gap positions
        apart.
    unmasked: np.ndarray
        1D array of unmasked position indices (0-indexed).
    real_jackpot_score: float
        Jackpotting score of the real dataset.
    confidence_level: float
        Confidence level for computing a confidence interval of the
        jackpotting quotient.
    max_jackpot_quotient: float
        Stop bootstrapping once the confidence interval lies entirely
        above or entirely below this threshold.
    max_jackpot_sims: int
        Maximum number of simulations to calculate a null jackpotting
        score.
    seed: int | None
        Seed for the random number generator stream.

    Returns
    -------
    np.ndarray
        1D array of jackpotting scores from the null bootstrap replicates.
    """
    find_dims([(UNIQUE_READS,),
               (UNIQUE_READS,),
               (UNIQUE_READS,),
               (POSITIONS, CLUSTERS),
               (POSITIONS, POSITIONS),
               (CLUSTERS,)],
              [uniq_end5s,
               uniq_end3s,
               counts_per_uniq,
               p_mut,
               p_ends,
               p_clust],
              ["uniq_end5s",
               "uniq_end3s",
               "counts_per_uniq",
               "p_mut",
               "p_ends",
               "p_clust"])
    n_reads = counts_per_uniq.sum()
    logger.routine(f"Began boostrapping null jackpotting scores for a dataset "
                   f"with {n_reads} reads and a real jackpotting score of "
                   f"{real_jackpot_score}")
    # Simulate observed and expected read counts.
    end5s = np.repeat(uniq_end5s, counts_per_uniq)
    end3s = np.repeat(uniq_end3s, counts_per_uniq)
    null_jackpotting_scores = list()
    obs_exps = sim_obs_exp(end5s,
                            end3s,
                            p_mut,
                            p_ends,
                            p_clust,
                            min_mut_gap,
                            mut_collisions,
                            unmasked,
                            seed=seed)
    jq_ci_lo = np.nan
    jq_ci_up = np.nan
    conf_pct = (f"{round(confidence_level * 100., 1)} % confidence interval "
                "for the mean jackpotting quotient")
    for _ in range(max_jackpot_sims):
        num_obs, log_exp = next(obs_exps)
        # Calculate this null model's jackpotting score.
        null_g_anomalies = calc_semi_g_anomaly(num_obs, log_exp)
        null_jackpotting_score = calc_jackpot_score(null_g_anomalies,
                                                    n_reads)
        null_jackpotting_scores.append(null_jackpotting_score)
        logger.detail(f"Null jackpotting score: {null_jackpotting_score}")
        # Calculate a confidence interval for the mean jackpotting score
        # of the null models simulated so far.
        js_ci_lo, js_ci_up = calc_jackpot_score_ci(null_jackpotting_scores,
                                                   confidence_level)
        # Determine a confidence interval for the jackpotting quotient.
        # Since the jackpotting quotient depends inversely on the null
        # jackpotting score, the lower and upper bounds are swapped.
        jq_ci_lo = calc_jackpot_quotient(real_jackpot_score, js_ci_up)
        jq_ci_up = calc_jackpot_quotient(real_jackpot_score, js_ci_lo)
        if not np.isnan(jq_ci_lo) and not np.isnan(jq_ci_up):
            logger.detail(f"{conf_pct}: {jq_ci_lo} - {jq_ci_up}")
        # Stop when the confidence interval lies entirely below or above
        # max_jackpot_quotient, so it's clear whether the jackpotting
        # quotient is less or greater than max_jackpot_quotient.
        # Avoid "if not jq_ci_lo <= max_jackpot_quotient <= jq_ci_up"
        # because this expression will evaluate to True after the first
        # iteration, when jq_ci_lo and jq_ci_up will both be NaN.
        if jq_ci_lo > max_jackpot_quotient or max_jackpot_quotient > jq_ci_up:
            break
    else:
        # The confidence interval still contains max_jackpot_quotient
        # after max_jackpot_sims simulations.
        logger.warning(
            f"After the maximum of {max_jackpot_sims} simulations, the "
            f"{conf_pct} is {jq_ci_lo} - {jq_ci_up}, which still contains "
            f"the maximum jackpotting quotient {max_jackpot_quotient}, "
            "making the data ambiguously jackpotted"
        )
    logger.routine("Ended boostrapping null jackpotting scores")
    return np.array(null_jackpotting_scores)
