"""Numba-jitted helpers for cluster.jackpot.

Kept separate so that importing ``jackpot`` (and thus ``cluster``) does
not import numba, which is slow.  ``jackpot`` references these via the
lazy ``_jit`` proxy (imported only when a simulation actually runs).
"""

import numpy as np
from numba import njit


@njit()
def assign_clusters_dropped(
    clusters: np.ndarray,
    p_clust_given_read: np.ndarray,
    n_reads_per_clust: np.ndarray,
    rng: np.random.Generator,
):
    """Assign one cluster to each read.

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
        p_clust_i_sum = 0.0
        for k in range(n_clust):
            p_clust_i_sum += p_clust_i[k]
            if rand_i < p_clust_i_sum:
                # Assign the read to the cluster.
                clusters[i] = k
                # Decrement the number of reads remaining in that cluster.
                n_reads_per_clust[k] -= 1
                break
    return clusters


@njit()
def sim_muts_dropped_jit(
    muts: np.ndarray,
    end5s: np.ndarray,
    end3s: np.ndarray,
    clusts: np.ndarray,
    p_mut_t: np.ndarray,
    min_mut_gap: int,
    max_attempts: int,
    rng: np.random.Generator,
):
    """Simulate the array of mutations. Drop reads with mutations too
    close.

    This function is meant to be called only after the arguments have
    been validated and thus makes the following assumptions:

    - muts is a C-contiguous NumPy boolean array with shape
      (reads x positions), initially all False.
    - end5s, end3s, and clusts are 1D NumPy integer arrays of the same
      length.
    - p_mut_t is a C-contiguous 2D NumPy float array with shape
      (clusters x positions) — the transpose of the usual p_mut so that
      each cluster's per-position rates are a contiguous 1D row.
    - p_mut_t has at least one row.
    - max_attempts ≥ 1
    """
    (n_reads,) = clusts.shape
    n_pos = muts.shape[1]
    # Reusable scratch buffer of positions placed during the current
    # attempt. Deferred commit (write to muts only on success) avoids
    # ever having to clear muts on a failed attempt.
    placed = np.empty(n_pos, dtype=np.int64)
    error = False
    for i in range(n_reads):
        k = clusts[i]
        end5 = end5s[i]
        end3 = end3s[i]
        p_mut_k = p_mut_t[k]
        read_valid = False
        n_placed = 0
        for _ in range(max_attempts):
            read_valid = True
            n_placed = 0
            # Sentinel chosen so the first candidate j cannot collide.
            last_mut = end5 - min_mut_gap - 1
            for j in range(end5, end3 + 1):
                if rng.random() < p_mut_k[j]:
                    if j - last_mut <= min_mut_gap:
                        read_valid = False
                        break
                    placed[n_placed] = j
                    n_placed += 1
                    last_mut = j
            if read_valid:
                break
        if not read_valid:
            error = True
            break
        muts_i = muts[i]
        for p in range(n_placed):
            muts_i[placed[p]] = True
    return error


@njit()
def sim_muts_merged_jit(
    muts: np.ndarray,
    end5s: np.ndarray,
    end3s: np.ndarray,
    clusts: np.ndarray,
    p_mut_t: np.ndarray,
    min_mut_gap: int,
    rng: np.random.Generator,
):
    """Simulate the array of mutations. Merge mutations too close.

    This function is meant to be called only after the arguments have
    been validated and thus makes the following assumptions:

    - muts is a C-contiguous NumPy boolean array with shape
      (reads x positions), initially all False.
    - end5s, end3s, and clusts are 1D NumPy integer arrays of the same
      length.
    - p_mut_t is a C-contiguous 2D NumPy float array with shape
      (clusters x positions) — the transpose of the usual p_mut so that
      each cluster's per-position rates are a contiguous 1D row.
    - p_mut_t has at least one row.
    - min_mut_gap ≥ 0
    """
    (n_reads,) = clusts.shape
    skip = min_mut_gap + 1
    for i in range(n_reads):
        k = clusts[i]
        end5 = end5s[i]
        end3 = end3s[i]
        p_mut_k = p_mut_t[k]
        muts_i = muts[i]
        j = end3
        while j >= end5:
            if rng.random() < p_mut_k[j]:
                muts_i[j] = True
                j -= skip
            else:
                j -= 1
