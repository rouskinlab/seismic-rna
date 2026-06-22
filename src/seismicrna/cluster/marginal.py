from __future__ import annotations

from ..core.arg.cli import MUT_COLLISIONS_DROP, MUT_COLLISIONS_MERGE
from ..core.array import find_dims
from ..core.unbias import (
    UNIQUE_READS,
    POSITIONS,
    CLUSTERS,
    calc_p_noclose_given_clust,
    calc_p_noclose_given_ends_auto,
    calc_p_ends_given_clust_noclose,
    calc_p_clust_given_noclose,
)

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np


def _zero_masked(p_mut: np.ndarray, unmasked: np.ndarray):
    """Set mutation rates of masked positions to zero."""
    import numpy as np

    p_mut_unmasked = np.zeros_like(p_mut)
    p_mut_unmasked[unmasked] = p_mut[unmasked]
    return p_mut_unmasked


def _calc_logp_joint(
    p_mut: np.ndarray,
    p_ends: np.ndarray,
    p_clust: np.ndarray,
    end5s: np.ndarray,
    end3s: np.ndarray,
    unmasked: np.ndarray,
    muts_per_pos: list[np.ndarray],
    min_mut_gap: int,
    mut_collisions: str,
):
    """Compute the log joint probability of each unique read and cluster.

    For each unique read i and cluster k, the joint probability is:
    P(read i, cluster k) = P(cluster k) * P(end coords) * P(mutations | ends, cluster k)

    adjusted for the mut_collisions model (drop or merge).

    Parameters
    ----------
    p_mut: np.ndarray
        2D (positions x clusters) array of per-position mutation rates.
    p_ends: np.ndarray
        2D (positions x positions) array of the probability distribution
        of 5'/3' end coordinates.
    p_clust: np.ndarray
        1D (clusters) array of cluster proportions.
    end5s: np.ndarray
        1D (unique reads) array of 5' end coordinates (0-indexed).
    end3s: np.ndarray
        1D (unique reads) array of 3' end coordinates (0-indexed).
    unmasked: np.ndarray
        1D array of unmasked position indices (0-indexed).
    muts_per_pos: list[np.ndarray]
        For each unmasked position, a 1D array of unique-read indices
        that are mutated at that position.
    min_mut_gap: int
        Minimum number of non-mutated positions between two mutations.
    mut_collisions: str
        How to handle mutations closer than `min_mut_gap` positions
        apart; either MUT_COLLISIONS_DROP or MUT_COLLISIONS_MERGE.

    Returns
    -------
    np.ndarray
        2D (unique reads x clusters) array of log joint probabilities.
    """
    import numpy as np

    # Validate the dimensions.
    find_dims(
        [
            (POSITIONS, CLUSTERS),
            (POSITIONS, POSITIONS),
            (CLUSTERS,),
            (UNIQUE_READS,),
            (UNIQUE_READS,),
        ],
        [p_mut, p_ends, p_clust, end5s, end3s],
        ["p_mut", "p_ends", "p_clust", "end5s", "end3s"],
        nonzero=True,
    )
    # Ensure the mutation rates of masked positions are 0.
    p_mut = _zero_masked(p_mut, unmasked)
    if mut_collisions == MUT_COLLISIONS_DROP:
        # Correct for reads with mutations too close dropping out.
        p_noclose_given_ends = calc_p_noclose_given_ends_auto(p_mut, min_mut_gap)
        # Adjust the cluster probabilities.
        p_clust = calc_p_clust_given_noclose(
            p_clust, calc_p_noclose_given_clust(p_ends, p_noclose_given_ends)
        )
        # Adjust the end coordinate probabilities.
        p_ends = calc_p_ends_given_clust_noclose(p_ends, p_noclose_given_ends)
    elif mut_collisions == MUT_COLLISIONS_MERGE:
        p_noclose_given_ends = None
        # Make p_ends 3D (positions x positions x clusters) to match the
        # output of the mut_collisions == MUT_COLLISIONS_DROP branch.
        p_ends = p_ends[:, :, np.newaxis]
    else:
        raise ValueError(f"Invalid value for mut_collisions: {mut_collisions}")
    # Calculate the log mutation rates and anti-mutation rates.
    with np.errstate(divide="ignore"):
        # Suppress warnings about taking the log of zero, which is a
        # valid mutation rate.
        logp_mut = np.log(p_mut)
    logp_not = np.log(1.0 - p_mut)
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
    # and have no mutations.
    # 2D (unique reads x clusters)
    logp_joint = (
        np.log(p_clust)[np.newaxis, :]
        + np.log(p_ends[end5s, end3s])
        + logp_nomut_given_clust
    )
    if p_noclose_given_ends is not None:
        # Normalize logp_joint by the fraction of all reads that have
        # no two mutations too close.
        logp_joint -= np.log(p_noclose_given_ends[end5s, end3s])
    # For each unique read, compute the likelihood of observing it
    # (including its mutations) by adjusting the above likelihood
    # of observing the end coordinates with no mutations.
    for j, mut_reads in zip(unmasked, muts_per_pos, strict=True):
        # Adjust the probability of observing reads with a mutation at
        # position j by multiplying by the probability that position j
        # is mutated and dividing by the probability that position j
        # is not mutated (which is the base probability that was already
        # accounted for by adding logp_nomut_given_clust to logp_joint).
        adjustment = (logp_mut[j] - logp_not[j])[np.newaxis, :]
        if mut_collisions == MUT_COLLISIONS_MERGE:
            # After merging mutations less than min_mut_gap positions
            # apart, the probability of observing a read equals that of
            # observing any read that would have resulted in the same
            # set of mutations post-merge. This means that before the
            # merge, all of the min_mut_gap positions 5' of each post-
            # merge mutation could have been mutated or not, since any
            # mutation in that region would have been removed during
            # the merge. Therefore, the min_mut_gap positions 5' of the
            # post-merge mutation do not affect the probability of the
            # read after the merge. Thus, the probability that those
            # positions are mutated must be removed from logp_joint.
            # First, determine which positions in each mutated read are
            # within min_mut_gap positions of position j on the 5' side.
            nomut_window_end5s = np.maximum(end5s[mut_reads], j - min_mut_gap)
            # Calculate the log probability that none of those positions
            # are mutated: sum logp_not from nomut_window_end5 to (j-1).
            logp_nomut_window = logp_nomut_excl[j] - logp_nomut_excl[nomut_window_end5s]
            # Adjust logp_joint by dividing by that probability.
            adjustment = adjustment - logp_nomut_window
        logp_joint[mut_reads] += adjustment
    return logp_joint


def _calc_logp_marginal(logp_joint: np.ndarray):
    """Compute the log marginal probability of each unique read.

    Marginalizes over clusters by summing the joint probabilities
    (in probability space) via log-sum-exp.

    Parameters
    ----------
    logp_joint: np.ndarray
        2D (unique reads x clusters) array of log joint probabilities,
        as returned by `_calc_logp_joint`.

    Returns
    -------
    np.ndarray
        1D (unique reads) array of log marginal probabilities.
    """
    import numpy as np

    # Sum the joint probability over all clusters (axis 1).
    return np.logaddexp.reduce(logp_joint, axis=1)


def calc_marginal(*args, **kwargs):
    """Compute the log marginal probability of each unique read.

    A convenience wrapper that chains `_calc_logp_joint` and
    `_calc_logp_marginal`.  Accepts the same arguments as
    `_calc_logp_joint`.

    Returns
    -------
    np.ndarray
        1D (unique reads) array of log marginal probabilities.
    """
    return _calc_logp_marginal(_calc_logp_joint(*args, **kwargs))


def calc_marginal_resps(*args, **kwargs):
    """Compute log marginal probabilities and cluster responsibilities.

    A convenience wrapper that chains `_calc_logp_joint`,
    `_calc_logp_marginal`, and the responsibility calculation.
    Accepts the same arguments as `_calc_logp_joint`.

    Returns
    -------
    logp_marginal: np.ndarray
        1D (unique reads) array of log marginal probabilities.
    resps: np.ndarray
        2D (unique reads x clusters) array of posterior probabilities
        that each read belongs to each cluster.
    """
    import numpy as np

    logp_joint = _calc_logp_joint(*args, **kwargs)
    logp_marginal = _calc_logp_marginal(logp_joint)
    # Calculate the posterior probability that each read came from
    # each cluster by dividing the joint probability (observing the
    # read and coming from the cluster) by the marginal probability
    # (observing the read in any cluster).
    # 2D (unique reads x clusters)
    resps = np.exp(logp_joint - logp_marginal[:, np.newaxis])
    return logp_marginal, resps
