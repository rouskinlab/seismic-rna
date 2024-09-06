from itertools import combinations
from logging import getLogger

import numpy as np

from .array import get_length

logger = getLogger(__name__)
rng = np.random.default_rng()


def stochastic_bool(p: np.ndarray, validate: bool = True):
    """ Return a boolean array of the same shape as p where each value
    is True with probability equal to the corresponding value in p. """
    if validate and p.size > 0:
        if p.min() < 0.:
            raise ValueError(f"All p must be ≥ 0, but got {p[p < 0.]}")
        if p.max() > 1.:
            raise ValueError(f"All p must be ≤ 1, but got {p[p > 1.]}")
    return rng.random(p.shape) < p


def combinations_array(n: int, m: int):
    """ Like itertools.combinations(), but returns a NumPy array.

    Parameters
    ----------
    n: int
        Number of items; must be ≥ 0.
    m: int
        Number of items to choose; must be ≥ 0 and ≤ n.

    Returns
    -------
    np.ndarray
        Combinations of items.
    """
    if n < 0:
        raise ValueError(f"n must be ≥ 0, but got {n}")
    if not 0 <= m <= n:
        raise ValueError(f"m must be ≥ 0 and ≤ {n}, but got {m}")
    return np.array(list(combinations(range(n), m)), dtype=int)


def combinations_bools(n: int, m: int):
    """ All combinations of m items from a collection of n items.

    Parameters
    ----------
    n: int
        Number of items in total; must be ≥ 0.
    m: int
        Number of items to choose; must be ≥ 0 and ≤ n.

    Returns
    -------
    np.ndarray
        Combinations of items.
    """
    # Determine which indexes should be set to True.
    indexes = combinations_array(n, m)
    num_combinations, _ = indexes.shape
    # Initialize the matrix of combinations as boolean values.
    bools = np.zeros((num_combinations, n), dtype=bool)
    # Set the indexes to True.
    rows = np.repeat(np.arange(num_combinations), m)
    cols = indexes.reshape(-1, order="C")
    bools[rows, cols] = True
    return bools


def _calc_p_roundup_given_m(p_roundup: np.ndarray, m: int):
    """ Calculate the probability of each value being rounded up given
    that exactly m values are rounded up. """
    # Generate a boolean matrix of all possible combinations of exactly
    # m values to round up.
    n = get_length(p_roundup, "p_roundup")
    roundup_combos = combinations_bools(n, m)
    # Calculate the joint probability of each possible combination of m
    # values and that m values are rounded up.
    p_roundup_combo_and_m = np.where(roundup_combos,
                                     p_roundup,
                                     1. - p_roundup).prod(axis=1)
    # Calculate the probability that exactly m values are rounded up.
    p_m = p_roundup_combo_and_m.sum()
    # Calculate the conditional probability of each possible combination
    # of m values given that m values are rounded up.
    if p_m > 0.:
        p_roundup_combo_given_m = p_roundup_combo_and_m / p_m
    else:
        p_roundup_combo_given_m = np.zeros_like(p_roundup_combo_and_m)
    # Calculate the probability that each value is rounded up given that
    # exactly m values are rounded up.
    return p_roundup_combo_given_m @ roundup_combos


def _calc_p_roundup_given_ms(p_roundup: np.ndarray,
                             nums_roundup: np.ndarray,
                             p_nums_roundup: np.ndarray):
    """ Calculate the probability of each value being rounded up given
    that the number of values rounded up can be any of nums_roundup,
    each with probability p_nums_roundup. """
    # For each possible total number of values to round up, m, calculate
    # the conditional probability that each value would be rounded up
    # given that at total of exactly m values are rounded up.
    p_roundup_given_m = np.stack([_calc_p_roundup_given_m(p_roundup, m)
                                  for m in nums_roundup],
                                 axis=1)
    # Calculate the conditional probability that each value would be
    # rounded up given that the total number of values that are rounded
    # up is any of the possibilities, weighted by the probability of
    # each total number of values to round up.
    return p_roundup_given_m @ p_nums_roundup


def _calc_p_simulate_given_ms(p_roundup: np.ndarray):
    """ Calculate the round-up probabilities to use for simulation such
    that when simulations in which the wrong total number of values were
    rounded up are dropped, the probabilities of rounding up among the
    remaining simulations (with the correct number of values rounded up)
    will equal p_roundup. """
    # List the possible total numbers of values to round up and the
    # probability of rounding up each total number of values.
    if p_roundup.size == 0:
        return p_roundup
    num_roundup = p_roundup.sum()
    nums_roundup = np.array(sorted({int(np.floor(num_roundup)),
                                    int(np.ceil(num_roundup))}))
    p_nums_roundup = 1. - np.abs(nums_roundup - num_roundup)

    def objective(p_simulate_given_ms: np.ndarray):
        return (_calc_p_roundup_given_ms(p_simulate_given_ms,
                                         nums_roundup,
                                         p_nums_roundup)
                - p_roundup)

    from scipy.optimize import newton_krylov, NoConvergence
    try:
        return newton_krylov(objective, p_roundup)
    except NoConvergence as error:
        logger.warning("Failed to calculate the round-up probabilities for "
                       f"simulation from {p_roundup}: {error}")
        return p_roundup


def stochastic_round(values: np.ndarray | list | float | int,
                     preserve_sum: bool = False):
    """ Round values to integers stochastically, so that the probability
    of rounding up equals the fractional part of the original value.

    Parameters
    ----------
    values: np.ndarray | list | float | int
        Values to round; if scalar, a 0D integer array will be returned.
    preserve_sum: bool
        Whether the rounded values should have the same sum as the input
        values; if the sum is not an integer, it will be rounded using
        this function.

    Returns
    -------
    np.ndarray
        Values rounded to integers.
    """
    values = np.asarray_chkfinite(values)
    # Break each value into integer and fractional parts.
    rounded = np.asarray(np.floor(values), dtype=int)
    fractionals = (values - rounded).reshape(-1)
    if preserve_sum:
        # Choose the number of values that must be rounded up.
        n_round_up = int(stochastic_round(fractionals.sum()))
        # Calculate the probability that each value is rounded up.
        p = _calc_p_simulate_given_ms(fractionals)
    else:
        n_round_up = None
        p = fractionals
    # Randomly choose values to round up.
    round_up = stochastic_bool(p, validate=False)
    if n_round_up is not None:
        # Make sure the number being rounded up equals n_round_up.
        while np.count_nonzero(round_up) != n_round_up:
            round_up = stochastic_bool(p, validate=False)
    # Round up the chosen values.
    rounded.reshape(-1)[round_up] += 1
    return rounded
