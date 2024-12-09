import numpy as np

from .array import calc_inverse

rng = np.random.default_rng()


def _stochastic_round(values: np.ndarray | list | float | int):
    """ Round values to integers stochastically, so that the probability
    of rounding up equals the fractional part of the original value.

    Parameters
    ----------
    values: np.ndarray | list | float | int
        Values to round; if scalar, a 0D integer array will be returned.

    Returns
    -------
    np.ndarray
        Values rounded to integers.
    """
    values = np.asarray_chkfinite(values)
    # Break each value into integer and fractional parts.
    rounded = np.asarray(np.floor(values), dtype=int)
    # Reshape the fractional parts into a 1D array.
    fractionals = np.reshape(values - rounded, -1)
    # Choose which values to round up.
    round_up = rng.random(fractionals.shape) < fractionals
    # Round up the chosen values.
    rounded.reshape(-1)[round_up] += 1
    return rounded


def _stochastic_round_sum(values: np.ndarray | list | float | int):
    """ Like stochastic_round, but guarantees that the sums before and
    after rounding are equal. If the former is not an integer, then the
    sum after rounding will be either the nearest integer greater than
    the sum (with probability equal to the fractional part of the sum)
    or the nearest integer less than the sum (with the complementary
    probability).

    Parameters
    ----------
    values: np.ndarray | list | float | int
        Values to round; if scalar, a 0D integer array will be returned.

    Returns
    -------
    np.ndarray
        Values rounded to integers, with the original sum preserved.
    """
    values = np.asarray_chkfinite(values)
    if values.size == 0:
        return np.zeros(values.shape, dtype=bool)
    # Shuffle the values so that the outcome of each will be independent
    # of the others.
    order = rng.permutation(values.size)
    shuffled = values.reshape(-1, order="C")[order]
    # Choose a random offset between 0 and 1 to remove bias against the
    # shuffled item that comes first.
    offset = rng.random()
    # Lay out the values along the real number line, such that they are
    # in shuffled order and abutting exactly (without gaps or overlaps),
    # each occupies a width of the number line equal to its value, and
    # the first item starts at the offset; calculate the start and end
    # coordinates of each item.
    ends = offset + np.cumsum(shuffled)
    starts = np.concatenate([[offset], ends[:-1]])
    # When the values are laid out in this way, the number of integers
    # that each value straddles must equal the value rounded either up
    # or down, and the probability that it equals the value rounded up
    # is the fractional part of the value; thus, count the number of
    # integers that each value straddles in order to round it up/down.
    rounded = np.asarray(np.floor(ends) - np.floor(starts), dtype=int)
    # Un-shuffle the values and reshape to the original shape.
    return rounded[calc_inverse(order)].reshape(values.shape, order="C")


def stochastic_round(values: np.ndarray | list | float | int,
                     preserve_sum: bool = False):
    """ Round values to integers stochastically, so that the probability
    of rounding up equals the fractional part of the original value.

    Parameters
    ----------
    values: np.ndarray | list | float | int
        Values to round; if scalar, a 0D integer array will be returned.
    preserve_sum: bool
        Whether to ensure that the sum of the rounded values equals the
        sum of the original values.

    Returns
    -------
    np.ndarray
        Values rounded to integers, with the original sum preserved.
    """
    if preserve_sum:
        return _stochastic_round_sum(values)
    return _stochastic_round(values)

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
