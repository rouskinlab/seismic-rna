import numpy as np

from ..arrays import np_internal


@np_internal
def get_ranks(mus: np.ndarray):
    """ Compute the ranks of the mutation rates.

    Parameters
    ----------
    mus: numpy.ndarray
        Mutation rates. Multiple sets of mutation rates can be given as
        columns of a multidimensional array or DataFrame.

    Returns
    -------
    numpy.ndarray
        Ranks of the mutation rates.
    """
    if mus.ndim == 0:
        raise ValueError(f"mus must have at least 1 dimension, but got {mus}")
    if np.any(np.isnan(mus)):
        raise ValueError(f"Cannot rank mutation rates with NaN values: {mus}")
    # Initialize all ranks to -1.
    ranks = np.full_like(mus, -1, dtype=int)
    # Make coordinate arrays to index every element of ranks in order.
    coords = np.unravel_index(np.arange(ranks.size), ranks.shape)
    # Replace the coordinate array of the first axis with the indexes
    # in order of their ranks.
    coords[0][:] = np.argsort(mus, axis=0).reshape(coords[0].shape)
    # Fill the ranks.
    ranks[coords] = np.ravel(
        np.broadcast_to(np.expand_dims(np.arange(mus.shape[0]),
                                       axis=tuple(range(1, mus.ndim))),
                        ranks.shape)
    )
    # Confirm all ranks were filled.
    if np.any(ranks < 0):
        # This error should be impossible, but checking just in case.
        raise ValueError(f"Ranks cannot be negative, but got {ranks}")
    return ranks
