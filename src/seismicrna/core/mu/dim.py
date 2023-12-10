import numpy as np
import pandas as pd


def get_common_num_pos(*mus: np.ndarray | pd.Series | pd.DataFrame):
    """ Determine the number of positions for the mutation rates.

    Parameters
    ----------
    *mus: numpy.ndarray | pandas.Series | pandas.DataFrame
        Groups of mutation rates; each can contain multiple sets as the
        columns of a multidimensional array.

    Returns
    -------
    int
        Number of positions for the mutation rates.
    """
    # Determine the number of positions in each array.
    try:
        nums_pos = list(set(mu.shape[0] for mu in mus))
    except IndexError:
        raise ValueError("Cannot count positions in 0-D arrays")
    if len(nums_pos) == 0:
        raise ValueError("No arrays were given to count positions")
    if len(nums_pos) > 1:
        raise ValueError(f"Got multiple numbers of positions: {nums_pos}")
    return nums_pos[0]
