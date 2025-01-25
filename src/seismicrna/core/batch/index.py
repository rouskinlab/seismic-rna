import numpy as np
import pandas as pd

from ..seq import BASE_NAME

# Indexes of read and batch numbers.
READ_NUM = "Read Number"
BATCH_NUM = "Batch Number"
RB_INDEX_NAMES = [BATCH_NUM, READ_NUM]


def list_batch_nums(num_batches: int):
    """ List the batch numbers. """
    return list(range(num_batches))


def count_base_types(base_pos_index: pd.Index):
    """ Return the number of each type of base in the index of positions
    and bases. """
    base_types, counts = np.unique(base_pos_index.get_level_values(BASE_NAME),
                                   return_counts=True)
    return pd.Series(counts, base_types)


def iter_base_types(base_pos_index: pd.Index):
    """ For each type of base in the index of positions and bases, yield
    the positions in the index with that type of base. """
    bases, inverse = np.unique(base_pos_index.get_level_values(BASE_NAME),
                               return_inverse=True)
    for i, base in enumerate(bases):
        yield base, base_pos_index[inverse == i]
