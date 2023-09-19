import unittest as ut

import numpy as np
import pandas as pd

from ..blank import blank_relvec
from ...core.rel import NOCOV
from ...core.sect import seq_pos_to_index
from ...core.seq import DNA


class TestBlankRelvec(ut.TestCase):
    """ Test function `blank_relvec`. """

    def compare_numpy(self, result, expect: np.ndarray):
        self.assertIsInstance(result, np.ndarray)
        self.assertIs(type(result), type(expect))
        self.assertIs(result.dtype, expect.dtype)
        self.assertTrue(np.array_equal(result, expect))

    def compare_pandas(self, result, expect: pd.Series | pd.DataFrame):
        self.assertIsInstance(result, (pd.Series, pd.DataFrame))
        self.assertIs(type(result), type(expect))
        self.assertTrue(expect.equals(result))

    def test_numpy_1d(self):
        """ Test returning a 1D NumPy array. """
        for length in range(10):
            self.compare_numpy(blank_relvec(length),
                               np.full(shape=(length,),
                                       fill_value=NOCOV,
                                       dtype=np.uint8))

    def test_numpy_2d_int(self):
        """ Test returning a 2D NumPy array with integer reads. """
        for length in range(10):
            for n_reads in range(10):
                self.compare_numpy(blank_relvec(length, n_reads),
                                   np.full(shape=(n_reads, length),
                                           fill_value=NOCOV,
                                           dtype=np.uint8))

    def test_numpy_2d_list(self):
        """ Test returning a 2D NumPy array with a list of reads. """
        for length in range(10):
            for n_reads in range(10):
                for list_func in [list, np.array, pd.Index]:
                    reads = list_func(list(map(str, range(n_reads))))
                    self.compare_numpy(blank_relvec(length, reads),
                                       np.full(shape=(n_reads, length),
                                               fill_value=NOCOV,
                                               dtype=np.uint8))

    def test_pandas_series(self):
        """ Test returning a Pandas Series. """
        for length in range(1, 10):
            seq = DNA.random(length)
            index = seq_pos_to_index(seq, list(range(1, length + 1)), 1)
            self.compare_pandas(blank_relvec(seq),
                                pd.Series(NOCOV, index=index, dtype=np.uint8))

    def test_pandas_dataframe_int(self):
        """ Test returning a Pandas DataFrame with integer reads. """
        for length in range(1, 10):
            seq = DNA.random(length)
            index = seq_pos_to_index(seq, list(range(1, length + 1)), 1)
            for n_reads in range(10):
                reads = [f"Read_{i + 1}" for i in range(n_reads)]
                self.compare_pandas(blank_relvec(seq, n_reads),
                                    pd.DataFrame(NOCOV,
                                                 index=reads,
                                                 columns=index,
                                                 dtype=np.uint8))

    def test_pandas_dataframe_list(self):
        """ Test returning a Pandas DataFrame with integer reads. """
        for length in range(1, 10):
            seq = DNA.random(length)
            index = seq_pos_to_index(seq, list(range(1, length + 1)), 1)
            for n_reads in range(10):
                for list_func in [list, np.array, pd.Index]:
                    reads = list_func(list(map(str, range(n_reads))))
                    self.compare_pandas(blank_relvec(seq, reads),
                                        pd.DataFrame(NOCOV,
                                                     index=reads,
                                                     columns=index,
                                                     dtype=np.uint8))
