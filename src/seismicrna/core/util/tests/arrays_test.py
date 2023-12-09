import unittest as ut
from itertools import permutations

import numpy as np
import pandas as pd

from seismicrna.core.util.arrays import (get_priority,
                                         get_shared_index,
                                         get_shared_indexes,
                                         np_internal,
                                         np2_internal,
                                         promote_arrays,
                                         rearray)


class TestGetPriority(ut.TestCase):

    def test_array1d(self):
        a1 = np.array([0.0, 0.1, 0.2])
        self.assertIs(get_priority([a1]), np.ndarray)

    def test_array2d(self):
        a2 = np.array([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]])
        self.assertIs(get_priority([a2]), np.ndarray)

    def test_series(self):
        s = pd.Series([0.0, 0.1, 0.2])
        self.assertIs(get_priority([s]), pd.Series)

    def test_dataframe(self):
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]])
        self.assertIs(get_priority([f]), pd.DataFrame)

    def test_array_series(self):
        a = np.array([0.0, 0.1, 0.2])
        s = pd.Series([0.0, 0.1, 0.2])
        self.assertIs(get_priority([a, s]), pd.Series)
        self.assertIs(get_priority([s, a]), pd.Series)

    def test_array_dataframe(self):
        a = np.array([0.0, 0.1, 0.2])
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]])
        self.assertIs(get_priority([a, f]), pd.DataFrame)
        self.assertIs(get_priority([f, a]), pd.DataFrame)

    def test_series_dataframe(self):
        s = pd.Series([0.0, 0.1, 0.2])
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]])
        self.assertIs(get_priority([s, f]), pd.DataFrame)
        self.assertIs(get_priority([f, s]), pd.DataFrame)

    def test_all(self):
        a = np.array([0.0, 0.1, 0.2])
        s = pd.Series([0.0, 0.1, 0.2])
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]])
        for p in permutations([a, s, f]):
            self.assertIs(get_priority(p), pd.DataFrame)

    def test_empty(self):
        self.assertRaisesRegex(ValueError,
                               "Cannot determine priority of zero arrays",
                               get_priority,
                               [])

    def test_non_array(self):
        self.assertRaisesRegex(TypeError,
                               "Cannot determine priority of array types",
                               get_priority,
                               [0.0, np.array([0.0, 0.1, 0.2])])


class TestGetSharedIndex(ut.TestCase):

    def test_empty(self):
        self.assertRaisesRegex(ValueError,
                               "No indexes were given",
                               get_shared_index,
                               [])

    def test_one_index(self):
        index = pd.Index([2, 5, 6, 8])
        self.assertIs(get_shared_index([index]), index)

    def test_one_multiindex(self):
        index = pd.MultiIndex.from_arrays([["A", "B", "G", "H"], [2, 5, 6, 8]],
                                          names=["letter", "number"])
        self.assertIs(get_shared_index([index]), index)

    def test_identical_indexes(self):
        a = pd.Index([2, 5, 6, 8])
        b = pd.Index([2, 5, 6, 8])
        self.assertIs(get_shared_index([a, b]), a)

    def test_equal_indexes(self):
        a = pd.Index([2, 5, 6, 8])
        b = pd.Index([2., 5., 6., 8.])
        self.assertIs(get_shared_index([a, b]), a)

    def test_equal_multiindexes(self):
        a = pd.MultiIndex.from_arrays([["A", "B", "G", "H"], [2, 5, 6, 8]],
                                      names=["letter", "number"])
        b = pd.MultiIndex.from_arrays([["A", "B", "G", "H"], [2., 5., 6., 8.]],
                                      names=["letter", "number"])
        self.assertIs(get_shared_index([a, b]), a)

    def test_unequal_values_indexes(self):
        a = pd.Index([2, 5, 6, 8])
        b = pd.Index([2, 5, 6, 9])
        self.assertRaisesRegex(ValueError,
                               "Indexes 0 and 1 differ",
                               get_shared_index,
                               [a, b])

    def test_unequal_names_indexes(self):
        a = pd.Index([2, 5, 6, 8], name="a")
        b = pd.Index([2, 5, 6, 8], name="b")
        self.assertRaisesRegex(ValueError,
                               "Indexes 0 and 1 differ",
                               get_shared_index,
                               [a, b])

    def test_unequal_values_multiindexes(self):
        a = pd.MultiIndex.from_arrays([["A", "B", "G", "H"], [2, 5, 6, 8]],
                                      names=["letter", "number"])
        b = pd.MultiIndex.from_arrays([["A", "B", "G", "H"], [2, 5, 6, 9]],
                                      names=["letter", "number"])
        self.assertRaisesRegex(ValueError,
                               "Indexes 0 and 1 differ",
                               get_shared_index,
                               [a, b])

    def test_unequal_names_multiindexes(self):
        a = pd.MultiIndex.from_arrays([["A", "B", "G", "H"], [2, 5, 6, 8]],
                                      names=["letter", "number"])
        b = pd.MultiIndex.from_arrays([["A", "B", "G", "H"], [2, 5, 6, 8]],
                                      names=["Letter", "Number"])
        self.assertRaisesRegex(ValueError,
                               "Indexes 0 and 1 differ",
                               get_shared_index,
                               [a, b])

    def test_mismatched_types(self):
        index = pd.Index([2, 5, 6, 8])
        multi = pd.MultiIndex.from_arrays([["A", "B", "G", "H"], [2, 5, 6, 8]],
                                          names=["letter", "number"])
        self.assertRaisesRegex(TypeError,
                               "Expected Index, but got MultiIndex",
                               get_shared_index,
                               [index, multi])
        self.assertRaisesRegex(TypeError,
                               "Expected MultiIndex, but got Index",
                               get_shared_index,
                               [multi, index])

    def test_not_index(self):
        index = pd.Index([2, 5, 6, 8])
        array = np.array([2, 5, 6, 8])
        self.assertRaisesRegex(TypeError,
                               "Expected Index, but got ndarray",
                               get_shared_index,
                               [index, array])
        self.assertRaisesRegex(TypeError,
                               "Expected Index, but got ndarray",
                               get_shared_index,
                               [array, index])


class TestGetSharedIndexes(ut.TestCase):

    def test_array(self):
        a = np.array([0.0, 0.1, 0.2])
        self.assertEqual(get_shared_indexes([a]), dict())

    def test_series(self):
        s = pd.Series([0.0, 0.1, 0.2], index=[4, 7, 9])
        indexes = get_shared_indexes([s])
        self.assertEqual(list(indexes), ["index"])
        self.assertTrue(indexes["index"].equals(pd.Index([4, 7, 9])))

    def test_dataframe(self):
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                         index=[3, 6],
                         columns=["x", "y", "z"])
        indexes = get_shared_indexes([f])
        self.assertEqual(list(indexes), ["index", "columns"])
        self.assertTrue(indexes["index"].equals(pd.Index([3, 6])))
        self.assertTrue(indexes["columns"].equals(pd.Index(["x", "y", "z"])))


    def test_empty(self):
        self.assertRaisesRegex(ValueError,
                               "Cannot determine priority of zero arrays",
                               get_shared_indexes,
                               [])

    def test_non_array(self):
        self.assertRaisesRegex(TypeError,
                               "Cannot determine priority of array types",
                               get_shared_indexes,
                               [0.0, np.array([0.0, 0.1, 0.2])])


if __name__ == "__main__":
    ut.main()
