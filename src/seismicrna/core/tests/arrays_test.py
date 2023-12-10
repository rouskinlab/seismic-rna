import unittest as ut
from itertools import combinations, permutations

import numpy as np
import pandas as pd

from seismicrna.core.arrays import (get_first_index,
                                    get_last_index,
                                    get_priority,
                                    get_shared_index,
                                    get_shared_indexes,
                                    get_shared_shape,
                                    np_internal,
                                    np_internal_broadcast,
                                    np2_internal,
                                    promote_arrays,
                                    rearray)

rng = np.random.default_rng()


class TestGetFirstIndex(ut.TestCase):

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "Cannot get first index of 0-D array",
                               get_first_index,
                               np.array(1))

    def test_array1d(self):
        for n in range(10):
            a = rng.random(n)
            self.assertTrue(np.array_equal(get_first_index(a),
                                           np.arange(n)))

    def test_array2d(self):
        for m in range(10):
            for n in range(5):
                a = rng.random((m, n))
                self.assertTrue(np.array_equal(get_first_index(a),
                                               np.arange(m)))

    def test_series(self):
        for n in range(10):
            a = pd.Series(rng.random(n))
            i = get_first_index(a)
            self.assertIsInstance(i, pd.Index)
            self.assertTrue(i.equals(pd.RangeIndex(n)))

    def test_dataframe(self):
        for m in range(10):
            for n in range(5):
                a = pd.DataFrame(rng.random((m, n)))
                i = get_first_index(a)
                self.assertTrue(i.equals(pd.RangeIndex(m)))


class TestGetLastIndex(ut.TestCase):

    def test_array0d(self):
        self.assertRaisesRegex(ValueError,
                               "Cannot get last index of 0-D array",
                               get_last_index,
                               np.array(1))

    def test_array1d(self):
        for n in range(10):
            a = rng.random(n)
            self.assertTrue(np.array_equal(get_last_index(a),
                                           np.arange(n)))

    def test_array2d(self):
        for m in range(10):
            for n in range(5):
                a = rng.random((m, n))
                self.assertTrue(np.array_equal(get_last_index(a),
                                               np.arange(n)))

    def test_series(self):
        for n in range(10):
            a = pd.Series(rng.random(n))
            i = get_last_index(a)
            self.assertIsInstance(i, pd.Index)
            self.assertTrue(i.equals(pd.RangeIndex(n)))

    def test_dataframe(self):
        for m in range(10):
            for n in range(5):
                a = pd.DataFrame(rng.random((m, n)))
                i = get_last_index(a)
                self.assertTrue(i.equals(pd.RangeIndex(n)))


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
        for arrays in ([a], [a, a]):
            self.assertEqual(get_shared_indexes(arrays), dict())

    def test_series(self):
        s = pd.Series([0.0, 0.1, 0.2], index=[4, 7, 9])
        for arrays in ([s], [s, s]):
            idxs = get_shared_indexes(arrays)
            self.assertEqual(list(idxs), ["index"])
            self.assertTrue(idxs["index"].equals(pd.Index([4, 7, 9])))

    def test_dataframe(self):
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                         index=[3, 6],
                         columns=["x", "y", "z"])
        for arrays in ([f], [f, f]):
            idxs = get_shared_indexes(arrays)
            self.assertEqual(list(idxs), ["index", "columns"])
            self.assertTrue(idxs["index"].equals(pd.Index([3, 6])))
            self.assertTrue(idxs["columns"].equals(pd.Index(["x", "y", "z"])))

    def test_mismatched(self):
        a = np.array([0.0, 0.1, 0.2])
        s = pd.Series([0.0, 0.1, 0.2], index=[4, 7, 9])
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                         index=[3, 6],
                         columns=["x", "y", "z"])
        for arrays in combinations([a, s, f], 2):
            self.assertRaisesRegex(TypeError,
                                   "Expected every array to be",
                                   get_shared_indexes,
                                   arrays)

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


class TestGetSharedShape(ut.TestCase):

    def test_empty(self):
        self.assertRaisesRegex(ValueError,
                               "Cannot determine shape of zero arrays",
                               get_shared_shape,
                               [])

    def test_array(self):
        a = np.array([0.0, 0.1, 0.2])
        for arrays in ([a], [a, a]):
            self.assertEqual(get_shared_shape(arrays), (3,))

    def test_series(self):
        s = pd.Series([0.0, 0.1, 0.2], index=[4, 7, 9])
        for arrays in ([s], [s, s]):
            self.assertEqual(get_shared_shape(arrays), (3,))

    def test_dataframe(self):
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                         index=[3, 6],
                         columns=["x", "y", "z"])
        for arrays in ([f], [f, f]):
            self.assertEqual(get_shared_shape(arrays), (2, 3))

    def test_mismatched(self):
        a = np.array([[0.0], [0.1], [0.2]])
        s = pd.Series([0.0, 0.1, 0.2], index=[4, 7, 9])
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                         index=[3, 5],
                         columns=["x", "y", "z"])
        for arrays in combinations([a, s, f], 2):
            self.assertRaisesRegex(ValueError,
                                   "Got more than one shape for arrays",
                                   get_shared_shape,
                                   arrays)

    def test_non_array(self):
        self.assertRaisesRegex(AttributeError,
                               "object has no attribute 'shape'",
                               get_shared_shape,
                               [0.0, np.array([0.0, 0.1, 0.2])])


class TestPromoteArrays(ut.TestCase):

    def test_array(self):
        a = np.array([0.0, 0.1, 0.2])
        p = promote_arrays(a)
        self.assertIsInstance(p, tuple)
        self.assertEqual(len(p), 1)
        self.assertIsInstance(p[0], np.ndarray)
        self.assertTrue(np.array_equal(p[0], a))

    def test_series(self):
        s = pd.Series([0.0, 0.1, 0.2], index=[4, 7, 9])
        p = promote_arrays(s)
        self.assertIsInstance(p, tuple)
        self.assertEqual(len(p), 1)
        self.assertIsInstance(p[0], pd.Series)
        self.assertTrue(p[0].equals(s))

    def test_dataframe(self):
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                         index=[3, 6],
                         columns=["x", "y", "z"])
        p = promote_arrays(f)
        self.assertIsInstance(p, tuple)
        self.assertEqual(len(p), 1)
        self.assertIsInstance(p[0], pd.DataFrame)
        self.assertTrue(p[0].equals(f))

    def test_array1d_series(self):
        a = np.array([0.0, 0.1, 0.2])
        s = pd.Series([0.3, 0.4, 0.5], index=[4, 7, 9])
        for args in [(a, s), (s, a)]:
            ps = promote_arrays(*args)
            self.assertIsInstance(ps, tuple)
            self.assertEqual(len(ps), 2)
            for p in ps:
                self.assertIsInstance(p, pd.Series)
                self.assertTrue(p.index.equals(s.index))

    def test_array1d_series_misdims(self):
        a = np.array([0.0, 0.1])
        s = pd.Series([0.3, 0.4, 0.5], index=[4, 7, 9])
        for args in [(a, s), (s, a)]:
            self.assertRaisesRegex(ValueError,
                                   "Got more than one shape for arrays",
                                   promote_arrays,
                                   *args)

    def test_array2d_series(self):
        a = np.array([[0.0, 0.1], [0.2, 0.3], [0.4, 0.5]])
        s = pd.Series([0.3, 0.4, 0.5], index=[4, 7, 9])
        for args in [(a, s), (s, a)]:
            self.assertRaisesRegex(ValueError,
                                   "Got more than one shape for arrays",
                                   promote_arrays,
                                   *args)

    def test_array1d_dataframe(self):
        a = np.array([0.0, 0.1])
        f = pd.DataFrame([[0.0], [0.3]],
                         index=[3, 6],
                         columns=["x"])
        for args in [(a, f), (f, a)]:
            self.assertRaisesRegex(ValueError,
                                   "Got more than one shape for arrays",
                                   promote_arrays,
                                   *args)

    def test_array2d_dataframe(self):
        a = np.array([[0.5, 0.3, 0.2], [0.4, 0.0, 0.1]])
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                         index=[3, 6],
                         columns=["x", "y", "z"])
        for args in [(a, f), (f, a)]:
            ps = promote_arrays(*args)
            self.assertIsInstance(ps, tuple)
            self.assertEqual(len(ps), 2)
            for p in ps:
                self.assertIsInstance(p, pd.DataFrame)
                self.assertTrue(p.index.equals(f.index))
                self.assertTrue(p.columns.equals(f.columns))

    def test_array2d_dataframe_misdims(self):
        a = np.array([[0.5, 0.3, 0.2]])
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                         index=[3, 6],
                         columns=["x", "y", "z"])
        for args in [(a, f), (f, a)]:
            self.assertRaisesRegex(ValueError,
                                   "Got more than one shape for arrays",
                                   promote_arrays,
                                   *args)

    def test_series_dataframe(self):
        s = pd.Series([0.3, 0.4, 0.5], index=[4, 7, 9])
        f = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                         index=[3, 6],
                         columns=["x", "y", "z"])
        for args in [(s, f), (f, s)]:
            self.assertRaisesRegex(ValueError,
                                   "Got more than one shape for arrays",
                                   promote_arrays,
                                   *args)


class TestRearray(ut.TestCase):

    def test_array1d_array0d(self):
        a1 = np.array([0.1, 0.7, 0.4])
        a2 = np.array(0.5)
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1))
        self.assertRaisesRegex(ValueError,
                               "cannot broadcast a non-scalar to a scalar",
                               rearray,
                               a1, a2, broadcast=True)
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, np.ndarray)
        self.assertTrue(np.array_equal(r2, a2))
        r2b = rearray(a2, a1, broadcast=True)
        self.assertIsInstance(r2b, np.ndarray)
        self.assertTrue(np.array_equal(r2b, np.array([0.5, 0.5, 0.5])))

    def test_array1d_array1d(self):
        a1 = np.array([0.1, 0.7, 0.4])
        a2 = np.array([0.5, 0.2, 0.8])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1))
        r1b = rearray(a1, a2, broadcast=True)
        self.assertIsInstance(r1b, np.ndarray)
        self.assertTrue(np.array_equal(r1b, a1))
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, np.ndarray)
        self.assertTrue(np.array_equal(r2, a2))
        r2b = rearray(a2, a1, broadcast=True)
        self.assertIsInstance(r2b, np.ndarray)
        self.assertTrue(np.array_equal(r2b, a2))

    def test_array1d_array1d_misdims(self):
        a1 = np.array([0.1, 0.7, 0.4])
        a2 = np.array([0.5, 0.2, 0.8, 0.0])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1))
        self.assertRaisesRegex(ValueError,
                               "operands could not be broadcast together",
                               rearray,
                               a1, a2, broadcast=True)
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, np.ndarray)
        self.assertTrue(np.array_equal(r2, a2))
        self.assertRaisesRegex(ValueError,
                               "operands could not be broadcast together",
                               rearray,
                               a2, a1, broadcast=True)

    def test_array1d_array2d(self):
        a1 = np.array([0.1, 0.7, 0.4])
        a2 = np.array([[0.5, 0.2, 0.8], [0.0, 0.6, 0.9]])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1))
        r1b = rearray(a1, a2, broadcast=True)
        self.assertIsInstance(r1b, np.ndarray)
        self.assertTrue(np.array_equal(r1b, np.array([[0.1, 0.7, 0.4],
                                                      [0.1, 0.7, 0.4]])))
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, np.ndarray)
        self.assertTrue(np.array_equal(r2, a2))
        self.assertRaisesRegex(ValueError,
                               "input operand has more dimensions than allowed",
                               rearray,
                               a2, a1, broadcast=True)

    def test_array1d_array2d_misdims(self):
        a1 = np.array([0.1, 0.7, 0.4])
        a2 = np.array([[0.5, 0.2], [0.8, 0.0], [0.6, 0.9]])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1))
        self.assertRaisesRegex(ValueError,
                               "operands could not be broadcast together",
                               rearray,
                               a1, a2, broadcast=True)
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, np.ndarray)
        self.assertTrue(np.array_equal(r2, a2))
        self.assertRaisesRegex(ValueError,
                               "input operand has more dimensions than allowed",
                               rearray,
                               a2, a1, broadcast=True)

    def test_series_array0d(self):
        a1 = pd.Series([0.5, 0.2, 0.8], [2, 7, 9])
        a2 = np.array(0.5)
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1.values))
        self.assertRaisesRegex(ValueError,
                               "cannot broadcast a non-scalar to a scalar",
                               rearray,
                               a1, a2, broadcast=True)
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, np.ndarray)
        self.assertTrue(np.array_equal(r2, a2))
        r2b = rearray(a2, a1, broadcast=True)
        self.assertIsInstance(r2b, pd.Series)
        self.assertTrue(r2b.equals(pd.Series([0.5, 0.5, 0.5], [2, 7, 9])))

    def test_series_array1d(self):
        a1 = pd.Series([0.5, 0.2, 0.8], [2, 7, 9])
        a2 = np.array([0.1, 0.7, 0.4])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1.values))
        r1b = rearray(a1, a2, broadcast=True)
        self.assertIsInstance(r1b, np.ndarray)
        self.assertTrue(np.array_equal(r1b, a1.values))
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, pd.Series)
        self.assertTrue(np.array_equal(r2.values, a2))
        self.assertTrue(r2.index.equals(a1.index))
        r2b = rearray(a2, a1, broadcast=True)
        self.assertIsInstance(r2b, pd.Series)
        self.assertTrue(np.array_equal(r2b.values, a2))
        self.assertTrue(r2b.index.equals(a1.index))

    def test_series_array2d(self):
        a1 = pd.Series([0.5, 0.2, 0.8], [2, 7, 9])
        a2 = np.array([[0.5, 0.2, 0.8], [0.0, 0.6, 0.9]])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1.values))
        r1b = rearray(a1, a2, broadcast=True)
        self.assertIsInstance(r1b, np.ndarray)
        self.assertTrue(np.array_equal(r1b, np.array([[0.5, 0.2, 0.8],
                                                      [0.5, 0.2, 0.8]])))
        self.assertRaisesRegex(ValueError,
                               "Cannot return a 2-dimensional array",
                               rearray,
                               a2, a1)
        self.assertRaisesRegex(ValueError,
                               "input operand has more dimensions than allowed",
                               rearray,
                               a2, a1, broadcast=True)

    def test_dataframe_array0d(self):
        a1 = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                          index=[3, 6],
                          columns=["x", "y", "z"])
        a2 = np.array(0.5)
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1.values))
        self.assertRaisesRegex(ValueError,
                               "cannot broadcast a non-scalar to a scalar",
                               rearray,
                               a1, a2, broadcast=True)
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, np.ndarray)
        self.assertTrue(np.array_equal(r2, a2))
        r2b = rearray(a2, a1, broadcast=True)
        self.assertIsInstance(r2b, pd.DataFrame)
        self.assertTrue(r2b.equals(pd.DataFrame([[0.5, 0.5, 0.5],
                                                 [0.5, 0.5, 0.5]],
                                                index=[3, 6],
                                                columns=["x", "y", "z"])))

    def test_dataframe_array1d(self):
        a1 = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                          index=[3, 6],
                          columns=["x", "y", "z"])
        a2 = np.array([0.1, 0.7, 0.4])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1.values))
        self.assertRaisesRegex(ValueError,
                               "input operand has more dimensions than allowed",
                               rearray,
                               a1, a2, broadcast=True)
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, pd.Series)
        self.assertTrue(r2.equals(pd.Series([0.1, 0.7, 0.4], ["x", "y", "z"])))
        r2b = rearray(a2, a1, broadcast=True)
        self.assertIsInstance(r2b, pd.DataFrame)
        self.assertTrue(r2b.equals(pd.DataFrame([[0.1, 0.7, 0.4],
                                                 [0.1, 0.7, 0.4]],
                                                index=[3, 6],
                                                columns=["x", "y", "z"])))

    def test_dataframe_array1d_misdims(self):
        a1 = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                          index=[3, 6],
                          columns=["x", "y", "z"])
        a2 = np.array([0.1, 0.7])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1.values))
        self.assertRaisesRegex(ValueError,
                               "input operand has more dimensions than allowed",
                               rearray,
                               a1, a2, broadcast=True)
        self.assertRaisesRegex(ValueError,
                               "does not match length of index",
                               rearray,
                               a2, a1)
        self.assertRaisesRegex(ValueError,
                               "operands could not be broadcast together",
                               rearray,
                               a2, a1, broadcast=True)

    def test_dataframe_array2d(self):
        a1 = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                          index=[3, 6],
                          columns=["x", "y", "z"])
        a2 = np.array([[0.1, 0.7, 0.4], [0.9, 0.3, 0.2]])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1.values))
        r1b = rearray(a1, a2, broadcast=True)
        self.assertIsInstance(r1b, np.ndarray)
        self.assertTrue(np.array_equal(r1b, a1.values))
        r2 = rearray(a2, a1)
        self.assertIsInstance(r2, pd.DataFrame)
        self.assertTrue(r2.equals(pd.DataFrame([[0.1, 0.7, 0.4],
                                                [0.9, 0.3, 0.2]],
                                               index=[3, 6],
                                               columns=["x", "y", "z"])))
        r2b = rearray(a2, a1, broadcast=True)
        self.assertIsInstance(r2b, pd.DataFrame)
        self.assertTrue(r2b.equals(pd.DataFrame([[0.1, 0.7, 0.4],
                                                 [0.9, 0.3, 0.2]],
                                                index=[3, 6],
                                                columns=["x", "y", "z"])))

    def test_dataframe_array2d_misdims(self):
        a1 = pd.DataFrame([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5]],
                          index=[3, 6],
                          columns=["x", "y", "z"])
        a2 = np.array([[0.1, 0.7], [0.4, 0.9], [0.3, 0.2]])
        r1 = rearray(a1, a2)
        self.assertIsInstance(r1, np.ndarray)
        self.assertTrue(np.array_equal(r1, a1.values))
        self.assertRaisesRegex(ValueError,
                               "operands could not be broadcast together",
                               rearray,
                               a1, a2, broadcast=True)
        self.assertRaisesRegex(ValueError,
                               "indices imply",
                               rearray,
                               a2, a1)
        self.assertRaisesRegex(ValueError,
                               "operands could not be broadcast together",
                               rearray,
                               a2, a1, broadcast=True)


class TestNpInternal(ut.TestCase):

    def test_np_internal(self):

        def _error_if_not_array(array: np.ndarray):
            if not isinstance(array, np.ndarray):
                raise TypeError(f"Not an array: {array}")
            return array

        wrapped = np_internal(_error_if_not_array)

        a = np.array([0.3, 0.6, 0.7])
        s = pd.Series([0.3, 0.6, 0.7], [3, 4, 6])
        f = pd.DataFrame([[0.3, 0.4], [0.6, 0.7], [0.7, 0.8]],
                         index=[3, 4, 6],
                         columns=[8, 9])
        self.assertTrue(isinstance(_error_if_not_array(a), np.ndarray))
        self.assertTrue(np.array_equal(_error_if_not_array(a), a))
        for x in [s, f]:
            self.assertRaisesRegex(TypeError,
                                   "Not an array",
                                   _error_if_not_array,
                                   x)
        self.assertTrue(isinstance(wrapped(a), np.ndarray))
        self.assertTrue(np.array_equal(wrapped(a), a))
        self.assertTrue(isinstance(wrapped(s), pd.Series))
        self.assertTrue(wrapped(s).equals(s))
        self.assertTrue(isinstance(wrapped(f), pd.DataFrame))
        self.assertTrue(wrapped(f).equals(f))

    def test_np_internal_broadcast(self):

        def _error_if_not_array(array: np.ndarray):
            if not isinstance(array, np.ndarray):
                raise TypeError(f"Not an array: {array}")
            return array.sum(axis=0)

        wrapped = np_internal_broadcast(_error_if_not_array)

        a = np.array([0., 1., 2.])
        s = pd.Series([1., 2., 3.],
                      index=[3, 4, 6])
        f = pd.DataFrame([[0., 2.],
                          [4., 6.],
                          [8., 0.]],
                         index=[3, 4, 6],
                         columns=[8, 9])
        self.assertTrue(isinstance(_error_if_not_array(a), float))
        self.assertTrue(np.array_equal(_error_if_not_array(a), a.sum()))
        for x in [s, f]:
            self.assertRaisesRegex(TypeError,
                                   "Not an array",
                                   _error_if_not_array,
                                   x)
        self.assertTrue(isinstance(wrapped(a), np.ndarray))
        self.assertTrue(np.array_equal(wrapped(a),
                                       np.array([3., 3., 3.])))
        self.assertTrue(isinstance(wrapped(s), pd.Series))
        self.assertTrue((wrapped(s).equals(pd.Series([6., 6., 6.],
                                                     index=[3, 4, 6]))))
        self.assertTrue(isinstance(wrapped(f), pd.DataFrame))
        self.assertTrue(wrapped(f).equals(pd.DataFrame([[12., 8.],
                                                        [12., 8.],
                                                        [12., 8.]],
                                                       index=[3, 4, 6],
                                                       columns=[8, 9])))


class TestNp2Internal(ut.TestCase):

    def test_np2_internal(self):
        def _error_if_not_array(array1: np.ndarray, array2: np.ndarray):
            for array in [array1, array2]:
                if not isinstance(array, np.ndarray):
                    raise TypeError(f"Not an array: {array}")
            return array1 + array2

        wrapped = np2_internal(_error_if_not_array)

        a1 = np.array([0., 1., 2.])
        s = pd.Series([1., 2., 3.],
                      index=[0, 2, 5])
        a2 = np.array([[9., 7.],
                       [3., 3.],
                       [1., 3.]])
        f = pd.DataFrame([[0., 2.],
                          [4., 6.],
                          [8., 0.]],
                         index=[3, 4, 6],
                         columns=[8, 9])
        for x in [a1, a2]:
            self.assertTrue(isinstance(_error_if_not_array(x, x), np.ndarray))
            self.assertTrue(np.array_equal(_error_if_not_array(x, x), x + x))
        for xs in [(a1, s), (s, a1), (a2, f), (f, a2)]:
            self.assertRaisesRegex(TypeError,
                                   "Not an array",
                                   _error_if_not_array,
                                   *xs)
        for x in [a1, a2]:
            self.assertTrue(isinstance(wrapped(x, x), np.ndarray))
            self.assertTrue(np.array_equal(wrapped(x, x), x + x))
        for xs in [(a1, s), (s, a1)]:
            self.assertTrue(isinstance(wrapped(*xs), pd.Series))
            self.assertTrue(wrapped(*xs).equals(pd.Series([1., 3., 5.],
                                                          index=[0, 2, 5])))
        for xs in [(a2, f), (f, a2)]:
            self.assertTrue(isinstance(wrapped(*xs), pd.DataFrame))
            self.assertTrue(wrapped(*xs).equals(pd.DataFrame([[9., 9.],
                                                              [7., 9.],
                                                              [9., 3.]],
                                                             index=[3, 4, 6],
                                                             columns=[8, 9])))


if __name__ == "__main__":
    ut.main()
