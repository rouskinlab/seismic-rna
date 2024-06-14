import unittest as ut

import numpy as np

from seismicrna.core.array import (calc_inverse,
                                   ensure_same_length,
                                   find_dims,
                                   get_length,
                                   locate_elements,
                                   find_true_dists,
                                   triangular)

rng = np.random.default_rng()


class TestGetLength(ut.TestCase):

    def test_1d(self):
        """ A 1D array has a length. """
        for length in range(10):
            array = rng.random(length)
            self.assertEqual(get_length(array), length)

    def test_other_dims(self):
        """ Non-1D arrays are not valid. """
        for d in range(4):
            if d == 1:
                continue
            for length in range(10):
                array = rng.random((length,) * d)
                self.assertRaisesRegex(ValueError,
                                       f"x must have 1 dimension, but got {d}",
                                       get_length,
                                       array, "x")


class TestEnsureSameLength(ut.TestCase):

    def test_same_length(self):
        """ The length is returned when the arrays share the length. """
        for length in range(10):
            x = rng.random(length)
            y = rng.random(length)
            self.assertEqual(ensure_same_length(x, y), length)

    def test_diff_lengths(self):
        """ Unequal lengths are not permitted. """
        for length1 in range(5):
            x = rng.random(length1)
            for length2 in range(5):
                if length1 == length2:
                    continue
                y = rng.random(length2)
                self.assertRaisesRegex(ValueError,
                                       "Lengths differ",
                                       ensure_same_length,
                                       x, y)

    def test_other_dims(self):
        """ Non-1D arrays are not valid. """
        for length in range(5):
            for d1 in range(4):
                x = rng.random((length,) * d1)
                for d2 in range(4):
                    if d1 == 1 == d2:
                        continue
                    y = rng.random((length,) * d2)
                    self.assertRaisesRegex(
                        ValueError,
                        f"{'x' if d1 != 1 else 'y'} must have 1 dimension",
                        ensure_same_length,
                        x, y, "x", "y"
                    )


class TestCalcInverse(ut.TestCase):

    def test_empty(self):
        target = np.array([])
        expect = np.array([])
        self.assertTrue(np.array_equal(calc_inverse(target), expect))

    def test_empty_max(self):
        target = np.array([])
        for maximum in range(5):
            expect = np.full(maximum + 1, -1)
            self.assertTrue(np.array_equal(
                calc_inverse(target, require=maximum),
                expect
            ))

    def test_calc_inverse(self):
        """ Invert an array with a known inverse. """
        target = np.array([4, 7, 3, 8, 1])
        expect = np.array([-1, 4, -1, 2, 0, -1, -1, 1, 3])
        self.assertTrue(np.array_equal(calc_inverse(target), expect))

    def test_calc_inverse_max(self):
        """ Invert an array with a known inverse. """
        target = np.array([4, 7, 3, 8, 1])
        expect = np.array([-1, 4, -1, 2, 0, -1, -1, 1, 3, -1])
        self.assertTrue(np.array_equal(calc_inverse(target, require=9),
                                       expect))

    def test_calc_inverse_fill_fwd(self):
        """ Invert an array with a known inverse; fill forward. """
        target = np.array([4, 7, 3, 8, 1])
        expect = np.array([-1, 4, 4, 2, 0, 0, 0, 1, 3])
        self.assertTrue(np.array_equal(calc_inverse(target, fill=True),
                                       expect))

    def test_calc_inverse_fill_fwd_max(self):
        target = np.array([4, 7, 3, 8, 1])
        expect = np.array([-1, 4, 4, 2, 0, 0, 0, 1, 3, 3])
        self.assertTrue(np.array_equal(
            calc_inverse(target, require=9, fill=True),
            expect
        ))

    def test_calc_inverse_fill_fwd_max_default(self):
        target = np.array([4, 7, 3, 8, 1])
        expect = np.array([23, 4, 4, 2, 0, 0, 0, 1, 3, 3])
        self.assertTrue(np.array_equal(
            calc_inverse(target, require=9, fill=True, fill_default=23),
            expect
        ))

    def test_calc_inverse_fill_rev(self):
        """ Invert an array with a known inverse; fill reverse. """
        target = np.array([4, 7, 3, 8, 1])
        expect = np.array([4, 4, 2, 2, 0, 1, 1, 1, 3])
        self.assertTrue(np.array_equal(
            calc_inverse(target, fill=True, fill_rev=True),
            expect
        ))

    def test_calc_inverse_fill_rev_max(self):
        target = np.array([4, 7, 3, 8, 1])
        expect = np.array([4, 4, 2, 2, 0, 1, 1, 1, 3, 5])
        self.assertTrue(np.array_equal(
            calc_inverse(target, require=9, fill=True, fill_rev=True),
            expect
        ))

    def test_calc_inverse_fill_rev_max_default(self):
        target = np.array([4, 7, 3, 8, 1])
        expect = np.array([4, 4, 2, 2, 0, 1, 1, 1, 3, 23])
        self.assertTrue(np.array_equal(
            calc_inverse(target,
                         require=9,
                         fill=True,
                         fill_rev=True,
                         fill_default=23),
            expect
        ))

    def test_is_inverse(self):
        """ Verify every position in the inverse of a random array. """
        target = rng.choice(16, 8, replace=False)
        for maximum in [-1, 0, 1, 20]:
            inverse = calc_inverse(target, require=maximum)
            for i, inv in enumerate(inverse):
                if inv == -1:
                    self.assertNotIn(i, target)
                else:
                    self.assertEqual(i, target[inv])

    def test_negative(self):
        """ Negative values in the target are not permitted. """
        target = np.array([4, 7, -3, 8, 1])
        self.assertRaisesRegex(ValueError,
                               "test_negative has negative values",
                               calc_inverse,
                               target,
                               what="test_negative")

    def test_repeated(self):
        """ Repeated values in the target are not permitted. """
        target = np.array([4, 7, 3, 8, 7, 1])
        self.assertRaisesRegex(ValueError,
                               "test_repeated has repeated values",
                               calc_inverse,
                               target,
                               what="test_repeated")


class TestLocateElements(ut.TestCase):

    def test_locate_0(self):
        collection = np.array([4, 1, 2, 7, 5, 3])
        self.assertEqual(locate_elements(collection), ())

    def test_locate_1(self):
        collection = np.array([4, 1, 2, 7, 5, 3])
        elements = np.array([5, 2, 5])
        expect = np.array([4, 2, 4])
        self.assertTrue(np.array_equal(locate_elements(collection, elements),
                                       expect))

    def test_locate_2(self):
        collection = np.array([4, 1, 2, 7, 5, 3])
        elements = np.array([3, 7, 4]), np.array([1, 5])
        expect = (np.array([5, 3, 0]), np.array([1, 4]))
        result = locate_elements(collection, *elements)
        self.assertTrue(all(np.array_equal(res, exp)
                            for res, exp in zip(result, expect, strict=True)))

    def test_locate_invalid(self):
        collection = np.array([4, 1, 2, 7, 5, 3])
        elements = np.array([5, 2, 6])
        self.assertRaisesRegex(ValueError,
                               r"Elements \[6\] are not in collection",
                               locate_elements,
                               collection,
                               elements)


class TestFindDims(ut.TestCase):

    def test_empty(self):
        arrays = []
        dims = []
        self.assertEqual(find_dims(dims, arrays), {})

    def test_0d(self):
        arrays = [np.zeros(())]
        dims = [()]
        self.assertEqual(find_dims(dims, arrays), {})

    def test_1d(self):
        for x in range(5):
            arrays = [np.zeros(x)]
            dims = [("x",)]
            self.assertEqual(find_dims(dims, arrays), {"x": x})

    def test_2d(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros((x, y))]
                dims = [("x", "y")]
                self.assertEqual(find_dims(dims, arrays), {"x": x, "y": y})

    def test_1d_1d_separate(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros(x), np.zeros(y)]
                dims = [("x",), ("y",)]
                self.assertEqual(find_dims(dims, arrays), {"x": x, "y": y})

    def test_1d_1d_crossed(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros(x), np.zeros(y)]
                dims = [("x",), ("x",)]
                if x == y:
                    self.assertEqual(find_dims(dims, arrays), {"x": x})
                else:
                    self.assertRaisesRegex(ValueError,
                                           "Got multiple sizes for dimension",
                                           find_dims,
                                           dims,
                                           arrays)

    def test_1d_2d_congruent(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros(x), np.zeros((y, x))]
                dims = [("x",), ("y", "x")]
                self.assertEqual(find_dims(dims, arrays), {"x": x, "y": y})

    def test_1d_2d_crossed(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros(y), np.zeros((y, x))]
                dims = [("x",), ("y", "x")]
                if x == y:
                    self.assertEqual(find_dims(dims, arrays), {"x": x, "y": y})
                else:
                    self.assertRaisesRegex(ValueError,
                                           "Got multiple sizes for dimension",
                                           find_dims,
                                           dims,
                                           arrays)

    def test_2d_2d_congruent(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros((x, y)), np.zeros((y, x))]
                dims = [("x", "y"), ("y", "x")]
                self.assertEqual(find_dims(dims, arrays), {"x": x, "y": y})

    def test_2d_2d_crossed(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros((x, y)), np.zeros((y, x))]
                dims = [("x", "y"), ("x", "y")]
                if x == y:
                    self.assertEqual(find_dims(dims, arrays), {"x": x, "y": y})
                else:
                    self.assertRaisesRegex(ValueError,
                                           "Got multiple sizes for dimension",
                                           find_dims,
                                           dims,
                                           arrays)

    def test_0d_none(self):
        arrays = [np.zeros(())]
        dims = [(None,)]
        self.assertEqual(find_dims(dims, arrays), {})

    def test_0d_1dim(self):
        arrays = [np.zeros(())]
        dims = [("x",)]
        self.assertRaisesRegex(ValueError,
                               "Array 'array0' must have 1",
                               find_dims,
                               dims,
                               arrays)

    def test_1d_0dim_none(self):
        for x in range(5):
            arrays = [np.zeros(x)]
            dims = [(None,)]
            self.assertEqual(find_dims(dims, arrays), {})

    def test_1d_1dim_none(self):
        for x in range(5):
            arrays = [np.zeros(x)]
            dims = [("x", None)]
            self.assertEqual(find_dims(dims, arrays), {"x": x})

    def test_1d_2dim(self):
        for x in range(5):
            arrays = [np.zeros(x)]
            dims = [("x", "y")]
            self.assertRaisesRegex(ValueError,
                                   "Array 'array0' must have 2",
                                   find_dims,
                                   dims,
                                   arrays)

    def test_1d_2dim_none(self):
        for x in range(5):
            arrays = [np.zeros(x)]
            dims = [("x", "y", None)]
            self.assertRaisesRegex(ValueError,
                                   "Array 'array0' must have ≥ 2",
                                   find_dims,
                                   dims,
                                   arrays)

    def test_2d_1dim_none(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros((x, y))]
                dims = [("x", None)]
                self.assertEqual(find_dims(dims, arrays), {"x": x})

    def test_none_2d(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros((x, y))]
                dims = [(None, "y")]
                self.assertRaisesRegex(TypeError,
                                       "The name of each dimension must be str",
                                       find_dims,
                                       dims,
                                       arrays)

    def test_0d_nonzero(self):
        arrays = [np.zeros(())]
        dims = [(None,)]
        self.assertEqual(find_dims(dims, arrays, nonzero=True), {})

    def test_0d_nonzero_extra(self):
        arrays = [np.zeros(())]
        dims = [(None,)]
        self.assertRaisesRegex(ValueError,
                               "Unknown dimensions for nonzero: {'x'}",
                               find_dims,
                               dims,
                               arrays,
                               nonzero="x")

    def test_1d_nonzero(self):
        for x in range(5):
            arrays = [np.zeros(x)]
            dims = [("x",)]
            if x:
                self.assertEqual(find_dims(dims, arrays, nonzero="x"), {"x": x})
            else:
                self.assertRaisesRegex(ValueError,
                                       "Size of dimension 'x' must be ≥ 1, "
                                       "but got 0",
                                       find_dims,
                                       dims,
                                       arrays,
                                       nonzero="x")

    def test_2d_nonzero(self):
        for x in range(5):
            for y in range(5):
                arrays = [np.zeros((x, y))]
                dims = [("x", "y")]
                if x:
                    self.assertEqual(find_dims(dims, arrays, nonzero="x"),
                                     {"x": x, "y": y})
                else:
                    self.assertRaisesRegex(ValueError,
                                           "Size of dimension 'x' must be ≥ 1, "
                                           "but got 0",
                                           find_dims,
                                           dims,
                                           arrays,
                                           nonzero="x")
                if y:
                    self.assertEqual(find_dims(dims, arrays, nonzero="y"),
                                     {"x": x, "y": y})
                else:
                    self.assertRaisesRegex(ValueError,
                                           "Size of dimension 'y' must be ≥ 1, "
                                           "but got 0",
                                           find_dims,
                                           dims,
                                           arrays,
                                           nonzero="y")
                self.assertRaisesRegex(ValueError,
                                       "Unknown dimensions for nonzero: {'z'}",
                                       find_dims,
                                       dims,
                                       arrays,
                                       nonzero="z")


class TestFindTrueDists(ut.TestCase):

    def test_all_true(self):
        for n in range(5):
            array = np.ones(n, dtype=bool)
            expect = np.zeros(n, dtype=int)
            self.assertTrue(np.array_equal(find_true_dists(array), expect))

    def test_all_false(self):
        for n in range(5):
            array = np.zeros(n, dtype=bool)
            expect = np.full(n, n)
            self.assertTrue(np.array_equal(find_true_dists(array), expect))

    def test_custom_1(self):
        array = np.array([1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0])
        expect = np.array([0, 0, 1, 0, 1, 1, 0, 1, 2, 1, 0, 1, 2, 3, 4])
        self.assertTrue(np.array_equal(find_true_dists(array), expect))

    def test_custom_2(self):
        array = np.array([0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1])
        expect = np.array([4, 3, 2, 1, 0, 1, 2, 1, 0, 1, 1, 0, 1, 0, 0])
        self.assertTrue(np.array_equal(find_true_dists(array), expect))

    def test_vs_slow(self):

        def slow_find_true_dists(x: np.ndarray):
            length = get_length(x)
            dists = np.full(length, length)
            for i in range(length):
                j = i
                while j >= 0:
                    if x[j]:
                        dists[i] = min(i - j, dists[i])
                        break
                    j -= 1
                j = i
                while j < length:
                    if x[j]:
                        dists[i] = min(j - i, dists[i])
                        break
                    j += 1
            return dists

        for n in range(15):
            array = rng.random(n) < 0.2
            expect = slow_find_true_dists(array)
            self.assertTrue(np.array_equal(find_true_dists(array), expect))


class TestTriangular(ut.TestCase):

    def test_whole_numbers(self):
        for n in range(10):
            self.assertEqual(triangular(n), sum(range(1, n + 1)))


if __name__ == "__main__":
    ut.main()
