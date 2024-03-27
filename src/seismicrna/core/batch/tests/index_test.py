import unittest as ut

import numpy as np

from seismicrna.core.batch.index import (contiguous_mates,
                                         ensure_same_length,
                                         get_inverse,
                                         get_length,
                                         stack_end_coords)

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


class TestGetInverse(ut.TestCase):

    def test_get_inverse(self):
        """ Invert an array with a known inverse. """
        target = np.array([4, 7, 3, 8, 1])
        expect = np.array([-1, 4, -1, 2, 0, -1, -1, 1, 3])
        self.assertTrue(np.array_equal(get_inverse(target), expect))

    def test_is_inverse(self):
        """ Verify every position in the inverse of a random array. """
        target = rng.choice(16, 8, replace=False)
        inverse = get_inverse(target)
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
                               get_inverse,
                               target,
                               "test_negative")

    def test_repeated(self):
        """ Repeated values in the target are not permitted. """
        target = np.array([4, 7, 3, 8, 7, 1])
        self.assertRaisesRegex(ValueError,
                               "test_repeated has repeated values",
                               get_inverse,
                               target,
                               "test_repeated")


class TestContiguousMates(ut.TestCase):

    def test_contiguous(self):
        """ Contiguous mates can be found from mid5s and mid3s. """
        mid5s = np.array([1, 1, 6, 6, 6])
        mid3s = np.array([1, 2, 3, 4, 5])
        expect = np.array([True, True, False, False, True])
        self.assertTrue(np.array_equal(contiguous_mates(mid5s, mid3s), expect))

    def test_unequal_lengths(self):
        """ Unequal lengths are not permitted. """
        mid5s = np.array([1, 1, 6, 6])
        mid3s = np.array([1, 2, 3, 4, 5])
        self.assertRaisesRegex(ValueError,
                               "Lengths differ",
                               contiguous_mates,
                               mid5s,
                               mid3s)


class TestStackEndCoords(ut.TestCase):

    def test_stack(self):
        for length in range(10):
            end5s = rng.integers(0, 10, length)
            end3s = rng.integers(0, 10, length)
            result = stack_end_coords(end5s, end3s)
            self.assertIsInstance(result, np.ndarray)
            self.assertEqual(result.shape, (length, 2))
            self.assertTrue(np.array_equal(result[:, 0], end5s))
            self.assertTrue(np.array_equal(result[:, 1], end3s))


if __name__ == "__main__":
    ut.main()
