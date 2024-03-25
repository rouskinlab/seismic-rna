import unittest as ut

import numpy as np

from seismicrna.core.dims import find_dims, triangular


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
                               "Unknown dimension for nonzero: 'x'",
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
                                       "Unknown dimension for nonzero: 'z'",
                                       find_dims,
                                       dims,
                                       arrays,
                                       nonzero="z")


class TestTriangular(ut.TestCase):

    def test_whole_numbers(self):
        for n in range(10):
            self.assertEqual(triangular(n), sum(range(1, n + 1)))

    def test_negative(self):
        self.assertRaisesRegex(ValueError,
                               "n must be ≥ 0, but got -1",
                               triangular,
                               -1)

    def test_float(self):
        self.assertRaisesRegex(TypeError,
                               "n must be int, but got float",
                               triangular,
                               1.)


if __name__ == "__main__":
    ut.main()
