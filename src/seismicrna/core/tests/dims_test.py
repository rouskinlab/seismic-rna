import unittest as ut

import numpy as np

from seismicrna.core.dims import find_dims


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

    def test_2d_none(self):
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

    def test_wrong_dims(self):
        for x in range(5):
            arrays = [np.zeros(x)]
            dims = [("x", "y")]
            self.assertRaisesRegex(ValueError,
                                   "Array 'array0' must have 2",
                                   find_dims,
                                   dims,
                                   arrays)

    def test_wrong_dims_none(self):
        for x in range(5):
            arrays = [np.zeros(x)]
            dims = [("x", "y", None)]
            self.assertRaisesRegex(ValueError,
                                   "Array 'array0' must have ≥ 2",
                                   find_dims,
                                   dims,
                                   arrays)


if __name__ == "__main__":
    ut.main()
