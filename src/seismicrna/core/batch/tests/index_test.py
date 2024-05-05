import unittest as ut

import numpy as np

from seismicrna.core.batch.index import stack_end_coords

rng = np.random.default_rng()


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
