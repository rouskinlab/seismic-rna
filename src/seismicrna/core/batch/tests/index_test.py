import unittest as ut

import numpy as np

from seismicrna.core.batch.index import (get_num_segments,
                                         split_ends,
                                         stack_end_coords)

rng = np.random.default_rng()


class TestGetNumSegments(ut.TestCase):

    def test_even(self):
        for nseg in range(10):
            ncol = nseg * 2
            ends = rng.integers(1, 11, (8, ncol))
            self.assertEqual(get_num_segments(ends), nseg)

    def test_odd(self):
        for ncol in range(1, 10, 2):
            ends = rng.integers(1, 11, (8, ncol))
            self.assertRaisesRegex(
                ValueError,
                f"Number of end coordinates must be even, but got {ncol}",
                get_num_segments,
                ends
            )

    def test_wrong_dim(self):
        for ndim in range(5):
            if ndim == 2:
                continue
            ends = rng.integers(1, 11, (2,) * ndim)
            self.assertRaisesRegex(
                ValueError,
                f"ends must have 2 dimensions, but got {ndim}",
                get_num_segments,
                ends
            )


class TestSplitEnds(ut.TestCase):

    def test_1_segment(self):
        for nreads in range(5):
            end5s = rng.integers(1, 11, nreads)
            end3s = rng.integers(1, 11, nreads)
            ends = np.stack([end5s, end3s], axis=1)
            self.assertEqual(ends.shape, (nreads, 2))
            for result, expect in zip(split_ends(ends),
                                      [end5s[: np.newaxis],
                                       end3s[: np.newaxis]],
                                      strict=True):
                self.assertTrue(np.array_equal(result, expect))

    def test_2_segments(self):
        for nreads in range(5):
            end5s1 = rng.integers(1, 11, nreads)
            end3s1 = rng.integers(1, 11, nreads)
            end5s2 = rng.integers(1, 11, nreads)
            end3s2 = rng.integers(1, 11, nreads)
            ends = np.stack([end5s1, end3s1, end5s2, end3s2], axis=1)
            self.assertEqual(ends.shape, (nreads, 4))
            for result, expect in zip(split_ends(ends),
                                      [np.stack([end5s1, end5s2], axis=1),
                                       np.stack([end3s1, end3s2], axis=1)],
                                      strict=True):
                self.assertTrue(np.array_equal(result, expect))

    def test_3_segments(self):
        for nreads in range(5):
            end5s1 = rng.integers(1, 11, nreads)
            end3s1 = rng.integers(1, 11, nreads)
            end5s2 = rng.integers(1, 11, nreads)
            end3s2 = rng.integers(1, 11, nreads)
            end5s3 = rng.integers(1, 11, nreads)
            end3s3 = rng.integers(1, 11, nreads)
            ends = np.stack([end5s1, end3s1, end5s2, end3s2, end5s3, end3s3],
                            axis=1)
            self.assertEqual(ends.shape, (nreads, 6))
            for result, expect in zip(split_ends(ends),
                                      [np.stack([end5s1, end5s2, end5s3],
                                                axis=1),
                                       np.stack([end3s1, end3s2, end3s3],
                                                axis=1)],
                                      strict=True):
                self.assertTrue(np.array_equal(result, expect))


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
