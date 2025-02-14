import unittest as ut
from itertools import product

import numpy as np

from seismicrna.core.batch.ends import (BadSegmentEndsError,
                                        count_reads_segments,
                                        find_contiguous_reads,
                                        find_read_end5s,
                                        find_read_end3s,
                                        mask_segment_ends,
                                        match_reads_segments,
                                        merge_read_ends,
                                        merge_segment_ends,
                                        sort_segment_ends)

rng = np.random.default_rng(0)


class TestCountReadsSegments(ut.TestCase):

    def test_valid(self):
        for num_reads in range(5):
            for num_segs in range(5):
                ends = rng.integers(1, 10, (num_reads, num_segs))
                self.assertEqual(count_reads_segments(ends),
                                 (num_reads, num_segs))

    def test_non_array(self):
        self.assertRaisesRegex(
            TypeError,
            ("xyz must be an instance of <class 'numpy.ndarray'>, "
             "but got 8 of type <class 'numpy.int64'>"),
            count_reads_segments,
            rng.integers(1, 10), "xyz"
        )

    def test_wrong_dims(self):
        for ndim in range(5):
            if ndim != 2:
                ends = rng.integers(1, 10, (3,) * ndim)
                self.assertRaisesRegex(
                    BadSegmentEndsError,
                    f"Must have xyz.ndim = 2, but got {ndim}",
                    count_reads_segments,
                    ends, "xyz"
                )


class TestMatchReadsSegments(ut.TestCase):

    def test_match(self):
        for r1, r2 in product(range(5), repeat=2):
            for s1, s2 in product(range(5), repeat=2):
                end5s = rng.integers(1, 10, (r1, s1))
                end3s = rng.integers(1, 10, (r2, s2))
                if r1 != r2:
                    self.assertRaisesRegex(
                        BadSegmentEndsError,
                        ("Must have num_reads_5 = num_reads_3, but got "
                         f"num_reads_5={r1} and num_reads_3={r2}"),
                        match_reads_segments,
                        end5s, end3s
                    )
                elif s1 != s2:
                    self.assertRaisesRegex(
                        BadSegmentEndsError,
                        ("Must have num_segs_5 = num_segs_3, but got "
                         f"num_segs_5={s1} and num_segs_3={s2}"),
                        match_reads_segments,
                        end5s, end3s
                    )
                else:
                    self.assertEqual(match_reads_segments(end5s, end3s),
                                     end5s.shape)


class TestMaskSegmentEnds(ut.TestCase):

    def test_nomask(self):
        for num_reads in range(5):
            for num_segs in range(5):
                end5s = rng.integers(1, 10, (num_reads, num_segs))
                end3s = end5s + rng.integers(1, 10, (num_reads, num_segs))
                mask5s, mask3s = mask_segment_ends(end5s, end3s)
                self.assertIs(mask5s, end5s)
                self.assertIs(mask3s, end3s)

    def test_mask(self):
        end5s = np.array([[1, 2],
                          [1, 3],
                          [2, 1],
                          [2, 3],
                          [3, 2],
                          [2, 3],
                          [3, 2]])
        end3s = np.array([[2, 3],
                          [1, 3],
                          [3, 2],
                          [1, 2],
                          [2, 1],
                          [2, 2],
                          [2, 2]])
        mask5s, mask3s = mask_segment_ends(end5s, end3s)
        self.assertIsNot(mask5s, end5s)
        self.assertIsNot(mask3s, end3s)
        self.assertIsInstance(mask5s, np.ma.MaskedArray)
        self.assertIsInstance(mask3s, np.ma.MaskedArray)
        masked = end5s > end3s
        self.assertTrue(np.array_equal(mask5s.mask, masked))
        self.assertTrue(np.array_equal(mask3s.mask, masked))


class TestMergeSegmentEnds(ut.TestCase):

    def test_nomask(self):
        end5s = np.array([[1, 2],
                          [1, 3],
                          [2, 1],
                          [2, 3],
                          [3, 2],
                          [2, 3],
                          [3, 2]])
        end3s = np.array([[2, 3],
                          [1, 3],
                          [3, 2],
                          [1, 2],
                          [2, 1],
                          [2, 2],
                          [2, 2]])
        expect = np.array([[1, 2, 2, 3],
                           [1, 3, 1, 3],
                           [2, 1, 3, 2],
                           [2, 3, 1, 2],
                           [3, 2, 2, 1],
                           [2, 3, 2, 2],
                           [3, 2, 2, 2]])
        result = merge_segment_ends(end5s, end3s)
        self.assertNotIsInstance(result, np.ma.MaskedArray)
        self.assertTrue(np.array_equal(result, expect))

    def test_mask(self):
        end5s = np.array([[1, 2],
                          [1, 3],
                          [2, 1],
                          [2, 3],
                          [3, 2],
                          [2, 3],
                          [3, 2]])
        end3s = np.array([[2, 3],
                          [1, 3],
                          [3, 2],
                          [1, 2],
                          [2, 1],
                          [2, 2],
                          [2, 2]])
        expect_data = np.array([[1, 2, 2, 3],
                                [1, 3, 1, 3],
                                [2, 1, 3, 2],
                                [2, 3, 1, 2],
                                [3, 2, 2, 1],
                                [2, 3, 2, 2],
                                [3, 2, 2, 2]])
        expect_mask = np.hstack([end5s > end3s] * 2)
        result = merge_segment_ends(*mask_segment_ends(end5s, end3s))
        self.assertIsInstance(result, np.ma.MaskedArray)
        self.assertTrue(np.array_equal(result.data, expect_data))
        self.assertTrue(np.array_equal(result.mask, expect_mask))


class TestFindReadEnd5s(ut.TestCase):

    def test_0_segments(self):
        end5s = np.array([[],
                          []])
        expect = np.array([1, 1])
        self.assertTrue(np.array_equal(find_read_end5s(end5s), expect))

    def test_1_segment_none_masked(self):
        end5s = np.array([[5],
                          [2]])
        expect = np.array([5, 2])
        self.assertTrue(np.array_equal(find_read_end5s(end5s), expect))

    def test_1_segment_masked(self):
        end5s = np.ma.masked_array([[5],
                                    [2]],
                                   [[True],
                                    [False]],
                                   fill_value=0)
        expect = np.array([1, 2])
        self.assertTrue(np.array_equal(find_read_end5s(end5s), expect))

    def test_2_segments_none_masked(self):
        end5s = np.array([[5, 4],
                          [2, 9]])
        expect = np.array([4, 2])
        self.assertTrue(np.array_equal(find_read_end5s(end5s), expect))

    def test_2_segments_masked(self):
        end5s = np.ma.masked_array([[5, 4],
                                    [2, 9]],
                                   [[False, True],
                                    [True, False]])
        expect = np.array([5, 9])
        self.assertTrue(np.array_equal(find_read_end5s(end5s), expect))

    def test_2_segments_both_masked(self):
        end5s = np.ma.masked_array([[5, 4],
                                    [2, 9]],
                                   [[False, True],
                                    [True, True]],
                                   fill_value=0)
        expect = np.array([5, 1])
        self.assertTrue(np.array_equal(find_read_end5s(end5s), expect))


class TestFindReadEnd3s(ut.TestCase):

    def test_0_segments(self):
        end3s = np.array([[],
                          []])
        expect = np.array([0, 0])
        self.assertTrue(np.array_equal(find_read_end3s(end3s), expect))

    def test_1_segment_none_masked(self):
        end3s = np.array([[5],
                          [2]])
        expect = np.array([5, 2])
        self.assertTrue(np.array_equal(find_read_end3s(end3s), expect))

    def test_1_segment_masked(self):
        end3s = np.ma.masked_array([[5],
                                    [2]],
                                   [[True],
                                    [False]],
                                   fill_value=0)
        expect = np.array([0, 2])
        self.assertTrue(np.array_equal(find_read_end3s(end3s), expect))

    def test_2_segments_none_masked(self):
        end3s = np.array([[5, 4],
                          [2, 9]])
        expect = np.array([5, 9])
        self.assertTrue(np.array_equal(find_read_end3s(end3s), expect))

    def test_2_segments_masked(self):
        end3s = np.ma.masked_array([[5, 4],
                                    [2, 9]],
                                   [[True, False],
                                    [False, True]],
                                   fill_value=0)
        expect = np.array([4, 2])
        self.assertTrue(np.array_equal(find_read_end3s(end3s), expect))

    def test_2_segments_both_masked(self):
        end3s = np.ma.masked_array([[5, 4],
                                    [2, 9]],
                                   [[False, True],
                                    [True, True]],
                                   fill_value=0)
        expect = np.array([5, 0])
        self.assertTrue(np.array_equal(find_read_end3s(end3s), expect))


class TestMergeReadEnds(ut.TestCase):

    def test_valid(self):
        for num_reads in range(5):
            end5s = rng.integers(1, 10, (num_reads,))
            end3s = rng.integers(1, 10, (num_reads,))
            expect = np.stack([end5s, end3s], axis=1)
            self.assertTrue(np.array_equal(merge_read_ends(end5s, end3s),
                                           expect))

    def test_mismatch(self):
        for r1, r2 in product(range(5), repeat=2):
            if r1 != r2:
                end5s = rng.integers(1, 10, (r1,))
                end3s = rng.integers(1, 10, (r2,))
                self.assertRaises(ValueError, merge_read_ends, end5s, end3s)

    def test_wrong_dims(self):
        for ndim in range(3):
            if ndim != 1:
                dims = (3,) * ndim
                end5s = rng.integers(1, 10, dims)
                end3s = rng.integers(1, 10, dims)
                self.assertRaises(ValueError, merge_read_ends, end5s, end3s)


class TestSortSegmentEnds(ut.TestCase):

    def test_0_segments(self):
        end5s = np.array([[],
                          []])
        end3s = np.array([[],
                          []])
        expect_ends = np.array([[],
                                []])
        expect_seg = np.array([[],
                               []])
        expect_contig = np.array([[],
                                  []])
        result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                   end3s)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_seg, expect_seg))
        self.assertTrue(np.array_equal(result_contig, expect_contig))

    def test_1_segment_nomask(self):
        end5s = np.array([[1],
                          [1]])
        end3s = np.array([[1],
                          [2]])
        expect_ends = np.array([[0, 1],
                                [0, 2]])
        expect_seg = np.array([[True, False],
                               [True, False]])
        expect_contig = np.array([[False, True],
                                  [False, True]])
        result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                   end3s)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_seg, expect_seg))
        self.assertTrue(np.array_equal(result_contig, expect_contig))

    def test_1_segment_mask(self):
        mask = np.array([[True],
                         [True]])
        end5s = np.ma.masked_array([[1],
                                    [1]],
                                   mask)
        end3s = np.ma.masked_array([[1],
                                    [2]],
                                   mask)
        mask = np.array([[True, True],
                         [True, True]])
        expect_ends = np.ma.masked_array([[0, 1],
                                          [0, 2]],
                                         mask)
        expect_seg = np.ma.masked_array([[True, False],
                                         [True, False]],
                                        mask)
        expect_contig = np.ma.masked_array([[False, True],
                                            [False, True]],
                                           mask)
        result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                   end3s)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_seg, expect_seg))
        self.assertTrue(np.array_equal(result_contig, expect_contig))

    def test_2_segments_nomask(self):
        end5s = np.array([[1, 4],
                          [1, 3],
                          [1, 2],
                          [1, 1],
                          [2, 1],
                          [3, 1],
                          [4, 1]])
        end3s = np.array([[2, 5],
                          [2, 4],
                          [2, 3],
                          [2, 2],
                          [3, 2],
                          [4, 2],
                          [5, 2]])
        expect_ends = np.array([[0, 2, 3, 5],
                                [0, 2, 2, 4],
                                [0, 1, 2, 3],
                                [0, 0, 2, 2],
                                [0, 1, 2, 3],
                                [0, 2, 2, 4],
                                [0, 2, 3, 5]])
        expect_seg = np.array([[True, False, True, False],
                               [True, True, False, False],
                               [True, True, False, False],
                               [True, True, False, False],
                               [True, True, False, False],
                               [True, True, False, False],
                               [True, False, True, False]])
        expect_contig = np.array([[False, True, False, True],
                                  [False, False, False, True],
                                  [False, False, False, True],
                                  [False, False, False, True],
                                  [False, False, False, True],
                                  [False, False, False, True],
                                  [False, True, False, True]])
        result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                   end3s)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_seg, expect_seg))
        self.assertTrue(np.array_equal(result_contig, expect_contig))

    def test_2_segments_mask(self):
        mask = np.array([[True, False],
                         [True, False],
                         [True, False],
                         [True, False],
                         [False, True],
                         [False, True],
                         [False, True],
                         [False, True]])
        end5s = np.ma.masked_array([[3, 4],
                                    [3, 3],
                                    [3, 2],
                                    [3, 1],
                                    [1, 3],
                                    [2, 3],
                                    [3, 3],
                                    [4, 3]],
                                   mask)
        end3s = np.ma.masked_array([[2, 5],
                                    [2, 4],
                                    [2, 3],
                                    [2, 2],
                                    [2, 2],
                                    [3, 2],
                                    [4, 2],
                                    [5, 2]],
                                   mask)
        mask = np.array([[True, True, False, False],
                         [True, True, False, False],
                         [True, True, False, False],
                         [True, False, True, False],
                         [False, True, True, False],
                         [True, True, False, False],
                         [True, True, False, False],
                         [True, True, False, False]])
        expect_ends = np.ma.masked_array([[2, 2, 3, 5],
                                          [2, 2, 2, 4],
                                          [2, 2, 1, 3],
                                          [2, 0, 2, 2],
                                          [0, 2, 2, 2],
                                          [2, 2, 1, 3],
                                          [2, 2, 2, 4],
                                          [2, 2, 3, 5]],
                                         mask)
        expect_seg = np.ma.masked_array([[True, False, True, False],
                                         [True, False, True, False],
                                         [True, False, True, False],
                                         [True, True, False, False],
                                         [True, True, False, False],
                                         [True, False, True, False],
                                         [True, False, True, False],
                                         [True, False, True, False]],
                                        mask)
        expect_contig = np.ma.masked_array([[False, True, False, True],
                                            [False, True, False, True],
                                            [False, True, False, True],
                                            [False, False, False, True],
                                            [False, False, False, True],
                                            [False, True, False, True],
                                            [False, True, False, True],
                                            [False, True, False, True]],
                                           mask)
        result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                   end3s)
        self.assertTrue(np.all(result_ends == expect_ends))
        self.assertTrue(np.array_equal(result_ends.mask,
                                       expect_ends.mask))
        self.assertTrue(np.all(result_seg == expect_seg))
        self.assertTrue(np.array_equal(result_seg.mask,
                                       expect_seg.mask))
        self.assertTrue(np.all(result_contig == expect_contig))
        self.assertTrue(np.array_equal(result_contig.mask,
                                       expect_contig.mask))

    def test_2_segments_fill_mask(self):
        mask = np.array([[True, False],
                         [True, False],
                         [True, False],
                         [True, False],
                         [False, True],
                         [False, True],
                         [False, True],
                         [False, True]])
        end5s = np.ma.masked_array([[3, 4],
                                    [3, 3],
                                    [3, 2],
                                    [3, 1],
                                    [1, 3],
                                    [2, 3],
                                    [3, 3],
                                    [4, 3]],
                                   mask)
        end3s = np.ma.masked_array([[2, 5],
                                    [2, 4],
                                    [2, 3],
                                    [2, 2],
                                    [2, 2],
                                    [3, 2],
                                    [4, 2],
                                    [5, 2]],
                                   mask)
        expect_ends = np.array([[0, 0, 3, 5],
                                [0, 0, 2, 4],
                                [0, 0, 1, 3],
                                [0, 0, 0, 2],
                                [0, 0, 0, 2],
                                [0, 0, 1, 3],
                                [0, 0, 2, 4],
                                [0, 0, 3, 5]])
        expect_seg = np.array([[True, False, True, False],
                               [True, False, True, False],
                               [True, False, True, False],
                               [True, True, False, False],
                               [True, True, False, False],
                               [True, False, True, False],
                               [True, False, True, False],
                               [True, False, True, False]])
        expect_contig = np.array([[False, True, False, True],
                                  [False, True, False, True],
                                  [False, True, False, True],
                                  [False, False, False, True],
                                  [False, False, False, True],
                                  [False, True, False, True],
                                  [False, True, False, True],
                                  [False, True, False, True]])
        result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                   end3s,
                                                                   True)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_seg, expect_seg))
        self.assertTrue(np.array_equal(result_contig, expect_contig))

    def test_n_segments_contig(self):
        for n in range(1, 10):
            ends = rng.permutation(np.arange(1, n + 1))[np.newaxis, :]
            expect_seg = np.array([True]
                                  + [True, False] * (n - 1)
                                  + [False])[np.newaxis, :]
            expect_contig = np.array([False] * (2 * n - 1)
                                     + [True])[np.newaxis, :]
            expect_ends = np.concatenate([[0]]
                                         + [[i, i] for i in range(1, n)]
                                         + [[n]])[np.newaxis, :]
            result_ends, result_seg, result_contig = sort_segment_ends(ends,
                                                                       ends)
            self.assertTrue(np.array_equal(result_ends, expect_ends))
            self.assertTrue(np.array_equal(result_seg, expect_seg))
            self.assertTrue(np.array_equal(result_contig, expect_contig))

    def test_n_segments_discontig(self):
        for n in range(1, 10):
            ends = rng.permutation(np.arange(1, n + 1) * 2)[np.newaxis, :]
            expect_seg = np.array([True, False] * n)[np.newaxis, :]
            expect_contig = np.logical_not(expect_seg)
            expect_ends = np.arange(1, 2 * n + 1)[np.newaxis, :]
            result_ends, result_seg, result_contig = sort_segment_ends(ends,
                                                                       ends)
            self.assertTrue(np.array_equal(result_ends, expect_ends))
            self.assertTrue(np.array_equal(result_seg, expect_seg))
            self.assertTrue(np.array_equal(result_contig, expect_contig))


class TestFindContiguousReads(ut.TestCase):

    def test_0_segments(self):
        for n_reads in range(10):
            end5s = np.zeros((n_reads, 0))
            end3s = np.zeros((n_reads, 0))
            self.assertTrue(np.array_equal(find_contiguous_reads(end5s, end3s),
                                           np.ones(n_reads, dtype=bool)))

    def test_1_segment(self):
        for n_reads in range(10):
            end5s = rng.integers(1, 10, (n_reads, 1))
            end3s = rng.integers(1, 10, (n_reads, 1)) + end5s
            self.assertTrue(np.array_equal(find_contiguous_reads(end5s, end3s),
                                           np.ones(n_reads, dtype=bool)))

    def test_2_segments_nomask(self):
        end5s = np.array([[1, 4],
                          [1, 3],
                          [1, 2],
                          [1, 1],
                          [2, 1],
                          [3, 1],
                          [4, 1]])
        end3s = np.array([[2, 5],
                          [2, 4],
                          [2, 3],
                          [2, 2],
                          [3, 2],
                          [4, 2],
                          [5, 2]])
        expect = np.array([False, True, True, True, True, True, False])
        self.assertTrue(np.array_equal(find_contiguous_reads(end5s, end3s),
                                       expect))

    def test_2_segments_mask(self):
        mask = np.array([[False, False],
                         [False, False],
                         [False, False],
                         [False, False],
                         [False, False],
                         [False, False],
                         [False, False],
                         [True, False],
                         [True, False],
                         [True, False],
                         [True, False],
                         [False, True],
                         [False, True],
                         [False, True],
                         [False, True]])
        end5s = np.ma.masked_array([[1, 4],
                                    [1, 3],
                                    [1, 2],
                                    [1, 1],
                                    [2, 1],
                                    [3, 1],
                                    [4, 1],
                                    [3, 4],
                                    [3, 3],
                                    [3, 2],
                                    [3, 1],
                                    [1, 3],
                                    [2, 3],
                                    [3, 3],
                                    [4, 3]],
                                   mask)
        end3s = np.ma.masked_array([[2, 5],
                                    [2, 4],
                                    [2, 3],
                                    [2, 2],
                                    [3, 2],
                                    [4, 2],
                                    [5, 2],
                                    [2, 5],
                                    [2, 4],
                                    [2, 3],
                                    [2, 2],
                                    [2, 2],
                                    [3, 2],
                                    [4, 2],
                                    [5, 2]],
                                   mask)
        expect = np.array([False, True, True, True, True, True, False,
                           True, True, True, True, True, True, True, True])
        self.assertTrue(np.array_equal(find_contiguous_reads(end5s, end3s),
                                       expect))


if __name__ == "__main__":
    ut.main(verbosity=2)
