import unittest as ut
from itertools import product

import numpy as np

from seismicrna.core.batch.ends import (BadSegmentEndsError,
                                        count_reads_segments,
                                        find_contiguous_reads,
                                        find_read_end5s,
                                        find_read_end3s,
                                        match_reads_segments,
                                        merge_read_ends,
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
             "but got np.int64[(]8[)] of type <class 'numpy.int64'>"),
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
        for r5, r3, rm in product(range(4), repeat=3):
            for s5, s3, sm in product(range(4), repeat=3):
                end5s = rng.integers(1, 10, (r5, s5))
                end3s = rng.integers(1, 10, (r3, s3))
                mask = np.zeros((rm, sm), dtype=bool)
                if r5 != r3:
                    self.assertRaisesRegex(
                        BadSegmentEndsError,
                        ("Must have num_reads_5 = num_reads_3, but got "
                         f"num_reads_5={r5} and num_reads_3={r3}"),
                        match_reads_segments,
                        end5s, end3s, mask
                    )
                elif s5 != s3:
                    self.assertRaisesRegex(
                        BadSegmentEndsError,
                        ("Must have num_segs_5 = num_segs_3, but got "
                         f"num_segs_5={s5} and num_segs_3={s3}"),
                        match_reads_segments,
                        end5s, end3s, mask
                    )
                else:
                    self.assertEqual(match_reads_segments(end5s, end3s, None),
                                     end5s.shape)
                    if r5 != rm:
                        self.assertRaisesRegex(
                            BadSegmentEndsError,
                            ("Must have num_reads_mask = num_reads_5, but got "
                             f"num_reads_mask={rm} and num_reads_5={r5}"),
                            match_reads_segments,
                            end5s, end3s, mask
                        )
                    elif s5 != sm:
                        self.assertRaisesRegex(
                            BadSegmentEndsError,
                            ("Must have num_segs_mask = num_segs_5, but got "
                             f"num_segs_mask={sm} and num_segs_5={s5}"),
                            match_reads_segments,
                            end5s, end3s, mask
                        )
                    else:
                        self.assertEqual(match_reads_segments(end5s,
                                                              end3s,
                                                              mask),
                                         end5s.shape)


class TestFindReadEnd5s(ut.TestCase):

    def test_0_segments(self):
        end5s = np.array([[],
                          []],
                         dtype=np.uint8)
        original = np.copy(end5s)
        expect = np.array([1, 1])
        self.assertTrue(np.array_equal(find_read_end5s(end5s, None), expect))
        self.assertTrue(np.array_equal(end5s, original))

    def test_1_segment_none_masked(self):
        end5s = np.array([[255],
                          [2]],
                         dtype=np.uint8)
        original = np.copy(end5s)
        expect = np.array([255, 2])
        self.assertTrue(np.array_equal(find_read_end5s(end5s, None), expect))
        self.assertTrue(np.array_equal(end5s, original))

    def test_1_segment_masked(self):
        end5s = np.array([[255],
                          [2]],
                         dtype=np.uint8)
        original = np.copy(end5s)
        mask = np.array([[True],
                         [False]])
        expect = np.array([1, 2])
        self.assertTrue(np.array_equal(find_read_end5s(end5s, mask), expect))
        self.assertTrue(np.array_equal(end5s, original))

    def test_2_segments_none_masked(self):
        end5s = np.array([[5, 4],
                          [2, 9]],
                         dtype=np.uint8)
        original = np.copy(end5s)
        expect = np.array([4, 2])
        self.assertTrue(np.array_equal(find_read_end5s(end5s, None), expect))
        self.assertTrue(np.array_equal(end5s, original))

    def test_2_segments_masked(self):
        end5s = np.array([[5, 4],
                          [2, 9]],
                         dtype=np.uint8)
        original = np.copy(end5s)
        mask = np.array([[False, True],
                         [True, False]])
        expect = np.array([5, 9])
        self.assertTrue(np.array_equal(find_read_end5s(end5s, mask), expect))
        self.assertTrue(np.array_equal(end5s, original))

    def test_2_segments_both_masked(self):
        end5s = np.array([[5, 4],
                          [2, 9]],
                         dtype=np.uint8)
        original = np.copy(end5s)
        mask = np.array([[False, True],
                         [True, True]])
        expect = np.array([5, 1])
        self.assertTrue(np.array_equal(find_read_end5s(end5s, mask), expect))
        self.assertTrue(np.array_equal(end5s, original))


class TestFindReadEnd3s(ut.TestCase):

    def test_0_segments(self):
        end3s = np.array([[],
                          []],
                         dtype=np.uint8)
        original = np.copy(end3s)
        expect = np.array([0, 0])
        self.assertTrue(np.array_equal(find_read_end3s(end3s, None), expect))
        self.assertTrue(np.array_equal(end3s, original))

    def test_1_segment_none_masked(self):
        end3s = np.array([[5],
                          [2]],
                         dtype=np.uint8)
        original = np.copy(end3s)
        expect = np.array([5, 2])
        self.assertTrue(np.array_equal(find_read_end3s(end3s, None), expect))
        self.assertTrue(np.array_equal(end3s, original))

    def test_1_segment_masked(self):
        end3s = np.array([[5],
                          [2]],
                         dtype=np.uint8)
        original = np.copy(end3s)
        mask = np.array([[True],
                         [False]])
        expect = np.array([0, 2])
        self.assertTrue(np.array_equal(find_read_end3s(end3s, mask), expect))
        self.assertTrue(np.array_equal(end3s, original))

    def test_2_segments_none_masked(self):
        end3s = np.array([[5, 4],
                          [2, 9]],
                         dtype=np.uint8)
        original = np.copy(end3s)
        expect = np.array([5, 9])
        self.assertTrue(np.array_equal(find_read_end3s(end3s, None), expect))
        self.assertTrue(np.array_equal(end3s, original))

    def test_2_segments_masked(self):
        end3s = np.array([[5, 4],
                          [2, 9]],
                         dtype=np.uint8)
        original = np.copy(end3s)
        mask = np.array([[True, False],
                         [False, True]])
        expect = np.array([4, 2])
        self.assertTrue(np.array_equal(find_read_end3s(end3s, mask), expect))
        self.assertTrue(np.array_equal(end3s, original))

    def test_2_segments_both_masked(self):
        end3s = np.array([[5, 4],
                          [2, 9]],
                         dtype=np.uint8)
        original = np.copy(end3s)
        mask = np.array([[False, True],
                         [True, True]])
        expect = np.array([5, 0])
        self.assertTrue(np.array_equal(find_read_end3s(end3s, mask), expect))
        self.assertTrue(np.array_equal(end3s, original))


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
        expect_contig5 = np.array([[],
                                   []])
        expect_contig3 = np.array([[],
                                   []])
        result_ends, result_contig5, result_contig3 = sort_segment_ends(end5s,
                                                                        end3s,
                                                                        None)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_contig5, expect_contig5))
        self.assertTrue(np.array_equal(result_contig3, expect_contig3))

    def test_1_segment_nomask(self):
        end5s = np.array([[1],
                          [2]])
        end3s = np.array([[3],
                          [4]])
        expect_ends = np.array([[0, 3],
                                [1, 4]])
        expect_contig5 = np.array([[True, False],
                                   [True, False]])
        expect_contig3 = np.array([[False, True],
                                   [False, True]])
        result_ends, result_contig5, result_contig3 = sort_segment_ends(end5s,
                                                                        end3s,
                                                                        None)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_contig5, expect_contig5))
        self.assertTrue(np.array_equal(result_contig3, expect_contig3))

    def test_1_segment_mask(self):
        mask = np.array([[False],
                         [True]])
        end5s = np.array([[1],
                          [2]])
        end3s = np.array([[3],
                          [0]])
        expect_ends = np.array([[0, 3],
                                [0, 0]])
        expect_contig5 = np.array([[True, False],
                                   [False, False]])
        expect_contig3 = np.array([[False, True],
                                   [False, False]])
        result_ends, result_contig5, result_contig3 = sort_segment_ends(end5s,
                                                                        end3s,
                                                                        mask)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_contig5, expect_contig5))
        self.assertTrue(np.array_equal(result_contig3, expect_contig3))

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
        expect_contig5 = np.array([[True, False, True, False],
                                   [True, False, False, False],
                                   [True, False, False, False],
                                   [True, False, False, False],
                                   [True, False, False, False],
                                   [True, False, False, False],
                                   [True, False, True, False]])
        expect_contig3 = np.array([[False, True, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, True, False, True]])
        result_ends, result_contig5, result_contig3 = sort_segment_ends(end5s,
                                                                        end3s,
                                                                        None)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_contig3, expect_contig3))
        self.assertTrue(np.array_equal(result_contig5, expect_contig5))

    def test_2_segments_mask(self):
        mask = np.array([[True, False],
                         [True, False],
                         [True, False],
                         [True, False],
                         [False, True],
                         [False, True],
                         [False, True],
                         [False, True],
                         [True, False],
                         [False, True]])
        end5s = np.array([[3, 4],
                          [3, 3],
                          [3, 2],
                          [3, 1],
                          [1, 3],
                          [2, 3],
                          [3, 3],
                          [4, 3],
                          [2, 3],
                          [3, 2]])
        end3s = np.array([[2, 5],
                          [2, 4],
                          [2, 3],
                          [2, 2],
                          [2, 2],
                          [3, 2],
                          [4, 2],
                          [5, 2],
                          [1, 4],
                          [4, 1]])
        expect_ends = np.array([[0, 0, 3, 5],
                                [0, 0, 2, 4],
                                [0, 0, 1, 3],
                                [0, 0, 0, 2],
                                [0, 0, 0, 2],
                                [0, 0, 1, 3],
                                [0, 0, 2, 4],
                                [0, 0, 3, 5],
                                [0, 0, 2, 4],
                                [0, 0, 2, 4]])
        expect_contig5 = np.array([[False, False, True, False],
                                   [False, False, True, False],
                                   [False, False, True, False],
                                   [False, True, False, False],
                                   [True, False, False, False],
                                   [False, False, True, False],
                                   [False, False, True, False],
                                   [False, False, True, False],
                                   [False, False, True, False],
                                   [False, False, True, False]])
        expect_contig3 = np.array([[False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True],
                                   [False, False, False, True]])
        result_ends, result_contig5, result_contig3 = sort_segment_ends(end5s,
                                                                        end3s,
                                                                        mask)
        self.assertTrue(np.array_equal(result_ends, expect_ends))
        self.assertTrue(np.array_equal(result_contig5, expect_contig5))
        self.assertTrue(np.array_equal(result_contig3, expect_contig3))

    def test_n_segments_contig(self):
        for n in range(1, 10):
            ends = rng.permutation(np.arange(1, n + 1))[np.newaxis, :]
            expect_contig5 = np.array([True]
                                      + [False] * (2 * n - 1))[np.newaxis, :]
            expect_contig3 = np.array([False] * (2 * n - 1)
                                      + [True])[np.newaxis, :]
            expect_ends = np.concatenate([[0]]
                                         + [[i, i] for i in range(1, n)]
                                         + [[n]])[np.newaxis, :]
            result_ends, result_contig5, result_contig3 = sort_segment_ends(ends,
                                                                            ends,
                                                                            None)
            self.assertTrue(np.array_equal(result_ends, expect_ends))
            self.assertTrue(np.array_equal(result_contig5, expect_contig5))
            self.assertTrue(np.array_equal(result_contig3, expect_contig3))

    def test_n_segments_discontig(self):
        for n in range(1, 10):
            ends = rng.permutation(np.arange(1, n + 1) * 2)[np.newaxis, :]
            expect_contig5 = np.array([True, False] * n)[np.newaxis, :]
            expect_contig3 = np.logical_not(expect_contig5)
            expect_ends = np.arange(1, 2 * n + 1)[np.newaxis, :]
            result_ends, result_contig5, result_contig3 = sort_segment_ends(ends,
                                                                            ends,
                                                                            None)
            self.assertTrue(np.array_equal(result_ends, expect_ends))
            self.assertTrue(np.array_equal(result_contig5, expect_contig5))
            self.assertTrue(np.array_equal(result_contig3, expect_contig3))


class TestFindContiguousReads(ut.TestCase):

    def test_0_segments(self):
        for n_reads in range(10):
            end5s = np.zeros((n_reads, 0))
            end3s = np.zeros((n_reads, 0))
            self.assertTrue(np.array_equal(find_contiguous_reads(end5s,
                                                                 end3s,
                                                                 None),
                                           np.ones(n_reads, dtype=bool)))

    def test_1_segment(self):
        for n_reads in range(10):
            end5s = rng.integers(1, 10, (n_reads, 1))
            end3s = rng.integers(1, 10, (n_reads, 1)) + end5s
            self.assertTrue(np.array_equal(find_contiguous_reads(end5s,
                                                                 end3s,
                                                                 None),
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
        self.assertTrue(np.array_equal(find_contiguous_reads(end5s,
                                                             end3s,
                                                             None),
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
                         [False, True],
                         [False, False],
                         [True, False],
                         [False, True],
                         [True, True]])
        end5s = np.array([[1, 4],
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
                          [4, 3],
                          [4, 1],
                          [4, 1],
                          [4, 1],
                          [4, 1]])
        end3s = np.array([[2, 5],
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
                          [5, 2],
                          [5, 2],
                          [5, 2],
                          [5, 2],
                          [5, 2]])
        expect = np.array([False, True, True, True, True, True, False,
                           True, True, True, True,
                           True, True, True, True,
                           False, True, True, True])
        self.assertTrue(np.array_equal(find_contiguous_reads(end5s,
                                                             end3s,
                                                             mask),
                                       expect))


if __name__ == "__main__":
    ut.main(verbosity=2)
