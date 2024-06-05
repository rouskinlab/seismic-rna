import unittest as ut
from itertools import product

import numpy as np

from seismicrna.core.batch.ends import (_check_no_coverage_reads,
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
            for num_segs in range(1, 5):
                ends = rng.integers(1, 10, (num_reads, num_segs))
                self.assertEqual(count_reads_segments(ends),
                                 (num_reads, num_segs))

    def test_0_segments(self):
        for num_reads in range(5):
            ends = rng.integers(1, 10, (num_reads, 0))
            self.assertRaisesRegex(ValueError,
                                   "xyz has 0 segments",
                                   count_reads_segments,
                                   ends, "xyz")

    def test_non_array(self):
        self.assertRaisesRegex(TypeError,
                               "xyz must be ndarray, but got int",
                               count_reads_segments,
                               rng.integers(1, 10), "xyz")

    def test_wrong_dims(self):
        for ndim in range(5):
            if ndim != 2:
                ends = rng.integers(1, 10, (3,) * ndim)
                self.assertRaisesRegex(ValueError,
                                       f"xyz must have 2 dimensions, "
                                       f"but got {ndim}",
                                       count_reads_segments,
                                       ends, "xyz")


class TestMatchReadsSegments(ut.TestCase):

    def test_match(self):
        for r1, r2 in product(range(5), repeat=2):
            for s1, s2 in product(range(1, 5), repeat=2):
                end5s = rng.integers(1, 10, (r1, s1))
                end3s = rng.integers(1, 10, (r2, s2))
                if r1 != r2:
                    self.assertRaisesRegex(ValueError,
                                           "Numbers of 5' and 3' reads must "
                                           f"equal, but got {r1} ≠ {r2}",
                                           match_reads_segments,
                                           end5s, end3s)
                elif s1 != s2:
                    self.assertRaisesRegex(ValueError,
                                           "Numbers of 5' and 3' segments must "
                                           f"equal, but got {s1} ≠ {s2}",
                                           match_reads_segments,
                                           end5s, end3s)
                else:
                    self.assertEqual(match_reads_segments(end5s, end3s),
                                     end5s.shape)


class TestMaskSegmentEnds(ut.TestCase):

    def test_nomask(self):
        for num_reads in range(5):
            for num_segs in range(1, 5):
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


class TestCheckNoCoverageReads(ut.TestCase):

    def test_1_segment_all_covered(self):
        ends = np.ma.masked_array([[1]], [[False]])
        self.assertEqual(_check_no_coverage_reads(ends),
                         count_reads_segments(ends))

    def test_1_segment_none_covered(self):
        ends = np.ma.masked_array([[1]], [[True]])
        self.assertRaisesRegex(ValueError,
                               "Got 1 read[(]s[)] with no coverage",
                               _check_no_coverage_reads,
                               ends)

    def test_2_segments_all_covered(self):
        ends = np.ma.masked_array([[1, 2]], [[False, False]])
        self.assertEqual(_check_no_coverage_reads(ends),
                         count_reads_segments(ends))

    def test_2_segments_half_covered(self):
        ends = np.ma.masked_array([[1, 2]], [[False, True]])
        self.assertEqual(_check_no_coverage_reads(ends),
                         count_reads_segments(ends))

    def test_2_segments_none_covered(self):
        ends = np.ma.masked_array([[1, 2]], [[True, True]])
        self.assertRaisesRegex(ValueError,
                               "Got 1 read[(]s[)] with no coverage",
                               _check_no_coverage_reads,
                               ends)


class TestFindReadEnd5s(ut.TestCase):

    def test_1_segment(self):
        end5s = np.ma.masked_array([[5],
                                    [2]],
                                   [[False],
                                    [False]])
        expect = np.array([5, 2])
        self.assertTrue(np.array_equal(find_read_end5s(end5s), expect))

    def test_2_segments_nomask(self):
        end5s = np.ma.masked_array([[5, 4],
                                    [2, 9]],
                                   [[False, False],
                                    [False, False]])
        expect = np.array([4, 2])
        self.assertTrue(np.array_equal(find_read_end5s(end5s), expect))

    def test_2_segments_mask(self):
        end5s = np.ma.masked_array([[5, 4],
                                    [2, 9]],
                                   [[False, True],
                                    [True, False]])
        expect = np.array([5, 9])
        self.assertTrue(np.array_equal(find_read_end5s(end5s), expect))


class TestFindReadEnd3s(ut.TestCase):

    def test_1_segment(self):
        end3s = np.ma.masked_array([[5],
                                    [2]],
                                   [[False],
                                    [False]])
        expect = np.array([5, 2])
        self.assertTrue(np.array_equal(find_read_end3s(end3s), expect))

    def test_2_segments_nomask(self):
        end3s = np.ma.masked_array([[5, 4],
                                    [2, 9]],
                                   [[False, False],
                                    [False, False]])
        expect = np.array([5, 9])
        self.assertTrue(np.array_equal(find_read_end3s(end3s), expect))

    def test_2_segments_mask(self):
        end3s = np.ma.masked_array([[5, 4],
                                    [2, 9]],
                                   [[True, False],
                                    [False, True]])
        expect = np.array([4, 2])
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

    def test_1_segment(self):
        end5s = np.array([[1],
                          [1]])
        end3s = np.array([[1],
                          [2]])
        expect_ends_0 = np.array([[0, 1],
                                  [0, 2]])
        expect_ends_1 = np.array([[1, 1],
                                  [1, 2]])
        expect_seg = np.array([[True, False],
                               [True, False]])
        expect_contig = np.array([[False, True],
                                  [False, True]])
        for idx0, expect_ends in [(True, expect_ends_0),
                                  (False, expect_ends_1)]:
            result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                       end3s,
                                                                       idx0)
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
        expect_ends_0 = np.array([[0, 2, 3, 5],
                                  [0, 2, 2, 4],
                                  [0, 1, 2, 3],
                                  [0, 0, 2, 2],
                                  [0, 1, 2, 3],
                                  [0, 2, 2, 4],
                                  [0, 2, 3, 5]])
        expect_ends_1 = np.array([[1, 2, 4, 5],
                                  [1, 3, 2, 4],
                                  [1, 2, 2, 3],
                                  [1, 1, 2, 2],
                                  [1, 2, 2, 3],
                                  [1, 3, 2, 4],
                                  [1, 2, 4, 5]])
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
        for idx0, expect_ends in [(True, expect_ends_0),
                                  (False, expect_ends_1)]:
            result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                       end3s,
                                                                       idx0)
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
        expect_ends_0 = np.ma.masked_array([[2, 2, 3, 5],
                                            [2, 2, 2, 4],
                                            [2, 2, 1, 3],
                                            [2, 0, 2, 2],
                                            [0, 2, 2, 2],
                                            [2, 2, 1, 3],
                                            [2, 2, 2, 4],
                                            [2, 2, 3, 5]],
                                           mask)
        expect_ends_1 = np.ma.masked_array([[3, 2, 4, 5],
                                            [3, 2, 3, 4],
                                            [3, 2, 2, 3],
                                            [3, 1, 2, 2],
                                            [1, 3, 2, 2],
                                            [3, 2, 2, 3],
                                            [3, 2, 3, 4],
                                            [3, 2, 4, 5]],
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
        for idx0, expect_ends in [(True, expect_ends_0),
                                  (False, expect_ends_1)]:
            result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                       end3s,
                                                                       idx0)
            self.assertTrue(np.all(result_ends == expect_ends))
            self.assertTrue(np.array_equal(result_ends.mask,
                                           expect_ends.mask))
            self.assertTrue(np.all(result_seg == expect_seg))
            self.assertTrue(np.array_equal(result_seg.mask,
                                           expect_seg.mask))
            self.assertTrue(np.all(result_contig == expect_contig))
            self.assertTrue(np.array_equal(result_contig.mask,
                                           expect_contig.mask))

    def test_2_segments_unmask(self):
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
        expect_ends_0 = np.array([[0, 0, 3, 5],
                                  [0, 0, 2, 4],
                                  [0, 0, 1, 3],
                                  [0, 0, 0, 2],
                                  [0, 0, 0, 2],
                                  [0, 0, 1, 3],
                                  [0, 0, 2, 4],
                                  [0, 0, 3, 5]])
        expect_ends_1 = np.array([[1, 0, 4, 5],
                                  [1, 0, 3, 4],
                                  [1, 0, 2, 3],
                                  [1, 1, 0, 2],
                                  [1, 1, 0, 2],
                                  [1, 0, 2, 3],
                                  [1, 0, 3, 4],
                                  [1, 0, 4, 5]])
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
        for idx0, expect_ends in [(True, expect_ends_0),
                                  (False, expect_ends_1)]:
            result_ends, result_seg, result_contig = sort_segment_ends(end5s,
                                                                       end3s,
                                                                       idx0,
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
            expect_ends_0 = np.concatenate([[0]]
                                           + [[i, i] for i in range(1, n)]
                                           + [[n]])[np.newaxis, :]
            expect_ends_1 = expect_ends_0 + expect_seg
            for idx0, expect_ends in [(True, expect_ends_0),
                                      (False, expect_ends_1)]:
                result_ends, result_seg, result_contig = sort_segment_ends(ends,
                                                                           ends,
                                                                           idx0)
                self.assertTrue(np.array_equal(result_ends, expect_ends))
                self.assertTrue(np.array_equal(result_seg, expect_seg))
                self.assertTrue(np.array_equal(result_contig, expect_contig))

    def test_n_segments_discontig(self):
        for n in range(1, 10):
            ends = rng.permutation(np.arange(1, n + 1) * 2)[np.newaxis, :]
            expect_seg = np.array([True, False] * n)[np.newaxis, :]
            expect_contig = np.logical_not(expect_seg)
            expect_ends_0 = np.arange(1, 2 * n + 1)[np.newaxis, :]
            expect_ends_1 = expect_ends_0 + expect_seg
            for idx0, expect_ends in [(True, expect_ends_0),
                                      (False, expect_ends_1)]:
                result_ends, result_seg, result_contig = sort_segment_ends(ends,
                                                                           ends,
                                                                           idx0)
                self.assertTrue(np.array_equal(result_ends, expect_ends))
                self.assertTrue(np.array_equal(result_seg, expect_seg))
                self.assertTrue(np.array_equal(result_contig, expect_contig))


class TestFindContiguousReads(ut.TestCase):

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
    ut.main()

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
