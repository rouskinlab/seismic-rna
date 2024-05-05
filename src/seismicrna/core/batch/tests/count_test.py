import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.batch.count import (count_end_coords,
                                         get_cover_per_read)
from seismicrna.core.batch.index import END5_COORD, END3_COORD
from seismicrna.core.seq.section import SEQ_INDEX_NAMES


class TestCountEndCoords(ut.TestCase):

    def test_length_0(self):
        end5s = np.array([], dtype=int)
        end3s = np.array([], dtype=int)
        result = count_end_coords(end5s, end3s)
        self.assertIsInstance(result, pd.Series)
        self.assertEqual(result.size, 0)
        self.assertIsInstance(result.index, pd.MultiIndex)
        self.assertEqual(result.index.names, [END5_COORD, END3_COORD])
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END5_COORD),
            np.array([], dtype=int)
        ))
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END3_COORD),
            np.array([], dtype=int)
        ))
        self.assertTrue(np.array_equal(
            result.values,
            np.array([])
        ))

    def test_length_1(self):
        end5s = np.array([4], dtype=int)
        end3s = np.array([6], dtype=int)
        result = count_end_coords(end5s, end3s)
        self.assertIsInstance(result, pd.Series)
        self.assertEqual(result.size, 1)
        self.assertIsInstance(result.index, pd.MultiIndex)
        self.assertEqual(result.index.names, [END5_COORD, END3_COORD])
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END5_COORD),
            np.array([4])
        ))
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END3_COORD),
            np.array([6])
        ))
        self.assertTrue(np.array_equal(
            result.values,
            np.array([1])
        ))

    def test_length_2_diff(self):
        end5s = np.array([4, 1], dtype=int)
        end3s = np.array([6, 8], dtype=int)
        result = count_end_coords(end5s, end3s)
        self.assertIsInstance(result, pd.Series)
        self.assertEqual(result.size, 2)
        self.assertIsInstance(result.index, pd.MultiIndex)
        self.assertEqual(result.index.names, [END5_COORD, END3_COORD])
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END5_COORD),
            np.array([1, 4])
        ))
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END3_COORD),
            np.array([8, 6])
        ))
        self.assertTrue(np.array_equal(
            result.values,
            np.array([1, 1])
        ))

    def test_length_2_same(self):
        end5s = np.array([4, 4], dtype=int)
        end3s = np.array([6, 6], dtype=int)
        result = count_end_coords(end5s, end3s)
        self.assertIsInstance(result, pd.Series)
        self.assertEqual(result.size, 1)
        self.assertIsInstance(result.index, pd.MultiIndex)
        self.assertEqual(result.index.names, [END5_COORD, END3_COORD])
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END5_COORD),
            np.array([4])
        ))
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END3_COORD),
            np.array([6])
        ))
        self.assertTrue(np.array_equal(
            result.values,
            np.array([2])
        ))

    def test_length_3_diff_weights(self):
        end5s = np.array([4, 1, 0], dtype=int)
        end3s = np.array([6, 8, 3], dtype=int)
        weights = pd.DataFrame(np.array([[2., 3.],
                                         [4., 5.],
                                         [6., 7.]]),
                               columns=["a", "b"])
        result = count_end_coords(end5s, end3s, weights)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape, weights.shape)
        self.assertIsInstance(result.index, pd.MultiIndex)
        self.assertEqual(result.index.names, [END5_COORD, END3_COORD])
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END5_COORD),
            np.array([0, 1, 4])
        ))
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END3_COORD),
            np.array([3, 8, 6])
        ))
        self.assertEqual(result.columns.to_list(),
                         weights.columns.to_list())
        self.assertTrue(np.array_equal(
            result.values,
            weights.values[[2, 1, 0]]
        ))

    def test_length_3_same_weights(self):
        end5s = np.array([4, 1, 4], dtype=int)
        end3s = np.array([6, 8, 6], dtype=int)
        weights = pd.DataFrame(np.array([[2., 3.],
                                         [4., 5.],
                                         [6., 7.]]),
                               columns=["a", "b"])
        result = count_end_coords(end5s, end3s, weights)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape, (2, 2))
        self.assertIsInstance(result.index, pd.MultiIndex)
        self.assertEqual(result.index.names, [END5_COORD, END3_COORD])
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END5_COORD),
            np.array([1, 4])
        ))
        self.assertTrue(np.array_equal(
            result.index.get_level_values(END3_COORD),
            np.array([8, 6])
        ))
        self.assertEqual(result.columns.to_list(),
                         weights.columns.to_list())
        self.assertTrue(np.array_equal(
            result.values,
            np.array([[4., 5.],
                      [8., 10.]])
        ))


class TestGetCoverPerRead(ut.TestCase):

    def test_1_segment(self):
        pos_index = pd.MultiIndex.from_tuples([(3, "G"),
                                               (4, "A"),
                                               (5, "C"),
                                               (6, "A"),
                                               (8, "T"),
                                               (9, "T"),
                                               (10, "G"),
                                               (11, "C")],
                                              names=SEQ_INDEX_NAMES)
        ends = np.array([[1, 1],
                         [1, 2],
                         [1, 3],
                         [2, 4],
                         [4, 4],
                         [4, 6],
                         [1, 7],
                         [1, 13],
                         [7, 13],
                         [8, 10],
                         [10, 10],
                         [10, 12],
                         [11, 13],
                         [12, 13],
                         [13, 13]])
        read_nums = np.arange(ends.shape[0])
        expect = {"A": pd.Series([0, 0, 0, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0],
                                 read_nums),
                  "C": pd.Series([0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 1, 1, 0, 0],
                                 read_nums),
                  "G": pd.Series([0, 0, 1, 1, 0, 0, 1, 2, 1, 1, 1, 1, 0, 0, 0],
                                 read_nums),
                  "T": pd.Series([0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0],
                                 read_nums),
                  "N": pd.Series([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 read_nums)}
        result = get_cover_per_read(pos_index, ends, read_nums)
        self.assertEqual(sorted(result), sorted(expect))
        for base in expect:
            self.assertIsInstance(result[base], pd.Series)
            self.assertTrue(np.array_equal(result[base].index.values,
                                           read_nums))
            self.assertTrue(result[base].equals(expect[base]))

    def test_2_segments(self):
        pos_index = pd.MultiIndex.from_tuples([(3, "G"),
                                               (4, "A"),
                                               (5, "C"),
                                               (6, "A"),
                                               (8, "T"),
                                               (9, "T"),
                                               (10, "G"),
                                               (11, "C")],
                                              names=SEQ_INDEX_NAMES)
        ends = np.array([[1, 5, 9, 13],
                         [9, 13, 1, 5],
                         [4, 8, 6, 10],
                         [6, 10, 4, 8],
                         [3, 6, 3, 6],
                         [9, 9, 9, 9],
                         [5, 8, 10, 13],
                         [10, 13, 5, 8],
                         [5, 9, 10, 13],
                         [10, 13, 5, 9],
                         [5, 10, 10, 13],
                         [10, 13, 5, 10],
                         [1, 13, 1, 13]])
        read_nums = np.arange(ends.shape[0])
        expect = {"A": pd.Series([1, 1, 2, 2, 2, 0, 1, 1, 1, 1, 1, 1, 2],
                                 read_nums),
                  "C": pd.Series([2, 2, 1, 1, 1, 0, 2, 2, 2, 2, 2, 2, 2],
                                 read_nums),
                  "G": pd.Series([2, 2, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 2],
                                 read_nums),
                  "T": pd.Series([1, 1, 2, 2, 0, 1, 1, 1, 2, 2, 2, 2, 2],
                                 read_nums),
                  "N": pd.Series([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 read_nums)}
        result = get_cover_per_read(pos_index, ends, read_nums)
        self.assertEqual(sorted(result), sorted(expect))
        for base in expect:
            self.assertIsInstance(result[base], pd.Series)
            self.assertTrue(np.array_equal(result[base].index.values,
                                           read_nums))
            self.assertTrue(result[base].equals(expect[base]))


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
