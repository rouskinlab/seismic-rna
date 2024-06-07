import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.batch.count import (calc_coverage,
                                         _calc_uniq_read_weights,
                                         count_end_coords)
from seismicrna.core.batch.ends import END5_COORD, END3_COORD
from seismicrna.core.seq.section import SEQ_INDEX_NAMES

rng = np.random.default_rng(0)


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


class TestCalcUniqReadWeights(ut.TestCase):

    def test_0_reads(self):
        for ncls in range(4):
            read_weights = np.ones((0, ncls), dtype=float)
            uniq_inverse = np.array([], dtype=int)
            num_uniq = 0
            expect = np.ones((0, ncls), dtype=float)
            result = _calc_uniq_read_weights(read_weights,
                                             uniq_inverse,
                                             num_uniq)
            self.assertTrue(np.array_equal(result, expect))

    def test_1_read(self):
        for ncls in range(4):
            read_weights = np.ones((1, ncls), dtype=float)
            uniq_inverse = np.array([0])
            num_uniq = 1
            expect = np.ones((1, ncls), dtype=float)
            result = _calc_uniq_read_weights(read_weights,
                                             uniq_inverse,
                                             num_uniq)
            self.assertTrue(np.array_equal(result, expect))

    def test_2_reads_1_uniq(self):
        read_weights = np.array([[0.1, 0.2],
                                 [0.3, 0.4]])
        uniq_inverse = np.array([0, 0])
        num_uniq = 1
        expect = np.array([[0.4, 0.6]])
        result = _calc_uniq_read_weights(read_weights,
                                         uniq_inverse,
                                         num_uniq)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))

    def test_2_reads_2_uniq(self):
        read_weights = np.array([[0.1, 0.2],
                                 [0.3, 0.4]])
        uniq_inverse = np.array([1, 0])
        num_uniq = 2
        expect = np.array([[0.3, 0.4],
                           [0.1, 0.2]])
        result = _calc_uniq_read_weights(read_weights,
                                         uniq_inverse,
                                         num_uniq)
        self.assertEqual(result.shape, expect.shape)
        self.assertTrue(np.allclose(result, expect))


class TestCalcCoverage(ut.TestCase):

    def test_0_positions(self):
        pos_index = pd.MultiIndex.from_tuples([], names=SEQ_INDEX_NAMES)
        end5s = np.array([[]], dtype=int)
        end3s = np.array([[]], dtype=int)
        read_nums = np.arange(0)
        exp_per_pos = pd.Series(0., index=pos_index)
        exp_per_read = {
            "A": pd.Series([], read_nums, dtype=int),
            "C": pd.Series([], read_nums, dtype=int),
            "G": pd.Series([], read_nums, dtype=int),
            "T": pd.Series([], read_nums, dtype=int),
            "N": pd.Series([], read_nums, dtype=int),
        }
        res_per_pos, res_per_read = calc_coverage(pos_index,
                                                  read_nums,
                                                  end5s,
                                                  end3s)
        self.assertIsInstance(res_per_pos, pd.Series)
        self.assertTrue(res_per_pos.equals(exp_per_pos))
        self.assertEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))

    def test_0_positions_weighted(self):
        pos_index = pd.MultiIndex.from_tuples([], names=SEQ_INDEX_NAMES)
        end5s = np.array([[]], dtype=int)
        end3s = np.array([[]], dtype=int)
        read_nums = np.arange(0)
        read_weights = pd.DataFrame.from_dict({
            "C1": pd.Series(np.array([])),
            "C2": pd.Series(np.array([])),
        })
        exp_per_pos = pd.DataFrame.from_dict({
            "C1": pd.Series(np.array([]), index=pos_index),
            "C2": pd.Series(np.array([]), index=pos_index),
        })
        exp_per_read = {
            "A": pd.Series([], read_nums, dtype=int),
            "C": pd.Series([], read_nums, dtype=int),
            "G": pd.Series([], read_nums, dtype=int),
            "T": pd.Series([], read_nums, dtype=int),
            "N": pd.Series([], read_nums, dtype=int),
        }
        res_per_pos, res_per_read = calc_coverage(pos_index,
                                                  read_nums,
                                                  end5s,
                                                  end3s,
                                                  read_weights)
        self.assertIsInstance(res_per_pos, pd.DataFrame)
        self.assertTrue(res_per_pos.equals(exp_per_pos))
        self.assertEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))

    def test_1_position(self):
        pos_index = pd.MultiIndex.from_tuples([(3, "G")],
                                              names=SEQ_INDEX_NAMES)
        end5s = np.array([[3],
                          [3],
                          [3],
                          [3]])
        end3s = np.array([[3],
                          [3],
                          [3],
                          [3]])
        read_nums = np.arange(4)
        exp_per_pos = pd.Series([4.], pos_index)
        exp_per_read = {
            "A": pd.Series([0, 0, 0, 0], read_nums),
            "C": pd.Series([0, 0, 0, 0], read_nums),
            "G": pd.Series([1, 1, 1, 1], read_nums),
            "T": pd.Series([0, 0, 0, 0], read_nums),
            "N": pd.Series([0, 0, 0, 0], read_nums)
        }
        res_per_pos, res_per_read = calc_coverage(pos_index,
                                                  read_nums,
                                                  end5s,
                                                  end3s)
        self.assertIsInstance(res_per_pos, pd.Series)
        self.assertTrue(res_per_pos.equals(exp_per_pos))
        self.assertEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))

    def test_1_segment(self):
        """
        1234567890123

        -
        --
        ---
         ---
           -
           ---
        -------
        -------------
              -------
               ---
                 -
                 ---
                  ---
                   --
                    -
        """
        pos_index = pd.MultiIndex.from_tuples([(3, "G"),
                                               (4, "A"),
                                               (5, "C"),
                                               (6, "A"),
                                               (8, "T"),
                                               (9, "T"),
                                               (10, "G"),
                                               (11, "C")],
                                              names=SEQ_INDEX_NAMES)
        end5s = np.array([[1],
                          [1],
                          [1],
                          [2],
                          [4],
                          [4],
                          [1],
                          [1],
                          [7],
                          [8],
                          [10],
                          [10],
                          [11],
                          [12],
                          [13]])
        end3s = np.array([[1],
                          [2],
                          [3],
                          [4],
                          [4],
                          [6],
                          [7],
                          [13],
                          [13],
                          [10],
                          [10],
                          [12],
                          [13],
                          [13],
                          [13]])
        read_nums = np.arange(15)
        exp_per_pos = pd.Series([4., 5., 3., 3., 3., 3., 5., 4.],
                                index=pos_index)
        exp_per_read = {
            "A": pd.Series([0, 0, 0, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0],
                           read_nums),
            "C": pd.Series([0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 1, 1, 0, 0],
                           read_nums),
            "G": pd.Series([0, 0, 1, 1, 0, 0, 1, 2, 1, 1, 1, 1, 0, 0, 0],
                           read_nums),
            "T": pd.Series([0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0],
                           read_nums),
            "N": pd.Series([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           read_nums)
        }
        res_per_pos, res_per_read = calc_coverage(pos_index,
                                                  read_nums,
                                                  end5s,
                                                  end3s)
        self.assertIsInstance(res_per_pos, pd.Series)
        self.assertTrue(res_per_pos.equals(exp_per_pos))
        self.assertEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))

    def test_1_segment_weighted(self):
        """
        1234567890123   C1  C2

        -               0.1 1.5
        --              0.2 1.4
        ---             0.3 1.3
         ---            0.4 1.2
           -            0.5 1.1
           ---          0.6 1.0
        -------         0.7 0.9
        -------------   0.8 0.8
              -------   0.9 0.7
               ---      1.0 0.6
                 -      1.1 0.5
                 ---    1.2 0.4
                  ---   1.3 0.3
                   --   1.4 0.2
                    -   1.5 0.1
        """
        pos_index = pd.MultiIndex.from_tuples([(3, "G"),
                                               (4, "A"),
                                               (5, "C"),
                                               (6, "A"),
                                               (8, "T"),
                                               (9, "T"),
                                               (10, "G"),
                                               (11, "C")],
                                              names=SEQ_INDEX_NAMES)
        end5s = np.array([[1],
                          [1],
                          [1],
                          [2],
                          [4],
                          [4],
                          [1],
                          [1],
                          [7],
                          [8],
                          [10],
                          [10],
                          [11],
                          [12],
                          [13]])
        end3s = np.array([[1],
                          [2],
                          [3],
                          [4],
                          [4],
                          [6],
                          [7],
                          [13],
                          [13],
                          [10],
                          [10],
                          [12],
                          [13],
                          [13],
                          [13]])
        read_nums = np.arange(15)
        read_weights = pd.DataFrame.from_dict({
            "C1": pd.Series(np.linspace(0.1, 1.5, 15)),
            "C2": pd.Series(np.linspace(1.5, 0.1, 15))
        })
        exp_per_pos = pd.DataFrame.from_dict({
            "C1": pd.Series([2.2, 3.0, 2.1, 2.1, 2.7, 2.7, 5.0, 4.2],
                            index=pos_index),
            "C2": pd.Series([4.2, 5.0, 2.7, 2.7, 2.1, 2.1, 3.0, 2.2],
                            index=pos_index),
        })
        exp_per_read = {
            "A": pd.Series([0, 0, 0, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0],
                           read_nums),
            "C": pd.Series([0, 0, 0, 0, 0, 1, 1, 2, 1, 0, 0, 1, 1, 0, 0],
                           read_nums),
            "G": pd.Series([0, 0, 1, 1, 0, 0, 1, 2, 1, 1, 1, 1, 0, 0, 0],
                           read_nums),
            "T": pd.Series([0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 0, 0, 0, 0],
                           read_nums),
            "N": pd.Series([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           read_nums)
        }
        res_per_pos, res_per_read = calc_coverage(pos_index,
                                                  read_nums,
                                                  end5s,
                                                  end3s,
                                                  read_weights)
        self.assertIsInstance(res_per_pos, pd.DataFrame)
        self.assertTrue(res_per_pos.round(6).equals(exp_per_pos.round(6)))
        self.assertEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))

    def test_1_segment_redundant(self):
        """
        1234567890123

        -------
           -------
              -------
           -------
        -------
        """
        pos_index = pd.MultiIndex.from_tuples([(3, "G"),
                                               (4, "A"),
                                               (5, "C"),
                                               (6, "A"),
                                               (8, "T"),
                                               (9, "T"),
                                               (10, "G"),
                                               (11, "C")],
                                              names=SEQ_INDEX_NAMES)
        end5s = np.array([[1],
                          [4],
                          [7],
                          [4],
                          [1]])
        end3s = np.array([[7],
                          [10],
                          [13],
                          [10],
                          [7]])
        read_nums = np.arange(5)
        exp_per_pos = pd.Series([2., 4., 4., 4., 3., 3., 3., 1.],
                                index=pos_index)
        exp_per_read = {
            "A": pd.Series([2, 2, 0, 2, 2], read_nums),
            "C": pd.Series([1, 1, 1, 1, 1], read_nums),
            "G": pd.Series([1, 1, 1, 1, 1], read_nums),
            "T": pd.Series([0, 2, 2, 2, 0], read_nums),
            "N": pd.Series([0, 0, 0, 0, 0], read_nums)
        }
        res_per_pos, res_per_read = calc_coverage(pos_index,
                                                  read_nums,
                                                  end5s,
                                                  end3s)
        self.assertIsInstance(res_per_pos, pd.Series)
        self.assertTrue(res_per_pos.equals(exp_per_pos))
        self.assertEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))

    def test_1_segment_redundant_weighted(self):
        """
        1234567890123   C1  C2

        -------         0.6 0.2
           -------      0.5 0.3
              -------   0.9 0.8
           -------      0.4 1.0
        -------         0.1 0.7
        """
        pos_index = pd.MultiIndex.from_tuples([(3, "G"),
                                               (4, "A"),
                                               (5, "C"),
                                               (6, "A"),
                                               (8, "T"),
                                               (9, "T"),
                                               (10, "G"),
                                               (11, "C")],
                                              names=SEQ_INDEX_NAMES)
        end5s = np.array([[1],
                          [4],
                          [7],
                          [4],
                          [1]])
        end3s = np.array([[7],
                          [10],
                          [13],
                          [10],
                          [7]])
        read_nums = np.arange(5)
        read_weights = pd.DataFrame.from_dict({
            "C1": pd.Series([0.6, 0.5, 0.9, 0.4, 0.1]),
            "C2": pd.Series([0.2, 0.3, 0.8, 1.0, 0.7])
        })
        exp_per_pos = pd.DataFrame.from_dict({
            "C1": pd.Series([0.7, 1.6, 1.6, 1.6, 1.8, 1.8, 1.8, 0.9],
                            index=pos_index),
            "C2": pd.Series([0.9, 2.2, 2.2, 2.2, 2.1, 2.1, 2.1, 0.8],
                            index=pos_index),
        })
        exp_per_read = {
            "A": pd.Series([2, 2, 0, 2, 2], read_nums),
            "C": pd.Series([1, 1, 1, 1, 1], read_nums),
            "G": pd.Series([1, 1, 1, 1, 1], read_nums),
            "T": pd.Series([0, 2, 2, 2, 0], read_nums),
            "N": pd.Series([0, 0, 0, 0, 0], read_nums)
        }
        res_per_pos, res_per_read = calc_coverage(pos_index,
                                                  read_nums,
                                                  end5s,
                                                  end3s,
                                                  read_weights)
        self.assertIsInstance(res_per_pos, pd.DataFrame)
        self.assertTrue(res_per_pos.round(6).equals(exp_per_pos.round(6)))
        self.assertEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))

    def test_2_segments_nomask(self):
        """
        1234567890123

        -----
                =====

           -----
             =====

          ----
          ====

                -
                =

            ----
                 ====

            -----
                 ====

            ------
                 ====

        -------------
        =============
        """
        pos_index = pd.MultiIndex.from_tuples([(3, "G"),
                                               (4, "A"),
                                               (5, "C"),
                                               (6, "A"),
                                               (8, "T"),
                                               (9, "T"),
                                               (10, "G"),
                                               (11, "C")],
                                              names=SEQ_INDEX_NAMES)
        end5s = np.array([[1, 9],
                          [9, 1],
                          [4, 6],
                          [6, 4],
                          [3, 3],
                          [9, 9],
                          [5, 10],
                          [10, 5],
                          [5, 10],
                          [10, 5],
                          [5, 10],
                          [10, 5],
                          [1, 1]])
        end3s = np.array([[5, 13],
                          [13, 5],
                          [8, 10],
                          [10, 8],
                          [6, 6],
                          [9, 9],
                          [8, 13],
                          [13, 8],
                          [9, 13],
                          [13, 9],
                          [10, 13],
                          [13, 10],
                          [13, 13]])
        read_nums = np.arange(13)
        exp_per_pos = pd.Series([4., 6., 12., 10., 9., 10., 11., 9.],
                                index=pos_index)
        exp_per_read = {
            "A": pd.Series([1, 1, 2, 2, 2, 0, 1, 1, 1, 1, 1, 1, 2],
                           read_nums),
            "C": pd.Series([2, 2, 1, 1, 1, 0, 2, 2, 2, 2, 2, 2, 2],
                           read_nums),
            "G": pd.Series([2, 2, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 2],
                           read_nums),
            "T": pd.Series([1, 1, 2, 2, 0, 1, 1, 1, 2, 2, 2, 2, 2],
                           read_nums),
            "N": pd.Series([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           read_nums)
        }
        res_per_pos, res_per_read = calc_coverage(pos_index,
                                                  read_nums,
                                                  end5s,
                                                  end3s)
        self.assertIsInstance(res_per_pos, pd.Series)
        self.assertTrue(res_per_pos.equals(exp_per_pos))
        self.assertEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))

    def test_2_segments_mask(self):
        """
        01234567890123

         [---]
        ][

         [---]
         ][

         [---]
            ][

         [---]
             ][

         [---]
              [===]

        ][
         [===]

         ][
         [===]

            ][
         [===]

             ][
         [===]
        """
        pos_index = pd.MultiIndex.from_tuples([(3, "G"),
                                               (4, "A"),
                                               (5, "C"),
                                               (6, "A"),
                                               (8, "T"),
                                               (9, "T"),
                                               (10, "G"),
                                               (11, "C")],
                                              names=SEQ_INDEX_NAMES)
        mask = np.array([[False, True],
                         [False, True],
                         [False, True],
                         [False, True],
                         [False, False],
                         [True, False],
                         [True, False],
                         [True, False],
                         [True, False]])
        end5s = np.ma.masked_array([[1, 1],
                                    [1, 2],
                                    [1, 5],
                                    [1, 6],
                                    [1, 6],
                                    [1, 1],
                                    [2, 1],
                                    [5, 1],
                                    [6, 1]],
                                   mask)
        end3s = np.ma.masked_array([[5, 0],
                                    [5, 1],
                                    [5, 4],
                                    [5, 5],
                                    [5, 10],
                                    [0, 5],
                                    [1, 5],
                                    [4, 5],
                                    [5, 5]],
                                   mask)
        read_nums = np.arange(9)
        read_weights = pd.DataFrame.from_dict({
            "C1": pd.Series([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
            "C2": pd.Series([0.2, 0.3, 0.8, 0.9, 0.7, 0.5, 0.6, 0.1, 0.4])
        })
        exp_per_pos = pd.DataFrame.from_dict({
            "C1": pd.Series([9.0, 9.0, 9.0, 1.0, 1.0, 1.0, 1.0, 0.0],
                            pos_index),
            "C2": pd.Series([4.5, 4.5, 4.5, 0.7, 0.7, 0.7, 0.7, 0.0],
                            pos_index),
        })
        exp_per_read = {
            "A": pd.Series([1, 1, 1, 1, 2, 1, 1, 1, 1],
                           read_nums),
            "C": pd.Series([1, 1, 1, 1, 1, 1, 1, 1, 1],
                           read_nums),
            "G": pd.Series([1, 1, 1, 1, 2, 1, 1, 1, 1],
                           read_nums),
            "T": pd.Series([0, 0, 0, 0, 2, 0, 0, 0, 0],
                           read_nums),
            "N": pd.Series([0, 0, 0, 0, 0, 0, 0, 0, 0],
                           read_nums)
        }
        res_per_pos, res_per_read = calc_coverage(pos_index,
                                                  read_nums,
                                                  end5s,
                                                  end3s,
                                                  read_weights)
        self.assertIsInstance(res_per_pos, pd.DataFrame)
        self.assertTrue(res_per_pos.round(6).equals(exp_per_pos.round(6)))
        self.assertEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))


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
