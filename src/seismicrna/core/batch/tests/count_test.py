import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.batch.count import (calc_coverage,
                                         _calc_uniq_read_weights,
                                         count_end_coords,
                                         calc_rels_per_pos,
                                         calc_rels_per_read,
                                         calc_count_per_pos,
                                         calc_count_per_read)
from seismicrna.core.batch.ends import END5_COORD, END3_COORD
from seismicrna.core.batch.read import calc_inverse
from seismicrna.core.rel import HalfRelPattern, RelPattern
from seismicrna.core.seq.section import SEQ_INDEX_NAMES, seq_pos_to_index
from seismicrna.core.seq.xna import DNA

rng = np.random.default_rng(0)


class TestCountEndCoords(ut.TestCase):

    def test_length_0(self):
        end5s = np.array([], dtype=int)
        end3s = np.array([], dtype=int)
        result = count_end_coords(end5s, end3s)
        expect = pd.Series(np.array([], dtype=int),
                           pd.MultiIndex.from_arrays([np.array([], dtype=int),
                                                      np.array([], dtype=int)],
                                                     names=[END5_COORD,
                                                            END3_COORD]))
        self.assertIsInstance(result, pd.Series)
        self.assertTrue(result.equals(expect))

    def test_length_1(self):
        end5s = np.array([4])
        end3s = np.array([6])
        result = count_end_coords(end5s, end3s)
        expect = pd.Series(np.array([1]),
                           pd.MultiIndex.from_arrays([np.array([4]),
                                                      np.array([6])],
                                                     names=[END5_COORD,
                                                            END3_COORD]))
        self.assertIsInstance(result, pd.Series)
        self.assertTrue(result.equals(expect))

    def test_length_2_diff(self):
        end5s = np.array([4, 1])
        end3s = np.array([6, 8])
        result = count_end_coords(end5s, end3s)
        expect = pd.Series(np.array([1, 1]),
                           pd.MultiIndex.from_arrays([np.array([1, 4]),
                                                      np.array([8, 6])],
                                                     names=[END5_COORD,
                                                            END3_COORD]))
        self.assertIsInstance(result, pd.Series)
        self.assertTrue(result.equals(expect))

    def test_length_2_same(self):
        end5s = np.array([4, 4])
        end3s = np.array([6, 6])
        result = count_end_coords(end5s, end3s)
        expect = pd.Series(np.array([2]),
                           pd.MultiIndex.from_arrays([np.array([4]),
                                                      np.array([6])],
                                                     names=[END5_COORD,
                                                            END3_COORD]))
        self.assertIsInstance(result, pd.Series)
        self.assertTrue(result.equals(expect))

    def test_length_3_diff_weights(self):
        end5s = np.array([4, 1, 0])
        end3s = np.array([6, 8, 3])
        weights = pd.DataFrame(np.array([[2., 3.],
                                         [4., 5.],
                                         [6., 7.]]),
                               columns=["a", "b"])
        result = count_end_coords(end5s, end3s, weights)
        expect = pd.DataFrame(weights.values[[2, 1, 0]],
                              pd.MultiIndex.from_arrays([np.array([0, 1, 4]),
                                                         np.array([3, 8, 6])],
                                                        names=[END5_COORD,
                                                               END3_COORD]),
                              ["a", "b"])
        self.assertIsInstance(result, pd.DataFrame)
        self.assertTrue(result.equals(expect))

    def test_length_3_same_weights(self):
        end5s = np.array([4, 1, 4])
        end3s = np.array([6, 8, 6])
        weights = pd.DataFrame(np.array([[2., 3.],
                                         [4., 5.],
                                         [6., 7.]]),
                               columns=["a", "b"])
        result = count_end_coords(end5s, end3s, weights)
        expect = pd.DataFrame(np.array([[4., 5.],
                                        [8., 10.]]),
                              pd.MultiIndex.from_arrays([np.array([1, 4]),
                                                         np.array([8, 6])],
                                                        names=[END5_COORD,
                                                               END3_COORD]),
                              ["a", "b"])
        self.assertIsInstance(result, pd.DataFrame)
        self.assertTrue(result.equals(expect))

    def test_length_3_same_weights_int(self):
        end5s = np.array([4, 1, 4])
        end3s = np.array([6, 8, 6])
        weights = pd.DataFrame(np.array([[2, 3],
                                         [4, 5],
                                         [6, 7]]),
                               columns=["a", "b"])
        result = count_end_coords(end5s, end3s, weights)
        expect = pd.DataFrame(np.array([[4, 5],
                                        [8, 10]]),
                              pd.MultiIndex.from_arrays([np.array([1, 4]),
                                                         np.array([8, 6])],
                                                        names=[END5_COORD,
                                                               END3_COORD]),
                              ["a", "b"])
        self.assertIsInstance(result, pd.DataFrame)
        self.assertTrue(result.equals(expect))


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
        self.assertListEqual(sorted(res_per_read), sorted(exp_per_read))
        for base in exp_per_read:
            self.assertIsInstance(res_per_read[base], pd.Series)
            self.assertTrue(res_per_read[base].equals(exp_per_read[base]))


class TestCalcRelsPerPos(ut.TestCase):

    def test_average(self):
        """
              11  13  14  15
        Read   A   C   G   T
        --------------------
           0   1  64   1 255
           1 255 255   1 255
           5 128   1   1  16
           6 255   1   1  17
           7   1  64   1 255
           9 128   1   1   1
          10   1   1   1 255
          11   1   1   1   1
          14  64   3   1   1
          16 255 255   1   1
        """
        positions = seq_pos_to_index(DNA("ANCGT"), [11, 13, 14, 15], 11)
        mutations = {11: {64: np.array([14]),
                          128: np.array([5, 9])},
                     13: {3: np.array([14]),
                          64: np.array([0, 7])},
                     14: {},
                     15: {16: np.array([5]),
                          17: np.array([6])}}
        num_reads = 10
        cover_per_pos = pd.Series([7, 8, 10, 6], positions)
        expect = {1: [4, 5, 10, 4],
                  3: [0, 1, 0, 0],
                  16: [0, 0, 0, 1],
                  17: [0, 0, 0, 1],
                  64: [1, 2, 0, 0],
                  128: [2, 0, 0, 0],
                  255: [3, 2, 0, 4]}
        rels_per_pos = calc_rels_per_pos(mutations,
                                         num_reads,
                                         cover_per_pos)
        self.assertIsInstance(rels_per_pos, dict)
        self.assertSetEqual(set(rels_per_pos), set(expect))
        for rel, rexp in expect.items():
            rres = rels_per_pos[rel]
            self.assertIsInstance(rres, pd.Series)
            self.assertTrue(rres.equals(pd.Series(rexp, positions)))

    def test_clusters(self):
        """
              11  13  14  15
        Read   A   C   G   T
        --------------------
           0   1  64   1 255
           1 255 255   1 255
           5 128   1   1  16
           6 255   1   1  17
           7   1  64   1 255
           9 128   1   1   1
          10   1   1   1 255
          11   1   1   1   1
          14  64   3   1   1
          16 255 255   1   1
        """
        positions = seq_pos_to_index(DNA("ANCGT"), [11, 13, 14, 15], 11)
        mutations = {11: {64: np.array([14]),
                          128: np.array([5, 9])},
                     13: {3: np.array([14]),
                          64: np.array([0, 7])},
                     14: {},
                     15: {16: np.array([5]),
                          17: np.array([6])}}
        cover_per_pos = pd.DataFrame([[3.9, 3.1],
                                      [4.3, 3.7],
                                      [5.5, 4.5],
                                      [4.0, 2.0]],
                                     positions,
                                     ["a", "b"])
        read_nums = np.array([0, 1, 5, 6, 7, 9, 10, 11, 14, 16])
        read_indexes = calc_inverse(read_nums)
        read_weights = pd.DataFrame([[0.1, 0.9],
                                     [0.2, 0.8],
                                     [0.3, 0.7],
                                     [0.4, 0.6],
                                     [0.5, 0.5],
                                     [0.6, 0.4],
                                     [0.7, 0.3],
                                     [0.8, 0.2],
                                     [0.9, 0.1],
                                     [1.0, 0.0]],
                                    index=read_nums,
                                    columns=["a", "b"])
        num_reads = read_weights.sum(axis=0)
        expect = {1: [[2.1, 2.8, 5.5, 3.3],
                      [1.9, 2.2, 4.5, 0.7]],
                  3: [[0.0, 0.9, 0.0, 0.0],
                      [0.0, 0.1, 0.0, 0.0]],
                  16: [[0.0, 0.0, 0.0, 0.3],
                       [0.0, 0.0, 0.0, 0.7]],
                  17: [[0.0, 0.0, 0.0, 0.4],
                       [0.0, 0.0, 0.0, 0.6]],
                  64: [[0.9, 0.6, 0.0, 0.0],
                       [0.1, 1.4, 0.0, 0.0]],
                  128: [[0.9, 0.0, 0.0, 0.0],
                        [1.1, 0.0, 0.0, 0.0]],
                  255: [[1.6, 1.2, 0.0, 1.5],
                        [1.4, 0.8, 0.0, 2.5]]}
        rels_per_pos = calc_rels_per_pos(mutations,
                                         num_reads,
                                         cover_per_pos,
                                         read_indexes,
                                         read_weights)
        self.assertIsInstance(rels_per_pos, dict)
        self.assertSetEqual(set(rels_per_pos), set(expect))
        for rel, rexp in expect.items():
            rres = rels_per_pos[rel]
            self.assertIsInstance(rres, pd.DataFrame)
            self.assertTrue(rres.equals(
                pd.DataFrame(rexp, ["a", "b"], positions).T
            ))


class TestCalcRelsPerRead(ut.TestCase):

    def test_average(self):
        """
              11  13  14  15  17
        Read   A   C   G   T   C
        ------------------------
           0   1  64   1 255 255
           1 255 255   1 255 255
           5 128   1   1  16   1
           6 255   1   1  17  64
           7   1  64   1 255 255
           9 128   1   1   1  16
          10   1   1   1 255 255
          11   1   1   1   1   1
          14  64   3   1   1 128
          16 255 255   1   1   2
        """
        positions = seq_pos_to_index(DNA("ANCGTNC"), [11, 13, 14, 15, 17], 11)
        mutations = {11: {64: np.array([14]),
                          128: np.array([5, 9])},
                     13: {3: np.array([14]),
                          64: np.array([0, 7])},
                     14: {},
                     15: {16: np.array([5]),
                          17: np.array([6])},
                     17: {2: np.array([16]),
                          16: np.array([9]),
                          64: np.array([6]),
                          128: np.array([14])}}
        read_nums = np.array([0, 1, 5, 6, 7, 9, 10, 11, 14, 16])
        read_indexes = calc_inverse(read_nums)
        cover_per_read = pd.DataFrame.from_dict({
            "A": pd.Series([1, 0, 1, 0, 1, 1, 1, 1, 1, 0], read_nums),
            "C": pd.Series([1, 0, 2, 2, 1, 2, 1, 2, 2, 1], read_nums),
            "G": pd.Series([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], read_nums),
            "T": pd.Series([0, 0, 1, 1, 0, 1, 0, 1, 1, 1], read_nums),
        })
        expect = {1: [[1, 0, 1, 0],
                      [0, 0, 1, 0],
                      [0, 2, 1, 0],
                      [0, 1, 1, 0],
                      [1, 0, 1, 0],
                      [0, 1, 1, 1],
                      [1, 1, 1, 0],
                      [1, 2, 1, 1],
                      [0, 0, 1, 1],
                      [0, 0, 1, 1]],
                  2: [[0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 1, 0, 0]],
                  3: [[0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 0, 0]],
                  16: [[0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0]],
                  17: [[0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 1],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0]],
                  64: [[0, 1, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 1, 0, 0],
                       [0, 1, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [1, 0, 0, 0],
                       [0, 0, 0, 0]],
                  128: [[0, 0, 0, 0],
                        [0, 0, 0, 0],
                        [1, 0, 0, 0],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0],
                        [1, 0, 0, 0],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, 0, 0]],
                  255: [[0, 1, 0, 1],
                        [1, 2, 0, 1],
                        [0, 0, 0, 0],
                        [1, 0, 0, 0],
                        [0, 1, 0, 1],
                        [0, 0, 0, 0],
                        [0, 1, 0, 1],
                        [0, 0, 0, 0],
                        [0, 0, 0, 0],
                        [1, 1, 0, 0]]}
        rels_per_pos = calc_rels_per_read(mutations,
                                          positions,
                                          cover_per_read,
                                          read_indexes)
        self.assertIsInstance(rels_per_pos, dict)
        self.assertSetEqual(set(rels_per_pos), set(expect))
        for rel, rexp in expect.items():
            rres = rels_per_pos[rel]
            self.assertIsInstance(rres, pd.DataFrame)
            self.assertTrue(rres.equals(
                pd.DataFrame(rexp, read_nums, ["A", "C", "G", "T"])
            ))


class TestCalcCountPerPos(ut.TestCase):

    def test_average(self):
        """
              11  13  14  15
        Read   A   C   G   T
        --------------------
           0   1  64   1 255
           1 255 255   1 255
           5 128   1   1  16
           6 255   1   1  17
           7   1  64   1 255
           9 128   1   1   1
          10   1   1   1 255
          11   1   1   1   1
          14  64   3   1   1
          16 255 255   1   1
        """
        positions = seq_pos_to_index(DNA("ANCGT"), [11, 13, 14, 15], 11)
        mutations = {11: {64: np.array([14]),
                          128: np.array([5, 9])},
                     13: {3: np.array([14]),
                          64: np.array([0, 7])},
                     14: {},
                     15: {16: np.array([5]),
                          17: np.array([6])}}
        num_reads = 10
        cover_per_pos = pd.Series([7, 8, 10, 6], positions)
        rels_per_pos = calc_rels_per_pos(mutations,
                                         num_reads,
                                         cover_per_pos)
        pattern = RelPattern.muts()
        iexp = pd.Series([7, 7, 10, 5], positions)
        fexp = pd.Series([3, 2, 0, 1], positions)
        info, fits = calc_count_per_pos(pattern, cover_per_pos, rels_per_pos)
        self.assertIsInstance(info, pd.Series)
        self.assertTrue(info.equals(iexp))
        self.assertIsInstance(fits, pd.Series)
        self.assertTrue(fits.equals(fexp))
        pattern = RelPattern.muts().invert()
        iexp = pd.Series([7, 7, 10, 5], positions)
        fexp = pd.Series([4, 5, 10, 4], positions)
        info, fits = calc_count_per_pos(pattern, cover_per_pos, rels_per_pos)
        self.assertIsInstance(info, pd.Series)
        self.assertTrue(info.equals(iexp))
        self.assertIsInstance(fits, pd.Series)
        self.assertTrue(fits.equals(fexp))

    def test_clusters(self):
        """
              11  13  14  15
        Read   A   C   G   T
        --------------------
           0   1  64   1 255
           1 255 255   1 255
           5 128   1   1  16
           6 255   1   1  17
           7   1  64   1 255
           9 128   1   1   1
          10   1   1   1 255
          11   1   1   1   1
          14  64   3   1   1
          16 255 255   1   1
        """
        positions = seq_pos_to_index(DNA("ANCGT"), [11, 13, 14, 15], 11)
        mutations = {11: {64: np.array([14]),
                          128: np.array([5, 9])},
                     13: {3: np.array([14]),
                          64: np.array([0, 7])},
                     14: {},
                     15: {16: np.array([5]),
                          17: np.array([6])}}
        cover_per_pos = pd.DataFrame([[3.9, 3.1],
                                      [4.3, 3.7],
                                      [5.5, 4.5],
                                      [4.0, 2.0]],
                                     positions,
                                     ["a", "b"])
        read_nums = np.array([0, 1, 5, 6, 7, 9, 10, 11, 14, 16])
        read_indexes = calc_inverse(read_nums)
        read_weights = pd.DataFrame([[0.1, 0.9],
                                     [0.2, 0.8],
                                     [0.3, 0.7],
                                     [0.4, 0.6],
                                     [0.5, 0.5],
                                     [0.6, 0.4],
                                     [0.7, 0.3],
                                     [0.8, 0.2],
                                     [0.9, 0.1],
                                     [1.0, 0.0]],
                                    index=read_nums,
                                    columns=["a", "b"])
        num_reads = read_weights.sum(axis=0)
        rels_per_pos = calc_rels_per_pos(mutations,
                                         num_reads,
                                         cover_per_pos,
                                         read_indexes,
                                         read_weights)
        pattern = RelPattern.muts()
        iexp = pd.DataFrame([[3.9, 3.1],
                             [3.4, 3.6],
                             [5.5, 4.5],
                             [3.6, 1.4]],
                            positions,
                            ["a", "b"])
        fexp = pd.DataFrame([[1.8, 1.2],
                             [0.6, 1.4],
                             [0.0, 0.0],
                             [0.3, 0.7]],
                            positions,
                            ["a", "b"])
        info, fits = calc_count_per_pos(pattern, cover_per_pos, rels_per_pos)
        self.assertIsInstance(info, pd.DataFrame)
        self.assertTrue(info.index.equals(iexp.index))
        self.assertTrue(info.columns.equals(iexp.columns))
        self.assertTrue(np.all(np.isclose(info.values, iexp.values)))
        self.assertIsInstance(fits, pd.DataFrame)
        self.assertTrue(fits.index.equals(fexp.index))
        self.assertTrue(fits.columns.equals(fexp.columns))
        self.assertTrue(np.all(np.isclose(fits.values, fexp.values)))
        pattern = RelPattern.muts().invert()
        iexp = pd.DataFrame([[3.9, 3.1],
                             [3.4, 3.6],
                             [5.5, 4.5],
                             [3.6, 1.4]],
                            positions,
                            ["a", "b"])
        fexp = pd.DataFrame([[2.1, 1.9],
                             [2.8, 2.2],
                             [5.5, 4.5],
                             [3.3, 0.7]],
                            positions,
                            ["a", "b"])
        info, fits = calc_count_per_pos(pattern, cover_per_pos, rels_per_pos)
        self.assertIsInstance(info, pd.DataFrame)
        self.assertTrue(info.index.equals(iexp.index))
        self.assertTrue(info.columns.equals(iexp.columns))
        self.assertTrue(np.all(np.isclose(info.values, iexp.values)))
        self.assertIsInstance(fits, pd.DataFrame)
        self.assertTrue(fits.index.equals(fexp.index))
        self.assertTrue(fits.columns.equals(fexp.columns))
        self.assertTrue(np.all(np.isclose(fits.values, fexp.values)))


class TestCalcCountPerRead(ut.TestCase):

    def test_average(self):
        """
              11  13  14  15  17
        Read   A   C   G   T   C
        ------------------------
           0   1  64   1 255 255
           1 255 255   1 255 255
           5 128   1   1  16   1
           6 255   1   1  17  64
           7   1  64   1 255 255
           9 128   1   1   1  16
          10   1   1   1 255 255
          11   1   1   1   1   1
          14  64   3   1   1 128
          16 255 255   1   1   2
        """
        positions = seq_pos_to_index(DNA("ANCGTNC"), [11, 13, 14, 15, 17], 11)
        mutations = {11: {64: np.array([14]),
                          128: np.array([5, 9])},
                     13: {3: np.array([14]),
                          64: np.array([0, 7])},
                     14: {},
                     15: {16: np.array([5]),
                          17: np.array([6])},
                     17: {2: np.array([16]),
                          16: np.array([9]),
                          64: np.array([6]),
                          128: np.array([14])}}
        read_nums = np.array([0, 1, 5, 6, 7, 9, 10, 11, 14, 16])
        read_indexes = calc_inverse(read_nums)
        cover_per_read = pd.DataFrame.from_dict({
            "A": pd.Series([1, 0, 1, 0, 1, 1, 1, 1, 1, 0], read_nums),
            "C": pd.Series([1, 0, 2, 2, 1, 2, 1, 2, 2, 1], read_nums),
            "G": pd.Series([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], read_nums),
            "T": pd.Series([0, 0, 1, 1, 0, 1, 0, 1, 1, 1], read_nums),
        })
        rels_per_read = calc_rels_per_read(mutations,
                                           positions,
                                           cover_per_read,
                                           read_indexes)
        pattern = RelPattern.muts()
        iexp = pd.Series([3, 1, 5, 3, 3, 5, 3, 5, 4, 3], read_nums)
        fexp = pd.Series([1, 0, 2, 1, 1, 2, 0, 0, 2, 1], read_nums)
        info, fits = calc_count_per_read(pattern,
                                         cover_per_read,
                                         rels_per_read)
        self.assertIsInstance(info, pd.Series)
        self.assertTrue(info.equals(iexp))
        self.assertIsInstance(fits, pd.Series)
        self.assertTrue(fits.equals(fexp))
        pattern = RelPattern(HalfRelPattern.from_counts(count_sub=True,
                                                        discount=["at", "cg"]),
                             HalfRelPattern.refs())
        iexp = pd.Series([2, 1, 4, 2, 2, 4, 3, 5, 4, 2], read_nums)
        fexp = pd.Series([0, 0, 1, 0, 0, 1, 0, 0, 2, 0], read_nums)
        info, fits = calc_count_per_read(pattern,
                                         cover_per_read,
                                         rels_per_read)
        self.assertIsInstance(info, pd.Series)
        self.assertTrue(info.equals(iexp))
        self.assertIsInstance(fits, pd.Series)
        self.assertTrue(fits.equals(fexp))


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
