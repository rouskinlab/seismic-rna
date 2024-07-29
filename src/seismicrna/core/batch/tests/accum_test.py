import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.cluster.batch import ClusterMutsBatch
from seismicrna.core.batch.accum import accumulate
from seismicrna.core.batch.index import RB_INDEX_NAMES
from seismicrna.core.batch.ends import END5_COORD, END3_COORD
from seismicrna.core.rel import RelPattern
from seismicrna.core.seq.section import Section
from seismicrna.core.seq.xna import DNA
from seismicrna.mask.batch import MaskMutsBatch
from seismicrna.relate.batch import RelateBatch


class TestAccumulate(ut.TestCase):

    def test_relate_1_batch(self):
        """
        . = Match
        ! = Mutation
        ? = Ambiguous
        _ = Not Covered

          1234
          ACGT
        ------
        0 _!?!
        1 .!..
        2 !!!_
        3 !!._
        4 ?.._
        """
        section = Section("myref", DNA("ACGT"))
        patterns = {"Matches": RelPattern.muts().intersect(None, invert=True),
                    "Mutations": RelPattern.muts()}
        batches = [
            RelateBatch(section=section,
                        batch=0,
                        muts={1: {128: np.array([2, 3]),
                                  129: np.array([4])},
                              2: {16: np.array([0, 1, 2, 3])},
                              3: {32: np.array([2]),
                                  33: np.array([0])},
                              4: {64: np.array([0])}},
                        seg_end5s=np.array([[2],
                                            [1],
                                            [1],
                                            [1],
                                            [1]]),
                        seg_end3s=np.array([[4],
                                            [4],
                                            [3],
                                            [3],
                                            [3]]))
        ]
        n, (fpp, ipp), (fpr, ipr), ends = accumulate(
            batches,
            section.seq,
            patterns,
            pos_nums=section.unmasked_int
        )
        self.assertEqual(n, 5)
        self.assertTrue(fpp.equals(pd.DataFrame(
            [[1, 2],
             [1, 4],
             [3, 1],
             [1, 1]],
            section.range,
            list(patterns)
        )))
        self.assertTrue(ipp.equals(pd.DataFrame(
            [[3, 3],
             [5, 5],
             [4, 4],
             [2, 2]],
            section.range,
            list(patterns)
        )))
        self.assertTrue(fpr.equals(pd.DataFrame(
            [[0, 2],
             [3, 1],
             [0, 3],
             [1, 2],
             [2, 0]],
            pd.MultiIndex.from_arrays([[0, 0, 0, 0, 0],
                                       [0, 1, 2, 3, 4]],
                                      names=RB_INDEX_NAMES),
            list(patterns)
        )))
        self.assertTrue(ipr.equals(pd.DataFrame(
            [[2, 2],
             [4, 4],
             [3, 3],
             [3, 3],
             [2, 2]],
            pd.MultiIndex.from_arrays([[0, 0, 0, 0, 0],
                                       [0, 1, 2, 3, 4]],
                                      names=RB_INDEX_NAMES),
            list(patterns)
        )))
        self.assertTrue(ends.equals(pd.Series(
            [3, 1, 1],
            pd.MultiIndex.from_tuples([(1, 3),
                                       (1, 4),
                                       (2, 4)],
                                      names=[END5_COORD, END3_COORD])
        )))

    def test_relate_2_batches(self):
        """
        . = Match
        ! = Mutation
        ? = Ambiguous
        _ = Not Covered

          1234
          ACGT
        ------
        0 _!?!
        1 .!..
        2 !!!_
        ------
        0 !!._
        1 ?.._
        """
        section = Section("myref", DNA("ACGT"))
        patterns = {"Matches": RelPattern.muts().intersect(None, invert=True),
                    "Mutations": RelPattern.muts()}
        batches = [
            RelateBatch(section=section,
                        batch=0,
                        muts={1: {128: np.array([2])},
                              2: {16: np.array([0, 1, 2])},
                              3: {32: np.array([2]),
                                  33: np.array([0])},
                              4: {64: np.array([0])}},
                        seg_end5s=np.array([[2],
                                            [1],
                                            [1]]),
                        seg_end3s=np.array([[4],
                                            [4],
                                            [3]])),
            RelateBatch(section=section,
                        batch=1,
                        muts={1: {128: np.array([0]),
                                  129: np.array([1])},
                              2: {16: np.array([0])},
                              3: {},
                              4: {}},
                        seg_end5s=np.array([[1],
                                            [1]]),
                        seg_end3s=np.array([[3],
                                            [3]]))
        ]
        n, (fpp, ipp), (fpr, ipr), ends = accumulate(
            batches,
            section.seq,
            patterns,
            pos_nums=section.unmasked_int
        )
        self.assertEqual(n, 5)
        self.assertTrue(fpp.equals(pd.DataFrame(
            [[1, 2],
             [1, 4],
             [3, 1],
             [1, 1]],
            section.range,
            list(patterns)
        )))
        self.assertTrue(ipp.equals(pd.DataFrame(
            [[3, 3],
             [5, 5],
             [4, 4],
             [2, 2]],
            section.range,
            list(patterns)
        )))
        self.assertTrue(fpr.equals(pd.DataFrame(
            [[0, 2],
             [3, 1],
             [0, 3],
             [1, 2],
             [2, 0]],
            pd.MultiIndex.from_arrays([[0, 0, 0, 1, 1],
                                       [0, 1, 2, 0, 1]],
                                      names=RB_INDEX_NAMES),
            list(patterns)
        )))
        self.assertTrue(ipr.equals(pd.DataFrame(
            [[2, 2],
             [4, 4],
             [3, 3],
             [3, 3],
             [2, 2]],
            pd.MultiIndex.from_arrays([[0, 0, 0, 1, 1],
                                       [0, 1, 2, 0, 1]],
                                      names=RB_INDEX_NAMES),
            list(patterns)
        )))
        self.assertTrue(ends.equals(pd.Series(
            [3, 1, 1],
            pd.MultiIndex.from_tuples([(1, 3),
                                       (1, 4),
                                       (2, 4)],
                                      names=[END5_COORD, END3_COORD])
        )))


if __name__ == "__main__":
    ut.main()
