import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.cluster.batch import ClusterMutsBatch
from seismicrna.core.batch.accum import accumulate
from seismicrna.core.batch.ends import END5_COORD, END3_COORD
from seismicrna.core.batch.index import RB_INDEX_NAMES
from seismicrna.core.header import RelClustHeader
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
        patterns = {"Matches": RelPattern.muts().invert(),
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
        patterns = {"Matches": RelPattern.muts().invert(),
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

    def test_mask_2_batches(self):
        """
        . = Match
        ! = Mutation
        ? = Ambiguous
        _ = Not Covered

          3468
          ACGT
        ------
        2 _!?!
        4 .!..
        7 !!!_
        ------
        0 !!._
        6 ?.._
        """
        section = Section("myref", DNA("NNACNGNT"))
        section.add_mask("mask", [3, 4, 6, 8], complement=True)
        patterns = {"Matches": RelPattern.muts().invert(),
                    "Mutations": RelPattern.muts()}
        batches = [
            MaskMutsBatch(section=section,
                          batch=0,
                          muts={3: {128: np.array([7])},
                                4: {16: np.array([2, 4, 7])},
                                6: {32: np.array([7]),
                                    33: np.array([2])},
                                8: {64: np.array([2])}},
                          seg_end5s=np.array([[4],
                                              [3],
                                              [2]]),
                          seg_end3s=np.array([[8],
                                              [8],
                                              [7]]),
                          read_nums=np.array([2, 4, 7])),
            MaskMutsBatch(section=section,
                          batch=1,
                          muts={3: {128: np.array([0]),
                                    129: np.array([6])},
                                4: {16: np.array([0])},
                                6: {},
                                8: {}},
                          seg_end5s=np.array([[2],
                                              [3]]),
                          seg_end3s=np.array([[7],
                                              [6]]),
                          read_nums=np.array([0, 6]))
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
            section.unmasked,
            list(patterns)
        )))
        self.assertTrue(ipp.equals(pd.DataFrame(
            [[3, 3],
             [5, 5],
             [4, 4],
             [2, 2]],
            section.unmasked,
            list(patterns)
        )))
        self.assertTrue(fpr.equals(pd.DataFrame(
            [[0, 2],
             [3, 1],
             [0, 3],
             [1, 2],
             [2, 0]],
            pd.MultiIndex.from_arrays([[0, 0, 0, 1, 1],
                                       [2, 4, 7, 0, 6]],
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
                                       [2, 4, 7, 0, 6]],
                                      names=RB_INDEX_NAMES),
            list(patterns)
        )))
        self.assertTrue(ends.equals(pd.Series(
            [2, 1, 1, 1],
            pd.MultiIndex.from_tuples([(2, 7),
                                       (3, 6),
                                       (3, 8),
                                       (4, 8)],
                                      names=[END5_COORD, END3_COORD])
        )))

    def test_cluster_2_batches(self):
        """
        . = Match
        ! = Mutation
        ? = Ambiguous
        _ = Not Covered

          3468
          ACGT  C1  C2
        --------------
        2 _!?! 0.1 0.9
        4 .!.. 0.3 0.7
        7 !!!_ 0.5 0.5
        --------------
        0 !!._ 0.6 0.4
        6 ?.._ 0.8 0.2
        """
        section = Section("myref", DNA("NNACNGNT"))
        section.add_mask("mask", [3, 4, 6, 8], complement=True)
        patterns = {"Matches": RelPattern.muts().invert(),
                    "Mutations": RelPattern.muts()}
        ks = [1, 2]
        rcheader = RelClustHeader(rels=list(patterns), ks=ks)
        rheader = rcheader.get_rel_header()
        cheader = rcheader.get_clust_header()
        batches = [
            ClusterMutsBatch(section=section,
                             batch=0,
                             muts={3: {128: np.array([7])},
                                   4: {16: np.array([2, 4, 7])},
                                   6: {32: np.array([7]),
                                       33: np.array([2])},
                                   8: {64: np.array([2])}},
                             seg_end5s=np.array([[4],
                                                 [3],
                                                 [2]]),
                             seg_end3s=np.array([[8],
                                                 [8],
                                                 [7]]),
                             resps=pd.DataFrame([[1.0, 0.1, 0.9],
                                                 [1.0, 0.3, 0.7],
                                                 [1.0, 0.5, 0.5]],
                                                [2, 4, 7],
                                                cheader.index)),
            ClusterMutsBatch(section=section,
                             batch=1,
                             muts={3: {128: np.array([0]),
                                       129: np.array([6])},
                                   4: {16: np.array([0])},
                                   6: {},
                                   8: {}},
                             seg_end5s=np.array([[2],
                                                 [3]]),
                             seg_end3s=np.array([[7],
                                                 [6]]),
                             resps=pd.DataFrame([[1.0, 0.6, 0.4],
                                                 [1.0, 0.8, 0.2]],
                                                [0, 6],
                                                cheader.index))
        ]
        n, (fpp, ipp), (fpr, ipr), ends = accumulate(
            batches,
            section.seq,
            patterns,
            ks=ks,
            pos_nums=section.unmasked_int
        )
        self.assertIsInstance(n, pd.Series)
        self.assertTrue(n.index.equals(cheader.index))
        self.assertTrue(np.allclose(n, [5.0, 2.3, 2.7]))
        self.assertIsInstance(fpp, pd.DataFrame)
        self.assertTrue(fpp.index.equals(section.unmasked))
        self.assertTrue(fpp.columns.equals(rcheader.index))
        self.assertTrue(np.allclose(fpp.values,
                                    [[1.0, 0.3, 0.7, 2.0, 1.1, 0.9],
                                     [1.0, 0.8, 0.2, 4.0, 1.5, 2.5],
                                     [3.0, 1.7, 1.3, 1.0, 0.5, 0.5],
                                     [1.0, 0.3, 0.7, 1.0, 0.1, 0.9]]))
        self.assertIsInstance(ipp, pd.DataFrame)
        self.assertTrue(ipp.index.equals(section.unmasked))
        self.assertTrue(ipp.columns.equals(rcheader.index))
        self.assertTrue(np.allclose(ipp.values,
                                    [[3.0, 1.4, 1.6, 3.0, 1.4, 1.6],
                                     [5.0, 2.3, 2.7, 5.0, 2.3, 2.7],
                                     [4.0, 2.2, 1.8, 4.0, 2.2, 1.8],
                                     [2.0, 0.4, 1.6, 2.0, 0.4, 1.6]]))
        read_nums = pd.MultiIndex.from_arrays([[0, 0, 0, 1, 1],
                                               [2, 4, 7, 0, 6]],
                                              names=RB_INDEX_NAMES)
        self.assertIsInstance(fpr, pd.DataFrame)
        self.assertTrue(fpr.equals(pd.DataFrame(
            [[0, 2],
             [3, 1],
             [0, 3],
             [1, 2],
             [2, 0]],
            read_nums,
            rheader.index
        )))
        self.assertIsInstance(ipr, pd.DataFrame)
        self.assertTrue(ipr.equals(pd.DataFrame(
            [[2, 2],
             [4, 4],
             [3, 3],
             [3, 3],
             [2, 2]],
            read_nums,
            rheader.index
        )))
        self.assertIsInstance(ends, pd.DataFrame)
        self.assertTrue(ends.index.equals(
            pd.MultiIndex.from_tuples([(2, 7),
                                       (3, 6),
                                       (3, 8),
                                       (4, 8)],
                                      names=[END5_COORD, END3_COORD])
        ))
        self.assertTrue(ends.columns.equals(cheader.index))
        self.assertTrue(np.allclose(ends.values,
                                    [[2.0, 1.1, 0.9],
                                     [1.0, 0.8, 0.2],
                                     [1.0, 0.3, 0.7],
                                     [1.0, 0.1, 0.9]]))


if __name__ == "__main__":
    ut.main()
