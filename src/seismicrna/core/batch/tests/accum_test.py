import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.cluster.batch import ClusterMutsBatch
from seismicrna.core.batch.accum import accumulate_batches
from seismicrna.core.batch.ends import END5_COORD, END3_COORD
from seismicrna.core.batch.index import RB_INDEX_NAMES
from seismicrna.core.batch.muts import RegionMutsBatch
from seismicrna.core.header import RelClustHeader
from seismicrna.core.rel import RelPattern
from seismicrna.core.seq.region import Region
from seismicrna.core.seq.xna import DNA
from seismicrna.mask.batch import MaskMutsBatch
from seismicrna.relate.batch import RelateRegionMutsBatch


def get_batch_count_all_func(batches: list[RegionMutsBatch]):
    def get_batch_count_all(batch_num: int, **kwargs):
        return batches[batch_num].count_all(**kwargs)

    return get_batch_count_all


class TestAccumulateBatches(ut.TestCase):

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
        region = Region("myref", DNA("ACGT"))
        patterns = {"Matches": RelPattern.muts().invert(),
                    "Mutations": RelPattern.muts()}
        batches = [
            RelateRegionMutsBatch(region=region,
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
        n, ends, fpp, fpr = accumulate_batches(
            get_batch_count_all_func(batches),
            len(batches),
            region.seq,
            region.unmasked_int,
            patterns,
        )
        self.assertEqual(n, 5)
        self.assertTrue(fpp.equals(pd.DataFrame(
            [[1, 2],
             [1, 4],
             [3, 1],
             [1, 1]],
            region.range,
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
        region = Region("myref", DNA("ACGT"))
        patterns = {"Matches": RelPattern.muts().invert(),
                    "Mutations": RelPattern.muts()}
        batches = [
            RelateRegionMutsBatch(region=region,
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
            RelateRegionMutsBatch(region=region,
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
        n, ends, fpp, fpr = accumulate_batches(
            get_batch_count_all_func(batches),
            len(batches),
            region.seq,
            region.unmasked_int,
            patterns,
        )
        self.assertEqual(n, 5)
        self.assertTrue(fpp.equals(pd.DataFrame(
            [[1, 2],
             [1, 4],
             [3, 1],
             [1, 1]],
            region.range,
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
        region = Region("myref", DNA("NNACNGNT"))
        region.add_mask("mask", [3, 4, 6, 8], complement=True)
        patterns = {"Matches": RelPattern.muts().invert(),
                    "Mutations": RelPattern.muts()}
        batches = [
            MaskMutsBatch(region=region,
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
            MaskMutsBatch(region=region,
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
        n, ends, fpp, fpr = accumulate_batches(
            get_batch_count_all_func(batches),
            len(batches),
            region.seq,
            region.unmasked_int,
            patterns,
        )
        self.assertEqual(n, 5)
        self.assertTrue(fpp.equals(pd.DataFrame(
            [[1, 2],
             [1, 4],
             [3, 1],
             [1, 1]],
            region.unmasked,
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
        region = Region("myref", DNA("NNACNGNT"))
        region.add_mask("mask", [3, 4, 6, 8], complement=True)
        patterns = {"Matches": RelPattern.muts().invert(),
                    "Mutations": RelPattern.muts()}
        ks = [1, 2]
        rcheader = RelClustHeader(rels=list(patterns), ks=ks)
        rheader = rcheader.get_rel_header()
        cheader = rcheader.get_clust_header()
        batches = [
            ClusterMutsBatch(region=region,
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
            ClusterMutsBatch(region=region,
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
        n, ends, fpp, fpr = accumulate_batches(
            get_batch_count_all_func(batches),
            len(batches),
            region.seq,
            region.unmasked_int,
            patterns,
            ks,
        )
        self.assertIsInstance(n, pd.Series)
        self.assertTrue(n.index.equals(cheader.index))
        self.assertTrue(np.allclose(n, [5.0, 2.3, 2.7]))
        self.assertIsInstance(fpp, pd.DataFrame)
        self.assertTrue(fpp.index.equals(region.unmasked))
        self.assertTrue(fpp.columns.equals(rcheader.index))
        self.assertTrue(np.allclose(fpp.values,
                                    [[1.0, 0.3, 0.7, 2.0, 1.1, 0.9],
                                     [1.0, 0.8, 0.2, 4.0, 1.5, 2.5],
                                     [3.0, 1.7, 1.3, 1.0, 0.5, 0.5],
                                     [1.0, 0.3, 0.7, 1.0, 0.1, 0.9]]))
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
    ut.main(verbosity=2)
