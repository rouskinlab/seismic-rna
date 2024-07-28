import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.cluster.batch import ClusterMutsBatch
from seismicrna.core.batch.accum import accumulate
from seismicrna.core.rel import RelPattern
from seismicrna.core.seq.section import Section
from seismicrna.core.seq.xna import DNA
from seismicrna.mask.batch import MaskMutsBatch
from seismicrna.relate.batch import RelateBatch


class TestAccumulate(ut.TestCase):

    def test_relate(self):
        """
        . = Match
        ! = Mutation

          1234
          ACGT
        0  !!!
        1 .!..
        2 !!!
        3 !!.
        4 !..
        """
        section = Section("myref", DNA("ACGT"))
        patterns = {"Matches": RelPattern.muts().intersect(None, invert=True),
                    "Mutations": RelPattern.muts()}
        batches = [
            RelateBatch(section=section,
                        batch=0,
                        muts={1: {128: np.array([2, 3, 4])},
                              2: {16: np.array([0, 1, 2, 3])},
                              3: {32: np.array([0, 2])},
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
        self.assertIsInstance(fpp, pd.DataFrame)
        self.assertTrue(fpp.equals(pd.DataFrame(
            [[1, 3],
             [1, 4],
             [3, 2],
             [1, 1]],
            section.range,
            list(patterns)
        )))


if __name__ == "__main__":
    ut.main()
