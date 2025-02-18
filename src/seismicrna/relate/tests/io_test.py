import unittest as ut

import numpy as np

from seismicrna.core.seq import DNA
from seismicrna.relate.io import from_reads


class TestFromReads(ut.TestCase):

    def test_from_0_reads(self):
        refseq = DNA("ACGT")
        reads = []
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True)
        self.assertDictEqual(
            batch.muts,
            {pos: dict() for pos in range(1, len(refseq) + 1)}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.zeros((0, 0), dtype=int)
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.zeros((0, 0), dtype=int)
        ))
        self.assertIsNone(batch.seg_ends_mask)
        self.assertTrue(np.array_equal(
            names.names,
            np.array([], dtype=str)
        ))

    def test_from_1_read_0_segs_drop_empty(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([], []), {}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=True)
        self.assertDictEqual(
            batch.muts,
            {pos: dict() for pos in range(1, len(refseq) + 1)}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.zeros((0, 0), dtype=int)
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.zeros((0, 0), dtype=int)
        ))
        self.assertIsNone(batch.seg_ends_mask)
        self.assertTrue(np.array_equal(
            names.names,
            np.array([], dtype=str)
        ))

    def test_from_1_read_0_segs_keep_empty(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([], []), {}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=False)
        self.assertDictEqual(
            batch.muts,
            {pos: dict() for pos in range(1, len(refseq) + 1)}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.zeros((1, 0), dtype=int)
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.zeros((1, 0), dtype=int)
        ))
        self.assertIsNone(batch.seg_ends_mask)
        self.assertTrue(np.array_equal(
            names.names,
            np.array(["Read1"])
        ))

    def test_from_1_read_1_segs_no_cover_drop_empty(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([3], [2]), {}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=True)
        self.assertDictEqual(
            batch.muts,
            {pos: dict() for pos in range(1, len(refseq) + 1)}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.zeros((0, 1), dtype=int)
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.zeros((0, 1), dtype=int)
        ))
        self.assertIsNone(batch.seg_ends_mask)
        self.assertTrue(np.array_equal(
            names.names,
            np.array([], dtype=str)
        ))

    def test_from_1_read_1_segs_no_cover_keep_empty(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([3], [2]), {}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=False)
        self.assertDictEqual(
            batch.muts,
            {pos: dict() for pos in range(1, len(refseq) + 1)}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.array([[3]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.array([[2]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_ends_mask,
            np.array([[True]])
        ))
        self.assertTrue(np.array_equal(
            names.names,
            np.array(["Read1"])
        ))

    def test_from_1_read_1_segs_no_muts(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([2], [2]), {}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=True)
        self.assertDictEqual(
            batch.muts,
            {pos: dict() for pos in range(1, len(refseq) + 1)}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.array([[2]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.array([[2]])
        ))
        self.assertIsNone(batch.seg_ends_mask)
        self.assertTrue(np.array_equal(
            names.names,
            np.array(["Read1"])
        ))

    def test_from_2_reads_1_segs(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([2], [2]), {2: 64})),
                 ("Read2", (([1], [4]), {1: 128, 3: 2}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=True)
        self.assertDictEqual(
            {pos: list(rels) for pos, rels in batch.muts.items()},
            {1: [128], 2: [64], 3: [2], 4: []}
        )
        self.assertTrue(np.array_equal(
            batch.muts[1][128],
            np.array([1])
        ))
        self.assertTrue(np.array_equal(
            batch.muts[2][64],
            np.array([0])
        ))
        self.assertTrue(np.array_equal(
            batch.muts[3][2],
            np.array([1])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.array([[2],
                      [1]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.array([[2],
                      [4]])
        ))
        self.assertIsNone(batch.seg_ends_mask)
        self.assertTrue(np.array_equal(
            names.names,
            np.array(["Read1", "Read2"])
        ))

    def test_from_2_reads_2_segs(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([2, 4], [2, 4]), {2: 64})),
                 ("Read2", (([1, 1], [4, 3]), {1: 128, 3: 2}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=True)
        self.assertDictEqual(
            {pos: {rel: reads.tolist() for rel, reads in rels.items()}
             for pos, rels in batch.muts.items()},
            {1: {128: [1]}, 2: {64: [0]}, 3: {2: [1]}, 4: {}}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.array([[2, 4],
                      [1, 1]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.array([[2, 4],
                      [4, 3]])
        ))
        self.assertIsNone(batch.seg_ends_mask)
        self.assertTrue(np.array_equal(
            names.names,
            np.array(["Read1", "Read2"])
        ))

    def test_from_2_reads_1_2_segs(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([2], [2]), {2: 64})),
                 ("Read2", (([1, 1], [4, 3]), {1: 128, 3: 2}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=True)
        self.assertDictEqual(
            {pos: {rel: reads.tolist() for rel, reads in rels.items()}
             for pos, rels in batch.muts.items()},
            {1: {128: [1]}, 2: {64: [0]}, 3: {2: [1]}, 4: {}}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.array([[2, 1],
                      [1, 1]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.array([[2, 0],
                      [4, 3]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_ends_mask,
            np.array([[False, True],
                      [False, False]])
        ))
        self.assertTrue(np.array_equal(
            names.names,
            np.array(["Read1", "Read2"])
        ))

    def test_from_2_reads_2_1_segs(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([2, 4], [2, 4]), {2: 64})),
                 ("Read2", (([1], [4]), {1: 128, 3: 2}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=True)
        self.assertDictEqual(
            {pos: {rel: reads.tolist() for rel, reads in rels.items()}
             for pos, rels in batch.muts.items()},
            {1: {128: [1]}, 2: {64: [0]}, 3: {2: [1]}, 4: {}}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.array([[2, 4],
                      [1, 1]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.array([[2, 4],
                      [4, 0]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_ends_mask,
            np.array([[False, False],
                      [False, True]])
        ))
        self.assertTrue(np.array_equal(
            names.names,
            np.array(["Read1", "Read2"])
        ))

    def test_from_4_reads_varied_segs_drop_empty(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([2, 4], [1, 3]), {})),
                 ("Read2", (([1, 3], [2, 4]), {1: 128, 2: 64})),
                 ("Read3", (([], []), {})),
                 ("Read4", (([1, 2, 3], [2, 3, 4]), {1: 128, 2: 16, 3: 2}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=True)
        self.assertDictEqual(
            {pos: {rel: reads.tolist() for rel, reads in rels.items()}
             for pos, rels in batch.muts.items()},
            {1: {128: [0, 1]}, 2: {16: [1], 64: [0]}, 3: {2: [1]}, 4: {}}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.array([[1, 3, 1],
                      [1, 2, 3]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.array([[2, 4, 0],
                      [2, 3, 4]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_ends_mask,
            np.array([[False, False, True],
                      [False, False, False]])
        ))
        self.assertTrue(np.array_equal(
            names.names,
            np.array(["Read2", "Read4"])
        ))

    def test_from_4_reads_varied_segs_keep_empty(self):
        refseq = DNA("ACGT")
        reads = [("Read1", (([2, 4], [1, 3]), {})),
                 ("Read2", (([1, 3], [2, 4]), {1: 128, 2: 64})),
                 ("Read3", (([], []), {})),
                 ("Read4", (([1, 2, 3], [2, 3, 4]), {1: 128, 2: 16, 3: 2}))]
        batch, names = from_reads(reads,
                                  sample="mysample",
                                  branches={},
                                  ref="myref",
                                  refseq=refseq,
                                  batch=0,
                                  write_read_names=True,
                                  drop_empty_reads=False)
        self.assertDictEqual(
            {pos: {rel: reads.tolist() for rel, reads in rels.items()}
             for pos, rels in batch.muts.items()},
            {1: {128: [1, 3]}, 2: {16: [3], 64: [1]}, 3: {2: [3]}, 4: {}}
        )
        self.assertTrue(np.array_equal(
            batch.seg_end5s,
            np.array([[2, 4, 1],
                      [1, 3, 1],
                      [1, 1, 1],
                      [1, 2, 3]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_end3s,
            np.array([[1, 3, 0],
                      [2, 4, 0],
                      [0, 0, 0],
                      [2, 3, 4]])
        ))
        self.assertTrue(np.array_equal(
            batch.seg_ends_mask,
            np.array([[True, True, True],
                      [False, False, True],
                      [True, True, True],
                      [False, False, False]])
        ))
        self.assertTrue(np.array_equal(
            names.names,
            np.array(["Read1", "Read2", "Read3", "Read4"])
        ))


if __name__ == "__main__":
    ut.main(verbosity=2)
