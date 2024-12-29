import unittest as ut
from abc import ABC, abstractmethod
from pathlib import Path
from shutil import rmtree
from typing import Iterable

from seismicrna.core import path
from seismicrna.core.batch.muts import RegionMutsBatch
from seismicrna.core.io.seq import RefseqIO
from seismicrna.core.logs import Level, set_config
from seismicrna.core.ngs.xam import SAM_DELIM
from seismicrna.core.seq.fasta import write_fasta
from seismicrna.core.seq.xna import DNA
from seismicrna.relate.data import RelateDataset
from seismicrna.relate.main import run

SAMPLE = "sample"
REF = "ref"
REF_SEQ = DNA("CGGCATATC")
REFS = "refs"

SAM_DATA_EMPTY = [
    ("@SQ", f"SN:{REF}", f"LN:{len(REF_SEQ)}"),
]

# Read CGGCATATC Qual
# -------------------
#    1 CG-C         6
#    2    CA--TC    5
#    3   gCaT       9
#    4    CtTt      7
# (Lower quality bases in lowercase.)
SAM_DATA_SINGLE = [
    ["@SQ", f"SN:{REF}", f"LN:{len(REF_SEQ)}"],
    ["Read1", 0, REF, 1, 6, "2M1D1M", "*", 0, 0, "CGC", "III"],
    ["Read2", 0, REF, 4, 5, "2M2D2M", "*", 0, 0, "CATC", "IIII"],
    ["Read3", 0, REF, 3, 9, "4M", "*", 0, 0, "GCAT", "HIHI"],
    ["Read4", 0, REF, 4, 7, "4M", "*", 0, 0, "CTTG", "IHIH"],
]

# Read CGGCATATC Qual
# -------------------
#  1R1 CG-C         6
#  1R2  gtCATG      6
#  2R2     At-TC    7
#  3R1    CAT--C    5
#  3R2   GCat       5
#  4R1  GGCTT       9
#  4R2    CaTATC    9
#  5R1 cGGC         7
#  6R1      TATG    8
#  6R2 cGGC         8
# (Lower quality bases in lowercase.)
SAM_DATA_PAIRED = [
    ["@SQ", f"SN:{REF}", f"LN:{len(REF_SEQ)}"],
    ["Read1", 83, REF, 1, 6, "2M1D1M", REF, 2, 7, "CGC", "III"],
    ["Read1", 163, REF, 2, 6, "6M", REF, 1, 7, "GTCATG", "HHIIII"],
    ["Read2", 145, REF, 5, 7, "2M1D2M", "*", 0, 5, "ATTC", "IHII"],
    ["Read3", 83, REF, 4, 5, "3M2D1M", REF, 3, 7, "CATC", "IIII"],
    ["Read3", 163, REF, 3, 5, "4M", REF, 4, 7, "GCAT", "IIHH"],
    ["Read4", 99, REF, 2, 9, "5M", REF, 4, 8, "GGCTT", "IIIII"],
    ["Read4", 147, REF, 4, 9, "6M", REF, 2, 8, "CATATC", "IHIIII"],
    ["Read5", 65, REF, 1, 7, "4M", "*", 0, 4, "CGGC", "HIII"],
    ["Read6", 99, REF, 6, 8, "4M", REF, 1, 9, "TATG", "IIII"],
    ["Read6", 147, REF, 1, 8, "4M", REF, 6, 9, "CGGC", "HIII"],
]


def write_sam_file(out_dir: Path, data: list[list]):
    sam_dir = out_dir.joinpath(SAMPLE, path.CMD_ALIGN_DIR)
    sam_dir.mkdir(parents=True)
    sam_file = sam_dir.joinpath(f"{REF}.sam")
    with open(sam_file, "x") as f:
        f.write("\n".join(SAM_DELIM.join(map(str, line)) for line in data))
    return sam_file


def write_fasta_file(out_dir: Path):
    fasta_file = out_dir.joinpath(f"{REFS}.fa")
    write_fasta(fasta_file, [(REF, REF_SEQ)])
    return fasta_file


def extract_batches(batches: Iterable[RegionMutsBatch]):
    batches = list(batches)
    read_nums = [list(batch.read_nums) for batch in batches]
    seg_end5s = [batch.seg_end5s.tolist() for batch in batches]
    seg_end3s = [batch.seg_end3s.tolist() for batch in batches]
    muts = [{pos: {mut: reads.tolist() for mut, reads in muts.items()}
             for pos, muts in batch.muts.items()}
            for batch in batches]
    return read_nums, seg_end5s, seg_end3s, muts


def load_refseq(out_dir: Path):
    return RefseqIO.load(RefseqIO.build_path(top=out_dir,
                                             sample=SAMPLE,
                                             cmd=path.CMD_REL_DIR,
                                             ref=REF),
                         checksum="").refseq


class TestRelate(ut.TestCase, ABC):

    @classmethod
    @abstractmethod
    def get_sam_data(cls) -> list:
        """ Data for the SAM file. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._out_dir: Path | None = None
        self._fasta_file: Path | None = None
        self._sam_file: Path | None = None

    def setUp(self):
        self._out_dir = path.randdir(prefix="out-")
        self._fasta_file = write_fasta_file(self._out_dir)
        self._sam_file = write_sam_file(self._out_dir, self.get_sam_data())
        set_config(verbosity=Level.SEVERE, raise_on_error=True)

    def tearDown(self):
        rmtree(self._out_dir)
        self._out_dir = None
        self._fasta_file = None
        self._sam_file = None
        set_config()

    def batches(self,
                min_reads: int = 0,
                min_mapq: int = 0,
                min_phred: int = 0,
                clip_end5: int = 0,
                clip_end3: int = 0,
                **kwargs):
        relate_dir, = run(str(self._fasta_file),
                          (str(self._sam_file),),
                          out_dir=str(self._out_dir),
                          min_reads=min_reads,
                          min_mapq=min_mapq,
                          min_phred=min_phred,
                          clip_end5=clip_end5,
                          clip_end3=clip_end3,
                          **kwargs)
        relate_report_file = relate_dir.joinpath("relate-report.json")
        return extract_batches(RelateDataset(relate_report_file).iter_batches())


class TestRelateEmpty(TestRelate):

    @classmethod
    def get_sam_data(cls):
        return SAM_DATA_EMPTY

    def test_noargs(self):
        read_nums, _, _, _ = self.batches()
        self.assertEqual(len(read_nums), 0)
        self.assertEqual(load_refseq(self._out_dir), REF_SEQ)

    def test_min_reads(self):
        self.assertRaisesRegex(ValueError,
                               "Insufficient reads in alignment map",
                               self.batches,
                               min_reads=1)


class TestRelateSingle(TestRelate):

    @classmethod
    def get_sam_data(cls):
        return SAM_DATA_SINGLE

    def test_noargs(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches()
        self.assertEqual(load_refseq(self._out_dir), REF_SEQ)
        self.assertListEqual(read_nums, [[0, 1, 2, 3]])
        self.assertListEqual(seg_end5s, [[[1], [4], [3], [4]]])
        self.assertListEqual(seg_end3s, [[[4], [9], [6], [7]]])
        self.assertListEqual(muts,
                             [{1: {},
                               2: {3: [0]},
                               3: {3: [0]},
                               4: {},
                               5: {3: [1], 128: [3]},
                               6: {3: [1]},
                               7: {3: [1], 64: [3]},
                               8: {3: [1]},
                               9: {}}])

    def test_min_reads(self):
        self.assertRaisesRegex(ValueError,
                               "Insufficient reads in alignment map",
                               self.batches,
                               min_reads=5)

    def test_min_phred(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(min_phred=40)
        self.assertListEqual(read_nums, [[0, 1, 2, 3]])
        self.assertListEqual(seg_end5s, [[[1], [4], [3], [4]]])
        self.assertListEqual(seg_end3s, [[[4], [9], [6], [7]]])
        self.assertListEqual(muts,
                             [{1: {},
                               2: {3: [0]},
                               3: {3: [0], 177: [2]},
                               4: {},
                               5: {3: [1], 225: [2, 3]},
                               6: {3: [1]},
                               7: {3: [1], 225: [3]},
                               8: {3: [1]},
                               9: {}}])

    def test_clip(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(clip_end5=1,
                                                             clip_end3=1)
        self.assertListEqual(read_nums, [[0, 1, 2, 3]])
        self.assertListEqual(seg_end5s, [[[2], [5], [4], [5]]])
        self.assertListEqual(seg_end3s, [[[3], [8], [5], [6]]])
        self.assertListEqual(muts,
                             [{1: {},
                               2: {3: [0]},
                               3: {3: [0]},
                               4: {},
                               5: {3: [1], 128: [3]},
                               6: {3: [1]},
                               7: {3: [1]},
                               8: {3: [1]},
                               9: {}}])

    def test_batch_size_1(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(batch_size=1)
        self.assertListEqual(read_nums, [[0], [0], [0], [0]])
        self.assertListEqual(seg_end5s, [[[1]], [[4]], [[3]], [[4]]])
        self.assertListEqual(seg_end3s, [[[4]], [[9]], [[6]], [[7]]])
        self.assertListEqual(muts,
                             [{1: {},
                               2: {3: [0]},
                               3: {3: [0]},
                               4: {},
                               5: {},
                               6: {},
                               7: {},
                               8: {},
                               9: {}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {3: [0]},
                               6: {3: [0]},
                               7: {3: [0]},
                               8: {3: [0]},
                               9: {}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {},
                               6: {},
                               7: {},
                               8: {},
                               9: {}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {128: [0]},
                               6: {},
                               7: {64: [0]},
                               8: {},
                               9: {}}])

    def test_batch_size_3(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(batch_size=3)
        self.assertListEqual(read_nums, [[0, 1, 2], [0]])
        self.assertListEqual(seg_end5s, [[[1], [4], [3]], [[4]]])
        self.assertListEqual(seg_end3s, [[[4], [9], [6]], [[7]]])
        self.assertListEqual(muts,
                             [{1: {},
                               2: {3: [0]},
                               3: {3: [0]},
                               4: {},
                               5: {3: [1]},
                               6: {3: [1]},
                               7: {3: [1]},
                               8: {3: [1]},
                               9: {}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {128: [0]},
                               6: {},
                               7: {64: [0]},
                               8: {},
                               9: {}}])


class TestRelatePaired(TestRelate):

    @classmethod
    def get_sam_data(cls):
        return SAM_DATA_PAIRED

    def test_noargs(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches()
        self.assertEqual(load_refseq(self._out_dir), REF_SEQ)
        self.assertListEqual(read_nums, [[0, 1, 2, 3, 4, 5]])
        self.assertListEqual(seg_end5s,
                             [[[2, 1], [3, 4], [2, 4], [6, 1], [5, 5], [1, 1]]])
        self.assertListEqual(seg_end3s,
                             [[[7, 4], [6, 9], [6, 9], [9, 4], [9, 9], [4, 4]]])
        self.assertListEqual(muts,
                             [{1: {},
                               2: {},
                               3: {0: [0]},
                               4: {},
                               5: {0: [2]},
                               6: {},
                               7: {2: [4], 3: [1], 64: [0]},
                               8: {3: [1]},
                               9: {64: [3]}}])

    def test_min_reads(self):
        self.assertRaisesRegex(ValueError,
                               "Insufficient reads in alignment map",
                               self.batches,
                               min_reads=7)

    def test_min_phred(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(min_phred=40)
        self.assertListEqual(read_nums, [[0, 1, 2, 3, 4, 5]])
        self.assertListEqual(seg_end5s,
                             [[[2, 1], [3, 4], [2, 4], [6, 1], [5, 5], [1, 1]]])
        self.assertListEqual(seg_end3s,
                             [[[7, 4], [6, 9], [6, 9], [9, 4], [9, 9], [4, 4]]])
        self.assertListEqual(muts,
                             [{1: {209: [3, 5]},
                               2: {},
                               3: {},
                               4: {},
                               5: {128: [2]},
                               6: {115: [4]},
                               7: {3: [1], 64: [0], 227: [4]},
                               8: {3: [1]},
                               9: {64: [3]}}])

    def test_clip(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(clip_end5=1,
                                                             clip_end3=1)
        self.assertListEqual(read_nums, [[0, 1, 2, 3, 4, 5]])
        self.assertListEqual(seg_end5s,
                             [[[3, 2], [4, 5], [3, 5], [7, 2], [6, 6], [2, 2]]])
        self.assertListEqual(seg_end3s,
                             [[[6, 3], [5, 8], [5, 8], [8, 3], [8, 8], [3, 3]]])
        self.assertListEqual(muts,
                             [{1: {},
                               2: {3: [0]},
                               3: {0: [0]},
                               4: {},
                               5: {0: [2]},
                               6: {3: [1]},
                               7: {2: [4], 3: [1]},
                               8: {3: [1]},
                               9: {}}])

    def test_batch_size_1(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(batch_size=1)
        self.assertListEqual(read_nums, [[0], [0], [0], [0], [0], [0]])
        self.assertListEqual(
            seg_end5s,
            [[[2, 1]], [[3, 4]], [[2, 4]], [[6, 1]], [[5, 5]], [[1, 1]]]
        )
        self.assertListEqual(
            seg_end3s,
            [[[7, 4]], [[6, 9]], [[6, 9]], [[9, 4]], [[9, 9]], [[4, 4]]]
        )
        self.assertListEqual(muts,
                             [{1: {},
                               2: {},
                               3: {0: [0]},
                               4: {},
                               5: {},
                               6: {},
                               7: {64: [0]},
                               8: {},
                               9: {}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {},
                               6: {},
                               7: {3: [0]},
                               8: {3: [0]},
                               9: {}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {0: [0]},
                               6: {},
                               7: {},
                               8: {},
                               9: {}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {},
                               6: {},
                               7: {},
                               8: {},
                               9: {64: [0]}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {},
                               6: {},
                               7: {2: [0]},
                               8: {},
                               9: {}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {},
                               6: {},
                               7: {},
                               8: {},
                               9: {}}])

    def test_batch_size_4(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(batch_size=4)
        self.assertListEqual(read_nums, [[0, 1, 2, 3], [0, 1]])
        self.assertListEqual(
            seg_end5s,
            [[[2, 1], [3, 4], [2, 4], [6, 1]], [[5, 5], [1, 1]]]
        )
        self.assertListEqual(
            seg_end3s,
            [[[7, 4], [6, 9], [6, 9], [9, 4]], [[9, 9], [4, 4]]]
        )
        self.assertListEqual(muts,
                             [{1: {},
                               2: {},
                               3: {0: [0]},
                               4: {},
                               5: {0: [2]},
                               6: {},
                               7: {3: [1], 64: [0]},
                               8: {3: [1]},
                               9: {64: [3]}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {},
                               6: {},
                               7: {2: [0]},
                               8: {},
                               9: {}}])

    def test_batch_size_5(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(batch_size=5)
        self.assertListEqual(read_nums, [[0, 1, 2, 3, 4], [0]])
        self.assertListEqual(
            seg_end5s,
            [[[2, 1], [3, 4], [2, 4], [6, 1], [5, 5]], [[1, 1]]]
        )
        self.assertListEqual(
            seg_end3s,
            [[[7, 4], [6, 9], [6, 9], [9, 4], [9, 9]], [[4, 4]]]
        )
        self.assertListEqual(muts,
                             [{1: {},
                               2: {},
                               3: {0: [0]},
                               4: {},
                               5: {0: [2]},
                               6: {},
                               7: {2: [4], 3: [1], 64: [0]},
                               8: {3: [1]},
                               9: {64: [3]}},
                              {1: {},
                               2: {},
                               3: {},
                               4: {},
                               5: {},
                               6: {},
                               7: {},
                               8: {},
                               9: {}}])

    def test_batch_size_6(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches(batch_size=6)
        self.assertListEqual(read_nums, [[0, 1, 2, 3, 4, 5]])
        self.assertListEqual(seg_end5s,
                             [[[2, 1], [3, 4], [2, 4], [6, 1], [5, 5], [1, 1]]])
        self.assertListEqual(seg_end3s,
                             [[[7, 4], [6, 9], [6, 9], [9, 4], [9, 9], [4, 4]]])
        self.assertListEqual(muts,
                             [{1: {},
                               2: {},
                               3: {0: [0]},
                               4: {},
                               5: {0: [2]},
                               6: {},
                               7: {2: [4], 3: [1], 64: [0]},
                               8: {3: [1]},
                               9: {64: [3]}}])


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
