import unittest as ut
from abc import ABC, abstractmethod
from pathlib import Path
from shutil import rmtree

from seismicrna.core import path
from seismicrna.core.batch.muts import SectionMutsBatch
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
SAM_DATA_SINGLE = [
    ("@SQ", f"SN:{REF}", f"LN:{len(REF_SEQ)}"),
    ("Read1", 0, REF, 1, 6, "2M1D1M", "*", 0, 0, "CGC", "III"),
    ("Read2", 0, REF, 4, 5, "2M2D2M", "*", 0, 0, "CATC", "IIII"),
    ("Read3", 0, REF, 3, 9, "4M", "*", 0, 0, "GCAT", "HIHI"),
    ("Read4", 0, REF, 4, 7, "4M", "*", 0, 0, "CGTA", "IHIH"),
]


def write_sam_file(tmp_dir: Path, data: list[tuple[[str | int], ...]]):
    sam_dir = tmp_dir.joinpath(SAMPLE, path.CMD_ALIGN_DIR)
    sam_dir.mkdir(parents=True)
    sam_file = sam_dir.joinpath(f"{REF}.sam")
    with open(sam_file, "x") as f:
        f.write("\n".join(SAM_DELIM.join(map(str, line)) for line in data))
    return sam_file


def write_fasta_file(tmp_dir: Path):
    fasta_file = tmp_dir.joinpath(f"{REFS}.fa")
    write_fasta(fasta_file, [(REF, REF_SEQ)])
    return fasta_file


def load_batches(tmp_dir: Path):
    report_file = tmp_dir.joinpath(SAMPLE,
                                   path.CMD_REL_DIR,
                                   REF,
                                   "relate-report.json")
    dataset = RelateDataset.load(report_file)
    return list(dataset.iter_batches())


def extract_batches(batches: list[SectionMutsBatch]):
    read_nums = [batch.read_nums.tolist() for batch in batches]
    seg_end5s = [batch.seg_end5s.tolist() for batch in batches]
    seg_end3s = [batch.seg_end3s.tolist() for batch in batches]
    muts = [{pos: {mut: reads.tolist() for mut, reads in muts.items()}
             for pos, muts in batch.muts.items()}
            for batch in batches]
    return read_nums, seg_end5s, seg_end3s, muts


class TestRelate(ut.TestCase, ABC):

    @classmethod
    @abstractmethod
    def get_sam_data(cls) -> list:
        """ Data for the SAM file. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._tmp_dir: Path | None = None
        self._fasta_file: Path | None = None
        self._sam_file: Path | None = None

    def setUp(self):
        self._tmp_dir = path.randdir(prefix="tmp-")
        self._fasta_file = write_fasta_file(self._tmp_dir)
        self._sam_file = write_sam_file(self._tmp_dir, self.get_sam_data())
        set_config(verbosity=Level.SEVERE, raise_on_error=True)

    def tearDown(self):
        rmtree(self._tmp_dir)
        self._tmp_dir = None
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
        run(str(self._fasta_file),
            (str(self._sam_file),),
            out_dir=str(self._tmp_dir),
            min_reads=min_reads,
            min_mapq=min_mapq,
            min_phred=min_phred,
            clip_end5=clip_end5,
            clip_end3=clip_end3,
            **kwargs)
        return load_batches(self._tmp_dir)


class TestRelateEmpty(TestRelate):

    @classmethod
    def get_sam_data(cls):
        return SAM_DATA_EMPTY

    def test_noargs(self):
        read_nums, seg_end5s, seg_end3s, muts = self.batches()
        self.assertEqual(len(read_nums), 0)

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
        read_nums, seg_end5s, seg_end3s, muts = extract_batches(self.batches())
        self.assertEqual(len(read_nums), 1)
        self.assertListEqual(read_nums[0], [0, 1, 2, 3])
        self.assertListEqual(seg_end5s[0], [[1], [4], [3], [4]])
        self.assertListEqual(seg_end3s[0], [[4], [9], [6], [7]])
        self.assertDictEqual(muts[0],
                             {1: {},
                              2: {3: [0]},
                              3: {3: [0]},
                              4: {},
                              5: {3: [1], 64: [3]},
                              6: {3: [1]},
                              7: {3: [1]},
                              8: {3: [1]},
                              9: {}})

    def test_min_reads(self):
        self.assertRaisesRegex(ValueError,
                               "Insufficient reads in alignment map",
                               self.batches,
                               min_reads=5)

    def test_min_mapq(self):
        self.assertRaisesRegex(ValueError,
                               "Read 'Read1' mapped with quality score 6, "
                               "less than the minimum of 7",
                               self.batches,
                               min_mapq=7)

    def test_min_phred(self):
        read_nums, seg_end5s, seg_end3s, muts = extract_batches(self.batches(
            min_phred=40
        ))
        self.assertEqual(len(read_nums), 1)
        self.assertListEqual(read_nums[0], [0, 1, 2, 3])
        self.assertListEqual(seg_end5s[0], [[1], [4], [3], [4]])
        self.assertListEqual(seg_end3s[0], [[4], [9], [6], [7]])
        self.assertDictEqual(muts[0],
                             {1: {},
                              2: {3: [0]},
                              3: {3: [0], 177: [2]},
                              4: {},
                              5: {3: [1], 225: [2, 3]},
                              6: {3: [1]},
                              7: {3: [1], 225: [3]},
                              8: {3: [1]},
                              9: {}})

    def test_clip(self):
        read_nums, seg_end5s, seg_end3s, muts = extract_batches(self.batches(
            clip_end5=1,
            clip_end3=1,
        ))
        self.assertEqual(len(read_nums), 1)
        self.assertListEqual(read_nums[0], [0, 1, 2, 3])
        self.assertListEqual(seg_end5s[0], [[2], [5], [4], [5]])
        self.assertListEqual(seg_end3s[0], [[3], [8], [5], [6]])
        self.assertDictEqual(muts[0],
                             {1: {},
                              2: {3: [0]},
                              3: {3: [0]},
                              4: {},
                              5: {3: [1], 64: [3]},
                              6: {3: [1]},
                              7: {3: [1]},
                              8: {3: [1]},
                              9: {}})


if __name__ == "__main__":
    ut.main()
