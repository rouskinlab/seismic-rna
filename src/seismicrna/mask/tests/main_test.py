import unittest as ut
from abc import ABC, abstractmethod
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from shutil import rmtree

import numpy as np

from seismicrna.core import path
from seismicrna.core.io.seq import RefseqIO
from seismicrna.core.logs import Level, set_config
from seismicrna.core.seq.region import Region
from seismicrna.core.seq.xna import DNA
from seismicrna.mask.data import MaskMutsDataset
from seismicrna.mask.main import run as run_mask
from seismicrna.pool import run as run_pool
from seismicrna.relate.io import RelateBatchIO, ReadNamesBatchIO
from seismicrna.relate.report import RelateReport

POOLED_SAMPLE = "pooled"
REF = "ref"
REF_SEQ = DNA("CGCAAATC")
REGION = Region(REF, REF_SEQ)

# Read CGCAAATC
#    0 ========
#    1 ?=======
#    2 ==?=?===
#    3 ==/=====
#    4 ====I===
#    5   ===G=G
#    6  ===C===
#    7 =/==T==
#    8  ===.===
READ_NAMES = [f"Read{i}" for i in range(9)]
SE5 = [[1], [1], [1], [1], [1], [3], [2], [1], [2]]
SE3 = [[8], [8], [8], [8], [8], [8], [8], [7], [8]]
PE5 = [[1, 2], [3, 1], [1, 4], [5, 1], [1, 8], [3, 6], [7, 2], [1, 1], [6, 2]]
PE3 = [[1, 8], [8, 2], [3, 8], [8, 4], [7, 8], [5, 8], [8, 6], [7, 7], [8, 4]]
MUTS = {1: {209: [1]},
        2: {2: [7]},
        3: {2: [3], 177: [2]},
        4: {},
        5: {8: [4], 32: [6], 128: [7], 225: [2]},
        6: {64: [5]},
        7: {},
        8: {64: [5]}}


def write_datasets(out_dir: Path,
                   reads_per_sample: dict[str, int],
                   reads_per_batch: int,
                   end5s: list[list[int]],
                   end3s: list[list[int]]):
    if sum(reads_per_sample.values()) != len(READ_NAMES):
        raise ValueError(f"reads_per_sample must sum to {len(READ_NAMES)}, "
                         f"but got {sum(reads_per_sample.values())}")
    report_files = list()
    read_num = 0
    for sample, num_reads in reads_per_sample.items():
        began = datetime.now()
        # Write the reference sequence.
        refseq = RefseqIO(sample=sample, ref=REF, refseq=REF_SEQ)
        _, refseq_checksum = refseq.save(out_dir)
        # Assign read numbers to each batch.
        read_nums = iter(range(read_num, read_num + num_reads))
        i_to_batches = list()
        # Map the read numbers (i) to the read numbers in the batch.
        while batch_read_nums := list(zip(range(reads_per_batch), read_nums)):
            i_to_batches.append({i: batch_i for batch_i, i in batch_read_nums})
        # Take the read numbers for this sample in groups equal to the
        # number of reads per batch.
        checksums = {ReadNamesBatchIO.btype(): list(),
                     RelateBatchIO.btype(): list()}
        for batch, i_to_batch in enumerate(i_to_batches):
            names = np.array([READ_NAMES[i] for i in i_to_batch])
            batch_end5s = np.array([end5s[i] for i in i_to_batch])
            batch_end3s = np.array([end3s[i] for i in i_to_batch])
            # When assembling the mutations for this batch, map the read
            # numbers to the corresponding numbers within the batch.
            muts = dict()
            for j, rels in MUTS.items():
                muts[j] = defaultdict(list)
                for rel, reads in rels.items():
                    for i in reads:
                        batch_i = i_to_batch.get(i)
                        if batch_i is not None:
                            muts[j][rel].append(batch_i)
            # Write the batches of relate data and read names.
            relate_batch = RelateBatchIO(
                sample=sample,
                region=REGION,
                batch=batch,
                seg_end5s=batch_end5s,
                seg_end3s=batch_end3s,
                muts={j: {rel: np.array(reads)
                          for rel, reads in rels.items()}
                      for j, rels in muts.items()}
            )
            checksums[RelateBatchIO.btype()].append(
                relate_batch.save(out_dir)[1]
            )
            name_batch = ReadNamesBatchIO(sample=sample,
                                          ref=REF,
                                          batch=batch,
                                          names=names)
            checksums[ReadNamesBatchIO.btype()].append(
                name_batch.save(out_dir)[1]
            )
        # Write the report file for this sample.
        report = RelateReport(sample=sample,
                              ref=REF,
                              min_mapq=0,
                              min_phred=0,
                              phred_enc=33,
                              overhangs=True,
                              insert3=True,
                              ambindel=False,
                              clip_end5=0,
                              clip_end3=0,
                              min_reads=0,
                              n_reads_xam=num_reads,
                              n_reads_rel=num_reads,
                              n_batches=len(i_to_batches),
                              checksums=checksums,
                              refseq_checksum=refseq_checksum,
                              began=began,
                              ended=datetime.now())
        report_files.append(report.save(out_dir))
    # Pool the samples if there is more than one.
    if len(report_files) == 1:
        return report_files[0]
    pooled_files = run_pool(report_files, pool=POOLED_SAMPLE)
    if len(pooled_files) != 1:
        raise ValueError("Expected exactly 1 pooled report file, "
                         f"but got {len(pooled_files)}")
    return pooled_files[0]


def extract_read_nums(dataset: MaskMutsDataset):
    return [list(batch.read_nums) for batch in dataset.iter_batches()]


def extract_positions(dataset: MaskMutsDataset):
    return list(dataset.region.unmasked_int)


class TestMask(ut.TestCase, ABC):

    @classmethod
    @abstractmethod
    def reads_per_sample(cls) -> dict[str, int]:
        """ Number of reads per sample. """

    @classmethod
    @abstractmethod
    def reads_per_batch(cls) -> int:
        """ Number of reads per batch. """

    @classmethod
    @abstractmethod
    def end5s(cls) -> list[list[int]]:
        """ 5' end coordinates. """

    @classmethod
    @abstractmethod
    def end3s(cls) -> list[list[int]]:
        """ 3' end coordinates. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._out_dir: Path | None = None
        self._report_file: Path | None = None

    def setUp(self):
        self._out_dir = path.randdir(prefix="out-")
        self._report_file = write_datasets(self._out_dir,
                                           self.reads_per_sample(),
                                           self.reads_per_batch(),
                                           self.end5s(),
                                           self.end3s())
        set_config(verbosity=Level.SEVERE, raise_on_error=True)

    def tearDown(self):
        rmtree(self._out_dir)
        self._out_dir = None
        self._report_file = None
        set_config()

    def dataset(self, *,
                mask_del: bool = False,
                mask_ins: bool = False,
                mask_polya: int = 0,
                mask_gu: bool = False,
                mask_discontig: bool = False,
                min_ncov_read: int = 1,
                min_finfo_read: float = 0.,
                max_fmut_read: float = 1.,
                min_mut_gap: int = 0,
                min_ninfo_pos: int = 1,
                max_fmut_pos: float = 1.,
                **kwargs):
        mask_dir, = run_mask((str(self._report_file),),
                             mask_del=mask_del,
                             mask_ins=mask_ins,
                             mask_polya=mask_polya,
                             mask_gu=mask_gu,
                             mask_discontig=mask_discontig,
                             min_ncov_read=min_ncov_read,
                             min_finfo_read=min_finfo_read,
                             max_fmut_read=max_fmut_read,
                             min_mut_gap=min_mut_gap,
                             min_ninfo_pos=min_ninfo_pos,
                             max_fmut_pos=max_fmut_pos,
                             **kwargs)
        return MaskMutsDataset(mask_dir.joinpath("mask-report.json"))


class TestMaskSingle(TestMask, ABC):

    @classmethod
    def end5s(cls):
        return SE5

    @classmethod
    def end3s(cls):
        return SE3


class TestMaskPaired(TestMask, ABC):

    @classmethod
    def end5s(cls):
        return PE5

    @classmethod
    def end3s(cls):
        return PE3


class TestMask1Sample(TestMask, ABC):

    @classmethod
    def reads_per_sample(cls):
        return {"sample": 9}


class TestMask2Samples(TestMask, ABC):

    @classmethod
    def reads_per_sample(cls):
        return {"sample1": 5, "sample2": 4}


class TestMask1Batch(TestMask, ABC):

    @classmethod
    def reads_per_batch(cls):
        return max(cls.reads_per_sample().values())


class TestMaskBatches(TestMask, ABC):

    @classmethod
    def reads_per_batch(cls):
        return 3


class TestMaskSingle1Sample1Batch(TestMaskSingle,
                                  TestMask1Sample,
                                  TestMask1Batch):

    def test_nomask(self):
        dataset = self.dataset()
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_gu(self):
        dataset = self.dataset(mask_gu=True)
        self.assertListEqual(extract_positions(dataset),
                             [1, 3, 4, 5, 6, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_polya_3(self):
        dataset = self.dataset(mask_polya=3)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_polya_4(self):
        dataset = self.dataset(mask_polya=4)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_pos(self):
        dataset = self.dataset(mask_pos=[(REF, 1), (REF + "2", 4), (REF, 5)])
        self.assertListEqual(extract_positions(dataset),
                             [2, 3, 4, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_pos_file(self):
        mask_pos_file = self._out_dir.joinpath("mask_pos.csv")
        with open(mask_pos_file, "x") as f:
            f.write("\n".join(["Reference,Position",
                               f"{REF},1",
                               f"{REF}2,4",
                               f"{REF},6"]))
        dataset = self.dataset(mask_pos_file=mask_pos_file)
        self.assertListEqual(extract_positions(dataset),
                             [2, 3, 4, 5, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_pos_and_mask_pos_file(self):
        mask_pos_file = self._out_dir.joinpath("mask_pos.csv")
        with open(mask_pos_file, "x") as f:
            f.write("\n".join(["Reference,Position",
                               f"{REF},1",
                               f"{REF}2,4",
                               f"{REF},6"]))
        dataset = self.dataset(mask_pos=[(REF, 1), (REF + "2", 4), (REF, 5)],
                               mask_pos_file=mask_pos_file)
        self.assertListEqual(extract_positions(dataset),
                             [2, 3, 4, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_pos_multiple(self):
        mask_pos_file = self._out_dir.joinpath("mask_pos.csv")
        with open(mask_pos_file, "x") as f:
            f.write("\n".join(["Reference,Position",
                               f"{REF},8"]))
        dataset = self.dataset(mask_gu=True,
                               mask_polya=3,
                               mask_pos=[(REF, 3)],
                               mask_pos_file=mask_pos_file)
        self.assertListEqual(extract_positions(dataset),
                             [1])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 7]])

    def test_mask_pos_all(self):
        dataset = self.dataset(mask_pos=[(REF, j + 1) for j in range(8)])
        self.assertListEqual(extract_positions(dataset),
                             [])
        self.assertListEqual(extract_read_nums(dataset),
                             [[]])

    def test_mask_read(self):
        dataset = self.dataset(mask_read=["Read8", "Read2"])
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 3, 4, 5, 6, 7]])

    def test_mask_read_file(self):
        mask_read_file = self._out_dir.joinpath("mask_read.txt")
        with open(mask_read_file, "x") as f:
            f.write("\n".join(["Read0", "Read2"]))
        dataset = self.dataset(mask_read_file=mask_read_file)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[1, 3, 4, 5, 6, 7, 8]])

    def test_mask_read_and_mask_read_file(self):
        mask_read_file = self._out_dir.joinpath("mask_read.txt")
        with open(mask_read_file, "x") as f:
            f.write("\n".join(["Read0", "Read2"]))
        dataset = self.dataset(mask_read=["Read8", "Read2"],
                               mask_read_file=mask_read_file)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[1, 3, 4, 5, 6, 7]])

    def test_mask_discontig(self):
        dataset = self.dataset()
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_min_ncov_read_6(self):
        dataset = self.dataset(min_ncov_read=6)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_min_ncov_read_7(self):
        dataset = self.dataset(min_ncov_read=7)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 6, 7, 8]])

    def test_min_ncov_read_8(self):
        dataset = self.dataset(min_ncov_read=8)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4]])

    def test_mask_gu_min_ncov_read_5(self):
        dataset = self.dataset(mask_gu=True,
                               min_ncov_read=5)
        self.assertListEqual(extract_positions(dataset),
                             [1, 3, 4, 5, 6, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_gu_min_ncov_read_6(self):
        dataset = self.dataset(mask_gu=True,
                               min_ncov_read=6)
        self.assertListEqual(extract_positions(dataset),
                             [1, 3, 4, 5, 6, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4]])

    def test_mask_pos_min_ncov_read_5(self):
        dataset = self.dataset(mask_pos=[(REF, 4), (REF, 8)],
                               min_ncov_read=5)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 5, 6, 7])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 6, 7, 8]])

    def test_mask_pos_min_ncov_read_6(self):
        dataset = self.dataset(mask_pos=[(REF, 4), (REF, 8)],
                               min_ncov_read=6)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 5, 6, 7])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 7]])

    def test_mask_all_muts_min_ncov_read_7(self):
        dataset = self.dataset(mask_del=True,
                               mask_ins=True,
                               mask_mut=["ac", "ag", "at",
                                         "ca", "cg", "ct",
                                         "ga", "gc", "gt",
                                         "ta", "tc", "tg"],
                               min_ncov_read=7)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 1, 2, 3, 4, 6, 7, 8]])

    def test_min_finfo_read_1(self):
        dataset = self.dataset(min_finfo_read=1.)
        self.assertListEqual(extract_positions(dataset),
                             [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset),
                             [[0, 3, 4, 5, 6, 7, 8]])


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
