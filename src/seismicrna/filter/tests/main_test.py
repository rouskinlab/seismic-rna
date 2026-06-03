import tempfile
import unittest as ut
from abc import ABC, abstractmethod
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np

from seismicrna.core.logs import Level, set_config
from seismicrna.core.seq.region import Region
from seismicrna.core.seq.xna import DNA
from seismicrna.filter.dataset import FilterMutsDataset
from seismicrna.filter.main import run as run_filter
from seismicrna.pool import run as run_pool
from seismicrna.idmut.io import IDmutBatchIO, ReadNamesBatchIO, RefseqIO
from seismicrna.idmut.report import IDmutReport


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
MUTS = {
    1: {209: [1]},
    2: {2: [7]},
    3: {2: [3], 177: [2]},
    4: {},
    5: {8: [4], 32: [6], 128: [7], 225: [2]},
    6: {64: [5]},
    7: {},
    8: {64: [5]},
}


def write_datasets(
    out_dir: Path,
    reads_per_sample: dict[str, int],
    reads_per_batch: int,
    end5s: list[list[int]],
    end3s: list[list[int]],
):
    if sum(reads_per_sample.values()) != len(READ_NAMES):
        raise ValueError(
            f"reads_per_sample must sum to {len(READ_NAMES)}, "
            f"but got {sum(reads_per_sample.values())}"
        )
    branches = dict()
    report_files = list()
    read_num = 0
    for sample, num_reads in reads_per_sample.items():
        began = datetime.now()
        # Write the reference sequence.
        refseq = RefseqIO(sample=sample, branches=branches, ref=REF, refseq=REF_SEQ)
        _, refseq_checksum = refseq.save(out_dir)
        # Assign read numbers to each batch.
        read_nums = iter(range(read_num, read_num + num_reads))
        i_to_batches = list()
        # Map the read numbers (i) to the read numbers in the batch.
        while batch_read_nums := list(zip(range(reads_per_batch), read_nums)):
            i_to_batches.append({i: batch_i for batch_i, i in batch_read_nums})
        # Take the read numbers for this sample in groups equal to the
        # number of reads per batch.
        checksums = {ReadNamesBatchIO.btype(): list(), IDmutBatchIO.btype(): list()}
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
            # Write the batches of IDmut data and read names.
            idmut_batch = IDmutBatchIO(
                sample=sample,
                branches=branches,
                region=REGION,
                batch=batch,
                seg_end5s=batch_end5s,
                seg_end3s=batch_end3s,
                muts={
                    j: {rel: np.array(reads) for rel, reads in rels.items()}
                    for j, rels in muts.items()
                },
            )
            checksums[IDmutBatchIO.btype()].append(idmut_batch.save(out_dir)[1])
            name_batch = ReadNamesBatchIO(
                sample=sample, branches=branches, ref=REF, batch=batch, names=names
            )
            checksums[ReadNamesBatchIO.btype()].append(name_batch.save(out_dir)[1])
        # Write the report file for this sample.
        report = IDmutReport(
            sample=sample,
            branches=branches,
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
            ended=datetime.now(),
        )
        report_files.append(report.save(out_dir))
    # Pool the samples if there is more than one.
    if len(report_files) == 1:
        return report_files[0]
    pooled_files = run_pool(POOLED_SAMPLE, report_files)
    if len(pooled_files) != 1:
        raise ValueError(
            f"Expected exactly 1 pooled report file, but got {len(pooled_files)}"
        )
    return pooled_files[0]


def extract_read_nums(dataset: FilterMutsDataset):
    return [list(map(int, batch.read_nums)) for batch in dataset.iter_batches()]


def extract_positions(dataset: FilterMutsDataset):
    return list(map(int, dataset.region.unmasked_int))


class TestFilter(ut.TestCase, ABC):
    @classmethod
    @abstractmethod
    def reads_per_sample(cls) -> dict[str, int]:
        """Number of reads per sample."""

    @classmethod
    @abstractmethod
    def reads_per_batch(cls) -> int:
        """Number of reads per batch."""

    @classmethod
    @abstractmethod
    def end5s(cls) -> list[list[int]]:
        """5' end coordinates."""

    @classmethod
    @abstractmethod
    def end3s(cls) -> list[list[int]]:
        """3' end coordinates."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._tmp: tempfile.TemporaryDirectory | None = None
        self._out_dir: Path | None = None
        self._report_file: Path | None = None

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self._out_dir = Path(self._tmp.name)
        self._report_file = write_datasets(
            self._out_dir,
            self.reads_per_sample(),
            self.reads_per_batch(),
            self.end5s(),
            self.end3s(),
        )
        set_config(verbosity=Level.ERROR, exit_on_error=True)

    def tearDown(self):
        self._tmp.cleanup()
        self._tmp = None
        self._out_dir = None
        self._report_file = None
        set_config()

    def dataset(
        self,
        *,
        count_del: bool = True,
        count_ins: bool = True,
        mask_polya: int = 0,
        probe: str = "none",
        mask_a: bool | None = None,
        mask_c: bool | None = None,
        mask_g: bool | None = None,
        mask_u: bool | None = None,
        drop_discontig: bool = False,
        min_ncov_read: int = 1,
        min_fcov_read: float = 0.0,
        min_finfo_read: float = 0.0,
        max_fmut_read: float = 1.0,
        min_mut_gap: int = 0,
        min_ninfo_pos: int = 1,
        max_fmut_pos: float = 1.0,
        **kwargs,
    ):
        (filter_dir,) = run_filter(
            [self._report_file],
            count_del=count_del,
            count_ins=count_ins,
            mask_polya=mask_polya,
            probe=probe,
            mask_a=mask_a,
            mask_c=mask_c,
            mask_g=mask_g,
            mask_u=mask_u,
            drop_discontig=drop_discontig,
            min_ncov_read=min_ncov_read,
            min_fcov_read=min_fcov_read,
            min_finfo_read=min_finfo_read,
            max_fmut_read=max_fmut_read,
            min_mut_gap=min_mut_gap,
            min_ninfo_pos=min_ninfo_pos,
            max_fmut_pos=max_fmut_pos,
            **kwargs,
        )
        return FilterMutsDataset(filter_dir.joinpath("filter-report.json"))


class TestFilterSingle(TestFilter, ABC):
    @classmethod
    def end5s(cls):
        return SE5

    @classmethod
    def end3s(cls):
        return SE3


class TestFilterPaired(TestFilter, ABC):
    @classmethod
    def end5s(cls):
        return PE5

    @classmethod
    def end3s(cls):
        return PE3


class TestFilter1Sample(TestFilter, ABC):
    @classmethod
    def reads_per_sample(cls):
        return {"sample": 9}


class TestFilter2Samples(TestFilter, ABC):
    @classmethod
    def reads_per_sample(cls):
        return {"sample1": 5, "sample2": 4}


class TestFilter1Batch(TestFilter, ABC):
    @classmethod
    def reads_per_batch(cls):
        return max(cls.reads_per_sample().values())


class TestFilterBatches(TestFilter, ABC):
    @classmethod
    def reads_per_batch(cls):
        return 3


class TestFilterSingle1Sample1Batch(
    TestFilterSingle, TestFilter1Sample, TestFilter1Batch
):
    def test_nomask(self):
        dataset = self.dataset()
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_a(self):
        dataset = self.dataset(mask_a=True)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_c(self):
        dataset = self.dataset(mask_c=True)
        self.assertListEqual(extract_positions(dataset), [2, 4, 5, 6, 7])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_g(self):
        dataset = self.dataset(mask_g=True)
        self.assertListEqual(extract_positions(dataset), [1, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_u(self):
        dataset = self.dataset(mask_u=True)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_shape(self):
        dataset = self.dataset(probe="SHAPE")
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_dms(self):
        dataset = self.dataset(probe="DMS")
        self.assertListEqual(extract_positions(dataset), [1, 3, 4, 5, 6, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_dms_mask_a(self):
        dataset = self.dataset(probe="DMS", mask_a=True)
        self.assertListEqual(extract_positions(dataset), [1, 3, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_dms_mask_c(self):
        dataset = self.dataset(probe="DMS", mask_c=True)
        self.assertListEqual(extract_positions(dataset), [4, 5, 6])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_dms_nomask_g(self):
        dataset = self.dataset(probe="DMS", mask_g=False)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_dms_nomask_u(self):
        dataset = self.dataset(probe="DMS", mask_u=False)
        self.assertListEqual(extract_positions(dataset), [1, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_etc(self):
        dataset = self.dataset(probe="ETC")
        self.assertListEqual(extract_positions(dataset), [2, 7])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_etc_nomask_a(self):
        dataset = self.dataset(probe="ETC", mask_a=False)
        self.assertListEqual(extract_positions(dataset), [2, 4, 5, 6, 7])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_etc_nomask_c(self):
        dataset = self.dataset(probe="ETC", mask_c=False)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_etc_mask_g(self):
        dataset = self.dataset(probe="ETC", mask_g=True)
        self.assertListEqual(extract_positions(dataset), [7])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_probe_etc_mask_u(self):
        dataset = self.dataset(probe="ETC", mask_u=True)
        self.assertListEqual(extract_positions(dataset), [2])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 6, 7, 8]])

    def test_mask_polya_3(self):
        dataset = self.dataset(mask_polya=3)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_polya_4(self):
        dataset = self.dataset(mask_polya=4)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_pos(self):
        dataset = self.dataset(mask_pos=[(REF, 1), (REF + "2", 4), (REF, 5)])
        self.assertListEqual(extract_positions(dataset), [2, 3, 4, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_pos_file(self):
        mask_pos_file = self._out_dir / "filter-position-list.csv"
        with open(mask_pos_file, "x") as f:
            f.write(
                "\n".join(["Reference,Position", f"{REF},1", f"{REF}2,4", f"{REF},6"])
            )
        dataset = self.dataset(mask_pos_file=[mask_pos_file])
        self.assertListEqual(extract_positions(dataset), [2, 3, 4, 5, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_pos_files(self):
        mask_pos_files = [
            self._out_dir / "1" / "filter-position-list.csv",
            self._out_dir / "2" / "filter-position-list.csv",
            self._out_dir / "3" / "filter-position-list.csv",
        ]
        for file in mask_pos_files:
            file.parent.mkdir()
        with open(mask_pos_files[0], "x") as f:
            f.write("\n".join(["Reference,Position", f"{REF}2,4", f"{REF},6"]))
        with open(mask_pos_files[1], "x") as f:
            f.write("\n".join(["Reference,Position", f"{REF}2,2", f"{REF}2,5"]))
        with open(mask_pos_files[2], "x") as f:
            f.write("\n".join(["Reference,Position", f"{REF},8", f"{REF},3"]))
        dataset = self.dataset(mask_pos_file=mask_pos_files)
        self.assertListEqual(extract_positions(dataset), [1, 2, 4, 5, 7])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_pos_and_mask_pos_file(self):
        mask_pos_file = self._out_dir / "filter-position-list.csv"
        with open(mask_pos_file, "x") as f:
            f.write(
                "\n".join(["Reference,Position", f"{REF},1", f"{REF}2,4", f"{REF},6"])
            )
        dataset = self.dataset(
            mask_pos=[(REF, 1), (REF + "2", 4), (REF, 5)], mask_pos_file=[mask_pos_file]
        )
        self.assertListEqual(extract_positions(dataset), [2, 3, 4, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_mask_pos_multiple(self):
        mask_pos_file = self._out_dir / "filter-position-list.csv"
        with open(mask_pos_file, "x") as f:
            f.write("\n".join(["Reference,Position", f"{REF},8"]))
        dataset = self.dataset(
            probe="DMS",
            mask_polya=3,
            mask_pos=[(REF, 3)],
            mask_pos_file=[mask_pos_file],
        )
        self.assertListEqual(extract_positions(dataset), [1])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 7]])

    def test_mask_pos_all(self):
        dataset = self.dataset(mask_pos=[(REF, j + 1) for j in range(8)])
        self.assertListEqual(extract_positions(dataset), [])
        self.assertListEqual(extract_read_nums(dataset), [[]])

    def test_drop_read(self):
        dataset = self.dataset(drop_read=["Read8", "Read2"])
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 3, 4, 5, 6, 7]])

    def test_drop_read_file(self):
        drop_read_file = self._out_dir / "filter-read-list.csv.gz"
        with open(drop_read_file, "x") as f:
            f.write("\n".join(["Read0", "Read2"]))
        dataset = self.dataset(drop_read_file=[drop_read_file])
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[1, 3, 4, 5, 6, 7, 8]])

    def test_drop_read_files(self):
        drop_read_files = [
            self._out_dir / "1" / "filter-read-list.csv.gz",
            self._out_dir / "2" / "filter-read-list.csv.gz",
            self._out_dir / "3" / "filter-read-list.csv.gz",
        ]
        for file in drop_read_files:
            file.parent.mkdir()
        with open(drop_read_files[0], "x") as f:
            f.write("\n".join(["Read0", "Read2"]))
        with open(drop_read_files[1], "x") as f:
            f.write("\n".join(["Read10"]))
        with open(drop_read_files[2], "x") as f:
            f.write("\n".join(["Read5", "Read2", "Read6"]))
        dataset = self.dataset(drop_read_file=drop_read_files)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[1, 3, 4, 7, 8]])

    def test_drop_read_and_drop_read_file(self):
        drop_read_file = self._out_dir / "filter-read-list.csv.gz"
        with open(drop_read_file, "x") as f:
            f.write("\n".join(["Read0", "Read2"]))
        dataset = self.dataset(
            drop_read=["Read8", "Read2"], drop_read_file=[drop_read_file]
        )
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[1, 3, 4, 5, 6, 7]])

    def test_drop_discontig(self):
        dataset = self.dataset()
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_min_ncov_read_6(self):
        dataset = self.dataset(min_ncov_read=6)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_min_ncov_read_7(self):
        dataset = self.dataset(min_ncov_read=7)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 6, 7, 8]])

    def test_min_ncov_read_8(self):
        dataset = self.dataset(min_ncov_read=8)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4]])

    def test_dms_min_ncov_read_5(self):
        dataset = self.dataset(probe="DMS", min_ncov_read=5)
        self.assertListEqual(extract_positions(dataset), [1, 3, 4, 5, 6, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_dms_min_ncov_read_6(self):
        dataset = self.dataset(probe="DMS", min_ncov_read=6)
        self.assertListEqual(extract_positions(dataset), [1, 3, 4, 5, 6, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4]])

    def test_mask_pos_min_ncov_read_5(self):
        dataset = self.dataset(mask_pos=[(REF, 4), (REF, 8)], min_ncov_read=5)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 5, 6, 7])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 6, 7, 8]])

    def test_mask_pos_min_ncov_read_6(self):
        dataset = self.dataset(mask_pos=[(REF, 4), (REF, 8)], min_ncov_read=6)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 5, 6, 7])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 7]])

    def test_all_no_mut_min_ncov_read_7(self):
        dataset = self.dataset(
            count_del=False,
            count_ins=False,
            no_mut=[
                "ac",
                "ag",
                "at",
                "ca",
                "cg",
                "ct",
                "ga",
                "gc",
                "gt",
                "ta",
                "tc",
                "tg",
            ],
            min_ncov_read=7,
        )
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 6, 7, 8]])

    def test_only_mut_ag_min_ninfo_pos_6(self):
        # only_mut=["ag"] counts only A->G substitutions (SUB_G=64) as mutations;
        # all reference bases count as non-mutations.
        # Informative reads per position:
        #   pos 1 (C): reads 0,2,3,4,7 have MATCH -> 5 informative
        #   pos 5 (A): reads 0,1,3,5,8 have MATCH -> 5 informative
        #   all other positions: >= 7 informative reads
        # At min_ninfo_pos=6, positions 1 and 5 are masked.
        dataset = self.dataset(only_mut=["ag"], min_ninfo_pos=6)
        self.assertListEqual(extract_positions(dataset), [2, 3, 4, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 5, 6, 7, 8]])

    def test_min_finfo_read_1(self):
        dataset = self.dataset(min_finfo_read=1.0)
        self.assertListEqual(extract_positions(dataset), [1, 2, 3, 4, 5, 6, 7, 8])
        self.assertListEqual(extract_read_nums(dataset), [[0, 3, 4, 5, 6, 7, 8]])

    def test_min_fcov_read_amplicons(self):
        """Two-PCR-amplicon scenario: filter reads that partially cover the useful
        region.

        Positions 1 and 8 simulate 5' and 3' primers and are masked.
        Kept region after primer masking: positions 2-7 (6 nt, region.size=6).

        Read5 (start=3, end=8) represents a read from an adjacent amplicon
        that begins inside amplicon A's primer region.  After masking position
        8 (3' primer), it covers positions 3-7 = 5/6 ≈ 0.83 → filtered.

        All other reads start at position 1 or 2 and end at position 7 or 8,
        so after primer masking they each cover all 6 kept positions
        (positions 2-7) → fcov = 1.0 → kept.
        """
        dataset = self.dataset(mask_pos=[(REF, 1), (REF, 8)], min_fcov_read=1.0)
        self.assertListEqual(extract_positions(dataset), [2, 3, 4, 5, 6, 7])
        self.assertListEqual(extract_read_nums(dataset), [[0, 1, 2, 3, 4, 6, 7, 8]])


if __name__ == "__main__":
    ut.main(verbosity=2)
