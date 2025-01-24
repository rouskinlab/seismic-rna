import os
import shutil
import unittest as ut
from itertools import chain, combinations
from pathlib import Path
from typing import Iterable

from seismicrna.align import run as run_align
from seismicrna.cluster import run as run_cluster
from seismicrna.core import path
from seismicrna.core.arg.cli import opt_out_dir, opt_sim_dir
from seismicrna.core.logs import Level, get_config, set_config
from seismicrna.core.ngs import DuplicateSampleReferenceError
from seismicrna.fold import run as run_fold
from seismicrna.join import run as run_join
from seismicrna.mask import run as run_mask
from seismicrna.pool import run as run_pool
from seismicrna.relate import run as run_relate
from seismicrna.sim.fastq import run as run_sim_fastq
from seismicrna.sim.fold import run as run_sim_fold
from seismicrna.sim.params import run as run_sim_params
from seismicrna.sim.ref import run as run_sim_ref
from seismicrna.sim.total import run as run_sim_total
from seismicrna.wf import run as run_wf

STEPS = ["align", "relate", "mask", "cluster", "fold", "graph"]


class TestWorkflow(ut.TestCase):
    OUT_DIR = Path(opt_out_dir.default).absolute()
    SIM_DIR = Path(opt_sim_dir.default).absolute()
    REFS = "test_refs"
    REF = "test_ref"
    SAMPLE = "test_sample"

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR,
                   log_file_path=None,
                   raise_on_error=True)
        self.SIM_DIR.mkdir()
        self.OUT_DIR.mkdir()

    def tearDown(self):
        shutil.rmtree(self.SIM_DIR)
        shutil.rmtree(self.OUT_DIR)
        set_config(*self._config)

    def test_wf_sim_20000reads_2clusts(self):
        # Simulate the data to be processed with wf.
        fastqs = run_sim_total(sim_dir=self.SIM_DIR,
                               sample=self.SAMPLE,
                               ref=self.REF,
                               refs=self.REFS,
                               fold_max=2,
                               num_reads=20000,
                               center_fvar=0.001,
                               length_fvar=0.001)
        sample_dir = self.SIM_DIR.joinpath("samples", self.SAMPLE)
        for fastq, mate in zip(fastqs, [1, 2], strict=True):
            self.assertEqual(fastq,
                             sample_dir.joinpath(f"{self.REF}_R{mate}.fq.gz"))
            self.assertTrue(os.path.isfile(fastq))
        fasta = self.SIM_DIR.joinpath("refs", f"{self.REFS}.fa")
        self.assertTrue(fasta.is_file())
        # Process the data with wf.
        run_wf(fasta=fasta,
               input_path=(),
               dmfastqx=fastqs,
               cluster=True,
               fold=True,
               quantile=0.95,
               export=True,
               out_dir=self.OUT_DIR)
        for step in STEPS:
            step_dir = self.OUT_DIR.joinpath(self.SAMPLE, step)
            self.assertTrue(step_dir.is_dir())

    def test_wf_sim_2samples_2refs_20000reads_2clusts(self):
        # Simulate the data to be processed with wf.
        samples = ["sample1", "sample2"]
        refs = ["refA", "refB"]
        samples_dir = self.SIM_DIR.joinpath("samples")
        all_fastas = list()
        for ref in refs:
            run_sim_ref(refs=ref, ref=ref)
            fasta = self.SIM_DIR.joinpath("refs", f"{ref}.fa")
            all_fastas.append(fasta)
            run_sim_fold(fasta, fold_max=2)
            param_dir = self.SIM_DIR.joinpath("params", ref, "full")
            ct_file = param_dir.joinpath("simulated.ct")
            run_sim_params(ct_file=(ct_file,))
            for sample in samples:
                fastqs = run_sim_fastq(input_path=(),
                                       param_dir=(param_dir,),
                                       sample=sample,
                                       num_reads=10000)
                sample_dir = samples_dir.joinpath(sample)
                for fastq, mate in zip(fastqs, [1, 2], strict=True):
                    self.assertEqual(
                        fastq,
                        sample_dir.joinpath(f"{ref}_R{mate}.fq.gz")
                    )
                    self.assertTrue(os.path.isfile(fastq))
        # Merge the FASTA files for all references.
        fasta = self.SIM_DIR.joinpath("refs", f"{self.REFS}.fa")
        with open(fasta, "x") as f:
            for ref in all_fastas:
                with open(ref) as r:
                    f.write(r.read())
        # Process the data with wf.
        run_wf(fasta=fasta,
               input_path=(),
               dmfastqx=(samples_dir,),
               cluster=True,
               fold=True,
               quantile=0.95,
               export=True,
               out_dir=self.OUT_DIR)
        for sample in samples:
            for step in STEPS:
                step_dir = self.OUT_DIR.joinpath(sample, step)
                self.assertTrue(step_dir.is_dir())


class TestWorkflowTwoOutDirs(ut.TestCase):
    NUMBERS = [1, 2]
    OUT_DIR = Path(opt_out_dir.default).absolute()
    SIM_DIR = Path(opt_sim_dir.default).absolute()
    SIM_DIRS = tuple(Path(f"sim{i}").absolute() for i in NUMBERS)
    OUT_DIRS = tuple(Path(f"out{i}").absolute() for i in NUMBERS)
    REFS = "test_refs"
    REF = "test_ref"
    SAMPLE = "test_sample"
    POOLED = "pooled_sample"
    MJOINED = "mjoined_region"
    CJOINED = "cjoined_region"

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR,
                   log_file_path=None,
                   raise_on_error=True)
        self.SIM_DIR.mkdir()
        self.OUT_DIR.mkdir()
        for sim_dir, out_dir in zip(self.SIM_DIRS, self.OUT_DIRS, strict=True):
            sim_dir.mkdir()
            out_dir.mkdir()

    def tearDown(self):
        shutil.rmtree(self.SIM_DIR)
        shutil.rmtree(self.OUT_DIR)
        for sim_dir, out_dir in zip(self.SIM_DIRS, self.OUT_DIRS, strict=True):
            shutil.rmtree(sim_dir)
            shutil.rmtree(out_dir)
        set_config(*self._config)

    def check_no_identical(self, files: Iterable[Path | str], binary: bool):
        """ Confirm no two files have identical contents. """
        mode = "rb" if binary else "rt"
        for file1, file2 in combinations(files, 2):
            with open(file1, mode) as f1, open(file2, mode) as f2:
                self.assertNotEqual(f1.read(), f2.read())

    def test_wf_two_out_dirs(self):
        fasta = run_sim_ref(sim_dir=str(self.SIM_DIR),
                            ref=self.REF,
                            refs=self.REFS,
                            reflen=60)
        ct_file, = run_sim_fold(fasta, fold_max=1)
        dmfastqxs = list()
        for sim_dir, out_dir in zip(self.SIM_DIRS, self.OUT_DIRS, strict=True):
            ct_file_copy = path.transpath(sim_dir,
                                          self.SIM_DIR,
                                          ct_file,
                                          strict=True)
            param_dir = ct_file_copy.parent
            param_dir.mkdir(parents=True)
            ct_file_copy.hardlink_to(ct_file)
            run_sim_params(ct_file=[ct_file_copy])
            dmfastqxs.append(run_sim_fastq(input_path=[],
                                           param_dir=[param_dir],
                                           sample=self.SAMPLE,
                                           read_length=30,
                                           num_reads=2000,
                                           fq_gzip=False))
        # Align FASTQ files in the same output directory.
        min_reads = 1
        min_mapq = 10
        align_kwargs = dict(min_reads=min_reads,
                            min_mapq=min_mapq,
                            bt2_score_min_loc="L,1,0.5",
                            fastp_poly_g="yes",
                            force=True)
        self.assertRaisesRegex(DuplicateSampleReferenceError,
                               str((self.SAMPLE, self.REF)),
                               run_align,
                               fasta,
                               dmfastqx=list(chain(*dmfastqxs)),
                               out_dir=str(self.OUT_DIR),
                               **align_kwargs)
        # Align FASTQ files in different output directories.
        bam_files = list()
        for dmfastqx, out_dir in zip(dmfastqxs, self.OUT_DIRS, strict=True):
            bam_dir, = run_align(fasta,
                                 dmfastqx=dmfastqx,
                                 out_dir=str(out_dir),
                                 **align_kwargs)
            bam_file = out_dir.joinpath(self.SAMPLE, "align", f"{self.REF}.bam")
            self.assertTrue(bam_file.is_file())
            self.assertEqual(bam_dir, bam_file.parent)
            bam_files.append(bam_file)
        self.check_no_identical(bam_files, True)
        # Relate BAM files in the same output directory.
        relate_kwargs = dict(min_reads=min_reads,
                             min_mapq=min_mapq)
        self.assertRaisesRegex(DuplicateSampleReferenceError,
                               str((self.SAMPLE, self.REF)),
                               run_relate,
                               fasta,
                               bam_files,
                               out_dir=str(self.OUT_DIR),
                               **relate_kwargs)
        # Relate BAM files in different output directories.
        relate_reports = list()
        for bam_file, out_dir in zip(bam_files, self.OUT_DIRS, strict=True):
            relate_dir, = run_relate(fasta,
                                     (str(bam_file),),
                                     out_dir=str(out_dir),
                                     **relate_kwargs)
            relate_report = out_dir.joinpath(self.SAMPLE,
                                             "relate",
                                             self.REF,
                                             "relate-report.json")
            self.assertTrue(relate_report.is_file())
            self.assertEqual(relate_dir, relate_report.parent)
            relate_reports.append(relate_report)
        self.check_no_identical(relate_reports, False)
        # Pool relate reports.
        pool_dirs = run_pool(relate_reports, pooled=self.POOLED)
        pool_reports = [out_dir.joinpath(self.POOLED,
                                         "relate",
                                         self.REF,
                                         "relate-report.json")
                        for out_dir in self.OUT_DIRS]
        for pool_report, pool_dir in zip(pool_reports, pool_dirs, strict=True):
            self.assertTrue(pool_report.is_file())
            self.assertEqual(pool_dir, pool_report.parent)
        # Mask relate reports.
        mask_dirs = run_mask(relate_reports,
                             mask_coords=[(self.REF, 5, 50)],
                             min_ninfo_pos=1)
        mask_reports = [out_dir.joinpath(self.SAMPLE,
                                         "mask",
                                         self.REF,
                                         "5-50",
                                         "mask-report.json")
                        for out_dir in self.OUT_DIRS]
        for mask_report, mask_dir in zip(mask_reports, mask_dirs, strict=True):
            self.assertTrue(mask_report.is_file())
            self.assertEqual(mask_dir, mask_report.parent)
        self.check_no_identical(mask_reports, False)
        # Join mask reports.
        mjoin_dirs = run_join(mask_reports, joined=self.MJOINED)
        mjoin_reports = [out_dir.joinpath(self.SAMPLE,
                                          "mask",
                                          self.REF,
                                          self.MJOINED,
                                          "mask-report.json")
                         for out_dir in self.OUT_DIRS]
        for mjoin_report, mjoin_dir in zip(mjoin_reports,
                                           mjoin_dirs,
                                           strict=True):
            self.assertTrue(mjoin_report.is_file())
            self.assertEqual(mjoin_dir, mjoin_report.parent)
        # Cluster mask reports.
        cluster_dirs = run_cluster(mask_dirs,
                                   max_clusters=1,
                                   jackpot=False)
        cluster_reports = [out_dir.joinpath(self.SAMPLE,
                                            "cluster",
                                            self.REF,
                                            "5-50",
                                            "cluster-report.json")
                           for out_dir in self.OUT_DIRS]
        for cluster_report, cluster_dir in zip(cluster_reports,
                                               cluster_dirs,
                                               strict=True):
            self.assertTrue(cluster_report.is_file())
            self.assertEqual(cluster_dir, cluster_report.parent)
        self.check_no_identical(cluster_reports, False)
        # Join cluster reports.
        cjoin_dirs = run_join(cluster_reports, joined=self.CJOINED)
        cjoin_reports = [out_dir.joinpath(self.SAMPLE,
                                          step,
                                          self.REF,
                                          self.CJOINED,
                                          f"{step}-report.json")
                         for step in ["mask", "cluster"]
                         for out_dir in self.OUT_DIRS]
        for cjoin_report, cjoin_dir in zip(cjoin_reports,
                                           cjoin_dirs,
                                           strict=True):
            self.assertTrue(cjoin_report.is_file())
            self.assertEqual(cjoin_dir, cjoin_report.parent)
        # Fold mask and cluster reports.
        fold_reports = run_fold(mask_dirs + cluster_dirs,
                                quantile=1.0)
        expect_fold_reports = list()
        for region in ["5-50"]:
            for profile in ["average", "cluster-1-1"]:
                fold_dirs = [out_dir.joinpath(self.SAMPLE,
                                              "fold",
                                              self.REF,
                                              "full")
                             for out_dir in self.OUT_DIRS]
                ct_name = f"{region}__{profile}.ct"
                self.check_no_identical([fold_dir.joinpath(ct_name)
                                         for fold_dir in fold_dirs],
                                        False)
                report_name = f"{region}__{profile}__fold-report.json"
                expect_fold_reports.extend([fold_dir.joinpath(report_name)
                                            for fold_dir in fold_dirs])
        self.assertListEqual(sorted(fold_reports), sorted(expect_fold_reports))
        for fold_report in fold_reports:
            self.assertTrue(fold_report.is_file())


if __name__ == "__main__":
    ut.main()

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
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
