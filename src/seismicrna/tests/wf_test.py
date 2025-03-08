import os
import shutil
import unittest as ut
from itertools import chain, combinations, product
from pathlib import Path
from typing import Iterable

from seismicrna.align import run as run_align
from seismicrna.cluster import run as run_cluster
from seismicrna.core import path
from seismicrna.core.arg.cli import GROUP_BY_K, GROUP_ALL, KEY_PEARSON
from seismicrna.core.header import list_ks_clusts, K_CLUST_KEY
from seismicrna.core.logs import Level, get_config, set_config
from seismicrna.core.ngs import DuplicateSampleReferenceError
from seismicrna.fold import run as run_fold
from seismicrna.graph.cgroup import make_tracks
from seismicrna.graph.profile import ProfileRunner
from seismicrna.graph.corroll import RollingCorrelationRunner
from seismicrna.graph.delprof import DeltaProfileRunner
from seismicrna.graph.histpos import PositionHistogramRunner
from seismicrna.graph.scatter import ScatterRunner
from seismicrna.graph.snrroll import RollingSNRRunner
from seismicrna.join import run as run_join
from seismicrna.mask import run as run_mask
from seismicrna.pool import run as run_pool
from seismicrna.relate import run as run_relate
from seismicrna.sim.fastq import run as run_sim_fastq
from seismicrna.sim.fold import run as run_sim_fold
from seismicrna.sim.params import run as run_sim_params
from seismicrna.sim.ref import run as run_sim_ref
from seismicrna.wf import run as run_wf
from seismicrna.cluster.data import ClusterMutsDataset


def list_step_dir_contents(parent_dir: Path,
                           step: str,
                           num_batches: int,
                           write_read_names: bool = True):
    files = [f"{step}-report.json", f"{step}-position-table.csv"]
    if step == "relate":
        files.append("refseq.brickle")
        if write_read_names:
            files.extend(f"names-batch-{i}.brickle" for i in range(num_batches))
    elif step == "mask":
        files.append(f"{step}-read-table.csv.gz")
    elif step == "cluster":
        files.extend([f"{step}-abundance-table.csv",
                      "parameters",
                      "read-counts",
                      "statistics"])
    files.extend(f"{step}-batch-{i}.brickle" for i in range(num_batches))
    return sorted(map(parent_dir.joinpath, files))


def list_actions(num_structs: int):
    profiles = ["filtered"]
    profiles.extend(f"clustered-{k}-x"
                    for k in range(1, num_structs + 1))
    return profiles


def list_profiles(num_structs: int):
    profiles = ["average"]
    profiles.extend(f"cluster-{k}-{clust}"
                    for k, clust in list_ks_clusts(range(1, num_structs + 1)))
    return profiles


class TestWorkflow(ut.TestCase):
    OUT_DIR = Path("wf-test-out").absolute()
    SIM_DIR = Path("wf-test-sim").absolute()

    def setUp(self):
        self.maxDiff = 10000
        self._config = get_config()
        set_config(verbosity=Level.ERROR,
                   log_file_path=None,
                   exit_on_error=True)
        self.SIM_DIR.mkdir()
        self.OUT_DIR.mkdir()

    def tearDown(self):
        shutil.rmtree(self.SIM_DIR)
        shutil.rmtree(self.OUT_DIR)
        set_config(*self._config)

    def test_wf_sim_2samples_2refs_20000reads_2clusts(self):
        # Simulate the data to be processed with wf.
        num_structs = 2
        num_reads = 10000
        batch_size = 2 ** 16
        num_batches = (num_reads + 1) // batch_size + 1
        samples = ["test-sample-1", "test-sample-2"]
        refs = {"ref-A": 100, "ref-B": 150}
        samples_dir = self.SIM_DIR.joinpath("samples")
        all_fastas = list()
        for ref, reflen in refs.items():
            run_sim_ref(refs=ref, ref=ref, reflen=reflen, sim_dir=self.SIM_DIR)
            fasta = self.SIM_DIR.joinpath("refs", f"{ref}.fa")
            all_fastas.append(fasta)
            run_sim_fold(fasta, fold_max=num_structs, sim_dir=self.SIM_DIR)
            param_dir = self.SIM_DIR.joinpath("params", ref, "full")
            ct_file = param_dir.joinpath("simulated.ct")
            run_sim_params(ct_file=[ct_file])
            for sample in samples:
                fastqs = run_sim_fastq(input_path=(),
                                       param_dir=(param_dir,),
                                       sample=sample,
                                       num_reads=num_reads)
                sample_dir = samples_dir.joinpath(sample)
                for fastq, mate in zip(fastqs, [1, 2], strict=True):
                    self.assertEqual(
                        fastq,
                        sample_dir.joinpath(f"{ref}_R{mate}.fq.gz")
                    )
                    self.assertTrue(os.path.isfile(fastq))
        # Merge the FASTA files for all references.
        fasta = self.SIM_DIR.joinpath("refs", "test-refs.fa")
        with open(fasta, "x") as f:
            for ref in all_fastas:
                with open(ref) as r:
                    f.write(r.read())
        # Process the data with wf.
        graph_kwargs = dict(cgroup=GROUP_BY_K,
                            csv=True,
                            html=True,
                            svg=True,
                            pdf=True,
                            png=True,
                            verify_times=True,
                            num_cpus=1,
                            force=False)
        rel_graph_kwargs = graph_kwargs | dict(rels=("m",),
                                               use_ratio=True,
                                               quantile=0.0)
        clust_rel_graph_kwargs = rel_graph_kwargs | {"cgroup": GROUP_ALL,
                                                     K_CLUST_KEY: [(1, 1),
                                                                   (2, 2)],
                                                     "k": None,
                                                     "clust": None}
        pair_graph_kwargs = rel_graph_kwargs | dict(out_dir=self.OUT_DIR,
                                                    comppair=True,
                                                    compself=False)
        refs_coords = {"ref-A": [(1, 70), (31, 90)],
                       "ref-B": [(31, 90), (101, 150)]}
        refs_regions = {ref: [f"{end5}-{end3}" for end5, end3 in ref_coords]
                        for ref, ref_coords in refs_coords.items()}
        run_wf(fasta=fasta,
               input_path=[],
               out_dir=self.OUT_DIR,
               dmfastqx=[samples_dir],
               batch_size=batch_size,
               brotli_level=0,
               write_read_names=True,
               mask_coords=[(ref, end5, end3)
                            for ref, ref_coords in refs_coords.items()
                            for end5, end3 in ref_coords],
               cluster=True,
               max_clusters=num_structs,
               try_all_ks=True,
               write_all_ks=True,
               em_thresh=1.,
               min_em_runs=2,
               jackpot=True,
               fold=True,
               quantile=0.95,
               export=True,
               graph_mprof=True,
               graph_tmprof=True,
               graph_ncov=True,
               graph_mhist=True,
               graph_abundance=True,
               graph_giniroll=True,
               graph_roc=True,
               graph_aucroll=True,
               graph_poscorr=True,
               graph_mutdist=True,
               mutdist_null=True,
               **graph_kwargs)
        # Create additional graphs that are not part of wf.
        RollingSNRRunner.run([self.OUT_DIR],
                             window=21,
                             winmin=7,
                             **rel_graph_kwargs)
        PositionHistogramRunner.run([self.OUT_DIR],
                                    hist_bins=37,
                                    hist_margin=0.01,
                                    **rel_graph_kwargs)
        RollingCorrelationRunner.run([self.OUT_DIR],
                                     metric=KEY_PEARSON,
                                     window=21,
                                     winmin=7,
                                     **pair_graph_kwargs)
        DeltaProfileRunner.run([self.OUT_DIR],
                               **pair_graph_kwargs)
        ScatterRunner.run([self.OUT_DIR],
                          metric=KEY_PEARSON,
                          **pair_graph_kwargs)
        # Re-run profile graph for clusters with arbitray (k, clust) list.
        ProfileRunner.run([self.OUT_DIR.joinpath(sample).joinpath("cluster")
                           for sample in samples],
                          **clust_rel_graph_kwargs)
        # Confirm that all expected output files exist.
        graph_formats = [".csv", ".html", ".svg", ".pdf", ".png"]
        for sample in samples:
            self.assertTrue(
                self.OUT_DIR.joinpath(f"{sample}__webapp.json").is_file()
            )
            sample_dir = self.OUT_DIR.joinpath(sample)
            align_dir = sample_dir.joinpath("align")
            self.assertListEqual(
                sorted(align_dir.iterdir()),
                sorted(align_dir.joinpath(file)
                       for ref in refs_coords
                       for file in [f"{ref}__align-report.json",
                                    f"{ref}.bam",
                                    f"{ref}__unaligned.fq.1.gz",
                                    f"{ref}__unaligned.fq.2.gz",
                                    f"{ref}__fastp.html",
                                    f"{ref}__fastp.json"])
            )
            for ref, ref_coords in refs_coords.items():
                relate_dir = sample_dir.joinpath("relate", ref)
                self.assertListEqual(sorted(relate_dir.iterdir()),
                                     sorted(list_step_dir_contents(
                                         relate_dir, "relate", num_batches
                                     )))
                graph_full_dir = sample_dir.joinpath("graph", ref, "full")
                for ext in graph_formats:
                    for name in ["giniroll_all_45-9_m-ratio-q0",
                                 "histpos_all_m-ratio-q0",
                                 "mutdist_all_m",
                                 "poscorr_all_m",
                                 "profile_all_acgtdi-ratio-q0",
                                 "profile_all_m-ratio-q0",
                                 "profile_all_n-count",
                                 "snrroll_all_21-7_m-ratio-q0"]:
                        file = graph_full_dir.joinpath(f"{name}{ext}")
                        with self.subTest(file=file):
                            self.assertTrue(file.is_file())
                    for reg in refs_regions[ref]:
                        for action in list_actions(num_structs):
                            for name in [
                                f"aucroll_{reg}__{action}_45-9_m-ratio-q0",
                                f"roc_{reg}__{action}_m-ratio-q0",
                            ]:
                                file = graph_full_dir.joinpath(f"{name}{ext}")
                                with self.subTest(file=file):
                                    self.assertTrue(file.is_file())
                for reg in refs_regions[ref]:
                    mask_dir = sample_dir.joinpath("mask", ref, reg)
                    self.assertListEqual(sorted(mask_dir.iterdir()),
                                         sorted(list_step_dir_contents(
                                             mask_dir, "mask", num_batches
                                         )))
                    cluster_dir = sample_dir.joinpath("cluster", ref, reg)
                    self.assertListEqual(sorted(cluster_dir.iterdir()),
                                         sorted(list_step_dir_contents(
                                             cluster_dir, "cluster", num_batches
                                         )))
                    graph_reg_dir = sample_dir.joinpath("graph", ref, reg)
                    for ext in graph_formats:
                        for name in ["histread_filtered_m-count"]:
                            file = graph_reg_dir.joinpath(f"{name}{ext}")
                            with self.subTest(file=file):
                                self.assertTrue(file.is_file())
                        for action in list_actions(num_structs):
                            for name in [f"giniroll_{action}_45-9_m-ratio-q0",
                                         f"histpos_{action}_m-ratio-q0",
                                         f"mutdist_{action}_m",
                                         f"poscorr_{action}_m",
                                         f"profile_{action}_m-ratio-q0",
                                         f"profile_{action}_acgtdi-ratio-q0",
                                         f"profile_{action}_n-count",
                                         f"snrroll_{action}_21-7_m-ratio-q0"]:
                                file = graph_reg_dir.joinpath(f"{name}{ext}")
                                with self.subTest(file=file):
                                    self.assertTrue(file.is_file())
                        for name in ["profile_clustered-x-x_m-ratio-q0",
                                     "abundance_clustered"]:
                            file = graph_reg_dir.joinpath(f"{name}{ext}")
                            with self.subTest(file=file):
                                self.assertTrue(file.is_file())
                fold_dir = sample_dir.joinpath("fold", ref, "full")
                self.assertListEqual(
                    sorted(fold_dir.iterdir()),
                    sorted([fold_dir.joinpath(file)
                            for reg in refs_regions[ref]
                            for profile in list_profiles(num_structs)
                            for file in [f"{reg}__{profile}.ct",
                                         f"{reg}__{profile}.db",
                                         f"{reg}__{profile}__fold-report.json",
                                         f"{reg}__{profile}__varna-color.txt"]
                            ])
                )
        cluster_report = cluster_dir.joinpath("cluster-report.json")
        cluster_dataset = ClusterMutsDataset(cluster_report)
        target_tracks = [(1, 1), (2, 2)]
        tracks = make_tracks(source=cluster_dataset,
                             k=None,
                             clust=None,
                             k_clust_list=target_tracks)
        self.assertListEqual(tracks, target_tracks)
        tracks = make_tracks(source=cluster_dataset,
                             k=1,
                             clust=None,
                             k_clust_list=[(2, 2)])
        self.assertListEqual(tracks, target_tracks)
        tracks = make_tracks(source=cluster_dataset,
                             k=None,
                             clust=2,
                             k_clust_list=[(1, 1)])
        self.assertListEqual(tracks, target_tracks)
        for sample1, sample2 in combinations(sorted(samples), 2):
            sample = f"{sample1}_VS_{sample2}"
            sample_dir = self.OUT_DIR.joinpath(sample)
            for ref, ref_regions in refs_regions.items():
                graph_full_dir = sample_dir.joinpath("graph", ref, "full")
                for ext in graph_formats:
                    for name in ["corroll_all_21-7_m-ratio-q0_pcc",
                                 "delprof_all_m-ratio-q0",
                                 "scatter_all_m-ratio-q0"]:
                        file = graph_full_dir.joinpath(f"{name}{ext}")
                        with self.subTest(file=file):
                            self.assertTrue(file.is_file())
                for reg in ref_regions:
                    graph_reg_dir = sample_dir.joinpath("graph", ref, reg)
                    for ext in graph_formats:
                        for action1, action2 in product(
                                list_actions(num_structs), repeat=2
                        ):
                            if action1 == action2:
                                action = action1
                            else:
                                action = f"{action1}_VS_{action2}"
                            for name in [
                                f"corroll_{action}_21-7_m-ratio-q0_pcc",
                                f"delprof_{action}_m-ratio-q0",
                                f"scatter_{action}_m-ratio-q0"
                            ]:
                                file = graph_reg_dir.joinpath(f"{name}{ext}")
                                with self.subTest(file=file):
                                    self.assertTrue(file.is_file())


class TestWorkflowTwoOutDirs(ut.TestCase):
    NUMBERS = [1, 2]
    SIM_DIR = Path("wf2-test-sim").absolute()
    OUT_DIR = Path("wf2-test-out").absolute()
    SIM_DIRS = tuple(Path(f"wf2-test-sim{i}").absolute() for i in NUMBERS)
    OUT_DIRS = tuple(Path(f"wf2-test-out{i}").absolute() for i in NUMBERS)
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
                   exit_on_error=True)
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
        fasta = run_sim_ref(sim_dir=self.SIM_DIR,
                            ref=self.REF,
                            refs=self.REFS,
                            reflen=60)
        ct_file, = run_sim_fold(fasta, sim_dir=self.SIM_DIR, fold_max=1)
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
                               out_dir=self.OUT_DIR,
                               **align_kwargs)
        # Align FASTQ files in different output directories.
        bam_files = list()
        for dmfastqx, out_dir in zip(dmfastqxs, self.OUT_DIRS, strict=True):
            bam_dir, = run_align(fasta,
                                 dmfastqx=dmfastqx,
                                 out_dir=out_dir,
                                 **align_kwargs)
            bam_file = out_dir.joinpath(self.SAMPLE, "align", f"{self.REF}.bam")
            self.assertTrue(bam_file.is_file())
            self.assertEqual(bam_dir, bam_file.parent)
            bam_files.append(bam_file)
        self.check_no_identical(bam_files, True)
        # Relate BAM files in the same output directory.
        relate_kwargs = dict(min_reads=min_reads,
                             min_mapq=min_mapq,
                             batch_size=256)
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
        pool_dirs = sorted(run_pool(relate_reports, pooled=self.POOLED))
        pool_reports = sorted(out_dir.joinpath(self.POOLED,
                                               "relate",
                                               self.REF,
                                               "relate-report.json")
                              for out_dir in self.OUT_DIRS)
        for pool_report, pool_dir in zip(pool_reports, pool_dirs, strict=True):
            self.assertTrue(pool_report.is_file())
            self.assertEqual(pool_dir, pool_report.parent)
        # Mask relate reports.
        mask_dirs = sorted(run_mask(relate_reports,
                                    mask_coords=[(self.REF, 5, 50)],
                                    min_ninfo_pos=1))
        mask_reports = sorted(out_dir.joinpath(self.SAMPLE,
                                               "mask",
                                               self.REF,
                                               "5-50",
                                               "mask-report.json")
                              for out_dir in self.OUT_DIRS)
        for mask_report, mask_dir in zip(mask_reports, mask_dirs, strict=True):
            self.assertTrue(mask_report.is_file())
            self.assertEqual(mask_dir, mask_report.parent)
        self.check_no_identical(mask_reports, False)
        # Join mask reports.
        mjoin_dirs = sorted(run_join(mask_reports, joined=self.MJOINED))
        mjoin_reports = sorted(out_dir.joinpath(self.SAMPLE,
                                                "mask",
                                                self.REF,
                                                self.MJOINED,
                                                "mask-report.json")
                               for out_dir in self.OUT_DIRS)
        for mjoin_report, mjoin_dir in zip(mjoin_reports,
                                           mjoin_dirs,
                                           strict=True):
            self.assertTrue(mjoin_report.is_file())
            self.assertEqual(mjoin_dir, mjoin_report.parent)
        # Cluster mask reports.
        cluster_dirs = sorted(run_cluster(mask_dirs,
                                          max_clusters=1,
                                          jackpot=False))
        cluster_reports = sorted(out_dir.joinpath(self.SAMPLE,
                                                  "cluster",
                                                  self.REF,
                                                  "5-50",
                                                  "cluster-report.json")
                                 for out_dir in self.OUT_DIRS)
        for cluster_report, cluster_dir in zip(cluster_reports,
                                               cluster_dirs,
                                               strict=True):
            self.assertTrue(cluster_report.is_file())
            self.assertEqual(cluster_dir, cluster_report.parent)
        self.check_no_identical(cluster_reports, False)
        # Join cluster reports.
        cjoin_dirs = sorted(run_join(cluster_reports, joined=self.CJOINED))
        cjoin_reports = sorted(out_dir.joinpath(self.SAMPLE,
                                                step,
                                                self.REF,
                                                self.CJOINED,
                                                f"{step}-report.json")
                               for step in ["mask", "cluster"]
                               for out_dir in self.OUT_DIRS)
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
    ut.main(verbosity=2)
