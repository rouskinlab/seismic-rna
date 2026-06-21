import os
import tempfile
import unittest as ut
from itertools import chain, combinations, product
from pathlib import Path
from typing import Iterable

from click.testing import CliRunner

import numpy as np

from seismicrna.align import run as run_align
from seismicrna.cluster import run as run_cluster
from seismicrna.core import path
from seismicrna.core.arg.cli import GROUP_BY_K, GROUP_ALL, KEY_PEARSON, PROBE_DMS
from seismicrna.core.header import list_ks_clusts
from seismicrna.core.logs import Level, get_config, set_config
from seismicrna.core.ngs.xam import DuplicateSampleReferenceError
from seismicrna.fold import run as run_fold
from seismicrna.graph.base import ACTION_IDMUT
from seismicrna.graph.cgroup import make_tracks
from seismicrna.graph.profile import ProfileRunner
from seismicrna.graph.corroll import RollingCorrelationRunner
from seismicrna.graph.delprof import DeltaProfileRunner
from seismicrna.graph.histpos import PositionHistogramRunner
from seismicrna.graph.scatter import ScatterRunner
from seismicrna.graph.snrroll import RollingSNRRunner
from seismicrna.join import run as run_join
from seismicrna.filter import run as run_filter
from seismicrna.filter.dataset import FilterMutsDataset
from seismicrna.pool import run as run_pool
from seismicrna.idmut import run as run_idmut
from seismicrna.sim.fastq import run as run_sim_fastq
from seismicrna.sim.fold import run as run_sim_fold
from seismicrna.sim.params import run as run_sim_params
from seismicrna.sim.ref import run as run_sim_ref
from seismicrna.wf import run as run_wf
from seismicrna.cluster.data import ClusterMutsDataset
from seismicrna.main import cli as seismic_cli


def list_step_dir_contents(
    parent_dir: Path, step: str, num_batches: int, write_read_names: bool = True
):
    files = [f"{step}-report.json", f"{step}-position-table.csv"]
    if step == "idmut":
        files.append("refseq.brickle")
        if write_read_names:
            files.extend(f"names-batch-{i}.brickle" for i in range(num_batches))
    elif step == "filter":
        files.append(f"{step}-read-table.csv.gz")
    elif step == "cluster":
        files.extend(
            [f"{step}-abundance-table.csv", "parameters", "read-counts", "statistics"]
        )
    files.extend(f"{step}-batch-{i}.brickle" for i in range(num_batches))
    return sorted(map(parent_dir.joinpath, files))


def list_actions(num_structs: int):
    profiles = ["filtered"]
    profiles.extend(f"clustered-{k}-x" for k in range(1, num_structs + 1))
    return profiles


def list_profiles(num_structs: int):
    profiles = ["average"]
    profiles.extend(
        f"cluster-{k}-{clust}" for k, clust in list_ks_clusts(range(1, num_structs + 1))
    )
    return profiles


class TestWorkflow(ut.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._config = None
        self._tmpdir = None

    @property
    def sim_dir(self):
        if self._tmpdir is None:
            return None
        return Path(self._tmpdir.name) / "sim"

    @property
    def out_dir(self):
        if self._tmpdir is None:
            return None
        return Path(self._tmpdir.name) / "out"

    def setUp(self):
        self.maxDiff = 10000
        self._config = get_config()
        set_config(verbosity=Level.ERROR, log_file_path=None, exit_on_error=True)
        self._tmpdir = tempfile.TemporaryDirectory()
        self.sim_dir.mkdir()
        self.out_dir.mkdir()

    def tearDown(self):
        self._tmpdir.cleanup()
        self._tmpdir = None
        set_config(*self._config)

    def test_wf_sim_2samples_2refs_20000reads_2clusts(self):
        # Simulate the data to be processed with wf.
        num_structs = 2
        num_reads = 10000
        batch_size = 2**16
        num_batches = (num_reads + 1) // batch_size + 1
        samples = ["test-sample-1", "test-sample-2"]
        refs = {"ref-A": 100, "ref-B": 150}
        samples_dir = self.sim_dir.joinpath("samples")
        all_fastas = list()
        for ref, reflen in refs.items():
            run_sim_ref(refs=ref, ref=ref, reflen=reflen, sim_dir=self.sim_dir, seed=0)
            fasta = self.sim_dir.joinpath("refs", f"{ref}.fa")
            all_fastas.append(fasta)
            run_sim_fold(
                fasta, probe=PROBE_DMS, fold_max=num_structs, sim_dir=self.sim_dir
            )
            param_dir = self.sim_dir.joinpath("params", ref, "full")
            ct_file = param_dir.joinpath("simulated.ct")
            run_sim_params(ct_file=[ct_file], probe=PROBE_DMS, seed=0)
            for sample in samples:
                fastqs = run_sim_fastq(
                    input_path=(),
                    param_dir=(param_dir,),
                    sample=sample,
                    num_reads=num_reads,
                    seed=0,
                )
                sample_dir = samples_dir.joinpath(sample)
                for fastq, mate in zip(fastqs, [1, 2], strict=True):
                    self.assertEqual(
                        fastq,
                        sample_dir.joinpath(path.DEMULT_STEP, f"{ref}_R{mate}.fq.gz"),
                    )
                    self.assertTrue(os.path.isfile(fastq))
        # Merge the FASTA files for all references.
        fasta = self.sim_dir.joinpath("refs", "test-refs.fa")
        with open(fasta, "x") as f:
            for ref in all_fastas:
                with open(ref) as r:
                    f.write(r.read())
        # Process the data with wf.
        graph_kwargs = dict(
            cgroup=GROUP_BY_K,
            csv=True,
            html=True,
            svg=True,
            pdf=True,
            png=True,
            verify_times=True,
            num_cpus=1,
            force=False,
        )
        rel_graph_kwargs = graph_kwargs | dict(
            rels=("m",), use_ratio=True, graph_quantile=0.0
        )
        clust_rel_graph_kwargs = rel_graph_kwargs | {
            "cgroup": GROUP_ALL,
            "k_clust_list": [(1, 1), (2, 2)],
            "k": None,
            "clust": None,
        }
        pair_graph_kwargs = rel_graph_kwargs | dict(
            out_dir=self.out_dir, comppair=True, compself=False
        )
        refs_coords = {"ref-A": [(1, 70), (31, 90)], "ref-B": [(31, 90), (101, 150)]}
        refs_regions = {
            ref: [f"{end5}-{end3}" for end5, end3 in ref_coords]
            for ref, ref_coords in refs_coords.items()
        }
        run_wf(
            fasta=fasta,
            input_path=[],
            out_dir=self.out_dir,
            dmfastqx=[samples_dir],
            seed=0,
            batch_size=batch_size,
            brotli_level=0,
            write_read_names=True,
            region_coords=[
                (ref, end5, end3)
                for ref, ref_coords in refs_coords.items()
                for end5, end3 in ref_coords
            ],
            cluster=True,
            max_clusters=num_structs,
            try_all_ks=True,
            write_all_ks=True,
            em_thresh=1.0,
            min_em_runs=2,
            jackpot=True,
            fold=True,
            fold_quantile=0.95,
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
            **graph_kwargs,
        )
        # Create additional graphs that are not part of wf.
        RollingSNRRunner.run([self.out_dir], window=21, winmin=7, **rel_graph_kwargs)
        PositionHistogramRunner.run(
            [self.out_dir], hist_bins=37, hist_margin=0.01, **rel_graph_kwargs
        )
        RollingCorrelationRunner.run(
            [self.out_dir], metric=KEY_PEARSON, window=21, winmin=7, **pair_graph_kwargs
        )
        DeltaProfileRunner.run([self.out_dir], **pair_graph_kwargs)
        ScatterRunner.run([self.out_dir], metric=KEY_PEARSON, **pair_graph_kwargs)
        # Re-run profile graph for clusters with arbitray (k, clust) list.
        ProfileRunner.run(
            [self.out_dir.joinpath(sample).joinpath("cluster") for sample in samples],
            **clust_rel_graph_kwargs,
        )
        # Confirm that all expected output files exist.
        graph_formats = [".csv", ".html", ".svg", ".pdf", ".png"]
        for sample in samples:
            self.assertTrue(self.out_dir.joinpath(f"{sample}__webapp.json").is_file())
            sample_dir = self.out_dir.joinpath(sample)
            align_dir = sample_dir.joinpath("align")
            self.assertListEqual(
                sorted(align_dir.iterdir()),
                sorted(
                    align_dir.joinpath(file)
                    for ref in refs_coords
                    for file in [
                        f"{ref}__align-report.json",
                        f"{ref}.bam",
                        f"{ref}__unaligned.fq.1.gz",
                        f"{ref}__unaligned.fq.2.gz",
                        f"{ref}__fastp.html",
                        f"{ref}__fastp.json",
                    ]
                ),
            )
            for ref, ref_coords in refs_coords.items():
                idmut_dir = sample_dir.joinpath("idmut", ref)
                self.assertListEqual(
                    sorted(idmut_dir.iterdir()),
                    sorted(list_step_dir_contents(idmut_dir, "idmut", num_batches)),
                )
                graph_full_dir = sample_dir.joinpath("graph", ref, "full")
                for ext in graph_formats:
                    for name in [
                        f"giniroll_{ACTION_IDMUT}_45-9_m-ratio-q0",
                        f"histpos_{ACTION_IDMUT}_m-ratio-q0",
                        f"mutdist_{ACTION_IDMUT}_m",
                        f"poscorr_{ACTION_IDMUT}_m",
                        f"profile_{ACTION_IDMUT}_acgtdi-ratio-q0",
                        f"profile_{ACTION_IDMUT}_m-ratio-q0",
                        f"profile_{ACTION_IDMUT}_n-count",
                        f"snrroll_{ACTION_IDMUT}_21-7_m-ratio-q0",
                    ]:
                        file = graph_full_dir.joinpath(f"{name}{ext}")
                        with self.subTest(file=file):
                            self.assertTrue(file.is_file())
                    for reg in refs_regions[ref]:
                        for action in list_actions(num_structs):
                            for name in [
                                f"aucroll_{reg}__{action}_45-9_m-ratio-q0_incl-term",
                                f"roc_{reg}__{action}_m-ratio-q0_incl-term",
                            ]:
                                file = graph_full_dir.joinpath(f"{name}{ext}")
                                with self.subTest(file=file):
                                    self.assertTrue(file.is_file())
                for reg in refs_regions[ref]:
                    filter_dir = sample_dir.joinpath("filter", ref, reg)
                    self.assertListEqual(
                        sorted(filter_dir.iterdir()),
                        sorted(
                            list_step_dir_contents(filter_dir, "filter", num_batches)
                        ),
                    )
                    cluster_dir = sample_dir.joinpath("cluster", ref, reg)
                    self.assertListEqual(
                        sorted(cluster_dir.iterdir()),
                        sorted(
                            list_step_dir_contents(cluster_dir, "cluster", num_batches)
                        ),
                    )
                    graph_reg_dir = sample_dir.joinpath("graph", ref, reg)
                    for ext in graph_formats:
                        for name in ["histread_filtered_m-count"]:
                            file = graph_reg_dir.joinpath(f"{name}{ext}")
                            with self.subTest(file=file):
                                self.assertTrue(file.is_file())
                        for action in list_actions(num_structs):
                            for name in [
                                f"giniroll_{action}_45-9_m-ratio-q0",
                                f"histpos_{action}_m-ratio-q0",
                                f"mutdist_{action}_m",
                                f"poscorr_{action}_m",
                                f"profile_{action}_m-ratio-q0",
                                f"profile_{action}_acgtdi-ratio-q0",
                                f"profile_{action}_n-count",
                                f"snrroll_{action}_21-7_m-ratio-q0",
                            ]:
                                file = graph_reg_dir.joinpath(f"{name}{ext}")
                                with self.subTest(file=file):
                                    self.assertTrue(file.is_file())
                        for name in [
                            "profile_clustered-x-x_m-ratio-q0",
                            "abundance_clustered",
                        ]:
                            file = graph_reg_dir.joinpath(f"{name}{ext}")
                            with self.subTest(file=file):
                                self.assertTrue(file.is_file())
                fold_dir = sample_dir.joinpath("fold", ref, "full")
                self.assertListEqual(
                    sorted(fold_dir.iterdir()),
                    sorted(
                        [
                            fold_dir.joinpath(file)
                            for reg in refs_regions[ref]
                            for profile in list_profiles(num_structs)
                            for file in [
                                f"{reg}__{profile}.ct",
                                f"{reg}__{profile}.db",
                                f"{reg}__{profile}__fold-report.json",
                                f"{reg}__{profile}__varna-color.txt",
                            ]
                        ]
                    ),
                )
        cluster_report = cluster_dir.joinpath("cluster-report.json")
        cluster_dataset = ClusterMutsDataset(cluster_report)
        target_tracks = [(1, 1), (2, 2)]
        tracks = make_tracks(
            source=cluster_dataset, k=None, clust=None, k_clust_list=target_tracks
        )
        self.assertListEqual(tracks, target_tracks)
        tracks = make_tracks(
            source=cluster_dataset, k=1, clust=None, k_clust_list=[(2, 2)]
        )
        self.assertListEqual(tracks, target_tracks)
        tracks = make_tracks(
            source=cluster_dataset, k=None, clust=2, k_clust_list=[(1, 1)]
        )
        self.assertListEqual(tracks, target_tracks)
        for sample1, sample2 in combinations(sorted(samples), 2):
            sample = f"{sample1}_VS_{sample2}"
            sample_dir = self.out_dir.joinpath(sample)
            for ref, ref_regions in refs_regions.items():
                graph_full_dir = sample_dir.joinpath("graph", ref, "full")
                for ext in graph_formats:
                    for name in [
                        f"corroll_{ACTION_IDMUT}_21-7_m-ratio-q0_pcc",
                        f"delprof_{ACTION_IDMUT}_m-ratio-q0",
                        f"scatter_{ACTION_IDMUT}_m-ratio-q0",
                    ]:
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
                                f"scatter_{action}_m-ratio-q0",
                            ]:
                                file = graph_reg_dir.joinpath(f"{name}{ext}")
                                with self.subTest(file=file):
                                    self.assertTrue(file.is_file())

    def test_wf_sim_2samples_2refs_20000reads_2clusts_cli(self):
        # Simulate the data to be processed with wf (same as API test).
        num_structs = 2
        num_reads = 10000
        batch_size = 2**16
        num_batches = (num_reads + 1) // batch_size + 1
        samples = ["test-sample-1", "test-sample-2"]
        refs = {"ref-A": 100, "ref-B": 150}
        samples_dir = self.sim_dir.joinpath("samples")
        all_fastas = list()
        for ref, reflen in refs.items():
            run_sim_ref(refs=ref, ref=ref, reflen=reflen, sim_dir=self.sim_dir, seed=0)
            fasta = self.sim_dir.joinpath("refs", f"{ref}.fa")
            all_fastas.append(fasta)
            run_sim_fold(
                fasta, probe=PROBE_DMS, fold_max=num_structs, sim_dir=self.sim_dir
            )
            param_dir = self.sim_dir.joinpath("params", ref, "full")
            ct_file = param_dir.joinpath("simulated.ct")
            run_sim_params(ct_file=[ct_file], probe=PROBE_DMS, seed=0)
            for sample in samples:
                fastqs = run_sim_fastq(
                    input_path=(),
                    param_dir=(param_dir,),
                    sample=sample,
                    num_reads=num_reads,
                    seed=0,
                )
                sample_dir = samples_dir.joinpath(sample)
                for fastq, mate in zip(fastqs, [1, 2], strict=True):
                    self.assertEqual(
                        fastq,
                        sample_dir.joinpath(path.DEMULT_STEP, f"{ref}_R{mate}.fq.gz"),
                    )
                    self.assertTrue(os.path.isfile(fastq))
        # Merge the FASTA files for all references.
        fasta = self.sim_dir.joinpath("refs", "test-refs.fa")
        with open(fasta, "x") as f:
            for ref in all_fastas:
                with open(ref) as r:
                    f.write(r.read())
        # Build region_coords for CLI: -c ref end5 end3 per region.
        refs_coords = {"ref-A": [(1, 70), (31, 90)], "ref-B": [(31, 90), (101, 150)]}
        refs_regions = {
            ref: [f"{end5}-{end3}" for end5, end3 in ref_coords]
            for ref, ref_coords in refs_coords.items()
        }
        region_coords_args = []
        for ref, ref_coords in refs_coords.items():
            for end5, end3 in ref_coords:
                region_coords_args += ["-c", ref, str(end5), str(end3)]
        # Invoke `seismic wf` via the CLI runner.
        runner = CliRunner()
        args = (
            [
                "-qq",
                "--exit-on-error",
                "wf",
                str(fasta),
                "-o",
                str(self.out_dir),
                "-X",
                str(samples_dir),
                "--seed",
                "0",
                "--batch-size",
                str(batch_size),
                "--brotli-level",
                "0",
                "--write-read-names",
            ]
            + region_coords_args
            + [
                "--cluster",
                "--max-clusters",
                str(num_structs),
                "--try-all-ks",
                "--write-all-ks",
                "--em-thresh",
                "1.0",
                "--min-em-runs",
                "2",
                "--jackpot",
                "--fold",
                "--fold-quantile",
                "0.95",
                "--export",
                "--graph-mprof",
                "--graph-tmprof",
                "--graph-ncov",
                "--graph-mhist",
                "--graph-abundance",
                "--graph-giniroll",
                "--graph-roc",
                "--graph-aucroll",
                "--graph-poscorr",
                "--graph-mutdist",
                "--mutdist-null",
                "--cgroup",
                GROUP_BY_K,
                "--csv",
                "--html",
                "--svg",
                "--pdf",
                "--png",
                "--verify-times",
                "--num-cpus",
                "1",
                "--no-force",
            ]
        )
        result = runner.invoke(seismic_cli, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, msg=result.output)
        # Restore the logging config changed by the CLI invocation so that
        # subsequent API calls behave the same as in the API test.
        set_config(verbosity=Level.ERROR, log_file_path=None, exit_on_error=True)
        # Create additional graphs that are not part of wf (same as API test).
        graph_kwargs = dict(
            cgroup=GROUP_BY_K,
            csv=True,
            html=True,
            svg=True,
            pdf=True,
            png=True,
            verify_times=True,
            num_cpus=1,
            force=False,
        )
        rel_graph_kwargs = graph_kwargs | dict(
            rels=("m",), use_ratio=True, graph_quantile=0.0
        )
        clust_rel_graph_kwargs = rel_graph_kwargs | {
            "cgroup": GROUP_ALL,
            "k_clust_list": [(1, 1), (2, 2)],
            "k": None,
            "clust": None,
        }
        PositionHistogramRunner.run(
            [self.out_dir], hist_bins=37, hist_margin=0.01, **rel_graph_kwargs
        )
        # Re-run profile graph for clusters with arbitrary (k, clust) list.
        ProfileRunner.run(
            [self.out_dir.joinpath(sample).joinpath("cluster") for sample in samples],
            **clust_rel_graph_kwargs,
        )
        graph_formats = [".csv", ".html", ".svg", ".pdf", ".png"]
        for sample in samples:
            self.assertTrue(self.out_dir.joinpath(f"{sample}__webapp.json").is_file())
            sample_dir = self.out_dir.joinpath(sample)
            align_dir = sample_dir.joinpath("align")
            self.assertListEqual(
                sorted(align_dir.iterdir()),
                sorted(
                    align_dir.joinpath(file)
                    for ref in refs_coords
                    for file in [
                        f"{ref}__align-report.json",
                        f"{ref}.bam",
                        f"{ref}__unaligned.fq.1.gz",
                        f"{ref}__unaligned.fq.2.gz",
                        f"{ref}__fastp.html",
                        f"{ref}__fastp.json",
                    ]
                ),
            )
            for ref, ref_coords in refs_coords.items():
                idmut_dir = sample_dir.joinpath("idmut", ref)
                self.assertListEqual(
                    sorted(idmut_dir.iterdir()),
                    sorted(list_step_dir_contents(idmut_dir, "idmut", num_batches)),
                )
                graph_full_dir = sample_dir.joinpath("graph", ref, "full")
                for ext in graph_formats:
                    for name in [
                        f"giniroll_{ACTION_IDMUT}_45-9_m-ratio-q0",
                        f"histpos_{ACTION_IDMUT}_m-ratio-q0",
                        f"mutdist_{ACTION_IDMUT}_m",
                        f"poscorr_{ACTION_IDMUT}_m",
                        f"profile_{ACTION_IDMUT}_acgtdi-ratio-q0",
                        f"profile_{ACTION_IDMUT}_m-ratio-q0",
                        f"profile_{ACTION_IDMUT}_n-count",
                    ]:
                        file = graph_full_dir.joinpath(f"{name}{ext}")
                        with self.subTest(file=file):
                            self.assertTrue(file.is_file())
                    for reg in refs_regions[ref]:
                        for action in list_actions(num_structs):
                            for name in [
                                f"aucroll_{reg}__{action}_45-9_m-ratio-q0_incl-term",
                                f"roc_{reg}__{action}_m-ratio-q0_incl-term",
                            ]:
                                file = graph_full_dir.joinpath(f"{name}{ext}")
                                with self.subTest(file=file):
                                    self.assertTrue(file.is_file())
                for reg in refs_regions[ref]:
                    filter_dir = sample_dir.joinpath("filter", ref, reg)
                    self.assertListEqual(
                        sorted(filter_dir.iterdir()),
                        sorted(
                            list_step_dir_contents(filter_dir, "filter", num_batches)
                        ),
                    )
                    cluster_dir = sample_dir.joinpath("cluster", ref, reg)
                    self.assertListEqual(
                        sorted(cluster_dir.iterdir()),
                        sorted(
                            list_step_dir_contents(cluster_dir, "cluster", num_batches)
                        ),
                    )
                    graph_reg_dir = sample_dir.joinpath("graph", ref, reg)
                    for ext in graph_formats:
                        for name in ["histread_filtered_m-count"]:
                            file = graph_reg_dir.joinpath(f"{name}{ext}")
                            with self.subTest(file=file):
                                self.assertTrue(file.is_file())
                        for action in list_actions(num_structs):
                            for name in [
                                f"giniroll_{action}_45-9_m-ratio-q0",
                                f"histpos_{action}_m-ratio-q0",
                                f"mutdist_{action}_m",
                                f"poscorr_{action}_m",
                                f"profile_{action}_m-ratio-q0",
                                f"profile_{action}_acgtdi-ratio-q0",
                                f"profile_{action}_n-count",
                            ]:
                                file = graph_reg_dir.joinpath(f"{name}{ext}")
                                with self.subTest(file=file):
                                    self.assertTrue(file.is_file())
                        for name in [
                            "profile_clustered-x-x_m-ratio-q0",
                            "abundance_clustered",
                        ]:
                            file = graph_reg_dir.joinpath(f"{name}{ext}")
                            with self.subTest(file=file):
                                self.assertTrue(file.is_file())
                fold_dir = sample_dir.joinpath("fold", ref, "full")
                self.assertListEqual(
                    sorted(fold_dir.iterdir()),
                    sorted(
                        [
                            fold_dir.joinpath(file)
                            for reg in refs_regions[ref]
                            for profile in list_profiles(num_structs)
                            for file in [
                                f"{reg}__{profile}.ct",
                                f"{reg}__{profile}.db",
                                f"{reg}__{profile}__fold-report.json",
                                f"{reg}__{profile}__varna-color.txt",
                            ]
                        ]
                    ),
                )
        cluster_report = cluster_dir.joinpath("cluster-report.json")
        cluster_dataset = ClusterMutsDataset(cluster_report)
        target_tracks = [(1, 1), (2, 2)]
        tracks = make_tracks(
            source=cluster_dataset, k=None, clust=None, k_clust_list=target_tracks
        )
        self.assertListEqual(tracks, target_tracks)
        tracks = make_tracks(
            source=cluster_dataset, k=1, clust=None, k_clust_list=[(2, 2)]
        )
        self.assertListEqual(tracks, target_tracks)
        tracks = make_tracks(
            source=cluster_dataset, k=None, clust=2, k_clust_list=[(1, 1)]
        )
        self.assertListEqual(tracks, target_tracks)


class TestWorkflowTwoOutDirs(ut.TestCase):
    NUMBERS = [1, 2]
    REFS = "test_refs"
    REF = "test_ref"
    SAMPLE = "test_sample"
    POOLED = "pooled_sample"
    MJOINED = "mjoined_region"
    CJOINED = "cjoined_region"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._config = None
        self._tmpdir = None

    @property
    def sim_dir(self):
        if self._tmpdir is None:
            return None
        return Path(self._tmpdir.name) / "sim"

    @property
    def out_dir(self):
        if self._tmpdir is None:
            return None
        return Path(self._tmpdir.name) / "out"

    @property
    def sim_dirs(self):
        if self._tmpdir is None:
            return None
        return tuple(Path(self._tmpdir.name) / f"sim{i}" for i in self.NUMBERS)

    @property
    def out_dirs(self):
        if self._tmpdir is None:
            return None
        return tuple(Path(self._tmpdir.name) / f"out{i}" for i in self.NUMBERS)

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR, log_file_path=None, exit_on_error=True)
        self._tmpdir = tempfile.TemporaryDirectory()
        self.sim_dir.mkdir()
        self.out_dir.mkdir()
        for sim_dir, out_dir in zip(self.sim_dirs, self.out_dirs, strict=True):
            sim_dir.mkdir()
            out_dir.mkdir()

    def tearDown(self):
        self._tmpdir.cleanup()
        self._tmpdir = None
        set_config(*self._config)

    def check_no_identical(self, files: Iterable[Path | str], binary: bool):
        """Confirm no two files have identical contents."""
        mode = "rb" if binary else "rt"
        for file1, file2 in combinations(files, 2):
            with open(file1, mode) as f1, open(file2, mode) as f2:
                self.assertNotEqual(f1.read(), f2.read())

    def test_wf_two_out_dirs(self):
        fasta = run_sim_ref(
            sim_dir=self.sim_dir, ref=self.REF, refs=self.REFS, reflen=60
        )
        (ct_file,) = run_sim_fold(
            fasta, probe=PROBE_DMS, sim_dir=self.sim_dir, fold_max=1
        )
        dmfastqxs = list()
        for sim_dir, out_dir in zip(self.sim_dirs, self.out_dirs, strict=True):
            ct_file_copy = path.transpath(sim_dir, self.sim_dir, ct_file, strict=True)
            param_dir = ct_file_copy.parent
            param_dir.mkdir(parents=True)
            ct_file_copy.hardlink_to(ct_file)
            run_sim_params(ct_file=[ct_file_copy], probe=PROBE_DMS)
            dmfastqxs.append(
                run_sim_fastq(
                    input_path=[],
                    param_dir=[param_dir],
                    sample=self.SAMPLE,
                    seed=0,
                    read_length=30,
                    num_reads=2000,
                    fq_gzip=False,
                )
            )
        # Align FASTQ files in the same output directory.
        min_reads = 1
        min_mapq = 10
        align_kwargs = dict(
            min_reads=min_reads,
            min_mapq=min_mapq,
            bt2_score_min_loc="L,1,0.5",
            fastp_poly_g="yes",
            seed=0,
            force=True,
        )
        self.assertRaisesRegex(
            DuplicateSampleReferenceError,
            str((self.SAMPLE, self.REF)),
            run_align,
            fasta,
            dmfastqx=list(chain(*dmfastqxs)),
            out_dir=self.out_dir,
            **align_kwargs,
        )
        # Align FASTQ files in different output directories.
        bam_files = list()
        for dmfastqx, out_dir in zip(dmfastqxs, self.out_dirs, strict=True):
            (bam_dir,) = run_align(
                fasta, dmfastqx=dmfastqx, out_dir=out_dir, **align_kwargs
            )
            bam_file = out_dir.joinpath(self.SAMPLE, "align", f"{self.REF}.bam")
            self.assertTrue(bam_file.is_file())
            self.assertEqual(bam_dir, bam_file.parent)
            bam_files.append(bam_file)
        self.check_no_identical(bam_files, True)
        # ID mutations within BAM files in the same output directory.
        idmut_kwargs = dict(min_reads=min_reads, min_mapq=min_mapq, batch_size=256)
        self.assertRaisesRegex(
            DuplicateSampleReferenceError,
            str((self.SAMPLE, self.REF)),
            run_idmut,
            fasta,
            bam_files,
            out_dir=str(self.out_dir),
            **idmut_kwargs,
        )
        # ID mutations in BAM files in different output directories.
        idmut_reports = list()
        for bam_file, out_dir in zip(bam_files, self.out_dirs, strict=True):
            (idmut_dir,) = run_idmut(
                fasta, (str(bam_file),), out_dir=str(out_dir), **idmut_kwargs
            )
            idmut_report = out_dir.joinpath(
                self.SAMPLE, "idmut", self.REF, "idmut-report.json"
            )
            self.assertTrue(idmut_report.is_file())
            self.assertEqual(idmut_dir, idmut_report.parent)
            idmut_reports.append(idmut_report)
        self.check_no_identical(idmut_reports, False)
        # Pool the IDmut datasets.
        pool_dirs = sorted(run_pool(self.POOLED, idmut_reports))
        pool_reports = sorted(
            out_dir.joinpath(self.POOLED, "idmut", self.REF, "idmut-report.json")
            for out_dir in self.out_dirs
        )
        for pool_report, pool_dir in zip(pool_reports, pool_dirs, strict=True):
            self.assertTrue(pool_report.is_file())
            self.assertEqual(pool_dir, pool_report.parent)
        # Filter the IDmut datasets.
        filter_dirs = sorted(
            run_filter(
                idmut_reports, region_coords=[(self.REF, 5, 50)], min_ninfo_pos=1
            )
        )
        filter_reports = sorted(
            out_dir.joinpath(
                self.SAMPLE, "filter", self.REF, "5-50", "filter-report.json"
            )
            for out_dir in self.out_dirs
        )
        for filter_report, filter_dir in zip(filter_reports, filter_dirs, strict=True):
            self.assertTrue(filter_report.is_file())
            self.assertEqual(filter_dir, filter_report.parent)
        self.check_no_identical(filter_reports, False)
        # Join the Filter datasets.
        sjoin_dirs = sorted(run_join(self.MJOINED, filter_reports))
        sjoin_reports = sorted(
            out_dir.joinpath(
                self.SAMPLE, "filter", self.REF, self.MJOINED, "filter-report.json"
            )
            for out_dir in self.out_dirs
        )
        for sjoin_report, sjoin_dir in zip(sjoin_reports, sjoin_dirs, strict=True):
            self.assertTrue(sjoin_report.is_file())
            self.assertEqual(sjoin_dir, sjoin_report.parent)
        # Cluster the Filter datasets.
        cluster_dirs = sorted(
            run_cluster(filter_dirs, max_clusters=1, seed=0, jackpot=False)
        )
        cluster_reports = sorted(
            out_dir.joinpath(
                self.SAMPLE, "cluster", self.REF, "5-50", "cluster-report.json"
            )
            for out_dir in self.out_dirs
        )
        for cluster_report, cluster_dir in zip(
            cluster_reports, cluster_dirs, strict=True
        ):
            self.assertTrue(cluster_report.is_file())
            self.assertEqual(cluster_dir, cluster_report.parent)
        self.check_no_identical(cluster_reports, False)
        # Join the Cluster datasets.
        cjoin_dirs = sorted(run_join(self.CJOINED, cluster_reports))
        cjoin_reports = sorted(
            out_dir.joinpath(
                self.SAMPLE, step, self.REF, self.CJOINED, f"{step}-report.json"
            )
            for step in ["filter", "cluster"]
            for out_dir in self.out_dirs
        )
        for cjoin_report, cjoin_dir in zip(cjoin_reports, cjoin_dirs, strict=True):
            self.assertTrue(cjoin_report.is_file())
            self.assertEqual(cjoin_dir, cjoin_report.parent)
        # Fold the Filter and Cluster datasets.
        fold_reports = run_fold(filter_dirs + cluster_dirs, fold_quantile=1.0)
        expect_fold_reports = list()
        for region in ["5-50"]:
            for profile in ["average", "cluster-1-1"]:
                fold_dirs = [
                    out_dir.joinpath(self.SAMPLE, "fold", self.REF, "full")
                    for out_dir in self.out_dirs
                ]
                ct_name = f"{region}__{profile}.ct"
                self.check_no_identical(
                    [fold_dir.joinpath(ct_name) for fold_dir in fold_dirs], False
                )
                report_name = f"{region}__{profile}__fold-report.json"
                expect_fold_reports.extend(
                    [fold_dir.joinpath(report_name) for fold_dir in fold_dirs]
                )
        self.assertListEqual(sorted(fold_reports), sorted(expect_fold_reports))
        for fold_report in fold_reports:
            self.assertTrue(fold_report.is_file())

    def test_self_contained_filter_and_cluster(self):
        """Filter and cluster batches written with self_contained=True carry full
        mutation data and load correctly from the self-contained path."""
        # Simulate data.
        fasta = run_sim_ref(
            sim_dir=self.sim_dir, ref=self.REF, refs=self.REFS, reflen=60
        )
        (ct_file,) = run_sim_fold(
            fasta, probe=PROBE_DMS, sim_dir=self.sim_dir, fold_max=1
        )
        param_dir = ct_file.parent
        run_sim_params(ct_file=[ct_file], probe=PROBE_DMS)
        dmfastqx = run_sim_fastq(
            input_path=[],
            param_dir=[param_dir],
            sample=self.SAMPLE,
            seed=0,
            read_length=30,
            num_reads=2000,
            fq_gzip=False,
        )
        # Align and ID mutations.
        (bam_dir,) = run_align(
            fasta,
            dmfastqx=dmfastqx,
            out_dir=self.out_dir,
            min_reads=1,
            min_mapq=10,
            bt2_score_min_loc="L,1,0.5",
            fastp_poly_g="yes",
            seed=0,
            force=True,
        )
        bam_file = self.out_dir.joinpath(self.SAMPLE, "align", f"{self.REF}.bam")
        (idmut_dir,) = run_idmut(
            fasta,
            (str(bam_file),),
            out_dir=str(self.out_dir),
            min_reads=1,
            min_mapq=10,
            batch_size=256,
        )
        idmut_report = idmut_dir.joinpath("idmut-report.json")
        # Run filter without self_contained to get a baseline.
        (filter_dir,) = run_filter(
            [idmut_report],
            region_coords=[(self.REF, 5, 50)],
            min_ninfo_pos=1,
            self_contained=False,
        )
        filter_report = filter_dir.joinpath("filter-report.json")
        filter_dataset_base = FilterMutsDataset(filter_report)
        base_filter_batch = filter_dataset_base.get_batch(0)
        # Re-run filter with self_contained=True (overwrite).
        run_filter(
            [idmut_report],
            region_coords=[(self.REF, 5, 50)],
            min_ninfo_pos=1,
            self_contained=True,
            force=True,
        )
        filter_dataset_sc = FilterMutsDataset(filter_report)
        raw_filter_batch = filter_dataset_sc.dataset2.get_batch(0)
        sc_filter_batch = filter_dataset_sc.get_batch(0)
        # The raw FilterBatchIO should be self-contained.
        self.assertTrue(raw_filter_batch.is_self_contained)
        self.assertIsNotNone(raw_filter_batch.seg_end5s)
        self.assertIsNotNone(raw_filter_batch.seg_end3s)
        self.assertIsNotNone(raw_filter_batch.muts)
        # The self-contained integrated batch must match the baseline.
        np.testing.assert_array_equal(
            sc_filter_batch.read_nums, base_filter_batch.read_nums
        )
        np.testing.assert_array_equal(
            sc_filter_batch.seg_end5s, base_filter_batch.seg_end5s
        )
        np.testing.assert_array_equal(
            sc_filter_batch.seg_end3s, base_filter_batch.seg_end3s
        )
        self.assertEqual(sc_filter_batch.muts.keys(), base_filter_batch.muts.keys())
        # Run cluster without self_contained to get a baseline.
        (cluster_dir,) = run_cluster(
            [filter_dir], max_clusters=1, seed=0, jackpot=False, self_contained=False
        )
        cluster_report = cluster_dir.joinpath("cluster-report.json")
        cluster_dataset_base = ClusterMutsDataset(cluster_report)
        base_cluster_batch = cluster_dataset_base.get_batch(0)
        # Re-run cluster with self_contained=True (overwrite).
        run_cluster(
            [filter_dir],
            max_clusters=1,
            seed=0,
            jackpot=False,
            self_contained=True,
            force=True,
        )
        cluster_dataset_sc = ClusterMutsDataset(cluster_report)
        raw_cluster_batch = cluster_dataset_sc.dataset2.get_batch(0)
        sc_cluster_batch = cluster_dataset_sc.get_batch(0)
        # The raw ClusterBatchIO should be self-contained.
        self.assertTrue(raw_cluster_batch.is_self_contained)
        self.assertIsNotNone(raw_cluster_batch.seg_end5s)
        self.assertIsNotNone(raw_cluster_batch.seg_end3s)
        self.assertIsNotNone(raw_cluster_batch.muts)
        # The self-contained integrated batch must match the baseline.
        np.testing.assert_array_equal(
            sc_cluster_batch.read_nums, base_cluster_batch.read_nums
        )
        np.testing.assert_array_equal(
            sc_cluster_batch.seg_end5s, base_cluster_batch.seg_end5s
        )
        np.testing.assert_array_equal(
            sc_cluster_batch.seg_end3s, base_cluster_batch.seg_end3s
        )
        self.assertEqual(sc_cluster_batch.muts.keys(), base_cluster_batch.muts.keys())


if __name__ == "__main__":
    ut.main(verbosity=2)
