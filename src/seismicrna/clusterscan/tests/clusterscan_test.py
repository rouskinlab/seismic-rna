import unittest as ut
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from click.testing import CliRunner

from seismicrna.core import path
from seismicrna.core.logs import Level, set_config
from seismicrna.core.report import DomainCoordsF, BestKsF, MergedDomainsF
from seismicrna.main import cli as seismic_cli
from seismicrna.clusterscan.report import ClusterScanReport
from seismicrna.clusterscan.gap import evaluate_gap
from seismicrna.core.seq.region import FULL_NAME
from seismicrna.filterscan.main import run as run_filterscan
from seismicrna.filterscan.tests.filterscan_test import ScanTestBase
from seismicrna.wf import run as run_wf


class _StubUniqReads(object):
    """Minimal stand-in for UniqReads exposing only what evaluate_gap uses."""

    def __init__(self, end5, positions, end5s, end3s, mut, cov):
        self.read_end5s = end5s
        self.read_end3s = end3s
        self._mut = mut
        self._cov = cov
        self.region = SimpleNamespace(
            end5=end5, unmasked_int=positions, length=mut.shape[1]
        )

    def get_mut_matrix(self):
        return self._mut

    def get_cov_matrix(self):
        return self._cov


def _simulate_spanning(assignments, mus_a, mus_b, seed=0):
    """Build a stub UniqReads of full-length spanning reads whose A/B halves
    are drawn from the given cluster profiles per read (a, b) assignment."""
    rng = np.random.default_rng(seed)
    a_pos = list(mus_a.index)
    b_pos = list(mus_b.index)
    positions = a_pos + b_pos
    end5 = positions[0]
    length = positions[-1] - end5 + 1
    n = len(assignments)
    mut = np.zeros((n, length), dtype=bool)
    cov = np.ones((n, length), dtype=bool)
    for i, (a, b) in enumerate(assignments):
        for p in a_pos:
            mut[i, p - end5] = rng.random() < mus_a.loc[p, a + 1]
        for p in b_pos:
            mut[i, p - end5] = rng.random() < mus_b.loc[p, b + 1]
    end5s = np.full(n, end5)
    end3s = np.full(n, positions[-1])
    return _StubUniqReads(end5, positions, end5s, end3s, mut, cov)


def _two_cluster_mus(positions, lo=0.02, hi=0.25):
    """Two clearly distinct cluster profiles over the given positions."""
    return pd.DataFrame(
        {1: [lo] * len(positions), 2: [hi] * len(positions)},
        index=pd.Index(positions, name="pos"),
    )


class TestEvaluateGap(ut.TestCase):
    """Unit-test the independence test on constructed spanning reads."""

    A_POS = list(range(1, 21))
    B_POS = list(range(21, 41))
    LEFT_END3 = 20
    RIGHT_END5 = 21

    def setUp(self):
        set_config(verbosity=Level.ERROR, log_file_path=None, exit_on_error=False)
        self.mus_a = _two_cluster_mus(self.A_POS)
        self.mus_b = _two_cluster_mus(self.B_POS)

    def test_independent_gap_kept(self):
        rng = np.random.default_rng(0)
        n = 3000
        assignments = [(int(rng.integers(2)), int(rng.integers(2))) for _ in range(n)]
        ur = _simulate_spanning(assignments, self.mus_a, self.mus_b)
        passed, info = evaluate_gap(
            ur, self.mus_a, self.mus_b, self.LEFT_END3, self.RIGHT_END5, 0.05
        )
        self.assertTrue(passed, msg=info)

    def test_coupled_gap_merged(self):
        # a=2 occurs only with b=2 (the (a2, b1) joint state is forbidden).
        rng = np.random.default_rng(1)
        n = 3000
        assignments = []
        for _ in range(n):
            a = int(rng.integers(2))
            b = int(rng.integers(2))
            if a == 1:
                b = 1  # force a2 -> b2
            assignments.append((a, b))
        ur = _simulate_spanning(assignments, self.mus_a, self.mus_b)
        passed, info = evaluate_gap(
            ur, self.mus_a, self.mus_b, self.LEFT_END3, self.RIGHT_END5, 0.05
        )
        self.assertFalse(passed, msg=info)

    def test_single_cluster_side_kept(self):
        rng = np.random.default_rng(2)
        n = 1000
        assignments = [(0, int(rng.integers(2))) for _ in range(n)]
        ur = _simulate_spanning(assignments, self.mus_a, self.mus_b)
        mus_a1 = self.mus_a[[1]]  # single A cluster
        passed, info = evaluate_gap(
            ur, mus_a1, self.mus_b, self.LEFT_END3, self.RIGHT_END5, 0.05
        )
        self.assertTrue(passed, msg=info)
        self.assertEqual(info.get("reason"), "single-cluster side")


class TestClusterScan(ScanTestBase):
    """Test the full filterscan -> clusterscan pipeline."""

    def load_report(self):
        top = self.sim_dir.joinpath(path.SIM_SAMPLES_DIR)
        filterscan_branches = path.add_branch(path.FILTERSCAN_STEP, "scan", {})
        clusterscan_branches = path.add_branch(
            path.CLUSTERSCAN_STEP, "", filterscan_branches
        )
        report_file = ClusterScanReport.build_path(
            {
                path.TOP: top,
                path.SAMPLE: self.SAMPLE,
                path.BRANCHES: clusterscan_branches,
                path.REF: self.REF,
                path.REG: FULL_NAME,
            }
        )
        return ClusterScanReport.load(report_file)

    def run_clusterscan_check(
        self,
        idmut_dirs: list[Path],
        expect_regions: dict[tuple[int, int], int],
        **kwargs,
    ):
        fasta = self.sim_dir.joinpath("refs", f"{self.REFS}.fa")
        run_wf(
            fasta=fasta,
            input_path=idmut_dirs,
            out_dir=self.sim_dir,
            demult=False,
            scan=True,
            wf_branch=[(path.FILTERSCAN_STEP, "scan")],
            cluster=True,
            filter_pos_table=False,
            filter_read_table=False,
            # Optimize for speed.
            min_em_runs=1,
            max_em_runs=1,
            jackpot=False,
            brotli_level=0,
            cluster_pos_table=False,
            cluster_abundance_table=False,
            **kwargs,
        )
        report = self.load_report()
        final_domains = [tuple(d) for d in report.get_field(DomainCoordsF)]
        best_ks = report.get_field(BestKsF)
        merged = report.get_field(MergedDomainsF)
        # Every final domain maps to at least one original filterscan domain,
        # and the mapping has one entry per final domain.
        self.assertEqual(len(merged), len(final_domains))
        domain_k = dict(zip(final_domains, best_ks))
        for (exp5, exp3), expect_k in expect_regions.items():
            expect_length = exp3 - exp5 + 1
            for (reg5, reg3), k in domain_k.items():
                overlap = max(min(reg3, exp3) - max(reg5, exp5), 0)
                if overlap >= expect_length / 2:
                    self.assertEqual(k, expect_k)
                    break
            else:
                raise ValueError(
                    f"Expected region {exp5, exp3} does not "
                    "overlap at least 50% of any final domain "
                    f"among {sorted(domain_k)}"
                )

    def test_domains012_read180(self):
        idmut_dirs = self.sim_data([0, 1, 2], 180, seed=0)
        self.run_clusterscan_check(idmut_dirs, {(1, 60): 2, (121, 180): 2}, seed=0)

    def test_domains012_read120(self):
        idmut_dirs = self.sim_data([0, 1, 2], 120, seed=0)
        self.run_clusterscan_check(idmut_dirs, {(1, 60): 2, (121, 180): 2}, seed=0)

    def test_domains012_read60(self):
        idmut_dirs = self.sim_data([0, 1, 2], 60, seed=0)
        self.run_clusterscan_check(idmut_dirs, {(1, 60): 2, (121, 180): 2}, seed=0)

    def test_domains02_read60(self):
        idmut_dirs = self.sim_data([0, 2], 60, seed=0)
        self.run_clusterscan_check(idmut_dirs, {(1, 60): 2, (61, 120): 2}, seed=0)

    def test_domains012_read180_cli(self):
        idmut_dirs = self.sim_data([0, 1, 2], 180, seed=0)
        filterscan_dirs = run_filterscan(
            idmut_dirs,
            branch="scan",
            brotli_level=0,
            filter_pos_table=False,
            filter_read_table=False,
        )
        runner = CliRunner()
        args = (
            ["-qq", "--exit-on-error", "clusterscan"]
            + [str(d) for d in filterscan_dirs]
            + [
                "--min-em-runs",
                "1",
                "--max-em-runs",
                "1",
                "--no-jackpot",
                "--brotli-level",
                "0",
                "--no-cluster-pos-table",
                "--no-cluster-abundance-table",
                "--seed",
                "0",
            ]
        )
        result = runner.invoke(seismic_cli, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, msg=result.output)
        set_config(verbosity=Level.ERROR, log_file_path=None, exit_on_error=True)


if __name__ == "__main__":
    ut.main(verbosity=2)
