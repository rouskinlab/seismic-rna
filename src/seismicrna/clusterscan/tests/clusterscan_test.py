import unittest as ut
from pathlib import Path

from click.testing import CliRunner

from seismicrna.cluster.data import ClusterMutsDataset
from seismicrna.core import path
from seismicrna.core.logs import Level, set_config
from seismicrna.core.report import ClusterDirsF
from seismicrna.main import cli as seismic_cli
from seismicrna.clusterscan.main import run as run_clusterscan
from seismicrna.clusterscan.report import ClusterScanReport
from seismicrna.filterscan.main import run as run_filterscan
from seismicrna.filterscan.tests.filterscan_test import ScanTestBase


class TestClusterScan(ScanTestBase):
    """Test the full filterscan -> clusterscan pipeline."""

    def run_clusterscan_check(
        self,
        idmut_dirs: list[Path],
        expect_regions: dict[tuple[int, int], int],
        **kwargs,
    ):
        filterscan_dirs = run_filterscan(
            idmut_dirs, brotli_level=0, filter_pos_table=False, filter_read_table=False
        )
        clusterscan_dirs = run_clusterscan(
            filterscan_dirs,
            # Optimize for speed.
            min_em_runs=1,
            max_em_runs=1,
            jackpot=False,
            brotli_level=0,
            cluster_pos_table=False,
            cluster_abundance_table=False,
            **kwargs,
        )
        cluster_dirs = {}
        for clusterscan_dir in clusterscan_dirs:
            for report_file in path.find_files_chain(
                [clusterscan_dir], ClusterScanReport.get_path_seg_types()
            ):
                # cluster_dirs are stored relative to the output directory.
                top, _ = ClusterScanReport.parse_path(report_file)
                report = ClusterScanReport.load(report_file)
                for rel in report.get_field(ClusterDirsF):
                    d = top.joinpath(rel)
                    reg5, reg3 = map(int, d.name.split("-"))
                    cluster_dirs[(reg5, reg3)] = d
        for (exp5, exp3), expect_k in expect_regions.items():
            expect_length = exp3 - exp5 + 1
            for (reg5, reg3), cluster_dir in cluster_dirs.items():
                overlap = max(min(reg3, exp3) - max(reg5, exp5), 0)
                if overlap >= expect_length / 2:
                    report_file = cluster_dir.joinpath("cluster-report.json")
                    dataset = ClusterMutsDataset(report_file)
                    self.assertListEqual(dataset.ks, [expect_k])
                    break
            else:
                raise ValueError(
                    f"Expected region {exp5, exp3} does not "
                    "overlap at least 50% of any region "
                    f"among {sorted(cluster_dirs)}"
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
            idmut_dirs, brotli_level=0, filter_pos_table=False, filter_read_table=False
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
