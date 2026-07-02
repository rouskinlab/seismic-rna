import tempfile
import unittest as ut
from itertools import product
from pathlib import Path

from click.testing import CliRunner

from seismicrna.core import path
from seismicrna.core.logs import Level, get_config, set_config
from seismicrna.core.report import DomainCoordsF
from seismicrna.main import cli as seismic_cli
from seismicrna.core.rna.convert import db_to_ct
from seismicrna.core.seq.fasta import write_fasta
from seismicrna.core.seq.region import FULL_NAME
from seismicrna.core.seq.xna import DNA
from seismicrna.filterscan.main import run as run_filterscan
from seismicrna.filterscan.report import FilterScanReport
from seismicrna.filterscan.write import (
    _calc_tiles,
    _aggregate_pairs,
    _select_pairs,
    _calc_domains_from_pairs,
    _insert_domains_into_gaps,
    _expand_domains_into_gaps,
    _filter_domains_length,
)
from seismicrna.sim.params import run as sim_params
from seismicrna.sim.idmut import run as sim_idmut


class TestCalcTiles(ut.TestCase):
    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR)

    def tearDown(self):
        set_config(*self._config)

    def test_tile_min_overlap_25(self):
        result = _calc_tiles(41, 145, 60, 0.25)
        expect = [(41, 100), (86, 145)]
        self.assertEqual(result, expect)

    def test_tile_min_overlap_75(self):
        result = _calc_tiles(41, 145, 60, 0.75)
        expect = [(41, 100), (56, 115), (71, 130), (86, 145)]
        self.assertEqual(result, expect)

    def test_region_not_divisible(self):
        result = _calc_tiles(41, 170, 60, 0.25)
        expect = [(41, 100), (76, 135), (111, 170)]
        self.assertEqual(result, expect)

    def test_tile_length_larger(self):
        result = _calc_tiles(41, 49, 52, 0.9)
        expect = [(41, 49)]
        self.assertEqual(result, expect)

    def test_total_tile_length_1(self):
        result = _calc_tiles(41, 41, 52, 0.9)
        expect = [(41, 41)]
        self.assertEqual(result, expect)


class TestAggregatePairs(ut.TestCase):
    def test_zero(self):
        result = _aggregate_pairs([])
        expect = []
        self.assertListEqual(result, expect)

    def test_one(self):
        result = _aggregate_pairs([(3, 9)])
        expect = [(3, 9)]
        self.assertListEqual(result, expect)

    def test_two_serial(self):
        result = _aggregate_pairs([(3, 10), (10, 15)])
        expect = [(3, 15)]
        self.assertListEqual(result, expect)

    def test_two_nested(self):
        result = _aggregate_pairs([(3, 15), (9, 10)])
        expect = [(3, 15)]
        self.assertListEqual(result, expect)

    def test_two_disjoint(self):
        result = _aggregate_pairs([(3, 9), (10, 15)])
        expect = [(3, 9), (10, 15)]
        self.assertListEqual(result, expect)

    def test_three_serial(self):
        result = _aggregate_pairs([(3, 10), (9, 15), (15, 20)])
        expect = [(3, 20)]
        self.assertListEqual(result, expect)

    def test_three_nested(self):
        result = _aggregate_pairs([(3, 20), (7, 11), (14, 16)])
        expect = [(3, 20)]
        self.assertListEqual(result, expect)

    def test_three_disjoint(self):
        result = _aggregate_pairs([(3, 6), (7, 11), (14, 16)])
        expect = [(3, 6), (7, 11), (14, 16)]
        self.assertListEqual(result, expect)


class TestSelectPairs(ut.TestCase):
    def test_select_pairs(self):
        pairs = [(1, 11), (5, 15), (9, 19)]
        result = _select_pairs(pairs, 1, 19)
        self.assertListEqual(result, pairs)
        result = _select_pairs(pairs, 2, 18)
        self.assertListEqual(result, [(5, 15)])


class TestFilterDomainsLength(ut.TestCase):
    def test_filter_default_length(self):
        domains = [(5, 20), (25, 30)]
        result = _filter_domains_length(domains)
        expect = domains
        self.assertListEqual(result, expect)

    def test_filter_min_length(self):
        domains = [(5, 20), (25, 30)]
        result = _filter_domains_length(domains, min_length=6)
        expect = domains
        self.assertListEqual(result, expect)
        result = _filter_domains_length(domains, min_length=7)
        expect = [(5, 20)]
        self.assertListEqual(result, expect)
        result = _filter_domains_length(domains, min_length=16)
        expect = [(5, 20)]
        self.assertListEqual(result, expect)
        result = _filter_domains_length(domains, min_length=17)
        expect = []
        self.assertListEqual(result, expect)

    def test_filter_max_length(self):
        domains = [(5, 20), (25, 30)]
        result = _filter_domains_length(domains, max_length=5)
        expect = []
        self.assertListEqual(result, expect)
        result = _filter_domains_length(domains, max_length=6)
        expect = [(25, 30)]
        self.assertListEqual(result, expect)
        result = _filter_domains_length(domains, max_length=15)
        expect = [(25, 30)]
        self.assertListEqual(result, expect)
        result = _filter_domains_length(domains, max_length=16)
        expect = domains
        self.assertListEqual(result, expect)


class TestInsertRegionsIntoGaps(ut.TestCase):
    def test_zero(self):
        result = _insert_domains_into_gaps([], 3, 9)
        expect = [(3, 9)]
        self.assertListEqual(result, expect)

    def test_one(self):
        result = _insert_domains_into_gaps([(3, 9)], 3, 9)
        expect = [(3, 9)]
        self.assertListEqual(result, expect)
        result = _insert_domains_into_gaps([(4, 9)], 3, 9)
        expect = [(3, 3), (4, 9)]
        self.assertListEqual(result, expect)
        result = _insert_domains_into_gaps([(3, 8)], 3, 9)
        expect = [(3, 8), (9, 9)]
        self.assertListEqual(result, expect)
        result = _insert_domains_into_gaps([(4, 8)], 3, 9)
        expect = [(3, 3), (4, 8), (9, 9)]
        self.assertListEqual(result, expect)

    def test_two(self):
        result = _insert_domains_into_gaps([(2, 10), (11, 20)], 2, 20)
        expect = [(2, 10), (11, 20)]
        self.assertListEqual(result, expect)
        result = _insert_domains_into_gaps([(3, 10), (11, 20)], 2, 20)
        expect = [(2, 2), (3, 10), (11, 20)]
        self.assertListEqual(result, expect)
        result = _insert_domains_into_gaps([(2, 10), (12, 20)], 2, 20)
        expect = [(2, 10), (11, 11), (12, 20)]
        self.assertListEqual(result, expect)
        result = _insert_domains_into_gaps([(2, 10), (11, 19)], 2, 20)
        expect = [(2, 10), (11, 19), (20, 20)]
        self.assertListEqual(result, expect)
        result = _insert_domains_into_gaps([(3, 9), (11, 19)], 2, 20)
        expect = [(2, 2), (3, 9), (10, 10), (11, 19), (20, 20)]
        self.assertListEqual(result, expect)


class TestExpandRegionsIntoGaps(ut.TestCase):
    def test_zero(self):
        result = _expand_domains_into_gaps([], 3, 9)
        expect = []
        self.assertListEqual(result, expect)

    def test_one(self):
        expect = [(3, 9)]
        result = _expand_domains_into_gaps([(3, 9)], 3, 9)
        self.assertListEqual(result, expect)
        result = _expand_domains_into_gaps([(4, 9)], 3, 9)
        self.assertListEqual(result, expect)
        result = _expand_domains_into_gaps([(3, 8)], 3, 9)
        self.assertListEqual(result, expect)
        result = _expand_domains_into_gaps([(5, 7)], 3, 9)
        self.assertListEqual(result, expect)

    def test_two(self):
        result = _expand_domains_into_gaps([(2, 10), (11, 20)], 2, 20)
        expect = [(2, 10), (11, 20)]
        self.assertListEqual(result, expect)
        result = _expand_domains_into_gaps([(4, 9), (11, 18)], 2, 20)
        expect = [(2, 9), (10, 20)]
        self.assertListEqual(result, expect)
        result = _expand_domains_into_gaps([(5, 9), (12, 17)], 2, 20)
        expect = [(2, 10), (11, 20)]
        self.assertListEqual(result, expect)
        result = _expand_domains_into_gaps([(6, 6), (10, 10)], 2, 20)
        expect = [(2, 7), (8, 20)]
        self.assertListEqual(result, expect)


class ScanTestBase(ut.TestCase):
    """Shared simulation infrastructure for the scan tests."""

    SAMPLE = "test_sample"
    REFS = "test_refs"
    REF = "test_ref"
    PROFILE = "scan"

    # Folding domains of the reference sequence (each 60 nt).
    DOMAINS = [
        (
            "TGACGAACAACGTGTTTGTGAACCATATAGGTAAACGCTGAATGCGTTCGCGCGGAGGGT",
            [
                "..(((((((...)))))))..(((......((.(((((.....))))).))......)))",
                "...(((...(((((((((((.(((.....)))...))).))))))))))).((.....))",
            ],
        ),
        (
            "TTTGCAGGAAGATGGTCAACTCTACACCTAGTTTTTACCAGTCCACAAGAGTTTGAACTG",
            [".(((..(((...((((.((..(((....)))..)).)))).))).))).(((....)))."],
        ),
        (
            "GTGCCTTAACCTGAGTACGCCCATATCATGGGAGACATTACAACTCAAATTCTAGGTGTG",
            [
                "..((((.....(((((...(((((...)))))..........)))))......))))...",
                "((((.(((...)))))))...((((((.((((..................))))))))))",
            ],
        ),
    ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._config = None
        self._tmpdir = None

    @property
    def sim_dir(self):
        if self._tmpdir is None:
            return None
        return Path(self._tmpdir.name)

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR, log_file_path=None, exit_on_error=True)
        self._tmpdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self._tmpdir.cleanup()
        self._tmpdir = None
        set_config(*self._config)

    def sim_data(self, domain_nums: list[int], read_length: int, seed: int):
        # Assemble and write the reference sequence.
        domains = dict(self.DOMAINS[m] for m in domain_nums)
        refseq = DNA("".join(domains.keys()))
        refs_dir = self.sim_dir.joinpath("refs")
        refs_dir.mkdir()
        fasta = refs_dir.joinpath(f"{self.REFS}.fa")
        write_fasta(fasta, [(self.REF, refseq)])
        # Assemble and write the secondary structures.
        structures = list(map("".join, product(*domains.values())))
        param_dir = self.sim_dir.joinpath("params", self.REF, FULL_NAME)
        param_dir.mkdir(parents=True)
        db_file = param_dir.joinpath(f"{self.PROFILE}.db")
        with open(db_file, "x") as f:
            for i, struct in enumerate(structures):
                if i == 0:
                    f.write(f">structure0\n{refseq.tr()}\n{struct}\n")
                else:
                    f.write(f">structure{i}\n{struct}\n")
        ct_file = db_to_ct(db_file)
        # Simulate data.
        sim_params(
            ct_file=[ct_file],
            # Make pmut_unpaired for A and C large (17%) so that
            # most reads get at least two mutations despite being
            # short and are thus useful for clustering.
            pmut_unpaired=[("am", 1 / 6), ("cm", 1 / 6)],
            # Make all reads the same length.
            length_fmean=(read_length / len(refseq)),
            length_fvar=0.0,
            # Make clust_conc very large so that the proportion
            # of each cluster is approximately equal, which makes
            # clustering easier.
            clust_conc=1000.0,
            seed=seed,
        )
        idmut_dirs = sim_idmut(
            param_dir=[param_dir],
            sample=self.SAMPLE,
            profile_name=self.PROFILE,
            num_reads=200000,
            paired_end=False,
            brotli_level=0,
            seed=seed,
        )
        return idmut_dirs


class TestCalcDomainsFromPairsSimulatedData(ut.TestCase):
    """Regression test using simulated DMS-MaPseq data of a 1200-nt RNA
    with a complex, multi-domain structure. Key expectations:

    - Pair (148, 196) must NOT appear in any returned domain: it was
      the original defect that motivated the algorithm overhaul.
    - Domains covering roughly 276–588 and 909–992 must be present,
      matching the known structural domains visible in the plot.
    """

    PAIRS_CSV = Path(__file__).parent / "long_rna_10_pairs.csv"
    PAIR_FDR = 0.05
    MIN_MUT_GAP = 4
    MIN_PAIRS = 2
    PAIR_DISTANCE_PERCENTILE = 95.0
    ENDPOINT_WINDOW = 2
    MIN_NEARBY_PAIRS = 1
    TOTAL_END5 = 1
    TOTAL_END3 = 1200

    @classmethod
    def setUpClass(cls):
        import csv

        pairs = []
        with open(cls.PAIRS_CSV) as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                pairs.append((int(row[0]), int(row[1])))
        cls._pairs = sorted(set(pairs))

    def _run(self):
        return _calc_domains_from_pairs(
            self._pairs,
            self.PAIR_FDR,
            self.MIN_MUT_GAP,
            self.MIN_PAIRS,
            self.TOTAL_END5,
            self.TOTAL_END3,
            self.PAIR_DISTANCE_PERCENTILE,
            self.ENDPOINT_WINDOW,
            self.MIN_NEARBY_PAIRS,
        )

    def test_expected_domains(self):
        result = self._run()
        self.assertListEqual(result, [(273, 591), (663, 796), (802, 896), (903, 992)])

    def test_spurious_pair_not_in_any_domain(self):
        # Pair (148, 196) lies in a sparse, weakly-supported transition
        # zone and must be excluded from every domain.
        result = self._run()
        self.assertFalse(
            any(a <= 148 and 196 <= b for a, b in result),
            f"Pair (148, 196) is unexpectedly inside a domain: {result}",
        )


class TestFilterScan(ScanTestBase):
    """Test that filterscan identifies domains without clustering."""

    def run_filterscan_check(
        self, idmut_dirs: list[Path], expect_regions: list[tuple[int, int]], **kwargs
    ):
        filterscan_dirs = run_filterscan(
            idmut_dirs,
            brotli_level=0,
            filter_pos_table=False,
            filter_read_table=False,
            **kwargs,
        )
        domains = list()
        for filterscan_dir in filterscan_dirs:
            for report_file in path.find_files_chain(
                [filterscan_dir], FilterScanReport.get_path_seg_types()
            ):
                report = FilterScanReport.load(report_file)
                for reg5, reg3 in report.get_field(DomainCoordsF):
                    domains.append((int(reg5), int(reg3)))
        for exp5, exp3 in expect_regions:
            expect_length = exp3 - exp5 + 1
            for reg5, reg3 in domains:
                overlap = max(min(reg3, exp3) - max(reg5, exp5), 0)
                if overlap >= expect_length / 2:
                    break
            else:
                raise ValueError(
                    f"Expected region {exp5, exp3} does not overlap at "
                    f"least 50% of any detected domain among {sorted(domains)}"
                )

    def test_domains012_read180(self):
        idmut_dirs = self.sim_data([0, 1, 2], 180, seed=0)
        self.run_filterscan_check(idmut_dirs, [(1, 60), (121, 180)])

    def test_domains012_read120(self):
        idmut_dirs = self.sim_data([0, 1, 2], 120, seed=0)
        self.run_filterscan_check(idmut_dirs, [(1, 60), (121, 180)])

    def test_domains012_read60(self):
        idmut_dirs = self.sim_data([0, 1, 2], 60, seed=0)
        self.run_filterscan_check(idmut_dirs, [(1, 60), (121, 180)])

    def test_domains02_read60(self):
        idmut_dirs = self.sim_data([0, 2], 60, seed=0)
        self.run_filterscan_check(idmut_dirs, [(1, 60), (61, 120)])

    def test_domains012_read180_cli(self):
        idmut_dirs = self.sim_data([0, 1, 2], 180, seed=0)
        runner = CliRunner()
        args = (
            ["-qq", "--exit-on-error", "filterscan"]
            + [str(d) for d in idmut_dirs]
            + ["--brotli-level", "0"]
        )
        result = runner.invoke(seismic_cli, args, catch_exceptions=False)
        self.assertEqual(result.exit_code, 0, msg=result.output)
        set_config(verbosity=Level.ERROR, log_file_path=None, exit_on_error=True)


if __name__ == "__main__":
    ut.main(verbosity=2)
