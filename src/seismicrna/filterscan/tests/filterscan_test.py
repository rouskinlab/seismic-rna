import tempfile
import unittest as ut
from itertools import product
from pathlib import Path

from click.testing import CliRunner

from seismicrna.core import path
from seismicrna.core.batch.confusion import POSITION_A, POSITION_B
from seismicrna.core.error import OutOfBoundsError
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
    N_COL,
    PHI_COL,
    CHI_SQUARE_COL,
    SIGNIFICANT_COL,
    _calc_tiles,
    _build_banded_table,
    _calc_domains_by_dp_segmentation,
    _insert_domains_into_gaps,
    _expand_domains_into_gaps,
    _filter_domains_length,
)
from seismicrna.sim.params import run as sim_params
from seismicrna.sim.idmut import run as sim_idmut


def _make_table(rows: list[tuple[int, int, float, float, bool]]):
    """Build a per-pair table like ``_find_correlated_pairs`` returns,
    from an iterable of (pos_a, pos_b, n, phi, significant)."""
    import pandas as pd

    if not rows:
        pos_a, pos_b, ns, phis, sig = (), (), (), (), ()
    else:
        pos_a, pos_b, ns, phis, sig = zip(*rows)
    index = pd.MultiIndex.from_arrays([pos_a, pos_b], names=[POSITION_A, POSITION_B])
    chi_square = [n * phi**2 for n, phi in zip(ns, phis)]
    return pd.DataFrame(
        {N_COL: ns, PHI_COL: phis, CHI_SQUARE_COL: chi_square, SIGNIFICANT_COL: sig},
        index=index,
    )


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


class TestBuildBandedTable(ut.TestCase):
    def test_empty(self):
        result = _build_banded_table([], band_width=0)
        self.assertEqual(len(result.index), 0)
        self.assertListEqual(
            list(result.columns), [N_COL, PHI_COL, CHI_SQUARE_COL, SIGNIFICANT_COL]
        )

    def test_dedup_keeps_max_n(self):
        # The same pair observed in two overlapping tiles: keep the
        # observation with the greater coverage (N).
        table1 = _make_table([(1, 5, 10.0, 0.5, False)])
        table2 = _make_table([(1, 5, 50.0, 0.2, True)])
        result = _build_banded_table([table1, table2], band_width=0)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result[N_COL].iloc[0], 50.0)
        self.assertEqual(bool(result[SIGNIFICANT_COL].iloc[0]), True)

    def test_band_filter(self):
        rows = [(1, 5, 10.0, 0.5, False), (1, 50, 10.0, 0.5, False)]
        table = _make_table(rows)
        result = _build_banded_table([table], band_width=10)
        self.assertListEqual(result.index.to_list(), [(1, 5)])

    def test_band_width_zero_applies_no_extra_cap(self):
        rows = [(1, 5, 10.0, 0.5, False), (1, 50, 10.0, 0.5, False)]
        table = _make_table(rows)
        result = _build_banded_table([table], band_width=0)
        self.assertListEqual(result.index.to_list(), [(1, 5), (1, 50)])

    def test_nan_chi_square_dropped(self):
        table = _make_table([(1, 5, 0.0, float("nan"), False)])
        result = _build_banded_table([table], band_width=0)
        self.assertEqual(len(result.index), 0)


def _rows_over_region(
    total_end5: int, total_end3: int, band: int, is_significant, n_cov: float = 20.0
) -> list[tuple[int, int, float, float, bool]]:
    """Generate (i, j, N, phi, significant) rows for every pair in
    [total_end5, total_end3] with ``j - i <= band``, with significance
    decided by ``is_significant(i, j)``. ``N``/``phi`` are placeholders
    (unused by the insulation functions, which read only ``N`` for the
    coverage filter and ``Significant`` for the density)."""
    rows = []
    for i in range(total_end5, total_end3 + 1):
        for j in range(i + 1, min(i + band, total_end3) + 1):
            rows.append((i, j, n_cov, 0.0, is_significant(i, j)))
    return rows


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


class TestCalcDomainsByDpSegmentation(ut.TestCase):
    """End-to-end tests of the DP block-diagonal domain caller."""

    TOTAL_END5, TOTAL_END3 = 1, 50
    DOMAIN_A = (1, 15)
    DOMAIN_B = (36, 50)
    BAND = 10
    N_COV = 2000.0
    BIC_MULTIPLIER = 1.0

    @classmethod
    def _in_domain(cls, i: int, j: int) -> bool:
        a_lo, a_hi = cls.DOMAIN_A
        b_lo, b_hi = cls.DOMAIN_B
        return (a_lo <= i and j <= a_hi) or (b_lo <= i and j <= b_hi)

    def _build_table(self):
        rows = _rows_over_region(
            self.TOTAL_END5,
            self.TOTAL_END3,
            self.BAND,
            self._in_domain,
            n_cov=self.N_COV,
        )
        return _make_table(rows)

    def test_domains_recovered_and_not_merged(self):
        table = self._build_table()
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=1,
        )
        self.assertTrue(domains, "Expected at least one domain to be called")
        for start, end in domains:
            self.assertFalse(
                start <= self.DOMAIN_A[0] and self.DOMAIN_B[1] <= end,
                f"Domain {(start, end)} spans both clusters",
            )
        for lo, hi in (self.DOMAIN_A, self.DOMAIN_B):
            expect_length = hi - lo + 1
            self.assertTrue(
                any(
                    max(min(end, hi) - max(start, lo) + 1, 0) >= expect_length / 2
                    for start, end in domains
                ),
                f"Domain {(lo, hi)} not recovered in {domains}",
            )

    def test_uniform_density_has_no_background_contrast(self):
        # Every pair in the whole scanned region is significant: with
        # no sparser background to contrast against, nothing can be
        # denser than background, so no domain is called (this is not
        # a false negative: a "domain" is meaningless without
        # surrounding background to distinguish it from).
        rows = _rows_over_region(
            self.TOTAL_END5, self.TOTAL_END3, self.BAND, lambda i, j: True
        )
        table = _make_table(rows)
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=1,
        )
        self.assertEqual(domains, [])

    def test_sparse_region_has_no_domain(self):
        rows = _rows_over_region(
            self.TOTAL_END5, self.TOTAL_END3, self.BAND, lambda i, j: False
        )
        table = _make_table(rows)
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=1,
        )
        self.assertEqual(domains, [])

    def test_min_pair_coverage_excludes_low_coverage_pairs(self):
        table = self._build_table()
        # Below the pairs' actual coverage (N_COV): pairs pass through
        # unfiltered, and the domains are detected as usual.
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=1,
        )
        self.assertTrue(domains)
        # Above the pairs' actual coverage: every pair is filtered out,
        # so no domain can be detected.
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=int(self.N_COV) + 1,
        )
        self.assertEqual(domains, [])

    def test_min_pair_coverage_must_be_at_least_one(self):
        table = self._build_table()
        with self.assertRaises(OutOfBoundsError):
            _calc_domains_by_dp_segmentation(
                table,
                self.TOTAL_END5,
                self.TOTAL_END3,
                bic_multiplier=self.BIC_MULTIPLIER,
                min_pair_coverage=0,
            )

    def test_bic_multiplier_must_be_nonnegative(self):
        table = self._build_table()
        with self.assertRaises(OutOfBoundsError):
            _calc_domains_by_dp_segmentation(
                table,
                self.TOTAL_END5,
                self.TOTAL_END3,
                bic_multiplier=-1.0,
                min_pair_coverage=1,
            )

    def test_masked_pairs_never_contribute(self):
        # A "masked" pocket in the gap between the two real domains:
        # pairs flagged significant, but at coverage far below
        # min_pair_coverage, as if these positions had too little read
        # depth to trust (even though, if not excluded, they would
        # look like a dense spurious third domain). The filter must
        # make them invisible, not merely down-weighted: calling with
        # and without the pocket present must give byte-identical
        # results.
        def build_rows(masked_pocket: bool):
            rows = []
            for i in range(self.TOTAL_END5, self.TOTAL_END3 + 1):
                for j in range(i, min(i + self.BAND, self.TOTAL_END3) + 1):
                    if masked_pocket and 16 <= i and j <= 35:
                        rows.append((i, j, 1.0, 0.0, True))
                    else:
                        rows.append((i, j, self.N_COV, 0.0, self._in_domain(i, j)))
            return rows

        table_baseline = _make_table(build_rows(masked_pocket=False))
        table_with_masked = _make_table(build_rows(masked_pocket=True))

        # Excludes the N=1 masked pocket while keeping the N=N_COV domains.
        min_pair_coverage = 100
        domains_baseline = _calc_domains_by_dp_segmentation(
            table_baseline,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=min_pair_coverage,
        )
        domains_with_masked = _calc_domains_by_dp_segmentation(
            table_with_masked,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=min_pair_coverage,
        )
        self.assertTrue(domains_baseline, "Expected at least one domain to be called")
        self.assertEqual(domains_baseline, domains_with_masked)

    def test_domain_longer_than_band_recovered(self):
        # Stress the far/near split in _triangle_sum_banded: a domain
        # longer than the band width, so most of a candidate block's
        # own triangle sum must come from the "far" (full-row) branch,
        # not just the "near" (partial-row) branch close to its end.
        total_end5, total_end3 = 1, 80
        band = 15
        domain = (20, 60)

        def is_significant(i: int, j: int) -> bool:
            return domain[0] <= i and j <= domain[1]

        rows = _rows_over_region(
            total_end5, total_end3, band, is_significant, n_cov=self.N_COV
        )
        table = _make_table(rows)
        domains = _calc_domains_by_dp_segmentation(
            table,
            total_end5,
            total_end3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=1,
        )
        expect_length = domain[1] - domain[0] + 1
        self.assertTrue(
            any(
                max(min(end, domain[1]) - max(start, domain[0]) + 1, 0)
                >= expect_length / 2
                for start, end in domains
            ),
            f"Domain {domain} not recovered in {domains}",
        )


class TestDpSegmentationAntiOverSplit(ut.TestCase):
    """The core design property of the DP block caller: splitting a
    domain forfeits the evidence of every significant pair that would
    have to move to the background. Two adjacent sub-domains [1,50]
    and [51,100], with a sparse background beyond [100,150] to give
    the model something to contrast against, are called as ONE domain
    when the pairs crossing between them (a in [1,50], b in [51,100])
    are also dense, but as TWO separate domains when that crossing
    region is sparse -- the exact mechanism discussed as the reason
    the DP model keeps a domain bridged by many crossing pairs whole
    instead of over-splitting it."""

    TOTAL_END5, TOTAL_END3 = 1, 150
    SUB_A = (1, 50)
    SUB_B = (51, 100)
    BAND = 149
    N_COV = 2000.0
    BIC_MULTIPLIER = 1.0

    @classmethod
    def _in_sub(cls, lo_hi, i: int, j: int) -> bool:
        lo, hi = lo_hi
        return lo <= i and j <= hi

    def _build_table(self, dense_bridge: bool):
        def is_significant(i: int, j: int) -> bool:
            if self._in_sub(self.SUB_A, i, j) or self._in_sub(self.SUB_B, i, j):
                return True
            if (
                dense_bridge
                and self.SUB_A[0] <= i <= self.SUB_A[1]
                and self.SUB_B[0] <= j <= self.SUB_B[1]
            ):
                return True
            return False

        rows = _rows_over_region(
            self.TOTAL_END5,
            self.TOTAL_END3,
            self.BAND,
            is_significant,
            n_cov=self.N_COV,
        )
        return _make_table(rows)

    def test_dense_bridge_keeps_domain_whole(self):
        table = self._build_table(dense_bridge=True)
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=1,
        )
        self.assertTrue(
            any(
                start <= self.SUB_A[0] and self.SUB_B[1] <= end
                for start, end in domains
            ),
            f"Expected one domain spanning both sub-clusters in {domains}",
        )

    def test_sparse_bridge_splits_domain(self):
        table = self._build_table(dense_bridge=False)
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=1,
        )
        self.assertFalse(
            any(
                start <= self.SUB_A[0] and self.SUB_B[1] <= end
                for start, end in domains
            ),
            f"Did not expect one domain spanning both sub-clusters in {domains}",
        )
        for lo, hi in (self.SUB_A, self.SUB_B):
            expect_length = hi - lo + 1
            self.assertTrue(
                any(
                    max(min(end, hi) - max(start, lo) + 1, 0) >= expect_length / 2
                    for start, end in domains
                ),
                f"Sub-domain {(lo, hi)} not recovered in {domains}",
            )


class TestPostHocMergeAdjacentBlocks(ut.TestCase):
    """The DP's exact partition optimum can decline to merge two
    blocks whose "in-between" material -- a gap plus any pairs
    crossing from one block into the other -- carries real signal but
    is weaker than either block's own density (see
    ``_merge_adjacent_blocks`` for why the DP itself never makes this
    trade: joining always scores lower there, by construction, since
    the DP already considered and rejected it). The post-hoc merge
    pass recovers a single flat domain in that case, matching what a
    person sees by eye as one region with a graded density change --
    while still leaving a genuinely silent bridge split, exactly as
    ``TestDpSegmentationAntiOverSplit.test_sparse_bridge_splits_domain``
    requires."""

    TOTAL_END5, TOTAL_END3 = 1, 60
    DOMAIN_A = (1, 15)
    DOMAIN_B = (31, 45)
    BAND = 8
    N_COV = 2000.0
    BIC_MULTIPLIER = 1.0

    def _build_table(self, bridge_frac: float, bg_frac: float):
        # A and B are fully dense. Every observable pair that touches
        # neither A nor B but still falls within domain B's end (the
        # gap between A and B, and short-range crossings into it --
        # BAND=8 is too narrow for any direct A-to-B pair to exist at
        # all) is the "bridge" tier. Everything past domain B sets the
        # outer background rate. Marking every 1/frac-th pair by
        # enumeration order (not a fixed (i, j) pattern) gives an
        # exact, robust rate, unlike an (i + j) % k rule, which can
        # alias unpredictably over a small, non-uniform (i, j) set.
        rows = []
        counters: dict[str, int] = {}

        def tier_significant(tier: str, frac: float) -> bool:
            count = counters.get(tier, 0)
            counters[tier] = count + 1
            period = max(1, round(1 / frac))
            return count % period == 0

        for i in range(self.TOTAL_END5, self.TOTAL_END3 + 1):
            for j in range(i + 1, min(i + self.BAND, self.TOTAL_END3) + 1):
                if self.DOMAIN_A[0] <= i and j <= self.DOMAIN_A[1]:
                    sig = True
                elif self.DOMAIN_B[0] <= i and j <= self.DOMAIN_B[1]:
                    sig = True
                elif j <= self.DOMAIN_B[1]:
                    sig = tier_significant("bridge", bridge_frac)
                else:
                    sig = tier_significant("background", bg_frac)
                rows.append((i, j, self.N_COV, 0.0, sig))
        return _make_table(rows)

    def test_weak_but_real_bridge_merges_split_domains(self):
        # The bridge (10% significant) is real signal, well above the
        # outer background (2%) -- but weaker than either domain (100%
        # significant), so the DP's exact optimum still declines to
        # fold it into either domain's own block, leaving it split as
        # [(1, 15), (31, 45)]. The merge pass should still join them
        # into one flat domain.
        table = self._build_table(bridge_frac=0.1, bg_frac=0.02)
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            bic_multiplier=self.BIC_MULTIPLIER,
            min_pair_coverage=1,
        )
        self.assertEqual(domains, [(self.DOMAIN_A[0], self.DOMAIN_B[1])])


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
