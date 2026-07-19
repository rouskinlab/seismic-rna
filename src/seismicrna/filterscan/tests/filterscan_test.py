import tempfile
import unittest as ut
from itertools import product
from pathlib import Path

import numpy as np
from click.testing import CliRunner

from seismicrna.core import path
from seismicrna.core.batch.confusion import (
    POSITION_A,
    POSITION_B,
    calc_bh_adjusted_pvals,
)
from seismicrna.core.error import IncompatibleValuesError, OutOfBoundsError
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
    CHI_SQUARE_COL,
    CONFUSION_COLS,
    NEITHER_COL,
    ONLY_A_COL,
    ONLY_B_COL,
    BOTH_COL,
    _calc_tiles,
    _build_banded_table,
    _write_pairs_with_confusion,
    _block_score,
    _calc_block_pvalue_cutoff,
    _calc_domains_by_dp_segmentation,
    _cut_crossing_scores,
    _merge_connected_blocks,
    _pair_band_row_cumsum,
    _triangle_sum_banded,
    _insert_domains_into_gaps,
    _expand_domains_into_gaps,
    _filter_domains_length,
)
from seismicrna.sim.params import run as sim_params
from seismicrna.sim.idmut import run as sim_idmut


def _make_table(rows: list[tuple[int, int, float, float]]):
    """Build a per-pair table like ``_confusion_to_table`` returns, from an
    iterable of (pos_a, pos_b, n, chi_square).

    The four confusion cells are synthesized as a simple, fixed split of
    ``n`` (not tied to ``chi_square``) purely so ``_analyzed_pairs_mask``
    (which needs them to compute the independence-expected both-mutated
    count) has columns to read; tests that care about the exact confusion
    counts (e.g. ``min_expect_both`` filtering) build their own table."""
    import pandas as pd

    if not rows:
        pos_a, pos_b, ns, chi_square = (), (), (), ()
    else:
        pos_a, pos_b, ns, chi_square = zip(*rows)
    index = pd.MultiIndex.from_arrays([pos_a, pos_b], names=[POSITION_A, POSITION_B])
    only_a = [n / 2 for n in ns]
    only_b = [n / 2 for n in ns]
    both = [n / 4 for n in ns]
    neither = [n - (a + b - ab) for n, a, b, ab in zip(ns, only_a, only_b, both)]
    return pd.DataFrame(
        {
            N_COL: ns,
            CHI_SQUARE_COL: chi_square,
            NEITHER_COL: neither,
            ONLY_A_COL: only_a,
            ONLY_B_COL: only_b,
            BOTH_COL: both,
        },
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
            list(result.columns), [N_COL, CHI_SQUARE_COL, *CONFUSION_COLS]
        )

    def test_dedup_keeps_max_n(self):
        # The same pair observed in two overlapping tiles: keep the
        # observation with the greater coverage (N).
        table1 = _make_table([(1, 5, 10.0, 2.5)])
        table2 = _make_table([(1, 5, 50.0, 2.0)])
        result = _build_banded_table([table1, table2], band_width=0)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result[N_COL].iloc[0], 50.0)
        self.assertEqual(result[CHI_SQUARE_COL].iloc[0], 2.0)

    def test_band_filter(self):
        rows = [(1, 5, 10.0, 2.5), (1, 50, 10.0, 2.5)]
        table = _make_table(rows)
        result = _build_banded_table([table], band_width=10)
        self.assertListEqual(result.index.to_list(), [(1, 5)])

    def test_band_width_zero_applies_no_extra_cap(self):
        rows = [(1, 5, 10.0, 2.5), (1, 50, 10.0, 2.5)]
        table = _make_table(rows)
        result = _build_banded_table([table], band_width=0)
        self.assertListEqual(result.index.to_list(), [(1, 5), (1, 50)])

    def test_nan_chi_square_dropped(self):
        table = _make_table([(1, 5, 0.0, float("nan"))])
        result = _build_banded_table([table], band_width=0)
        self.assertEqual(len(result.index), 0)


class TestCalcBlockPvalueCutoff(ut.TestCase):
    """The analytic Benjamini-Hochberg p-value cutoff (``chi2.sf(chi2_sum,
    df=n)`` per candidate block, no null replicates needed) that gates
    ``_dp_segment_blocks`` alongside ``_block_score``."""

    @staticmethod
    def _row_cums(rows, total_end5, total_end3):
        table = _make_table(rows)
        return _pair_band_row_cumsum(
            table, total_end5, total_end3, value_col=CHI_SQUARE_COL
        )

    def test_empty_table_returns_no_survivor(self):
        import numpy as np

        # No observable pairs at all: an all-zero grid of the shape
        # _pair_band_row_cumsum would produce (an empty table's index has no
        # dtype to infer a real max_gap from, so build the zero grid
        # directly rather than route an empty table through it).
        n_positions, max_gap = 60, 0
        obs_row_cum = np.zeros((n_positions, max_gap + 2), dtype=np.int64)
        chi2_row_cum = np.zeros((n_positions, max_gap + 2), dtype=float)
        cutoff = _calc_block_pvalue_cutoff(
            obs_row_cum, chi2_row_cum, n_positions, max_gap, detect_fdr=0.1
        )
        self.assertEqual(cutoff, -1.0)

    def test_all_null_like_has_no_survivor(self):
        rows = _rows_over_region(1, 60, 10, lambda i, j: False)
        obs_row_cum, chi2_row_cum, n_positions, max_gap = self._row_cums(rows, 1, 60)
        cutoff = _calc_block_pvalue_cutoff(
            obs_row_cum, chi2_row_cum, n_positions, max_gap, detect_fdr=0.1
        )
        self.assertEqual(cutoff, -1.0)

    def test_strong_domain_survives_against_a_mostly_null_pool(self):
        # A real domain amid a much larger null-like background: the domain's
        # own candidate block must clear BH correction against every other
        # (mostly null) candidate the scan considers.
        rows = _rows_over_region(1, 300, 10, lambda i, j: 20 <= i and j <= 40)
        obs_row_cum, chi2_row_cum, n_positions, max_gap = self._row_cums(rows, 1, 300)
        cutoff = _calc_block_pvalue_cutoff(
            obs_row_cum, chi2_row_cum, n_positions, max_gap, detect_fdr=0.1
        )
        self.assertGreater(cutoff, -1.0)
        # The domain's own exact block (0-indexed [19, 39]) must itself be
        # at or below the returned cutoff.
        n_se = float(_triangle_sum_banded(obs_row_cum, 39, max_gap)[19])
        chi2_se = float(_triangle_sum_banded(chi2_row_cum, 39, max_gap)[19])
        from scipy.stats import chi2 as chi2dist

        domain_pvalue = chi2dist.sf(chi2_se, df=n_se)
        self.assertLessEqual(domain_pvalue, cutoff)


class TestWritePairsWithConfusion(ut.TestCase):
    """pairs.csv must carry each written pair's 2x2 confusion-matrix counts
    (chi-square, and every other derived quantity, is omitted since it is
    recomputable from the counts alone)."""

    @staticmethod
    def _table(rows):
        # rows: (pos_a, pos_b, neither, only_a, only_b, both)
        import pandas as pd

        pa, pb, ne, oa, ob, bo = zip(*rows)
        index = pd.MultiIndex.from_arrays([pa, pb], names=[POSITION_A, POSITION_B])
        return pd.DataFrame(
            {
                NEITHER_COL: list(ne),
                ONLY_A_COL: list(oa),
                ONLY_B_COL: list(ob),
                BOTH_COL: list(bo),
            },
            index=index,
        )

    def test_writes_confusion_counts(self):
        import pandas as pd

        pos_table = self._table([(1, 5, 900, 30, 20, 50), (2, 9, 800, 60, 40, 100)])
        with tempfile.TemporaryDirectory() as tmp:
            csv_file = Path(tmp) / "pairs.csv"
            _write_pairs_with_confusion(pos_table, csv_file)
            out = pd.read_csv(csv_file)
        # Exactly the position columns and the four confusion cells;
        # chi-square (recomputable from the counts) is omitted.
        self.assertListEqual(
            list(out.columns), [POSITION_A, POSITION_B, *CONFUSION_COLS]
        )
        self.assertNotIn(CHI_SQUARE_COL, out.columns)
        self.assertListEqual(
            list(zip(out[POSITION_A], out[POSITION_B])), [(1, 5), (2, 9)]
        )
        self.assertListEqual(list(out[BOTH_COL]), [50, 100])
        self.assertListEqual(list(out[ONLY_A_COL]), [30, 60])
        # N is recoverable from the written cells.
        recovered_n = (
            out[NEITHER_COL] + out[ONLY_A_COL] + out[ONLY_B_COL] + out[BOTH_COL]
        )
        self.assertListEqual(list(recovered_n), [1000, 1000])


def _rows_over_region(
    total_end5: int,
    total_end3: int,
    band: int,
    is_domain,
    n_cov: float = 2000.0,
    chi2_domain: float = 30.0,
    chi2_bg: float = 1.0,
) -> list[tuple[int, int, float, float]]:
    """Generate (i, j, N, chi_square) rows for every pair in
    [total_end5, total_end3] with ``j - i <= band``. A pair inside a domain
    (``is_domain(i, j)`` true) gets a strongly correlated chi-square
    (``chi2_domain``); the rest get a null-like chi-square (``chi2_bg``, the
    ~1 mean of ``chi^2(1)``), so the domain caller's variance-inflation
    score sees real, coverage-scaled correlation signal."""
    rows = []
    for i in range(total_end5, total_end3 + 1):
        for j in range(i + 1, min(i + band, total_end3) + 1):
            chi2 = chi2_domain if is_domain(i, j) else chi2_bg
            rows.append((i, j, n_cov, chi2))
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
    DOMAIN_FDR = 0.1

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
            min_pair_coverage=1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
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

    def test_uniform_correlation_is_one_domain(self):
        # Every pair in the whole scanned region is strongly correlated:
        # unlike the old density caller (which needed a sparser background to
        # contrast against and called nothing), the distributional score sees
        # the entire region as drawn from a wider-than-null distribution, so
        # it is one domain spanning the whole region.
        rows = _rows_over_region(
            self.TOTAL_END5, self.TOTAL_END3, self.BAND, lambda i, j: True
        )
        table = _make_table(rows)
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
        )
        self.assertEqual(domains, [(self.TOTAL_END5, self.TOTAL_END3)])

    def test_sparse_region_has_no_domain(self):
        rows = _rows_over_region(
            self.TOTAL_END5, self.TOTAL_END3, self.BAND, lambda i, j: False
        )
        table = _make_table(rows)
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
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
            min_pair_coverage=1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
        )
        self.assertTrue(domains)
        # Above the pairs' actual coverage: every pair is filtered out,
        # so no domain can be detected.
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=int(self.N_COV) + 1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
        )
        self.assertEqual(domains, [])

    def test_min_expect_both_excludes_pairs_with_low_expected_count(self):
        # Same domain/background chi-squares as _build_table, but with the
        # confusion cells built so the independence-expected both-mutated
        # count (a * b / N) is far below 5 despite ample coverage (N_COV):
        # standard chi-square practice requires expected cell counts >= ~5.
        import pandas as pd

        def build_table(low_expected: bool):
            rows = _rows_over_region(
                self.TOTAL_END5,
                self.TOTAL_END3,
                self.BAND,
                self._in_domain,
                n_cov=self.N_COV,
            )
            pos_a, pos_b, ns, chi2 = zip(*rows)
            if low_expected:
                only_a = [2.0] * len(rows)
                only_b = [2.0] * len(rows)
                both = [1.0] * len(rows)
            else:
                only_a = [n / 2 for n in ns]
                only_b = [n / 2 for n in ns]
                both = [n / 4 for n in ns]
            neither = [
                n - (a + b - ab) for n, a, b, ab in zip(ns, only_a, only_b, both)
            ]
            index = pd.MultiIndex.from_arrays(
                [pos_a, pos_b], names=[POSITION_A, POSITION_B]
            )
            return pd.DataFrame(
                {
                    N_COL: ns,
                    CHI_SQUARE_COL: chi2,
                    NEITHER_COL: neither,
                    ONLY_A_COL: only_a,
                    ONLY_B_COL: only_b,
                    BOTH_COL: both,
                },
                index=index,
            )

        # Low expected-both: every pair is excluded, so no domain is called.
        table = build_table(low_expected=True)
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            min_expect_both=5,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
        )
        self.assertEqual(domains, [], "Pairs with expected-both < 5 must be excluded")
        # The same chi-squares with ample expected-both: recovered as usual.
        table = build_table(low_expected=False)
        domains = _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            min_expect_both=5,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
        )
        self.assertTrue(domains)

    def test_min_pair_coverage_must_be_at_least_one(self):
        table = self._build_table()
        with self.assertRaises(OutOfBoundsError):
            _calc_domains_by_dp_segmentation(
                table,
                self.TOTAL_END5,
                self.TOTAL_END3,
                min_pair_coverage=0,
                min_expect_both=0,
                detect_fdr=self.DOMAIN_FDR,
                merge_fdr=self.DOMAIN_FDR,
            )

    def test_masked_pairs_never_contribute(self):
        # A "masked" pocket in the gap between the two real domains:
        # strongly correlated pairs, but at coverage far below
        # min_pair_coverage, as if these positions had too little read
        # depth to trust (even though, if not excluded, they would look
        # like a spurious third domain). The filter must make them
        # invisible, not merely down-weighted: calling with and without
        # the pocket present must give byte-identical results -- including
        # every candidate block's own locally-fit variance-inflation shape
        # r, so the pocket is excluded from scoring too (at min_pair_coverage).
        min_pair_coverage = 100

        def build_rows(masked_pocket: bool):
            rows = []
            for i in range(self.TOTAL_END5, self.TOTAL_END3 + 1):
                for j in range(i + 1, min(i + self.BAND, self.TOTAL_END3) + 1):
                    if masked_pocket and 16 <= i and j <= 35:
                        # N=1 (below min_pair_coverage), spuriously strong.
                        rows.append((i, j, 1.0, 30.0))
                    else:
                        chi2 = 30.0 if self._in_domain(i, j) else 1.0
                        rows.append((i, j, self.N_COV, chi2))
            return rows

        table_baseline = _make_table(build_rows(masked_pocket=False))
        table_with_masked = _make_table(build_rows(masked_pocket=True))
        domains_baseline = _calc_domains_by_dp_segmentation(
            table_baseline,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=min_pair_coverage,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
        )
        domains_with_masked = _calc_domains_by_dp_segmentation(
            table_with_masked,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=min_pair_coverage,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
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
            min_pair_coverage=1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
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

    def test_out_of_order_ends_raise(self):
        table = self._build_table()
        with self.assertRaises(IncompatibleValuesError):
            _calc_domains_by_dp_segmentation(
                table,
                self.TOTAL_END3,
                self.TOTAL_END5,
                min_pair_coverage=1,
                min_expect_both=0,
                detect_fdr=self.DOMAIN_FDR,
                merge_fdr=self.DOMAIN_FDR,
            )

    def test_pair_outside_region_raises(self):
        # A pair whose position lies beyond total_end3 would otherwise be
        # scattered out of bounds (IndexError) or wrap a negative index;
        # it must raise an explicit, interpretable error instead.
        rows = _rows_over_region(
            self.TOTAL_END5, self.TOTAL_END3, self.BAND, self._in_domain
        )
        rows.append((self.TOTAL_END3, self.TOTAL_END3 + 5, self.N_COV, 0.0))
        table = _make_table(rows)
        with self.assertRaises(OutOfBoundsError):
            _calc_domains_by_dp_segmentation(
                table,
                self.TOTAL_END5,
                self.TOTAL_END3,
                min_pair_coverage=1,
                min_expect_both=0,
                detect_fdr=self.DOMAIN_FDR,
                merge_fdr=self.DOMAIN_FDR,
            )

    def test_realistic_noise_does_not_fragment_the_domain(self):
        # _block_score alone (no p-value gate) can prefer atomizing a real
        # domain into single-pair fragments over reporting it whole, once
        # per-pair chi-square carries realistic sampling noise rather than
        # an idealized constant (verified this session: ~3.6% margin for a
        # 150-pair domain). The Benjamini-Hochberg p-value gate closes this
        # gap by rejecting most such "lucky" fragments against the full
        # candidate pool. Build one real domain (every pair elevated, but
        # each pair's own chi-square independently noisy) and confirm it
        # comes back as one contiguous domain, not shattered into pieces
        # far smaller than its true extent.
        import numpy as np

        total_end5, total_end3, band = 1, 30, 10
        rng = np.random.default_rng(0)
        true_r = 29.0
        rows = []
        for i in range(total_end5, total_end3 + 1):
            for j in range(i + 1, min(i + band, total_end3) + 1):
                chi2 = float((1 + true_r) * rng.chisquare(df=1))
                rows.append((i, j, self.N_COV, chi2))
        table = _make_table(rows)
        domains = _calc_domains_by_dp_segmentation(
            table,
            total_end5,
            total_end3,
            min_pair_coverage=1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
        )
        self.assertTrue(domains, "Expected the domain to be recovered")
        expect_length = total_end3 - total_end5 + 1
        self.assertTrue(
            any((end - start + 1) >= expect_length / 2 for start, end in domains),
            f"Domain was fragmented into pieces far smaller than its true "
            f"extent instead of reported (near-)whole: {domains}",
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
    DOMAIN_FDR = 0.1

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
            min_pair_coverage=1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
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
            min_pair_coverage=1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
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


class TestCutCrossingScores(ut.TestCase):
    """``_cut_crossing_scores`` sums, at every cut ``m``, exactly the
    observable pairs ``(a, b)`` with ``a <= m < b`` -- the pairs that
    straddle that specific cut."""

    def test_no_pairs_straddle_any_cut(self):
        # A pair entirely on one side of every cut (or, trivially, an
        # empty table) contributes to no cut at all.
        table = _make_table([])
        n_cross, chi2_cross = _cut_crossing_scores(table, 1, 10)
        self.assertEqual(n_cross.shape, (9,))
        self.assertTrue((n_cross == 0).all())
        self.assertTrue((chi2_cross == 0).all())

    def test_single_pair_straddles_exactly_the_cuts_between_its_ends(self):
        # Pair (3, 7) straddles cuts 3, 4, 5, 6 (between position m and
        # m + 1, for m = 3..6) and no others.
        table = _make_table([(3, 7, 2000.0, 12.0)])
        n_cross, chi2_cross = _cut_crossing_scores(table, 1, 10)
        expect_n = [0, 0, 1, 1, 1, 1, 0, 0, 0]
        expect_chi2 = [0, 0, 12.0, 12.0, 12.0, 12.0, 0, 0, 0]
        self.assertEqual(n_cross.tolist(), expect_n)
        self.assertEqual(chi2_cross.tolist(), expect_chi2)

    def test_multiple_pairs_accumulate_at_a_shared_cut(self):
        table = _make_table(
            [(1, 5, 2000.0, 10.0), (2, 6, 2000.0, 20.0), (4, 5, 2000.0, 5.0)]
        )
        n_cross, chi2_cross = _cut_crossing_scores(table, 1, 6)
        # Cut m=4 (between positions 4 and 5) is straddled by all three
        # pairs: (1,5), (2,6), (4,5).
        self.assertEqual(n_cross[3], 3)
        self.assertEqual(chi2_cross[3], 35.0)


class TestMergeConnectedBlocks(ut.TestCase):
    """``_merge_connected_blocks`` merges adjacent blocks whenever *every*
    cut in the gap between them (including the block-edge cut itself) is
    connected: its crossing pairs clear both the effect-size bar
    (``_block_score >= 0``) and Benjamini-Hochberg significance. Both gates
    are required -- a cut's crossing-pair count can be in the thousands (see
    ``_cut_crossing_scores``), so significance alone flags a
    biologically-trivial elevation as "real" (verified below)."""

    def _score_and_pvalue(self, n: float, chi2_sum: float):
        from scipy.stats import chi2 as chi2dist

        gain = float(_block_score(np.array([n]), np.array([chi2_sum]))[0])
        pval = float(chi2dist.sf(chi2_sum, df=n))
        return gain, pval

    def test_transitive_merge_across_a_bridged_gap(self):
        # Three blocks; every cut in both gaps is strongly connected, so
        # all three should merge transitively into one domain.
        blocks = [(0, 9), (20, 29), (40, 49)]
        n_cross = np.full(49, 100.0)
        chi2_cross = np.full(49, 100.0 * 4.0)  # mean chi2 = 4: strong signal
        merged = _merge_connected_blocks(blocks, n_cross, chi2_cross, merge_fdr=0.1)
        self.assertEqual(merged, [(0, 49)])

    def test_one_null_cut_in_a_gap_keeps_that_gap_split(self):
        # Gap 1 (cuts 9..19) is strong everywhere except cut 14, which is
        # null; gap 2 (cuts 29..39) is strong everywhere. Only gap 2 should
        # merge.
        blocks = [(0, 9), (20, 29), (40, 49)]
        n_cross = np.full(49, 100.0)
        chi2_cross = np.full(49, 100.0 * 4.0)
        n_cross[14] = 100.0
        chi2_cross[14] = 100.0  # mean chi2 = 1: null
        merged = _merge_connected_blocks(blocks, n_cross, chi2_cross, merge_fdr=0.1)
        self.assertEqual(merged, [(0, 9), (20, 49)])

    def test_effect_size_gate_rejects_a_faint_but_bh_significant_crossing(self):
        # A gap whose crossing pairs are elevated only trivially above the
        # null (mean chi2 = 1.05) but with a huge n (5000, the scale of an
        # actual cut's crossing-pair count): BH-significance alone would
        # call this "connected" (p ~ 7e-3), but the exact block score is
        # negative (BIC-charged effect size below the bar), so the gap
        # must stay split.
        n, mean = 5000.0, 1.05
        gain, pval = self._score_and_pvalue(n, n * mean)
        self.assertLess(gain, 0.0)
        self.assertLess(pval, 0.1)  # BH-significant on its own
        blocks = [(0, 9), (20, 29)]
        n_cross = np.full(19, n)
        chi2_cross = np.full(19, n * mean)
        merged = _merge_connected_blocks(blocks, n_cross, chi2_cross, merge_fdr=0.1)
        self.assertEqual(merged, blocks)

    def test_merge_never_extends_past_the_given_blocks(self):
        blocks = [(0, 9), (20, 29)]
        n_cross = np.full(19, 100.0)
        chi2_cross = np.full(19, 100.0 * 4.0)
        merged = _merge_connected_blocks(blocks, n_cross, chi2_cross, merge_fdr=0.1)
        self.assertEqual(merged, [(0, 29)])

    def test_grouping_is_independent_of_sweep_direction(self):
        # Four blocks, three gaps: gap 1 and gap 3 strongly connected,
        # gap 2 null. Whether grouped left-to-right (as implemented) or
        # right-to-left, the connectivity of each gap is fixed in advance
        # (it depends only on that gap's own cut scores, never on merge
        # decisions elsewhere), so both sweep directions must produce the
        # identical partition.
        blocks = [(0, 9), (20, 29), (40, 49), (60, 69)]
        n_cross = np.full(69, 100.0)
        chi2_cross = np.full(69, 100.0 * 4.0)
        n_cross[29:40] = 100.0
        chi2_cross[29:40] = 100.0  # gap 2 (cuts 29..39) is null

        def merge_right_to_left(blocks, n_cross, chi2_cross, merge_fdr):
            from scipy.stats import chi2 as chi2dist

            valid = n_cross > 0
            gains = np.where(
                valid, _block_score(np.maximum(n_cross, 1), chi2_cross), -1.0
            )
            with np.errstate(divide="ignore", invalid="ignore"):
                pvals = np.where(
                    valid, chi2dist.sf(chi2_cross, df=np.maximum(n_cross, 1)), 1.0
                )
            qvals = np.ones_like(pvals)
            if valid.any():
                qvals[valid] = calc_bh_adjusted_pvals(pvals[valid])
            connected = valid & (gains >= 0.0) & (qvals <= merge_fdr)
            gap_connected = [
                bool(connected[blocks[i][1] : blocks[i + 1][0]].all())
                for i in range(len(blocks) - 1)
            ]
            merged = [blocks[-1]]
            for i in range(len(blocks) - 2, -1, -1):
                if gap_connected[i]:
                    _, e1 = merged[0]
                    s0, _ = blocks[i]
                    merged[0] = (s0, e1)
                else:
                    merged.insert(0, blocks[i])
            return merged

        forward = _merge_connected_blocks(blocks, n_cross, chi2_cross, merge_fdr=0.1)
        backward = merge_right_to_left(blocks, n_cross, chi2_cross, merge_fdr=0.1)
        self.assertEqual(forward, backward)
        self.assertEqual(forward, [(0, 29), (40, 69)])

    def test_fewer_than_two_blocks_returned_unchanged(self):
        self.assertEqual(
            _merge_connected_blocks([], np.zeros(0), np.zeros(0), merge_fdr=0.1), []
        )
        self.assertEqual(
            _merge_connected_blocks([(0, 9)], np.zeros(0), np.zeros(0), merge_fdr=0.1),
            [(0, 9)],
        )


class TestMergeConnectedBlocksIntegration(ut.TestCase):
    """The scenario that motivated ``_merge_connected_blocks``: two dense
    blocks (DOMAIN_A, DOMAIN_B) separated by a genuinely null gap --
    every gap-internal and block-to-gap pair is pure null, so the DP's
    own whole-bridge objective (union minus each block's own triangle)
    is necessarily very negative and the DP splits them (see
    ``_merge_connected_blocks``'s docstring for why a post-hoc re-test of
    that same whole-bridge quantity could never find anything to join).
    But when there is real, direct evidence connecting the two blocks
    -- pairs with one endpoint in each, independent of the gap between
    them -- the merge should still report one domain spanning both; when
    that direct evidence is absent, it should correctly leave them split.
    This is the real-world case discovered in a 1M-read dataset: two
    pieces of what was structurally one long-range-paired domain,
    separated by an unpaired linker with no correlations of its own,
    wrongly reported as two domains until this mechanism was added."""

    DOMAIN_A = (1, 10)
    GAP = (11, 30)
    DOMAIN_B = (31, 40)
    TOTAL_END5, TOTAL_END3 = DOMAIN_A[0], DOMAIN_B[1]
    BAND = TOTAL_END3 - TOTAL_END5
    N_COV = 2000.0
    DOMAIN_FDR = 0.1
    # Below this many elevated direct A-to-B crossing pairs (out of the
    # 10*10=100 possible), the crossing evidence itself nets negative;
    # comfortably above it, the join should fire (verified offline).
    N_CROSS_ELEVATED_JOINS = 10

    def _build_table(self, n_cross_elevated: int):
        a_lo, a_hi = self.DOMAIN_A
        b_lo, b_hi = self.DOMAIN_B
        cross_pairs = [
            (i, j) for i in range(a_lo, a_hi + 1) for j in range(b_lo, b_hi + 1)
        ]
        elevated = set(cross_pairs[:n_cross_elevated])
        rows = []
        for i in range(self.TOTAL_END5, self.TOTAL_END3 + 1):
            for j in range(i + 1, min(i + self.BAND, self.TOTAL_END3) + 1):
                if a_lo <= i <= a_hi and a_lo <= j <= a_hi:
                    chi2 = 30.0
                elif b_lo <= i <= b_hi and b_lo <= j <= b_hi:
                    chi2 = 30.0
                elif (i, j) in elevated:
                    chi2 = 30.0
                else:
                    chi2 = 1.0
                rows.append((i, j, self.N_COV, chi2))
        return _make_table(rows)

    def _call(self, n_cross_elevated: int):
        table = self._build_table(n_cross_elevated)
        return _calc_domains_by_dp_segmentation(
            table,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            min_expect_both=0,
            detect_fdr=self.DOMAIN_FDR,
            merge_fdr=self.DOMAIN_FDR,
        )

    def test_dp_alone_splits_across_the_null_gap(self):
        # With no direct crossing evidence at all, the DP itself must
        # split -- the gap is genuinely null -- and the merge must not
        # manufacture a connection that isn't there.
        domains = self._call(n_cross_elevated=0)
        self.assertFalse(
            any(
                start <= self.DOMAIN_A[0] and end >= self.DOMAIN_B[1]
                for start, end in domains
            ),
            f"Did not expect A and B joined with no crossing evidence, got {domains}",
        )
        for lo, hi in (self.DOMAIN_A, self.DOMAIN_B):
            self.assertTrue(
                any(start <= lo and hi <= end for start, end in domains),
                f"Sub-domain {(lo, hi)} not recovered whole in {domains}",
            )

    def test_real_crossing_evidence_joins_across_the_null_gap(self):
        # The DP still splits them (the gap alone makes the whole-bridge
        # objective very negative), but the direct A-to-B crossing pairs
        # alone are real signal against the null, so the merge fires.
        domains = self._call(n_cross_elevated=self.N_CROSS_ELEVATED_JOINS)
        self.assertTrue(
            any(
                start <= self.DOMAIN_A[0] and end >= self.DOMAIN_B[1]
                for start, end in domains
            ),
            f"Expected the real crossing evidence to join A and B, got {domains}",
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
