import tempfile
import unittest as ut
from itertools import product
from pathlib import Path

from click.testing import CliRunner

from seismicrna.core import path
from seismicrna.core.batch.confusion import POSITION_A, POSITION_B
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
    SCORE_COL,
    CONFUSION_COLS,
    NEITHER_COL,
    ONLY_A_COL,
    ONLY_B_COL,
    BOTH_COL,
    _calc_tiles,
    _build_banded_table,
    _write_pairs_with_confusion,
    _calc_pair_scores,
    _scan_block_gains,
    _calibrate_penalty_from_null,
    _calc_null_penalty,
    _calc_domains_by_dp_segmentation,
    _pair_band_row_cumsum,
    _insert_domains_into_gaps,
    _expand_domains_into_gaps,
    _filter_domains_length,
)
from seismicrna.sim.params import run as sim_params
from seismicrna.sim.idmut import run as sim_idmut


def _make_table(rows: list[tuple[int, int, float, float]]):
    """Build a per-pair table like ``_confusion_to_table`` returns, from an
    iterable of (pos_a, pos_b, n, chi_square)."""
    import pandas as pd

    if not rows:
        pos_a, pos_b, ns, chi_square = (), (), (), ()
    else:
        pos_a, pos_b, ns, chi_square = zip(*rows)
    index = pd.MultiIndex.from_arrays([pos_a, pos_b], names=[POSITION_A, POSITION_B])
    return pd.DataFrame(
        {N_COL: ns, CHI_SQUARE_COL: chi_square},
        index=index,
    )


def _scored(table, min_pair_coverage: int = 1):
    """Attach the per-pair block score (``SCORE_COL``) the DP consumes, the
    way ``_calc_cluster_domains`` does before segmentation (no null replicates:
    the variance-inflation shape ``r`` is estimated from the observed data)."""
    scored, _ = _calc_pair_scores(table, [], min_pair_coverage)
    return scored


def _null_tables(table, n_replicates: int = 10):
    """Build ``n_replicates`` synthetic null-replicate tables (unscored) for
    a real ``table``: same index/N, chi-square fixed at 1.0 (the chi^2(1)
    null's mean) everywhere -- a deterministic stand-in for 'no real
    structure' null replicates, for tests exercising the null-requiring
    paths (``_calc_null_penalty``, ``_join_bridged_blocks``) without going
    through the actual per-position-independence resampling."""
    import pandas as pd

    null = pd.DataFrame({N_COL: table[N_COL], CHI_SQUARE_COL: 1.0}, index=table.index)
    return [null.copy() for _ in range(n_replicates)]


def _scored_with_nulls(table, min_pair_coverage: int = 1, n_null_replicates: int = 10):
    """Like ``_scored``, but also returns ``n_null_replicates`` scored null
    tables (see ``_null_tables``), scored with the same real-table-derived
    shape ``r`` via ``_calc_pair_scores``, matching production."""
    nulls = _null_tables(table, n_null_replicates)
    scored, scored_nulls = _calc_pair_scores(table, nulls, min_pair_coverage)
    return scored, scored_nulls


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
            list(result.columns),
            [N_COL, CHI_SQUARE_COL, *CONFUSION_COLS],
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


class TestCalibratePenaltyFromNull(ut.TestCase):
    """BH-style FDR comparison of real candidate gains to the pooled null."""

    def test_no_real_gains_calls_none(self):
        import numpy as np

        self.assertEqual(
            _calibrate_penalty_from_null(np.array([]), np.array([1.0]), 10, 0.05),
            float("inf"),
        )

    def test_clear_signal_admitted_below_its_gain(self):
        import numpy as np

        real = np.array([100.0])
        null = np.arange(1.0, 11.0)  # all far below the real block
        penalty = _calibrate_penalty_from_null(real, null, 10, 0.05)
        self.assertTrue(0.0 < penalty < 100.0)  # admits the block

    def test_signal_buried_in_null_calls_none(self):
        import numpy as np

        real = np.array([5.0])
        null = np.array([6.0, 7.0, 8.0])  # 3 null blocks above 5 over R=10
        # FDR = (3/10)/1 = 0.3 > 0.05 -> reject.
        self.assertEqual(
            _calibrate_penalty_from_null(real, null, 10, 0.05), float("inf")
        )

    def test_fdr_step_up_selects_expected_cutoff(self):
        import numpy as np

        real = np.array([10.0, 9.0, 8.0, 7.0, 6.0])
        null = np.array([9.0, 7.5, 6.5, 6.2])
        # k=4 at T=7: (2/10)/4 = 0.05 <= 0.05; k=5 at T=6: (4/10)/5 = 0.08 > 0.05.
        # best_k=4 -> penalty between real[3]=7 and real[4]=6 -> 6.5.
        self.assertEqual(_calibrate_penalty_from_null(real, null, 10, 0.05), 6.5)


class TestCalcNullPenalty(ut.TestCase):
    """The block-level null distribution (which captures the cross-pair
    transitivity structure) has no analytical substitute, so at least one
    null replicate is required -- there is no BIC-style fallback."""

    def test_requires_at_least_one_null_replicate(self):
        with self.assertRaises(OutOfBoundsError):
            _calc_null_penalty(
                _make_table([]),
                [],
                None,
                min_pair_coverage=1,
                domain_fdr=0.05,
            )


class TestScanBlockGains(ut.TestCase):
    """The penalty-free block-gain scan (summed per-pair score) used for both
    the observed and the null tables."""

    def test_planted_domain_has_positive_gain(self):
        rows = _rows_over_region(1, 60, 10, lambda i, j: 20 <= i and j <= 40)
        gains = _scan_block_gains(_scored(_make_table(rows)), 1, 60, min_pair_coverage=1)
        self.assertGreater(gains.size, 0)
        self.assertGreater(gains.max(), 0.0)

    def test_no_correlated_pairs_has_no_blocks(self):
        # Every pair null-like (chi^2 ~ 1): the variance-inflation shape r is
        # floored, every score is ~<= 0, so the penalty-free scan carves out
        # no positive-score block.
        rows = _rows_over_region(1, 60, 10, lambda i, j: False)
        gains = _scan_block_gains(_scored(_make_table(rows)), 1, 60, min_pair_coverage=1)
        self.assertEqual(gains.size, 0)


class TestWritePairsWithConfusion(ut.TestCase):
    """pairs.csv must carry each positive-score pair's 2x2 confusion-matrix
    counts and its score (chi-square is omitted since it is recomputable
    from the counts, but the score also depends on the table-wide shape
    r, so it is written explicitly)."""

    @staticmethod
    def _table(rows):
        # rows: (pos_a, pos_b, neither, only_a, only_b, both, chi_square)
        import pandas as pd

        pa, pb, ne, oa, ob, bo, chi2 = zip(*rows)
        n = [ne_ + oa_ + ob_ + bo_ for ne_, oa_, ob_, bo_ in zip(ne, oa, ob, bo)]
        index = pd.MultiIndex.from_arrays([pa, pb], names=[POSITION_A, POSITION_B])
        return pd.DataFrame(
            {
                N_COL: n,
                CHI_SQUARE_COL: list(chi2),
                NEITHER_COL: list(ne),
                ONLY_A_COL: list(oa),
                ONLY_B_COL: list(ob),
                BOTH_COL: list(bo),
            },
            index=index,
        )

    def test_writes_counts_for_positive_score_pairs_only(self):
        import pandas as pd

        table = self._table(
            [
                (1, 5, 900, 30, 20, 50, 30.0),
                (1, 8, 995, 3, 2, 0, 1.0),  # null-like -> excluded (score <= 0)
                (2, 9, 800, 60, 40, 100, 30.0),
            ]
        )
        pos_table = _scored(table)
        pos_table = pos_table[pos_table[SCORE_COL] > 0]
        with tempfile.TemporaryDirectory() as tmp:
            csv_file = Path(tmp) / "pairs.csv"
            _write_pairs_with_confusion(pos_table, csv_file)
            out = pd.read_csv(csv_file)
        # Exactly the position columns, the four confusion cells, and the
        # score; chi-square (recomputable from the counts) is omitted.
        self.assertListEqual(
            list(out.columns), [POSITION_A, POSITION_B, *CONFUSION_COLS, SCORE_COL]
        )
        self.assertNotIn(CHI_SQUARE_COL, out.columns)
        # Only the two positive-score pairs, with their exact counts.
        self.assertListEqual(
            list(zip(out[POSITION_A], out[POSITION_B])), [(1, 5), (2, 9)]
        )
        self.assertListEqual(list(out[BOTH_COL]), [50, 100])
        self.assertListEqual(list(out[ONLY_A_COL]), [30, 60])
        self.assertTrue((out[SCORE_COL] > 0).all())
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
    # An explicit per-block penalty on the score scale, in the regime that
    # admits a real domain's block gain but rejects null-like background and
    # keeps two separate domains from merging across their gap (in production
    # this is calibrated from the null; see _calc_null_penalty).
    PENALTY = 10.0
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
        return _scored_with_nulls(_make_table(rows))

    def test_domains_recovered_and_not_merged(self):
        table, nulls = self._build_table()
        domains = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
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

    def test_explicit_penalty(self):
        # The null-calibrated path supplies an explicit per-block penalty.
        table, nulls = self._build_table()
        domains = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
        )
        self.assertTrue(domains, "A moderate penalty should recover the domains")
        # A penalty larger than any block's gain rejects everything.
        none = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            penalty=1e12,
            domain_fdr=self.DOMAIN_FDR,
        )
        self.assertEqual(none, [])

    def test_uniform_correlation_is_one_domain(self):
        # Every pair in the whole scanned region is strongly correlated:
        # unlike the old density caller (which needed a sparser background to
        # contrast against and called nothing), the distributional score sees
        # the entire region as drawn from a wider-than-null distribution, so
        # it is one domain spanning the whole region.
        rows = _rows_over_region(
            self.TOTAL_END5, self.TOTAL_END3, self.BAND, lambda i, j: True
        )
        table, nulls = _scored_with_nulls(_make_table(rows))
        domains = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
        )
        self.assertEqual(domains, [(self.TOTAL_END5, self.TOTAL_END3)])

    def test_sparse_region_has_no_domain(self):
        rows = _rows_over_region(
            self.TOTAL_END5, self.TOTAL_END3, self.BAND, lambda i, j: False
        )
        table, nulls = _scored_with_nulls(_make_table(rows))
        domains = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
        )
        self.assertEqual(domains, [])

    def test_min_pair_coverage_excludes_low_coverage_pairs(self):
        table, nulls = self._build_table()
        # Below the pairs' actual coverage (N_COV): pairs pass through
        # unfiltered, and the domains are detected as usual.
        domains = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
        )
        self.assertTrue(domains)
        # Above the pairs' actual coverage: every pair is filtered out,
        # so no domain can be detected.
        domains = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=int(self.N_COV) + 1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
        )
        self.assertEqual(domains, [])

    def test_min_pair_coverage_must_be_at_least_one(self):
        table, nulls = self._build_table()
        with self.assertRaises(OutOfBoundsError):
            _calc_domains_by_dp_segmentation(
                table,
                nulls,
                self.TOTAL_END5,
                self.TOTAL_END3,
                min_pair_coverage=0,
                penalty=self.PENALTY,
                domain_fdr=self.DOMAIN_FDR,
            )

    def test_masked_pairs_never_contribute(self):
        # A "masked" pocket in the gap between the two real domains:
        # strongly correlated pairs, but at coverage far below
        # min_pair_coverage, as if these positions had too little read
        # depth to trust (even though, if not excluded, they would look
        # like a spurious third domain). The filter must make them
        # invisible, not merely down-weighted: calling with and without
        # the pocket present must give byte-identical results -- including
        # the estimated variance-inflation shape r, so the pocket is
        # excluded from scoring too (score at min_pair_coverage).
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

        table_baseline, nulls_baseline = _scored_with_nulls(
            _make_table(build_rows(masked_pocket=False)), min_pair_coverage
        )
        table_with_masked, nulls_with_masked = _scored_with_nulls(
            _make_table(build_rows(masked_pocket=True)), min_pair_coverage
        )
        domains_baseline = _calc_domains_by_dp_segmentation(
            table_baseline,
            nulls_baseline,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=min_pair_coverage,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
        )
        domains_with_masked = _calc_domains_by_dp_segmentation(
            table_with_masked,
            nulls_with_masked,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=min_pair_coverage,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
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
        table, nulls = _scored_with_nulls(_make_table(rows))
        domains = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            total_end5,
            total_end3,
            min_pair_coverage=1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
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
        table, nulls = self._build_table()
        with self.assertRaises(IncompatibleValuesError):
            _calc_domains_by_dp_segmentation(
                table,
                nulls,
                self.TOTAL_END3,
                self.TOTAL_END5,
                min_pair_coverage=1,
                penalty=self.PENALTY,
                domain_fdr=self.DOMAIN_FDR,
            )

    def test_pair_outside_region_raises(self):
        # A pair whose position lies beyond total_end3 would otherwise be
        # scattered out of bounds (IndexError) or wrap a negative index;
        # it must raise an explicit, interpretable error instead.
        rows = _rows_over_region(
            self.TOTAL_END5, self.TOTAL_END3, self.BAND, self._in_domain
        )
        rows.append((self.TOTAL_END3, self.TOTAL_END3 + 5, self.N_COV, 0.0))
        table, nulls = _scored_with_nulls(_make_table(rows))
        with self.assertRaises(OutOfBoundsError):
            _calc_domains_by_dp_segmentation(
                table,
                nulls,
                self.TOTAL_END5,
                self.TOTAL_END3,
                min_pair_coverage=1,
                penalty=self.PENALTY,
                domain_fdr=self.DOMAIN_FDR,
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
    PENALTY = 10.0
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
        return _scored_with_nulls(_make_table(rows))

    def test_dense_bridge_keeps_domain_whole(self):
        table, nulls = self._build_table(dense_bridge=True)
        domains = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
        )
        self.assertTrue(
            any(
                start <= self.SUB_A[0] and self.SUB_B[1] <= end
                for start, end in domains
            ),
            f"Expected one domain spanning both sub-clusters in {domains}",
        )

    def test_sparse_bridge_splits_domain(self):
        table, nulls = self._build_table(dense_bridge=False)
        domains = _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
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


class TestJoinBridgedBlocks(ut.TestCase):
    """The scenario that motivated ``_join_bridged_blocks``: two dense
    blocks (DOMAIN_A, DOMAIN_B) separated by a genuinely null gap --
    every gap-internal and block-to-gap pair is pure null, so the DP's
    own whole-bridge objective (union minus each block's own triangle)
    is necessarily very negative and the DP splits them (see
    ``_join_bridged_blocks``'s docstring for why a post-hoc re-test of
    that same whole-bridge quantity could never find anything to join).
    But when there is real, direct evidence connecting the two blocks
    -- pairs with one endpoint in each, independent of the gap between
    them -- ``_join_bridged_blocks`` should still report one domain
    spanning both; when that direct evidence is absent, it should
    correctly leave them split. This is the real-world case discovered
    in a 1M-read dataset: two pieces of what was structurally one
    long-range-paired domain, separated by an unpaired linker with no
    correlations of its own, wrongly reported as two domains until this
    mechanism was added."""

    DOMAIN_A = (1, 10)
    GAP = (11, 30)
    DOMAIN_B = (31, 40)
    TOTAL_END5, TOTAL_END3 = DOMAIN_A[0], DOMAIN_B[1]
    BAND = TOTAL_END3 - TOTAL_END5
    N_COV = 2000.0
    PENALTY = 10.0
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
        return _scored_with_nulls(_make_table(rows))

    def _call(self, n_cross_elevated: int):
        table, nulls = self._build_table(n_cross_elevated)
        return _calc_domains_by_dp_segmentation(
            table,
            nulls,
            self.TOTAL_END5,
            self.TOTAL_END3,
            min_pair_coverage=1,
            penalty=self.PENALTY,
            domain_fdr=self.DOMAIN_FDR,
        )

    def test_dp_alone_splits_across_the_null_gap(self):
        # With no direct crossing evidence at all, the DP itself must
        # split -- the gap is genuinely null -- and the join must not
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
        # alone are real signal against the null, so the join fires.
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
            seed=0,
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
