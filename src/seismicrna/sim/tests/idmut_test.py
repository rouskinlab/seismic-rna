import unittest as ut
from itertools import product

import numpy as np
import pandas as pd

from seismicrna.core.arg import MUT_COLLISIONS_DROP, MUT_COLLISIONS_MERGE
from seismicrna.core.header import CLUST_NAME, NUM_CLUSTS_NAME, REL_NAME
from seismicrna.core.rel import MATCH, SUB_T, RelPattern
from seismicrna.core.seq.region import BASE_NAME, POS_NAME
from seismicrna.core.unbias import (
    calc_p_clust_given_noclose,
    calc_p_noclose_given_clust,
    calc_p_noclose_given_ends_auto,
)
from seismicrna.idmut.sim import simulate_batches, make_p_ends_2d, calc_pmut_pattern
from seismicrna.sim.idmut import parse_min_mut_gap_weights


class TestParseMinMutGapWeights(ut.TestCase):
    def test_empty_string(self):
        self.assertEqual(parse_min_mut_gap_weights(""), {})

    def test_single_pair_gap0(self):
        self.assertEqual(parse_min_mut_gap_weights("0:1.0"), {0: 1.0})

    def test_single_pair_nonzero_gap(self):
        self.assertEqual(parse_min_mut_gap_weights("3:1.0"), {3: 1.0})

    def test_multiple_pairs(self):
        self.assertEqual(parse_min_mut_gap_weights("0:0.3,2:0.7"), {0: 0.3, 2: 0.7})

    def test_output_sorted_by_gap(self):
        result = parse_min_mut_gap_weights("5:0.5,0:0.5")
        self.assertEqual(list(result.keys()), [0, 5])

    def test_zero_weight_excluded(self):
        self.assertEqual(parse_min_mut_gap_weights("0:0.0,3:1.0"), {3: 1.0})

    def test_multiple_zero_weights_excluded(self):
        self.assertEqual(
            parse_min_mut_gap_weights("0:0.5,2:0.0,4:0.5"), {0: 0.5, 4: 0.5}
        )

    def test_weights_sum_to_1_within_float_tolerance(self):
        # 0.1 + 0.2 + 0.7 may not be exactly 1.0 in floating point.
        result = parse_min_mut_gap_weights("0:0.1,2:0.2,4:0.7")
        self.assertAlmostEqual(sum(result.values()), 1.0)

    def test_single_pair_weight_zero_invalid(self):
        # A single pair with weight 0 sums to 0, not 1.
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("0:0.0")

    def test_pair_missing_colon(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("01.0")

    def test_pair_extra_colon(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("0:0.5:extra,2:0.5")

    def test_gap_not_integer(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("a:1.0")

    def test_gap_float_string(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("1.5:1.0")

    def test_weight_not_float(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("0:b")

    def test_negative_gap(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("-1:1.0")

    def test_negative_weight(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("0:-0.5,2:1.5")

    def test_weight_above_1(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("0:1.5")

    def test_duplicate_gap(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("0:0.5,0:0.5")

    def test_weights_do_not_sum_to_1(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("0:0.3,2:0.3")

    def test_weights_exceed_1_in_total(self):
        with self.assertRaises(ValueError):
            parse_min_mut_gap_weights("0:0.7,2:0.7")


class TestSimulateBatches(ut.TestCase):
    NUM_READS = 10_000_000
    PCLUST = [0.2, 0.5, 0.3]
    P_GAP = {0: 0.2, 5: 0.3, 15: 0.5}
    TOLERANCE = 0.005  # 0.5%

    # MATCH|SUB_T represents an ambiguous/low-quality base call that is
    # neither a definitive match nor a mutation.  RelPattern.muts() excludes
    # it (it has the MATCH bit set), while 1-pmut[MATCH] would incorrectly
    # count it as a mutation.  simulate_muts accepts it because it is not
    # NOCOV (255).  Including it here ensures the test catches any regression
    # where calc_pmut_pattern is replaced with 1-pmut[MATCH].
    AMB = MATCH | SUB_T  # = 129

    def _make_pmut(self, seed: int):
        rng = np.random.default_rng(seed)
        row_index = pd.MultiIndex.from_tuples(
            [(pos, "A") for pos in range(1, 31)], names=[POS_NAME, BASE_NAME]
        )
        col_index = pd.MultiIndex.from_tuples(
            [
                (rel, 3, clust)
                for rel, clust in product([MATCH, self.AMB, SUB_T], [1, 2, 3])
            ],
            names=[REL_NAME, NUM_CLUSTS_NAME, CLUST_NAME],
        )
        p_amb = 0.10
        p_mut = rng.beta(a=2, b=8, size=(row_index.size, 3)) * (1.0 - p_amb)
        p_match = 1.0 - p_amb - p_mut
        pmut = pd.DataFrame(0.0, index=row_index, columns=col_index)
        pmut[MATCH] = p_match
        pmut[self.AMB] = p_amb
        pmut[SUB_T] = p_mut
        return pmut

    def _make_pclust(self):
        index = pd.MultiIndex.from_tuples(
            [(3, 1), (3, 2), (3, 3)], names=[NUM_CLUSTS_NAME, CLUST_NAME]
        )
        return pd.Series(self.PCLUST, index=index)

    def _run(
        self,
        mut_collisions: str,
        seed: int,
        min_mut_gap: int = 0,
        min_mut_gap_weights: dict | None = None,
    ):
        return list(
            simulate_batches(
                batch_size=self.NUM_READS,
                pmut=self._make_pmut(seed),
                pclust=self._make_pclust(),
                pends=np.array([1.0]),
                num_reads=self.NUM_READS,
                seed=0,
                min_mut_gap=min_mut_gap,
                min_mut_gap_weights=(
                    self.P_GAP if min_mut_gap_weights is None else min_mut_gap_weights
                ),
                mut_collisions=mut_collisions,
                sample="test",
                branches={},
                ref="testref",
                write_read_names=False,
                uniq_end5s=np.array([1]),
                uniq_end3s=np.array([10]),
                paired=False,
                read_length=0,
                p_rev=0.5,
                injected_mut_probs=None,
            )
        )

    def _expected_drop_with_weights(self, seed: int):
        pmut = self._make_pmut(seed)
        pmut_given_clust = calc_pmut_pattern(
            pmut, RelPattern.muts(), normalize=False
        ).values
        p_ends_2d = make_p_ends_2d(
            np.array([1.0]), np.array([1]), np.array([10]), pmut.index
        )
        # The weights are the probability that a read came from the
        # distibution with the minimum gap given that it has no two
        # mutations too close.
        p_gap_given_noclose = np.array(list(self.P_GAP.values()))
        # Calculate the probability that a read has no mutations too
        # close for each cluster and mutation gap.
        # Shape: (num_gaps, num_clusters)
        p_noclose_given_gap_clust = np.stack(
            [
                calc_p_noclose_given_clust(
                    p_ends_2d, calc_p_noclose_given_ends_auto(pmut_given_clust, gap)
                )
                for gap in self.P_GAP.keys()
            ],
            axis=0,
        )
        # Calculate the unconditional proportion of each gap.
        # P(gap | noclose) = P(gap) / P(noclose) * P(noclose | gap)
        # P(gap) = p(noclose) * P(gap | noclose) / P(noclose | gap)
        pclust = np.array(self.PCLUST)
        p_noclose_given_gap = p_noclose_given_gap_clust @ pclust
        p_gap = p_gap_given_noclose / p_noclose_given_gap
        p_gap /= p_gap.sum()
        # Calculate the proportion of each cluster after dropping
        # reads with two mutations too close.
        p_noclose_given_clust = (
            p_gap[np.newaxis, :] @ p_noclose_given_gap_clust
        ).reshape(-1)
        p_clust_given_noclose = calc_p_clust_given_noclose(
            pclust, p_noclose_given_clust
        )
        # Calculate the number of reads expected after dropping those
        # with close mutations for each combination of number of
        # clusters and minimum gap between mutations.
        return np.array(
            [
                self.NUM_READS * p_clust_val * p_gap_val
                for p_clust_val, p_gap_val in product(
                    p_clust_given_noclose, p_gap_given_noclose
                )
            ]
        )

    def test_read_counts_with_gap_weights(self):
        rng = np.random.default_rng(seed=0)
        for mut_collisions in [MUT_COLLISIONS_DROP, MUT_COLLISIONS_MERGE]:
            with self.subTest(mut_collisions=mut_collisions):
                seed = int(rng.integers(2**32))
                batches = self._run(mut_collisions, seed)
                self.assertEqual(len(batches), len(self.PCLUST) * len(self.P_GAP))
                p_gaps = list(self.P_GAP.values())
                if mut_collisions == MUT_COLLISIONS_DROP:
                    expected = self._expected_drop_with_weights(seed)
                else:
                    expected = np.array(
                        [
                            self.NUM_READS * p_clust * p_gap
                            for p_clust, p_gap in product(self.PCLUST, p_gaps)
                        ]
                    )
                actual = np.array([batch.num_reads for batch, _ in batches])
                for act, exp in zip(actual, expected, strict=True):
                    self.assertAlmostEqual(
                        act / exp,
                        1.0,
                        delta=self.TOLERANCE,
                        msg=f"Expected ~{exp:.0f} reads, got {act}",
                    )

    def _expected_drop_no_weights(self, seed: int, min_mut_gap: int):
        """Expected per-cluster read counts for the single-gap DROP path."""
        pmut = self._make_pmut(seed)
        pmut_given_clust = calc_pmut_pattern(
            pmut, RelPattern.muts(), normalize=False
        ).values
        p_ends_2d = make_p_ends_2d(
            np.array([1.0]), np.array([1]), np.array([10]), pmut.index
        )
        p_noclose_given_clust = calc_p_noclose_given_clust(
            p_ends_2d, calc_p_noclose_given_ends_auto(pmut_given_clust, min_mut_gap)
        )
        pclust = np.array(self.PCLUST)
        p_clust_given_noclose = calc_p_clust_given_noclose(
            pclust, p_noclose_given_clust
        )
        return self.NUM_READS * p_clust_given_noclose

    def test_drop_read_counts_without_gap_weights(self):
        """Single-gap DROP path must yield ~num_reads final reads."""
        min_mut_gap = 1
        rng = np.random.default_rng(seed=0)
        seed = int(rng.integers(2**32))
        batches = self._run(
            MUT_COLLISIONS_DROP, seed, min_mut_gap=min_mut_gap, min_mut_gap_weights={}
        )
        # batch_size == num_reads, so each cluster produces exactly one batch.
        self.assertEqual(len(batches), len(self.PCLUST))
        expected = self._expected_drop_no_weights(seed, min_mut_gap)
        actual = np.array([batch.num_reads for batch, _ in batches])
        # After inflating simulated count by 1/P(noclose), total should be ~num_reads.
        self.assertAlmostEqual(
            actual.sum() / self.NUM_READS,
            1.0,
            delta=self.TOLERANCE,
            msg=f"Total reads: expected ~{self.NUM_READS}, got {actual.sum()}",
        )
        for act, exp in zip(actual, expected, strict=True):
            self.assertAlmostEqual(
                act / exp,
                1.0,
                delta=self.TOLERANCE,
                msg=f"Expected ~{exp:.0f} reads, got {act}",
            )


if __name__ == "__main__":
    ut.main(verbosity=2)
