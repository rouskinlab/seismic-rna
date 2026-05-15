import unittest as ut

import numpy as np
import pandas as pd

from seismicrna.core.batch.muts import calc_muts_matrix
from seismicrna.core.rel import RelPattern
from seismicrna.core.rel.code import (DELET,
                                      MATCH,
                                      NOCOV,
                                      SUB_A,
                                      SUB_C,
                                      SUB_G,
                                      SUB_T)
from seismicrna.core.seq import DNA, Region
from seismicrna.mask.batch import MaskMutsBatch


class TestCalcMutsMatrix(ut.TestCase):

    def test_full_reads_no_muts(self):
        for length in range(10):
            region = Region("myref", DNA.random(length))
            muts = dict()
            for num_reads in range(10):
                read_nums = np.arange(num_reads)
                seg_end5s = np.full((num_reads, 1), region.end5)
                seg_end3s = np.full((num_reads, 1), region.end3)
                mask = seg_end5s > seg_end3s
                with self.subTest(length=length, num_reads=num_reads):
                    result = calc_muts_matrix(region,
                                              read_nums,
                                              seg_end5s,
                                              seg_end3s,
                                              mask,
                                              muts)
                    expect = pd.DataFrame(MATCH, read_nums, region.unmasked)
                    self.assertTrue(expect.equals(result))

    def test_full_reads_no_muts_some_masked(self):
        region = Region("myref", DNA("GTACTCAG"))
        region.mask_g()
        region.mask_t()
        muts = dict()
        for num_reads in range(10):
            read_nums = np.arange(num_reads)
            seg_end5s = np.full((num_reads, 1), region.end5)
            seg_end3s = np.full((num_reads, 1), region.end3)
            mask = seg_end5s > seg_end3s
            with self.subTest(num_reads=num_reads):
                result = calc_muts_matrix(region,
                                          read_nums,
                                          seg_end5s,
                                          seg_end3s,
                                          mask,
                                          muts)
                expect = pd.DataFrame(MATCH, read_nums, region.unmasked)
                self.assertTrue(expect.equals(result))

    def test_partial_reads_no_muts(self):
        region = Region("myref", DNA.random(5))
        muts = dict()
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16, 19, 20])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6]]).T
        seg_end3s = np.array([[0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5]]).T
        mask = seg_end5s > seg_end3s
        result = calc_muts_matrix(region,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  mask,
                                  muts)
        expect = pd.DataFrame([[NOCOV, NOCOV, NOCOV, NOCOV, NOCOV],
                               [MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [MATCH, MATCH, NOCOV, NOCOV, NOCOV],
                               [MATCH, MATCH, MATCH, NOCOV, NOCOV],
                               [MATCH, MATCH, MATCH, MATCH, NOCOV],
                               [MATCH, MATCH, MATCH, MATCH, MATCH],
                               [NOCOV, MATCH, MATCH, MATCH, MATCH],
                               [NOCOV, NOCOV, MATCH, MATCH, MATCH],
                               [NOCOV, NOCOV, NOCOV, MATCH, MATCH],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH],
                               [NOCOV, NOCOV, NOCOV, NOCOV, NOCOV]],
                              read_nums,
                              region.unmasked)
        self.assertTrue(expect.equals(result))

    def test_partial_reads_no_muts_some_masked(self):
        region = Region("myref", DNA("TAGCT"))
        region.mask_g()
        region.mask_t()
        muts = dict()
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16, 19, 20])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6]]).T
        seg_end3s = np.array([[0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5]]).T
        mask = seg_end5s > seg_end3s
        result = calc_muts_matrix(region,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  mask,
                                  muts)
        expect = pd.DataFrame([[NOCOV, NOCOV],
                               [NOCOV, NOCOV],
                               [MATCH, NOCOV],
                               [MATCH, NOCOV],
                               [MATCH, MATCH],
                               [MATCH, MATCH],
                               [MATCH, MATCH],
                               [NOCOV, MATCH],
                               [NOCOV, MATCH],
                               [NOCOV, NOCOV],
                               [NOCOV, NOCOV]],
                              read_nums,
                              region.unmasked)
        self.assertTrue(expect.equals(result))

    def test_paired_reads_no_muts(self):
        region = Region("myref", DNA.random(5))
        muts = dict()
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 2, 3, 4, 5],
                              [1, 2, 3, 4, 5, 5, 5, 5, 5]]).T
        seg_end3s = np.array([[1, 1, 1, 1, 1, 2, 3, 4, 5],
                              [1, 2, 3, 4, 5, 5, 5, 5, 5]]).T
        mask = seg_end5s > seg_end3s
        result = calc_muts_matrix(region,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  mask,
                                  muts)
        expect = pd.DataFrame([[MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [MATCH, MATCH, NOCOV, NOCOV, NOCOV],
                               [MATCH, NOCOV, MATCH, NOCOV, NOCOV],
                               [MATCH, NOCOV, NOCOV, MATCH, NOCOV],
                               [MATCH, NOCOV, NOCOV, NOCOV, MATCH],
                               [NOCOV, MATCH, NOCOV, NOCOV, MATCH],
                               [NOCOV, NOCOV, MATCH, NOCOV, MATCH],
                               [NOCOV, NOCOV, NOCOV, MATCH, MATCH],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH]],
                              read_nums,
                              region.unmasked)
        self.assertTrue(expect.equals(result))

    def test_paired_reads_masked_segments(self):
        region = Region("myref", DNA.random(5))
        muts = dict()
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 2, 3, 4, 5],
                              [1, 2, 3, 4, 5, 5, 5, 5, 5]]).T
        seg_end3s = np.array([[1, 0, 1, 0, 1, 0, 3, 0, 5],
                              [0, 2, 0, 4, 0, 5, 0, 5, 0]]).T
        mask = seg_end5s > seg_end3s
        result = calc_muts_matrix(region,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  mask,
                                  muts)
        expect = pd.DataFrame([[MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [NOCOV, MATCH, NOCOV, NOCOV, NOCOV],
                               [MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [NOCOV, NOCOV, NOCOV, MATCH, NOCOV],
                               [MATCH, NOCOV, NOCOV, NOCOV, NOCOV],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH],
                               [NOCOV, NOCOV, MATCH, NOCOV, NOCOV],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH],
                               [NOCOV, NOCOV, NOCOV, NOCOV, MATCH]],
                              read_nums,
                              region.unmasked)
        self.assertTrue(expect.equals(result))

    def test_partial_reads_muts(self):
        region = Region("myref", DNA.random(5))
        muts = {1: {DELET: np.array([3]),
                    SUB_A: np.array([5]),
                    SUB_C: np.array([7]),
                    SUB_G: np.array([8]),
                    SUB_T: np.array([9])},
                2: {DELET: np.array([5]),
                    SUB_A: np.array([7]),
                    SUB_C: np.array([8]),
                    SUB_G: np.array([9]),
                    SUB_T: np.array([12])},
                3: {DELET: np.array([7]),
                    SUB_A: np.array([8]),
                    SUB_C: np.array([9]),
                    SUB_G: np.array([12]),
                    SUB_T: np.array([13])},
                4: {DELET: np.array([8]),
                    SUB_A: np.array([9]),
                    SUB_C: np.array([12]),
                    SUB_G: np.array([13]),
                    SUB_T: np.array([16])},
                5: {DELET: np.array([9]),
                    SUB_A: np.array([12]),
                    SUB_C: np.array([13]),
                    SUB_G: np.array([16]),
                    SUB_T: np.array([19])}}
        read_nums = np.array([2, 3, 5, 7, 8, 9, 12, 13, 16, 19, 20])
        seg_end5s = np.array([[1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6]]).T
        seg_end3s = np.array([[0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5]]).T
        mask = seg_end5s > seg_end3s
        result = calc_muts_matrix(region,
                                  read_nums,
                                  seg_end5s,
                                  seg_end3s,
                                  mask,
                                  muts)
        expect = pd.DataFrame([[NOCOV, NOCOV, NOCOV, NOCOV, NOCOV],
                               [DELET, NOCOV, NOCOV, NOCOV, NOCOV],
                               [SUB_A, DELET, NOCOV, NOCOV, NOCOV],
                               [SUB_C, SUB_A, DELET, NOCOV, NOCOV],
                               [SUB_G, SUB_C, SUB_A, DELET, NOCOV],
                               [SUB_T, SUB_G, SUB_C, SUB_A, DELET],
                               [NOCOV, SUB_T, SUB_G, SUB_C, SUB_A],
                               [NOCOV, NOCOV, SUB_T, SUB_G, SUB_C],
                               [NOCOV, NOCOV, NOCOV, SUB_T, SUB_G],
                               [NOCOV, NOCOV, NOCOV, NOCOV, SUB_T],
                               [NOCOV, NOCOV, NOCOV, NOCOV, NOCOV]],
                              read_nums,
                              region.unmasked)
        self.assertTrue(expect.equals(result))


class TestMergeCloseMuts(ut.TestCase):

    region = Region("r", DNA("AAAAA"))  # positions 1-5, all A

    def _make_batch(self, num_reads, sparse_muts):
        read_nums = np.arange(num_reads)
        seg_end5s = np.full((num_reads, 1), self.region.end5, dtype=int)
        seg_end3s = np.full((num_reads, 1), self.region.end3, dtype=int)
        full_muts = {pos: {} for pos in self.region.unmasked_int}
        for pos, rels in sparse_muts.items():
            full_muts[pos] = {rel: np.array(reads) for rel, reads in rels.items()}
        return MaskMutsBatch(batch=0,
                             region=self.region,
                             read_nums=read_nums,
                             seg_end5s=seg_end5s,
                             seg_end3s=seg_end3s,
                             muts=full_muts)

    def _assert_muts_equal(self, result, expected):
        self.assertEqual(set(result.keys()), set(expected.keys()))
        for pos in expected:
            self.assertEqual(set(result[pos].keys()), set(expected[pos].keys()),
                             f"Position {pos}: relationship keys differ")
            for rel in expected[pos]:
                np.testing.assert_array_equal(
                    np.sort(result[pos][rel]),
                    np.sort(np.asarray(expected[pos][rel])),
                    err_msg=f"Position {pos}, rel {rel}: read arrays differ"
                )

    def test_no_mutations(self):
        batch = self._make_batch(3, {})
        result = batch.merge_close_muts(RelPattern.muts(), min_gap=2)
        expected = {pos: {} for pos in self.region.unmasked_int}
        self._assert_muts_equal(result, expected)

    def test_single_mutation_survives(self):
        batch = self._make_batch(3, {3: {SUB_T: [0]}})
        result = batch.merge_close_muts(RelPattern.muts(), min_gap=2)
        expected = {1: {}, 2: {}, 3: {SUB_T: [0]}, 4: {}, 5: {}}
        self._assert_muts_equal(result, expected)

    def test_two_muts_gap_exceeds_min_gap(self):
        # gap = 4-1 = 3 > min_gap=2; both survive
        batch = self._make_batch(3, {1: {SUB_T: [0]}, 4: {SUB_T: [0]}})
        result = batch.merge_close_muts(RelPattern.muts(), min_gap=2)
        expected = {1: {SUB_T: [0]}, 2: {}, 3: {}, 4: {SUB_T: [0]}, 5: {}}
        self._assert_muts_equal(result, expected)

    def test_two_muts_gap_equals_min_gap(self):
        # gap = 4-1 = 3 == min_gap=3; pos 1 is an artifact (condition: 4>1+3=4 is False)
        batch = self._make_batch(3, {1: {SUB_T: [0]}, 4: {SUB_T: [0]}})
        result = batch.merge_close_muts(RelPattern.muts(), min_gap=3)
        expected = {1: {}, 2: {}, 3: {}, 4: {SUB_T: [0]}, 5: {}}
        self._assert_muts_equal(result, expected)

    def test_two_muts_gap_less_than_min_gap(self):
        # gap = 4-2 = 2 < min_gap=3; pos 2 is an artifact
        batch = self._make_batch(3, {2: {SUB_T: [0]}, 4: {SUB_T: [0]}})
        result = batch.merge_close_muts(RelPattern.muts(), min_gap=3)
        expected = {1: {}, 2: {}, 3: {}, 4: {SUB_T: [0]}, 5: {}}
        self._assert_muts_equal(result, expected)

    def test_three_muts_middle_is_artifact(self):
        # Positions 1, 3, 5 with min_gap=3:
        # pos 5 → mod (last=5); pos 3: 5>3+3=6? No → artifact; pos 1: 5>1+3=4? Yes → mod
        batch = self._make_batch(3, {1: {SUB_T: [0]}, 3: {SUB_T: [0]}, 5: {SUB_T: [0]}})
        result = batch.merge_close_muts(RelPattern.muts(), min_gap=3)
        expected = {1: {SUB_T: [0]}, 2: {}, 3: {}, 4: {}, 5: {SUB_T: [0]}}
        self._assert_muts_equal(result, expected)

    def test_reads_processed_independently(self):
        # read 0: muts at pos 3, 5 (gap=2, min_gap=2: 5>3+2=5? No → pos 3 artifact)
        # read 1: muts at pos 1, 5 (gap=4, min_gap=2: 5>1+2=3? Yes → pos 1 survives)
        batch = self._make_batch(3,
                                 {1: {SUB_T: [1]},
                                  3: {SUB_T: [0]},
                                  5: {SUB_T: [0, 1]}})
        result = batch.merge_close_muts(RelPattern.muts(), min_gap=2)
        expected = {1: {SUB_T: [1]}, 2: {}, 3: {}, 4: {}, 5: {SUB_T: [0, 1]}}
        self._assert_muts_equal(result, expected)

    def test_min_gap_zero_keeps_all(self):
        # With min_gap=0, last_mod_pos > pos+0 = pos is always True once set, so nothing filtered
        batch = self._make_batch(3, {1: {SUB_T: [0]}, 2: {SUB_T: [0]}, 3: {SUB_T: [0]}})
        result = batch.merge_close_muts(RelPattern.muts(), min_gap=0)
        expected = {1: {SUB_T: [0]}, 2: {SUB_T: [0]}, 3: {SUB_T: [0]}, 4: {}, 5: {}}
        self._assert_muts_equal(result, expected)

    def test_non_pattern_relationships_preserved(self):
        # RelPattern.from_counts() matches substitutions only (not DELET) at 'A' positions.
        # pos 5: SUB_T for reads [0,1] → both true mods
        # pos 3: SUB_T for read [0] (close artifact: 5>3+3=6? No), DELET for reads [2,3]
        # DELET is not matched by the sub-only pattern, so reads [2,3] are outside the
        # artifact logic and must be preserved unchanged.
        batch = self._make_batch(4,
                                 {3: {DELET: [2, 3], SUB_T: [0]},
                                  5: {SUB_T: [0, 1]}})
        result = batch.merge_close_muts(RelPattern.from_counts(), min_gap=3)
        expected = {1: {}, 2: {}, 3: {DELET: [2, 3]}, 4: {}, 5: {SUB_T: [0, 1]}}
        self._assert_muts_equal(result, expected)


class TestInjectCloseMuts(ut.TestCase):

    region = Region("r", DNA("AAAAA"))  # positions 1-5, all A

    def _make_batch(self, num_reads, sparse_muts):
        read_nums = np.arange(num_reads)
        seg_end5s = np.full((num_reads, 1), self.region.end5, dtype=int)
        seg_end3s = np.full((num_reads, 1), self.region.end3, dtype=int)
        full_muts = {pos: {} for pos in self.region.unmasked_int}
        for pos, rels in sparse_muts.items():
            full_muts[pos] = {rel: np.array(reads) for rel, reads in rels.items()}
        return MaskMutsBatch(batch=0,
                             region=self.region,
                             read_nums=read_nums,
                             seg_end5s=seg_end5s,
                             seg_end3s=seg_end3s,
                             muts=full_muts)

    def _reads_with_any_mut(self, muts, pos):
        arrays = list(muts[pos].values())
        if not arrays:
            return np.array([], dtype=int)
        return np.sort(np.unique(np.concatenate(arrays)))

    def _assert_muts_equal(self, result, expected):
        self.assertEqual(set(result.keys()), set(expected.keys()))
        for pos in expected:
            self.assertEqual(set(result[pos].keys()), set(expected[pos].keys()),
                             f"Position {pos}: relationship keys differ")
            for rel in expected[pos]:
                np.testing.assert_array_equal(
                    np.sort(result[pos][rel]),
                    np.sort(np.asarray(expected[pos][rel])),
                    err_msg=f"Position {pos}, rel {rel}: read arrays differ"
                )

    def test_empty_mut_probs(self):
        batch = self._make_batch(3, {3: {SUB_T: [0, 1, 2]}})
        result = batch.inject_close_muts(RelPattern.muts(), [], seed=0)
        expected = {1: {}, 2: {}, 3: {SUB_T: [0, 1, 2]}, 4: {}, 5: {}}
        self._assert_muts_equal(result, expected)

    def test_all_zero_probs(self):
        batch = self._make_batch(3, {3: {SUB_T: [0, 1, 2]}})
        result = batch.inject_close_muts(RelPattern.muts(), [0.0], seed=0)
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 2), [])

    def test_invalid_ndim(self):
        batch = self._make_batch(2, {})
        with self.assertRaises(ValueError):
            batch.inject_close_muts(RelPattern.muts(), np.array([[0.5]]), seed=0)

    def test_invalid_negative(self):
        batch = self._make_batch(2, {})
        with self.assertRaises(ValueError):
            batch.inject_close_muts(RelPattern.muts(), [-0.1], seed=0)

    def test_invalid_above_one(self):
        batch = self._make_batch(2, {})
        with self.assertRaises(ValueError):
            batch.inject_close_muts(RelPattern.muts(), [1.1], seed=0)

    def test_prob_one_injects_all(self):
        # All 3 reads have a mutation at pos 3; with prob=1.0 all should get one at pos 2.
        batch = self._make_batch(3, {3: {SUB_T: [0, 1, 2]}})
        result = batch.inject_close_muts(RelPattern.muts(), [1.0], seed=0)
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 2), [0, 1, 2])

    def test_original_mutation_preserved(self):
        batch = self._make_batch(3, {3: {SUB_T: [0, 1, 2]}})
        result = batch.inject_close_muts(RelPattern.muts(), [1.0], seed=0)
        np.testing.assert_array_equal(np.sort(result[3][SUB_T]), [0, 1, 2])

    def test_only_mutated_reads_injected(self):
        # Only read 0 is mutated at pos 3; reads 1-2 should not appear at pos 2.
        batch = self._make_batch(3, {3: {SUB_T: [0]}})
        result = batch.inject_close_muts(RelPattern.muts(), [1.0], seed=0)
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 2), [0])
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 1), [])

    def test_window_of_two(self):
        # Mutation at pos 4; window of 2 covers pos 3 and pos 2.
        batch = self._make_batch(2, {4: {SUB_T: [0, 1]}})
        result = batch.inject_close_muts(RelPattern.muts(), [1.0, 1.0], seed=0)
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 3), [0, 1])
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 2), [0, 1])
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 1), [])

    def test_read_not_covering_skipped(self):
        # Read 0 spans positions 4-5; read 1 spans 1-5.
        # Both mutated at pos 5 with prob=1.0 window of 4.
        # pos 4: both reads cover it (4<=4 for read 0, 1<=4 for read 1).
        # pos 3: read 0 does NOT cover it (4<=3 is False); read 1 does.
        read_nums = np.array([0, 1])
        seg_end5s = np.array([[4], [1]])
        seg_end3s = np.array([[5], [5]])
        full_muts = {1: {}, 2: {}, 3: {}, 4: {},
                     5: {SUB_T: np.array([0, 1])}}
        batch = MaskMutsBatch(batch=0,
                              region=self.region,
                              read_nums=read_nums,
                              seg_end5s=seg_end5s,
                              seg_end3s=seg_end3s,
                              muts=full_muts)
        result = batch.inject_close_muts(RelPattern.muts(), [1.0, 1.0, 1.0, 1.0], seed=0)
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 4), [0, 1])
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 3), [1])
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 2), [1])
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 1), [1])

    def test_masked_position_skipped(self):
        # Mask pos 2; mutation at pos 3 with window of 2 (covers pos 2 and pos 1).
        # pos 2 is masked so injection skips it; pos 1 is unmasked and receives injection.
        region = Region("r", DNA("AAAAA"))
        region.mask_list([2])
        num_reads = 2
        read_nums = np.arange(num_reads)
        seg_end5s = np.full((num_reads, 1), region.end5, dtype=int)
        seg_end3s = np.full((num_reads, 1), region.end3, dtype=int)
        full_muts = {pos: {} for pos in region.unmasked_int}
        full_muts[3] = {SUB_T: np.array([0, 1])}
        batch = MaskMutsBatch(batch=0,
                              region=region,
                              read_nums=read_nums,
                              seg_end5s=seg_end5s,
                              seg_end3s=seg_end3s,
                              muts=full_muts)
        result = batch.inject_close_muts(RelPattern.muts(), [1.0, 1.0], seed=0)
        self.assertNotIn(2, result)
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 1), [0, 1])

    def test_existing_mutation_not_duplicated(self):
        # Read 0 already has SUB_T at pos 2 and pos 3; read 1 only at pos 3.
        # Injection skips read 0 at pos 2 (already mutated there), so read 0
        # stays in SUB_T and read 1 gets a new injection — total == 2.
        batch = self._make_batch(2, {2: {SUB_T: [0]}, 3: {SUB_T: [0, 1]}})
        result = batch.inject_close_muts(RelPattern.muts(), [1.0], seed=0)
        np.testing.assert_array_equal(self._reads_with_any_mut(result, 2), [0, 1])
        total = sum(len(v) for v in result[2].values())
        self.assertEqual(total, 2)


if __name__ == "__main__":
    ut.main(verbosity=2)
