import shutil
import unittest as ut
from functools import reduce
from itertools import product
from operator import add
from pathlib import Path

import numpy as np
import pandas as pd

from seismicrna.cluster.data import ClusterMutsDataset
from seismicrna.core.arg.cli import opt_sim_dir
from seismicrna.core.logs import Level, get_config, set_config
from seismicrna.core.rna.convert import db_to_ct
from seismicrna.core.seq.fasta import write_fasta
from seismicrna.core.seq.region import FULL_NAME
from seismicrna.core.seq.xna import DNA
from seismicrna.ensembles import (_calc_tiles,
                                  _aggregate_pairs,
                                  _select_pairs,
                                  _calc_span_per_pos,
                                  _calc_null_span_per_pos_keep_dists,
                                  _calc_null_span_per_pos_rand_dists,
                                  _calc_modules_from_pairs,
                                  _insert_modules_into_gaps,
                                  _expand_modules_into_gaps,
                                  _filter_modules_length,
                                  run as run_ensembles)
from seismicrna.sim.params import run as sim_params
from seismicrna.sim.relate import run as sim_relate

rng = np.random.default_rng()


class TestCalcTiles(ut.TestCase):

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR)

    def tearDown(self):
        set_config(self._config)

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


class TestCalcSpanPerPos(ut.TestCase):

    def test_calc_span_per_pos(self):
        pairs = [(5, 10), (9, 9), (3, 11)]
        result = _calc_span_per_pos(pairs, 3, 11)
        expect = pd.Series({3: 1,
                            4: 1,
                            5: 2,
                            6: 2,
                            7: 2,
                            8: 2,
                            9: 3,
                            10: 2,
                            11: 1})
        self.assertIsInstance(result, pd.Series)
        self.assertTrue(result.equals(expect))
        self.assertTrue(result.index.equals(expect.index))


class TestCalcNullSpanPerPosKeepDists(ut.TestCase):

    def _run_test(self,
                  pairs: list[tuple[int, int]],
                  end5: int,
                  end3: int,
                  expect: pd.Series):
        result = _calc_null_span_per_pos_keep_dists(pairs, end5, end3)
        self.assertIsInstance(result, pd.Series)
        self.assertTrue(np.allclose(result, expect))
        self.assertTrue(result.index.equals(expect.index))

    def test_no_pairs(self):
        expect = pd.Series({4: 0., 5: 0., 6: 0., 7: 0.})
        self._run_test([], 4, 7, expect)

    def test_ref_1_read_1(self):
        expect = pd.Series({5: 1.})
        self._run_test([(5, 5)], 5, 5, expect)

    def test_ref_4_read_1(self):
        expect = pd.Series({4: 0.25, 5: 0.25, 6: 0.25, 7: 0.25})
        self._run_test([(5, 5)], 4, 7, expect)

    def test_ref_4_read_2(self):
        expect = pd.Series({4: 1 / 3, 5: 2 / 3, 6: 2 / 3, 7: 1 / 3})
        self._run_test([(6, 7)], 4, 7, expect)

    def test_ref_4_read_3(self):
        expect = pd.Series({4: 0.5, 5: 1., 6: 1., 7: 0.5})
        self._run_test([(4, 6)], 4, 7, expect)

    def test_ref_4_read_4(self):
        expect = pd.Series({4: 1., 5: 1., 6: 1., 7: 1.})
        self._run_test([(4, 7)], 4, 7, expect)

    def test_ref_5_read_1(self):
        expect = pd.Series({4: 0.2, 5: 0.2, 6: 0.2, 7: 0.2, 8: 0.2})
        self._run_test([(7, 7)], 4, 8, expect)

    def test_ref_5_read_2(self):
        expect = pd.Series({4: 0.25, 5: 0.5, 6: 0.5, 7: 0.5, 8: 0.25})
        self._run_test([(5, 6)], 4, 8, expect)

    def test_ref_5_read_3(self):
        expect = pd.Series({4: 1 / 3, 5: 2 / 3, 6: 1., 7: 2 / 3, 8: 1 / 3})
        self._run_test([(5, 7)], 4, 8, expect)

    def test_ref_5_read_4(self):
        expect = pd.Series({4: 0.5, 5: 1., 6: 1., 7: 1., 8: 0.5})
        self._run_test([(4, 7)], 4, 8, expect)

    def test_ref_5_read_5(self):
        expect = pd.Series({4: 1., 5: 1., 6: 1., 7: 1., 8: 1.})
        self._run_test([(4, 8)], 4, 8, expect)

    def test_multiple_reads(self):
        pairs = [(6, 7), (5, 5), (4, 8)]
        expect = reduce(add,
                        [_calc_null_span_per_pos_keep_dists([pair], 4, 8)
                         for pair in pairs])
        self._run_test(pairs, 4, 8, expect)


class TestCalcNullSpanPerPosRandDists(ut.TestCase):

    @staticmethod
    def _calc_expect(pairs: list[tuple[int, int]],
                     end5: int,
                     end3: int,
                     min_mut_gap: int):
        positions = np.arange(end5, end3 + 1)
        counts = pd.Series(0, index=positions)
        num_intervals = 0
        for a in positions:
            for b in positions:
                gap = b - a - 1
                if gap >= min_mut_gap:
                    counts.loc[a: b] += 1
                    num_intervals += 1
        if num_intervals > 0:
            factor = len(pairs) / num_intervals
        else:
            factor = 0.
        return factor * counts

    def _compare(self, result: pd.Series, expect: pd.Series):
        self.assertIsInstance(result, pd.Series)
        self.assertTrue(np.allclose(result, expect))
        self.assertTrue(result.index.equals(expect.index))

    def _run_test(self,
                  pairs: list[tuple[int, int]],
                  end5: int,
                  end3: int,
                  min_mut_gap: int):
        result = _calc_null_span_per_pos_rand_dists(pairs,
                                                    end5,
                                                    end3,
                                                    min_mut_gap)
        expect = self._calc_expect(pairs,
                                   end5,
                                   end3,
                                   min_mut_gap)
        self._compare(result, expect)

    def test_no_pairs(self):
        self._run_test([], 4, 7, 0)

    def test_1_pair(self):
        for end5 in range(1, 4):
            for end3 in range(end5, 7):
                for min_mut_gap in range(5):
                    self._run_test([(end5, end3)],
                                   end5,
                                   end3,
                                   min_mut_gap)

    def test_2_pairs(self):
        end5, end3 = 2, 7
        min_mut_gap = 0
        expect = (2 / 15) * pd.Series([5, 9, 11, 11, 9, 5],
                                      index=range(end5, end3 + 1))
        result = _calc_null_span_per_pos_rand_dists([(end5, end3),
                                                     (end5, end3)],
                                                    end5, end3, min_mut_gap)
        self._compare(result, expect)
        min_mut_gap = 1
        expect = (2 / 10) * pd.Series([4, 7, 9, 9, 7, 4],
                                      index=range(end5, end3 + 1))
        result = _calc_null_span_per_pos_rand_dists([(end5, end3),
                                                     (end5, end3)],
                                                    end5, end3, min_mut_gap)
        self._compare(result, expect)
        min_mut_gap = 2
        expect = (2 / 6) * pd.Series([3, 5, 6, 6, 5, 3],
                                     index=range(end5, end3 + 1))
        result = _calc_null_span_per_pos_rand_dists([(end5, end3),
                                                     (end5, end3)],
                                                    end5, end3, min_mut_gap)
        self._compare(result, expect)
        min_mut_gap = 3
        expect = (2 / 3) * pd.Series([2, 3, 3, 3, 3, 2],
                                     index=range(end5, end3 + 1))
        result = _calc_null_span_per_pos_rand_dists([(end5, end3),
                                                     (end5, end3)],
                                                    end5, end3, min_mut_gap)
        self._compare(result, expect)
        min_mut_gap = 4
        expect = pd.Series(2, index=range(end5, end3 + 1))
        result = _calc_null_span_per_pos_rand_dists([(end5, end3),
                                                     (end5, end3)],
                                                    end5, end3, min_mut_gap)
        self._compare(result, expect)
        min_mut_gap = 5
        expect = pd.Series(0, index=range(end5, end3 + 1))
        result = _calc_null_span_per_pos_rand_dists([(end5, end3),
                                                     (end5, end3)],
                                                    end5, end3, min_mut_gap)
        self._compare(result, expect)


class TestCalcModulesFromPairs(ut.TestCase):

    def test_1_module(self):
        pairs = [(11, 21),
                 (11, 31),
                 (11, 71),
                 (21, 31),
                 (51, 61),
                 (51, 71),
                 (61, 71)]
        result = _calc_modules_from_pairs(pairs, 0.05, 4)
        expect = [(11, 71)]
        self.assertListEqual(result, expect)

    def test_2_modules(self):
        pairs = [(11, 21),
                 (11, 31),
                 (21, 31),
                 (51, 61),
                 (51, 71),
                 (61, 71)]
        result = _calc_modules_from_pairs(pairs, 0.05, 4)
        expect = [(11, 31), (51, 71)]
        self.assertListEqual(result, expect)

    def test_split_modules(self):
        pairs = [(11, 21),
                 (11, 31),
                 (11, 71),
                 (21, 31),
                 (51, 61),
                 (51, 71),
                 (61, 71)]
        result = _calc_modules_from_pairs(pairs, 0.35, 4)
        expect = [(11, 31), (51, 71)]
        self.assertListEqual(result, expect)


class TestFilterModulesLength(ut.TestCase):

    def test_filter_default_length(self):
        modules = [(5, 20), (25, 30)]
        result = _filter_modules_length(modules)
        expect = modules
        self.assertListEqual(result, expect)

    def test_filter_min_length(self):
        modules = [(5, 20), (25, 30)]
        result = _filter_modules_length(modules, min_length=6)
        expect = modules
        self.assertListEqual(result, expect)
        result = _filter_modules_length(modules, min_length=7)
        expect = [(5, 20)]
        self.assertListEqual(result, expect)
        result = _filter_modules_length(modules, min_length=16)
        expect = [(5, 20)]
        self.assertListEqual(result, expect)
        result = _filter_modules_length(modules, min_length=17)
        expect = []
        self.assertListEqual(result, expect)

    def test_filter_max_length(self):
        modules = [(5, 20), (25, 30)]
        result = _filter_modules_length(modules, max_length=5)
        expect = []
        self.assertListEqual(result, expect)
        result = _filter_modules_length(modules, max_length=6)
        expect = [(25, 30)]
        self.assertListEqual(result, expect)
        result = _filter_modules_length(modules, max_length=15)
        expect = [(25, 30)]
        self.assertListEqual(result, expect)
        result = _filter_modules_length(modules, max_length=16)
        expect = modules
        self.assertListEqual(result, expect)


class TestInsertRegionsIntoGaps(ut.TestCase):

    def test_zero(self):
        result = _insert_modules_into_gaps([], 3, 9)
        expect = [(3, 9)]
        self.assertListEqual(result, expect)

    def test_one(self):
        result = _insert_modules_into_gaps([(3, 9)], 3, 9)
        expect = [(3, 9)]
        self.assertListEqual(result, expect)
        result = _insert_modules_into_gaps([(4, 9)], 3, 9)
        expect = [(3, 3), (4, 9)]
        self.assertListEqual(result, expect)
        result = _insert_modules_into_gaps([(3, 8)], 3, 9)
        expect = [(3, 8), (9, 9)]
        self.assertListEqual(result, expect)
        result = _insert_modules_into_gaps([(4, 8)], 3, 9)
        expect = [(3, 3), (4, 8), (9, 9)]
        self.assertListEqual(result, expect)

    def test_two(self):
        result = _insert_modules_into_gaps([(2, 10), (11, 20)], 2, 20)
        expect = [(2, 10), (11, 20)]
        self.assertListEqual(result, expect)
        result = _insert_modules_into_gaps([(3, 10), (11, 20)], 2, 20)
        expect = [(2, 2), (3, 10), (11, 20)]
        self.assertListEqual(result, expect)
        result = _insert_modules_into_gaps([(2, 10), (12, 20)], 2, 20)
        expect = [(2, 10), (11, 11), (12, 20)]
        self.assertListEqual(result, expect)
        result = _insert_modules_into_gaps([(2, 10), (11, 19)], 2, 20)
        expect = [(2, 10), (11, 19), (20, 20)]
        self.assertListEqual(result, expect)
        result = _insert_modules_into_gaps([(3, 9), (11, 19)], 2, 20)
        expect = [(2, 2), (3, 9), (10, 10), (11, 19), (20, 20)]
        self.assertListEqual(result, expect)


class TestExpandRegionsIntoGaps(ut.TestCase):

    def test_zero(self):
        result = _expand_modules_into_gaps([], 3, 9)
        expect = []
        self.assertListEqual(result, expect)

    def test_one(self):
        expect = [(3, 9)]
        result = _expand_modules_into_gaps([(3, 9)], 3, 9)
        self.assertListEqual(result, expect)
        result = _expand_modules_into_gaps([(4, 9)], 3, 9)
        self.assertListEqual(result, expect)
        result = _expand_modules_into_gaps([(3, 8)], 3, 9)
        self.assertListEqual(result, expect)
        result = _expand_modules_into_gaps([(5, 7)], 3, 9)
        self.assertListEqual(result, expect)

    def test_two(self):
        result = _expand_modules_into_gaps([(2, 10), (11, 20)], 2, 20)
        expect = [(2, 10), (11, 20)]
        self.assertListEqual(result, expect)
        result = _expand_modules_into_gaps([(4, 9), (11, 18)], 2, 20)
        expect = [(2, 9), (10, 20)]
        self.assertListEqual(result, expect)
        result = _expand_modules_into_gaps([(5, 9), (12, 17)], 2, 20)
        expect = [(2, 10), (11, 20)]
        self.assertListEqual(result, expect)
        result = _expand_modules_into_gaps([(6, 6), (10, 10)], 2, 20)
        expect = [(2, 7), (8, 20)]
        self.assertListEqual(result, expect)


class TestEnsembles(ut.TestCase):
    SIM_DIR = Path(opt_sim_dir.default).absolute()
    SAMPLE = "test_sample"
    REFS = "test_refs"
    REF = "test_ref"
    PROFILE = "ensembles"

    # Folding modules of the reference sequence (each 60 nt).
    MODULES = [
        ("TGACGAACAACGTGTTTGTGAACCATATAGGTAAACGCTGAATGCGTTCGCGCGGAGGGT",
         ["..(((((((...)))))))..(((......((.(((((.....))))).))......)))",
          "...(((...(((((((((((.(((.....)))...))).))))))))))).((.....))"]),
        ("TTTGCAGGAAGATGGTCAACTCTACACCTAGTTTTTACCAGTCCACAAGAGTTTGAACTG",
         [".(((..(((...((((.((..(((....)))..)).)))).))).))).(((....)))."]),
        ("GTGCCTTAACCTGAGTACGCCCATATCATGGGAGACATTACAACTCAAATTCTAGGTGTG",
         ["..((((.....(((((...(((((...)))))..........)))))......))))...",
          "((((.(((...)))))))...((((((.((((..................))))))))))"]),
    ]

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR,
                   log_file_path=None,
                   exit_on_error=True)
        self.SIM_DIR.mkdir()

    def tearDown(self):
        shutil.rmtree(self.SIM_DIR)
        set_config(*self._config)

    @classmethod
    def sim_data(cls, module_nums: list[int], read_length: int):
        # Assemble and write the reference sequence.
        modules = dict(cls.MODULES[m] for m in module_nums)
        refseq = DNA("".join(modules.keys()))
        refs_dir = cls.SIM_DIR.joinpath("refs")
        refs_dir.mkdir()
        fasta = refs_dir.joinpath(f"{cls.REFS}.fa")
        write_fasta(fasta, [(cls.REF, refseq)])
        # Assemble and write the secondary structures.
        structures = list(map("".join, product(*modules.values())))
        param_dir = cls.SIM_DIR.joinpath("params", cls.REF, FULL_NAME)
        param_dir.mkdir(parents=True)
        db_file = param_dir.joinpath(f"{cls.PROFILE}.db")
        with open(db_file, "x") as f:
            for i, struct in enumerate(structures):
                if i == 0:
                    f.write(f">structure0\n{refseq.tr()}\n{structures[0]}\n")
                else:
                    f.write(f">structure{i}\n{structures[i]}\n")
        ct_file = db_to_ct(db_file)
        # Simulate data.
        sim_params(ct_file=[ct_file],
                   # Make pmut_unpaired for A and C large (17%) so that
                   # most reads get at least two mutations despite being
                   # short and are thus useful for clustering.
                   pmut_unpaired=[("am", 1 / 6), ("cm", 1 / 6)],
                   # Make all reads the same length.
                   length_fmean=(read_length / len(refseq)),
                   length_fvar=0.,
                   # Make clust_conc very large so that the proportion
                   # of each cluster is approximately equal, which makes
                   # clustering easier.
                   clust_conc=1000.)
        relate_dirs = sim_relate(param_dir=[param_dir],
                                 sample=cls.SAMPLE,
                                 profile_name=cls.PROFILE,
                                 num_reads=200000,
                                 paired_end=False,
                                 brotli_level=0)
        return relate_dirs

    def run_ensembles(self,
                      relate_dirs: list[Path],
                      expect_regions: dict[tuple[int, int], int],
                      tolerance: int = 0,
                      **kwargs):
        cluster_dirs = {tuple(map(int, d.name.split("-"))): d
                        for d in run_ensembles(relate_dirs,
                                               # Optimize for speed.
                                               min_em_runs=1,
                                               max_em_runs=1,
                                               jackpot=False,
                                               brotli_level=0,
                                               mask_pos_table=False,
                                               mask_read_table=False,
                                               cluster_pos_table=False,
                                               cluster_abundance_table=False,
                                               **kwargs)}
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
                raise ValueError(f"Expected region {exp5, exp3} does not "
                                 "overlap at least 50% of any region "
                                 f"among {sorted(cluster_dirs)}")

    def test_modules012_read180(self):
        relate_dirs = self.sim_data([0, 1, 2], 180)
        self.run_ensembles(relate_dirs,
                           {(1, 60): 2,
                            (121, 180): 2},
                           tolerance=60)

    def test_modules012_read120(self):
        relate_dirs = self.sim_data([0, 1, 2], 120)
        self.run_ensembles(relate_dirs,
                           {(1, 60): 2,
                            (121, 180): 2},
                           tolerance=60)

    def test_modules012_read60(self):
        relate_dirs = self.sim_data([0, 1, 2], 60)
        self.run_ensembles(relate_dirs,
                           {(1, 60): 2,
                            (121, 180): 2},
                           tolerance=60)

    def test_modules02_read60(self):
        relate_dirs = self.sim_data([0, 2], 60)
        self.run_ensembles(relate_dirs,
                           {(1, 60): 2,
                            (61, 120): 2},
                           tolerance=30)


if __name__ == "__main__":
    ut.main(verbosity=2)
