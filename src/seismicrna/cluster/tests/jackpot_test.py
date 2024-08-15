import shutil
import unittest as ut
from pathlib import Path

import numpy as np
from scipy.stats import t as studentt

from seismicrna.cluster.em import EMRun
from seismicrna.cluster.jackpot import calc_jackpot_g_stat
from seismicrna.cluster.uniq import UniqReads
from seismicrna.core.arg.cli import (opt_sim_dir,
                                     opt_max_em_iter,
                                     opt_em_thresh,
                                     opt_jackpot_conf_level,
                                     opt_max_jackpot_index)
from seismicrna.core.logs import get_config, set_config
from seismicrna.mask import run as run_mask
from seismicrna.mask.data import MaskMutsDataset
from seismicrna.mask.report import MaskReport, NumReadsInitF, NumReadsKeptF
from seismicrna.sim.fold import run as run_sim_fold
from seismicrna.sim.params import run as run_sim_params
from seismicrna.sim.ref import run as run_sim_ref
from seismicrna.sim.relate import run as run_sim_relate

rng = np.random.default_rng()


class TestCalcJackpotGStat(ut.TestCase):

    def test_diff_lengths(self):
        num_obs = np.array([1])
        num_exp = np.array([1., 2.])
        self.assertRaisesRegex(ValueError,
                               (r"Lengths differ between observed \(1\) "
                                r"and expected \(2\)"),
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_neg_num_obs(self):
        num_obs = np.array([-1, 2])
        num_exp = np.array([1., 2.])
        self.assertRaisesRegex(ValueError,
                               r"All num_obs must be ≥ 0, but got \[-1\]",
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_neg_min_exp(self):
        num_obs = np.array([], dtype=int)
        num_exp = np.array([])
        min_exp = -0.1
        self.assertRaisesRegex(ValueError,
                               r"min_exp must be ≥ 0, but got -0.1",
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp),
                               min_exp)

    def test_each_exp_below_min_exp_all_exp_equal_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([2., 3., 1.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 0)

    def test_each_exp_below_min_exp_sum_exp_equal_sum_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([1., 2., 3.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 0)

    def test_each_exp_below_min_exp_sum_exp_less_than_sum_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([1., 2., 2.5])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 1)

    def test_each_exp_below_min_exp_sum_exp_greater_than_sum_obs(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([1., 2., 3.5])
        self.assertRaisesRegex(ValueError,
                               (r"Total observed reads \(6\) is less than "
                                r"total expected reads \(6.5\)"),
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_each_exp_at_least_min_exp_all_exp_equal_obs(self):
        min_obs = 4
        for n in range(5):
            num_obs = np.arange(min_obs, min_obs + n)
            num_exp = np.asarray(num_obs, dtype=float)
            g_stat, df = calc_jackpot_g_stat(num_obs,
                                             np.log(num_exp),
                                             float(min_obs))
            self.assertIsInstance(g_stat, float)
            self.assertEqual(g_stat, 0.)
            self.assertIsInstance(df, int)
            self.assertEqual(df, max(n - 1, 0))

    def test_each_exp_at_least_min_exp_sum_exp_equal_sum_obs(self):
        num_obs = np.array([5, 4, 10])
        num_exp = np.array([4., 5., 10.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   sum([10. * np.log(1.25),
                                        8. * np.log(0.8)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 2)

    def test_each_exp_at_least_min_exp_sum_exp_less_than_sum_obs(self):
        num_obs = np.array([5, 4, 11])
        num_exp = np.array([4., 5., 10.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   2. * sum([5 * np.log(5 / 4),
                                             4 * np.log(4 / 5),
                                             11 * np.log(11 / 10)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 3)

    def test_each_exp_at_least_min_exp_sum_exp_greater_than_sum_obs(self):
        num_obs = np.array([5, 4, 10])
        num_exp = np.array([4., 5., 11.])
        self.assertRaisesRegex(ValueError,
                               (r"Total observed reads \(19\) is less than "
                                r"total expected reads \(20.0\)"),
                               calc_jackpot_g_stat,
                               num_obs,
                               np.log(num_exp))

    def test_mixed_all_exp_equal_obs(self):
        num_obs = np.array([2, 4, 5, 8, 3])
        num_exp = np.asarray(num_obs, dtype=float)
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertEqual(g_stat, 0.)
        self.assertIsInstance(df, int)
        self.assertEqual(df, 3)

    def test_mixed_sum_exp_equal_sum_obs(self):
        num_obs = np.array([0, 3, 7, 6, 6])
        num_exp = np.array([2., 4., 5., 8., 3.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   2. * sum([3 * np.log(3 / 4),
                                             7 * np.log(7 / 5),
                                             6 * np.log(6 / 8),
                                             6 * np.log(6 / 5)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 3)

    def test_mixed_sum_exp_under_min_exp_less_than_sum_obs(self):
        num_obs = np.array([0, 3, 7, 6, 6])
        num_exp = np.array([1., 4., 5., 8., 3.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   2. * sum([3 * np.log(3 / 4),
                                             7 * np.log(7 / 5),
                                             6 * np.log(6 / 8),
                                             6 * np.log(6 / 5)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 4)

    def test_mixed_sum_exp_over_min_exp_less_than_sum_obs(self):
        num_obs = np.array([0, 3, 5, 6, 8])
        num_exp = np.array([2., 4., 5., 7., 3.])
        min_exp = 4.
        g_stat, df = calc_jackpot_g_stat(num_obs, np.log(num_exp), min_exp)
        self.assertIsInstance(g_stat, float)
        self.assertTrue(np.isclose(g_stat,
                                   2. * sum([3 * np.log(3 / 4),
                                             5 * np.log(5 / 5),
                                             6 * np.log(6 / 7),
                                             8 * np.log(8 / 6)])))
        self.assertIsInstance(df, int)
        self.assertEqual(df, 4)


class TestBootstrapNormGStats(ut.TestCase):
    SIM_DIR = Path(opt_sim_dir.default).absolute()
    REFS = "test_refs"
    REF = "test_ref"
    SAMPLE = "test_sample"
    CONF_LEVEL = 0.95
    CONF_WIDTH = 0.02

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._config = None

    def setUp(self):
        self.SIM_DIR.mkdir()
        self._config = get_config()
        set_config(verbose=2, quiet=0, log_file=None, raise_on_error=True)

    def tearDown(self):
        if self.SIM_DIR.exists():
            shutil.rmtree(self.SIM_DIR)
        set_config(**self._config._asdict())

    def simulate_jackpot_index(self):
        """ Simulate a dataset and calculate its jackpotting index. """
        n_pos = 120
        n_reads = 10000
        n_clusts = 2
        min_mut_gap = 3
        # Simulate "ideal" data with no low-quality bases or deletions.
        run_sim_ref(refs=self.REFS, ref=self.REF, reflen=n_pos)
        fasta = self.SIM_DIR.joinpath("refs", f"{self.REFS}.fa")
        run_sim_fold(fasta, fold_max=n_clusts)
        param_dir = self.SIM_DIR.joinpath("params", self.REF, "full")
        ct_file = param_dir.joinpath("simulated.ct")
        pmut = [("loq", 0.),
                ("ac", 0.25),
                ("ag", 0.25),
                ("at", 0.50),
                ("ca", 0.25),
                ("cg", 0.50),
                ("ct", 0.25),
                ("ga", 0.25),
                ("gc", 0.50),
                ("gt", 0.25),
                ("ta", 0.50),
                ("tc", 0.25),
                ("tg", 0.25)]
        run_sim_params(ct_file=(ct_file,),
                       pmut_paired=pmut,
                       pmut_unpaired=pmut,
                       insert_fmean=1.,
                       end3_fmean=1.,
                       clust_conc=2.)
        relate_report_file, = run_sim_relate(param_dir=(param_dir,),
                                             sample=self.SAMPLE,
                                             min_mut_gap=min_mut_gap,
                                             num_reads=n_reads)
        # Mask the data.
        mask_report_file, = run_mask((relate_report_file,),
                                     mask_polya=0,
                                     mask_del=False,
                                     mask_ins=False,
                                     min_mut_gap=min_mut_gap,
                                     quick_unbias_thresh=0.)
        mask_report = MaskReport.load(mask_report_file)
        self.assertEqual(mask_report.get_field(NumReadsInitF),
                         mask_report.get_field(NumReadsKeptF))
        # Cluster the data and calculate the jackpotting index.
        mask_dataset = MaskMutsDataset.load(mask_report_file)
        uniq_reads = UniqReads.from_dataset_contig(mask_dataset)
        em_run = EMRun(uniq_reads,
                       k=n_clusts,
                       seed=rng.integers(2 ** 32),
                       min_iter=2,
                       max_iter=opt_max_em_iter.default,
                       em_thresh=opt_em_thresh.default,
                       jackpot=True,
                       jackpot_conf_level=opt_jackpot_conf_level.default,
                       max_jackpot_index=opt_max_jackpot_index.default)
        # Delete the simulated files so that this function can run again
        # if necessary.
        shutil.rmtree(self.SIM_DIR)
        return em_run.jackpot_index

    def calc_confidence_interval(self, log_jackpot_indexes: list[float]):
        n = len(log_jackpot_indexes)
        if n <= 1:
            return np.nan, np.nan
        mean = np.mean(log_jackpot_indexes)
        std_err = np.std(log_jackpot_indexes) / np.sqrt(n)
        t_lo, t_up = studentt.interval(self.CONF_LEVEL, n - 1)
        j_lo = mean + std_err * t_lo
        j_up = mean + std_err * t_up
        return j_lo, j_up

    def test_ideal_jackpot(self):
        """ Test that bootstrapping "perfect" data correctly returns a
        jackpotting index that is expected to be 0. """
        log_jackpot_indexes = list()
        while True:
            log_jackpot_indexes.append(np.log(self.simulate_jackpot_index()))
            j_lo, j_up = self.calc_confidence_interval(log_jackpot_indexes)
            width = j_up - j_lo
            if width <= self.CONF_WIDTH:
                # Verify that the confidence interval contains 0.
                self.assertLessEqual(j_lo, 0.)
                self.assertGreaterEqual(j_up, 0.)
                break


if __name__ == "__main__":
    ut.main()
