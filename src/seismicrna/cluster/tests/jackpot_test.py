import shutil
import unittest as ut
from pathlib import Path

import numpy as np
from scipy.stats import binom, chi2, t as studentt

from seismicrna.cluster.em import EMRun
from seismicrna.cluster.jackpot import (calc_semi_g_anomaly,
                                        linearize_ends_matrix,
                                        _sim_reads)
from seismicrna.cluster.uniq import UniqReads
from seismicrna.core.arg.cli import (opt_sim_dir,
                                     opt_max_em_iter,
                                     opt_em_thresh,
                                     opt_jackpot_conf_level,
                                     opt_max_jackpot_quotient)
from seismicrna.core.array import find_dims
from seismicrna.core.logs import get_config, set_config
from seismicrna.core.unbias import (READS,
                                    calc_p_ends_given_clust_noclose,
                                    calc_p_noclose_given_clust,
                                    calc_p_clust_given_noclose,
                                    calc_p_clust_given_ends_noclose,
                                    calc_p_nomut_window,
                                    calc_p_noclose_given_ends,
                                    calc_p_ends_given_noclose,
                                    calc_p_mut_given_span_noclose,
                                    calc_rectangular_sum)
from seismicrna.mask import run as run_mask
from seismicrna.mask.data import MaskMutsDataset
from seismicrna.sim.fold import run as run_sim_fold
from seismicrna.sim.params import run as run_sim_params
from seismicrna.sim.ref import run as run_sim_ref
from seismicrna.sim.relate import run as run_sim_relate

rng = np.random.default_rng()


def g_test(obs: np.ndarray, exp: np.ndarray):
    """ Perform a G-test of the observed values. """
    if exp.shape != obs.shape:
        raise ValueError(
            f"obs and exp have different dimensions: {obs.shape} ≠ {exp.shape}"
        )
    if not np.isclose(obs.sum(), exp.sum()):
        raise ValueError(
            f"obs and exp have different sums: {obs.sum()} ≠ {exp.sum()}"
        )
    n, = obs.shape
    if n >= 2:
        # Calculate the G-test statistic and P-value.
        g_stat = 2. * np.sum(obs * np.log(obs / exp))
        p_value = 1. - chi2.cdf(g_stat, n - 1)
    else:
        # For fewer than 2 items, there can be no difference between
        # observed and expected counts.
        g_stat = 0.
        p_value = 1.
    return g_stat, p_value


class TestSimReads(ut.TestCase):

    @staticmethod
    def count_reads(reads: np.ndarray, clusts: np.ndarray, n_clust: int):
        """ Count reads for each cluster and position. """
        dims = find_dims([(READS, "positions + 2"), (READS,)],
                         [reads, clusts],
                         ["reads", "clusts"])
        n_pos = dims["positions + 2"] - 2
        clust_counts = np.zeros(n_clust, dtype=int)
        span_counts = np.zeros((n_pos, n_clust), dtype=int)
        mut_counts = np.zeros((n_pos, n_clust), dtype=int)
        for k in range(n_clust):
            reads_k = reads[clusts == k]
            clust_counts[k], _ = reads_k.shape
            for read in reads_k:
                span_counts[read[-2]: read[-1] + 1, k] += 1
            mut_counts[:, k] = reads_k[:, :-2].sum(axis=0)
        return clust_counts, span_counts, mut_counts

    def test_sim_reads(self):
        n_pos = 10
        n_reads = 100000
        n_clust = 2
        min_mut_gap = 3
        beta_a = 5.
        beta_b = 15.
        cluster_alpha = 2.
        confidence = 0.999
        p_mut = rng.beta(beta_a, beta_b, (n_pos, n_clust))
        p_ends = np.triu(rng.random((n_pos, n_pos)))
        p_ends /= p_ends.sum()
        p_clust = rng.dirichlet(np.full(n_clust, cluster_alpha))
        # Calculate the parameters with no two mutations too close.
        p_nomut_window = calc_p_nomut_window(
            p_mut, min_mut_gap
        )
        p_noclose_given_ends = calc_p_noclose_given_ends(
            p_mut, p_nomut_window
        )
        p_mut_given_noclose = calc_p_mut_given_span_noclose(
            p_mut, p_ends, p_noclose_given_ends, p_nomut_window
        )
        p_noclose_given_clust = calc_p_noclose_given_clust(
            p_ends, p_noclose_given_ends
        )
        p_ends_given_clust_noclose = calc_p_ends_given_clust_noclose(
            p_ends, p_noclose_given_ends
        )
        p_clust_given_noclose = calc_p_clust_given_noclose(
            p_clust, p_noclose_given_clust
        )
        p_clust_given_ends_noclose = calc_p_clust_given_ends_noclose(
            p_ends_given_clust_noclose, p_clust_given_noclose
        )
        uniq_end5s, uniq_end3s, p_ends_given_noclose = linearize_ends_matrix(
            calc_p_ends_given_noclose(
                p_ends_given_clust_noclose, p_clust_given_noclose)
        )
        # Choose 5'/3' end coordinates.
        ends = rng.choice(p_ends_given_noclose.size,
                          n_reads,
                          p=p_ends_given_noclose,
                          replace=True)
        end5s = uniq_end5s[ends]
        end3s = uniq_end3s[ends]
        # Simulate reads and clusters.
        reads, clusts = _sim_reads(end5s,
                                   end3s,
                                   p_clust_given_ends_noclose,
                                   p_mut,
                                   min_mut_gap)
        clust_counts, span_counts, mut_counts = self.count_reads(reads,
                                                                 clusts,
                                                                 n_clust)
        # Confirm the 5'/3' coordinates match.
        self.assertTupleEqual(reads.shape, (n_reads, n_pos + 2))
        self.assertTupleEqual(clusts.shape, (n_reads,))
        self.assertTrue(np.all(reads[:, -2] == end5s))
        self.assertTrue(np.all(reads[:, -1] == end3s))
        # G-test of the number of reads in each cluster.
        self.assertEqual(clust_counts.sum(), n_reads)
        _, p_value = g_test(clust_counts, n_reads * p_clust_given_noclose)
        self.assertGreaterEqual(p_value, 1. - confidence)
        # Binomial test of the reads covering each position.
        p_span = calc_rectangular_sum(p_ends_given_clust_noclose
                                      * p_clust_given_noclose)
        span_lo, span_up = binom.interval(confidence, n_reads, p_span)
        self.assertTrue(np.all(span_counts >= span_lo))
        self.assertTrue(np.all(span_counts <= span_up))
        # Binomial test of the mutations at each position.
        mut_lo, mut_up = binom.interval(confidence,
                                        span_counts,
                                        p_mut_given_noclose)
        self.assertTrue(np.all(mut_counts >= mut_lo))
        self.assertTrue(np.all(mut_counts <= mut_up))


class TestCalcSemiGAnomaly(ut.TestCase):

    def test_float_equal(self):
        num_obs = 8
        num_exp = 8.
        expect = 0.
        result = calc_semi_g_anomaly(num_obs, np.log(num_exp))
        self.assertIsInstance(result, float)
        self.assertTrue(np.isclose(result, expect))

    def test_float_unequal(self):
        num_obs = 8
        num_exp = 7.
        expect = 8 * np.log(8 / 7.)
        result = calc_semi_g_anomaly(num_obs, np.log(num_exp))
        self.assertIsInstance(result, float)
        self.assertTrue(np.isclose(result, expect))

    def test_array_all_equal(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([2., 3., 1.])
        expect = np.array([0., 0., 0.])
        result = calc_semi_g_anomaly(num_obs, np.log(num_exp))
        self.assertIsInstance(result, np.ndarray)
        self.assertTrue(np.allclose(result, expect))

    def test_array_all_unequal(self):
        num_obs = np.array([2, 3, 1])
        num_exp = np.array([1., 2., 3.])
        expect = np.array([2 * np.log(2 / 1.),
                           3 * np.log(3 / 2.),
                           1 * np.log(1 / 3.)])
        result = calc_semi_g_anomaly(num_obs, np.log(num_exp))
        self.assertIsInstance(result, np.ndarray)
        self.assertTrue(np.allclose(result, expect))


class TestBootstrapJackpotScores(ut.TestCase):
    SIM_DIR = Path(opt_sim_dir.default).absolute()
    REFS = "test_refs"
    REF = "test_ref"
    SAMPLE = "test_sample"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._config = None

    def setUp(self):
        self.SIM_DIR.mkdir()
        self._config = get_config()
        set_config(verbose=0, quiet=1, log_file=None, raise_on_error=True)

    def tearDown(self):
        if self.SIM_DIR.exists():
            shutil.rmtree(self.SIM_DIR)
        set_config(**self._config._asdict())

    def sim_jackpot_quotient(self):
        """ Simulate a dataset and return its jackpotting quotient. """
        n_pos = 60
        n_reads = 50000
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
        pmut_paired = pmut + [("am", 0.005),
                              ("cm", 0.005),
                              ("gm", 0.005),
                              ("tm", 0.005)]
        pmut_unpaired = pmut + [("am", 0.1),
                                ("cm", 0.1),
                                ("gm", 0.005),
                                ("tm", 0.005)]
        run_sim_params(ct_file=(ct_file,),
                       pmut_paired=pmut_paired,
                       pmut_unpaired=pmut_unpaired,
                       insert_fmean=0.5,
                       end3_fmean=0.75,
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
                                     min_finfo_read=1.,
                                     min_ninfo_pos=1,
                                     quick_unbias_thresh=0.)
        # Cluster the data and calculate the jackpotting quotient.
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
                       max_jackpot_quotient=opt_max_jackpot_quotient.default)
        # Delete the simulated files so that this function can run again
        # if necessary.
        shutil.rmtree(self.SIM_DIR)
        return em_run.jackpot_quotient

    @staticmethod
    def calc_confidence_interval(log_jackpot_quotients: list[float],
                                 confidence_level: float):
        n = len(log_jackpot_quotients)
        if n <= 1:
            return np.nan, np.nan
        mean = np.mean(log_jackpot_quotients)
        std_err = np.std(log_jackpot_quotients) / np.sqrt(n)
        t_lo, t_up = studentt.interval(confidence_level, n - 1)
        j_lo = mean + std_err * t_lo
        j_up = mean + std_err * t_up
        return j_lo, j_up

    def test_ideal_jackpot(self):
        """ Test that bootstrapping "perfect" data correctly returns a
        jackpotting quotient that is expected to be 1. """
        confidence_level = 0.99
        confidence_width = 0.02
        log_jackpot_quotients = list()
        while True:
            log_jackpot_quotients.append(np.log(self.sim_jackpot_quotient()))
            jq_lo, jq_up = self.calc_confidence_interval(log_jackpot_quotients,
                                                         confidence_level)
            if not np.isnan(jq_lo) and not np.isnan(jq_up):
                # Verify that the confidence interval contains 0.
                self.assertLessEqual(jq_lo, 0.)
                self.assertGreaterEqual(jq_up, 0.)
                if jq_up - jq_lo < confidence_width:
                    # The confidence interval has converged around 0.
                    break


if __name__ == "__main__":
    ut.main()
