import unittest as ut
from pathlib import Path
from string import ascii_letters

import numpy as np
import pandas as pd
from scipy.stats import binom

from seismicrna.core.batch import END5_COORD, END3_COORD
from seismicrna.core.seq import DNA
from seismicrna.relate.sim import (_list_naturals,
                                   choose_clusters,
                                   simulate_p_mut,
                                   simulate_p_ends,
                                   simulate_p_clust,
                                   simulate_relate_batch,
                                   simulate_relate,
                                   sub_options)

rng = np.random.default_rng()


class TestSubOptions(ut.TestCase):

    def test_bases(self):
        self.assertEqual(sub_options("A"), [2 ** 5, 2 ** 6, 2 ** 7])
        self.assertEqual(sub_options("C"), [2 ** 4, 2 ** 6, 2 ** 7])
        self.assertEqual(sub_options("G"), [2 ** 4, 2 ** 5, 2 ** 7])
        self.assertEqual(sub_options("T"), [2 ** 4, 2 ** 5, 2 ** 6])

    def test_invalid(self):
        for letter in ascii_letters:
            if letter not in "ACGT":
                self.assertRaisesRegex(ValueError,
                                       f"Invalid base: '{letter}'",
                                       sub_options,
                                       letter)


class TestChooseClusters(ut.TestCase):

    def test_1_cluster(self):
        """ With 1 cluster, all reads should be from that cluster. """
        p_clust = simulate_p_clust(1)
        for n_reads in range(5):
            with self.subTest(n_reads=n_reads):
                expect = np.ones(n_reads, dtype=int)
                result = choose_clusters(p_clust, n_reads)
                self.assertIsInstance(result, np.ndarray)
                self.assertIs(result.dtype, expect.dtype)
                self.assertTrue(np.array_equal(result, expect))

    def test_2_clusters(self):
        """ With 2 clusters. """
        confidence = 0.9995
        # Vary the proportion of cluster 1.
        for p1 in np.linspace(0., 1., 5):
            p_clust = pd.Series([p1, 1. - p1], [1, 2])
            for n_reads in [0, 1, 100, 10000]:
                with self.subTest(p1=p1, n_reads=n_reads):
                    result = choose_clusters(p_clust, n_reads)
                    self.assertIsInstance(result, np.ndarray)
                    n1 = np.count_nonzero(result == 1)
                    n2 = np.count_nonzero(result == 2)
                    self.assertEqual(n1 + n2, n_reads)
                    # Verify the number of reads assigned to cluster 1
                    # follows a binomial distribution.
                    lo, up = binom.interval(confidence, n_reads, p1)
                    self.assertGreaterEqual(n1, lo)
                    self.assertLessEqual(n1, up)


class TestSimulateRelateBatch(ut.TestCase):

    def test_nmut(self):
        """ Test number of mutated reads per position. """
        confidence = 0.9999
        n_reads = 20_000
        # Vary the number of positions.
        for npos in range(1, 10):
            # Define the positions and sequence.
            refseq = DNA.random(npos)
            # Simulate proportions of end coordinates.
            p_ends = simulate_p_ends(npos)
            # Estimate the fraction of reads covering each position.
            positions = _list_naturals(npos)
            coverage = pd.Series(
                [p_ends.values[np.logical_and(
                    pos >= p_ends.index.get_level_values(END5_COORD),
                    pos <= p_ends.index.get_level_values(END3_COORD)
                )].sum() for pos in positions],
                index=positions
            )
            # Vary the number of clusters.
            for ncls in range(1, 4):
                # Simulate mutation rates for each position.
                p_mut = simulate_p_mut(refseq, ncls)
                # Choose a cluster for each read.
                p_clust = simulate_p_clust(ncls)
                cluster_choices = choose_clusters(p_clust, n_reads)
                # Simulate the batch.
                batch = simulate_relate_batch(
                    "sample", "ref", 0, n_reads, p_mut, p_ends, cluster_choices
                )
                # Count the reads with mutations at each position.
                for pos, muts in batch.muts.items():
                    base = refseq[pos - 1]
                    # Check the types of substitutions.
                    self.assertEqual(sorted(muts), sub_options(base))
                    # Count all mutated reads at this position.
                    n_mut_reads = sum(map(len, muts.values()))
                    # Compute a confidence interval.
                    p_mut_reads = coverage[pos] * np.vdot(
                        p_mut.loc[pos].values, p_clust.values
                    )
                    lo, up = binom.interval(confidence, n_reads, p_mut_reads)
                    self.assertGreaterEqual(n_mut_reads, lo)
                    self.assertLessEqual(n_mut_reads, up)


class TestSimulateRelate(ut.TestCase):

    def test_simulate(self):
        # Define parameters.
        npos = 300
        nreads = 100_000
        ncls = 2
        batch_size = 16.
        # Simulate the reference sequence.
        refseq = DNA.random(npos)
        # Simulate proportions of end coordinates.
        p_ends = simulate_p_ends(npos)
        # Choose a cluster for each read.
        p_clust = simulate_p_clust(ncls)
        cluster_choices = choose_clusters(p_clust, nreads)
        # Simulate mutation rates for each position.
        p_mut = simulate_p_mut(refseq, ncls)
        # Simulate the relate step.
        report_file = simulate_relate(out_dir=Path.cwd(),
                                      sample="sample",
                                      ref="ref",
                                      refseq=refseq,
                                      batch_size=batch_size,
                                      n_reads=nreads,
                                      p_mut=p_mut,
                                      p_ends=p_ends,
                                      cluster_choices=cluster_choices,
                                      brotli_level=5,
                                      force=False)
        '''
        # Count the reads with mutations at each position.
        for pos, muts in batch.muts.items():
            base = refseq[pos - 1]
            # Check the types of substitutions.
            self.assertEqual(sorted(muts), sub_options(base))
            # Count all mutated reads at this position.
            n_mut_reads = sum(map(len, muts.values()))
            # Compute a confidence interval.
            p_mut_reads = coverage[pos] * np.vdot(
                p_mut.loc[pos].values, p_clust.values
            )
            lo, up = binom.interval(confidence, n_reads, p_mut_reads)
            self.assertGreaterEqual(n_mut_reads, lo)
            self.assertLessEqual(n_mut_reads, up)
        '''


if __name__ == "__main__":
    ut.main()
