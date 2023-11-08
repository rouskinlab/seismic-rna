import unittest as ut

import numpy as np
import pandas as pd

from ..frame import calc_f_obs_frame, calc_mu_adj_frame, calc_prop_adj_frame
from ..unbias import calc_f_obs_numpy, calc_mu_adj_numpy
from ...seq import DNA, Section, seq_pos_to_index
from ...rand import rng


class TestCalcDataFrame(ut.TestCase):
    """ Test functions `mu.calc_mu_adj_df` and mu.calc_f_obs_df. """

    def test_equals_numpy(self):
        """ Check if the output of `calc_mu_adj_df` equals that of
        `calc_mu_adj_numpy` """
        max_mu = 0.1
        start = 1
        gaps = [0, 3]
        for length in range(1, 10):
            # Generate a random reference sequence.
            refseq = DNA.random(length)
            # Make a section for the sequence.
            section = Section("myref", refseq)
            for n_pos in range(length):
                # Choose a random set of positions, and sort them.
                pos = np.sort(rng.choice(length, n_pos, replace=False))
                # Make an index from those positions.
                index = seq_pos_to_index(refseq, pos + start, start)
                for n_clust in range(5):
                    clusters = pd.Index([f"Cluster-{i}"
                                         for i in range(1, n_clust + 1)])
                    # Generate random mutation rates.
                    mus_obs_values = max_mu * rng.random((n_pos, n_clust))
                    mus_obs_df = pd.DataFrame(mus_obs_values,
                                              index=index,
                                              columns=clusters)
                    # To run calc_mu_adj_numpy, create an array of the
                    # mutation rates where values in missing positions
                    # are set to 0.
                    mus_obs_np = np.zeros((length, n_clust))
                    for i_value, i_numpy in enumerate(pos):
                        mus_obs_np[i_numpy] = mus_obs_values[i_value]
                    for gap in gaps:
                        # Run calc_mu_adj_df.
                        mus_adj_df = calc_mu_adj_frame(mus_obs_df, section, gap)
                        # Run calc_mu_adj_numpy.
                        mus_adj_np = calc_mu_adj_numpy(mus_obs_np, gap)
                        # Compare the results.
                        self.assertIsInstance(mus_adj_df, pd.DataFrame)
                        self.assertTrue(np.allclose(mus_adj_df.values,
                                                    mus_adj_np[pos]))
                        self.assertTrue(index.equals(mus_adj_df.index))
                        self.assertTrue(clusters.equals(mus_adj_df.columns))
                        # Run calc_f_obs_df.
                        f_obs_df = calc_f_obs_frame(mus_adj_df, section, gap)
                        # Run calc_f_obs_numpy.
                        f_obs_np = calc_f_obs_numpy(mus_adj_np, gap)
                        # Compare the results.
                        self.assertIsInstance(f_obs_df, pd.Series)
                        self.assertTrue(np.allclose(f_obs_df.values,
                                                    f_obs_np))
                        self.assertTrue(clusters.equals(f_obs_df.index))


class TestCalcSeries(ut.TestCase):
    """ Test `mu.calc_mu_adj_series` and mu.calc_f_obs_series. """

    def test_equals_numpy(self):
        """ Check if the output of `calc_mu_adj_df` equals that of
        `calc_mu_adj_numpy` """
        max_mu = 0.1
        start = 1
        gaps = [0, 3]
        for length in range(1, 10):
            # Generate a random reference sequence.
            refseq = DNA.random(length)
            # Make a section for the sequence.
            section = Section("myref", refseq)
            for n_pos in range(length):
                # Choose a random set of positions, and sort them.
                pos = np.sort(rng.choice(length, n_pos, replace=False))
                # Make an index from those positions.
                index = seq_pos_to_index(refseq, pos + start, start)
                # Generate random mutation rates.
                mus_obs_values = max_mu * rng.random(n_pos)
                mus_obs_series = pd.Series(mus_obs_values, index=index)
                # To run calc_mu_adj_numpy, create an array of the
                # mutation rates where values in missing positions
                # are set to 0.
                mus_obs_np = np.zeros(length)
                for i_value, i_numpy in enumerate(pos):
                    mus_obs_np[i_numpy] = mus_obs_values[i_value]
                for gap in gaps:
                    # Run calc_mu_adj_series.
                    mus_adj_series = calc_mu_adj_frame(mus_obs_series,
                                                       section,
                                                       gap)
                    # Run calc_mu_adj_numpy.
                    mus_adj_np = calc_mu_adj_numpy(mus_obs_np, gap)
                    # Compare the results.
                    self.assertIsInstance(mus_adj_series, pd.Series)
                    self.assertTrue(np.array_equal(mus_adj_series.values,
                                                   mus_adj_np[pos]))
                    self.assertTrue(index.equals(mus_adj_series.index))
                    # Run calc_f_obs_series.
                    f_obs_series = calc_f_obs_frame(mus_adj_series,
                                                    section,
                                                    gap)
                    # Run calc_f_obs_numpy.
                    f_obs_np = calc_f_obs_numpy(mus_adj_np, gap)
                    # Compare the results.
                    self.assertIsInstance(f_obs_series, float)
                    self.assertIsInstance(f_obs_np, float)
                    self.assertEqual(f_obs_series, f_obs_np)

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
