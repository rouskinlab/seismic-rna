"""

Tests for Alignment FASTQ Utilities Module

========================================================================

©2023, the Rouskin Lab.

This file is part of SEISMIC-RNA.

SEISMIC-RNA is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

SEISMIC-RNA is distributed in the hope that it will be useful, but WITH
NO WARRANTY; not even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along
with SEISMIC-RNA. If not, see https://www.gnu.org/licenses/.

========================================================================

"""

import unittest as ut

from .fqutil import (MATCH_BONUS, MISMATCH_PENALTY, N_PENALTY,
                     REF_GAP_PENALTY, READ_GAP_PENALTY)


class TestAlignmentScoreParams(ut.TestCase):
    """ Test that the bonuses and penalties satisfy both inequalities
    described in the docstring of `fqutil.py`. """

    @staticmethod
    def parse_scores(scores: str):
        """ Parse Bowtie2 scores (strings of comma-separated integers)
        into tuples of integers. """
        return tuple(map(int, scores.split(',')))

    @classmethod
    def parse_all_scores(cls):
        """ Parse all scores. """
        # Bonus for matches.
        mat, = cls.parse_scores(MATCH_BONUS)
        # Min and max penalties for substitutions, depending on quality.
        sun, sux = cls.parse_scores(MISMATCH_PENALTY)
        # Penalty for ambiguous (N) bases.
        amb, = cls.parse_scores(N_PENALTY)
        # Penalties for opening and extending gaps in the read.
        opd, exd = cls.parse_scores(READ_GAP_PENALTY)
        # Penalties for opening and extending gaps in the reference.
        opf, exf = cls.parse_scores(REF_GAP_PENALTY)
        return mat, sun, sux, amb, opd, exd, opf, exf

    def test_score_consistency(self):
        """ Test the self-consistency of the scores. """
        mat, sun, sux, amb, opd, exd, opf, exf = self.parse_all_scores()
        # The mismatch penalty should not depend on the quality score.
        self.assertEqual(sun, sux)
        # The gap opening penalties should equal for the read and ref.
        self.assertEqual(opd, opf)
        # The gap extension penalties should equal for the read and ref.
        self.assertEqual(exd, exf)

    def test_inequality_better_indel(self):
        """ Pass """
