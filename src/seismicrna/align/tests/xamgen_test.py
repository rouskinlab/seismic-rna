"""

Tests for Alignment XAM Generation Module

========================================================================

"""

import unittest as ut

from ..xamops import (MATCH_BONUS, MISMATCH_PENALTY, N_PENALTY,
                      REF_GAP_PENALTY, READ_GAP_PENALTY)


class TestAlignmentScoreParams(ut.TestCase):
    """ Test that the bonuses and penalties satisfy both inequalities
    described in the docstring of `fqutil.py` and repeated below. """

    @staticmethod
    def parse_scores(scores: str):
        """ Parse Bowtie2 scores (strings of comma-separated integers)
        into tuples of integers. """
        return tuple(map(int, scores.split(',')))

    def parse_all_scores(self):
        """ Parse all scores. """
        # Bonus for matches.
        match, = self.parse_scores(MATCH_BONUS)
        # Penalty for substitutions.
        subst, _ = self.parse_scores(MISMATCH_PENALTY)
        # Penalty for ambiguous (N) bases.
        ambig, = self.parse_scores(N_PENALTY)
        # Penalties for opening and extending gaps in the read.
        gapop, gapex = self.parse_scores(READ_GAP_PENALTY)
        return match, subst, ambig, gapop, gapex

    def test_score_consistency(self):
        """ Test the self-consistency of the scores. """
        match, subst, ambig, gapop, gapex = self.parse_all_scores()
        # Min and max penalties for substitutions, depending on quality.
        _, subst2 = self.parse_scores(MISMATCH_PENALTY)
        # Penalties for opening and extending gaps in the reference.
        gapop2, gapex2 = self.parse_scores(REF_GAP_PENALTY)
        # The mismatch penalty should not depend on the quality score.
        self.assertEqual(subst, subst2)
        # The gap opening penalties should equal for the read and ref.
        self.assertEqual(gapop, gapop2)
        # The gap extension penalties should equal for the read and ref.
        self.assertEqual(gapex, gapex2)

    def test_inequality_min_num_edits(self):
        """
        Test that alignments with two closely separated deletions would
        score better than alignments with two consecutive deletions and
        one substitution.

        Consider this example: Ref = ACGT, Read = AG

        Assume that we want to minimize the number of edits needed to
        convert the reference into the read sequence. The smallest
        number of edits is two, specifically these two deletions (/)
        from the reference: [A/G/] which gets an alignment score of
        (2 * match - 2 * gap_open - 2 * gap_extend).

        But there are two alternative alignments, each with 3 edits:
        [Ag//] and [A//g] (substitutions in lowercase). Each scores
        (match - substitution - gap_open - 2 * gap_extend).

        In order to favor the simpler alignment with two edits,
        (2 * match - 2 * gap_open - 2 * gap_extend) must be greater than
        (match - substitution - gap_open - 2 * gap_extend). Simplifying:
        (substitution > gap_open - match).
        """
        match, subst, ambig, gapop, gapex = self.parse_all_scores()
        self.assertGreater(subst, gapop - match)

    def test_inequality_subst_vs_indel(self):
        """
        Test that alignments with two adjacent substitutions would score
        better than alignments with an insertion and a deletion flanking
        a match.

        Consider this example: Ref = ATAT, Read = ACTT

        The alignment ActT has two mutations (two substitutions) and an
        alignment score of (2 * match - 2 * substitution). The alignment
        A{C}T/T (where {C} means a C was inserted into the and / means
        an A deleted from the read). This alignment would get a score of
        (3 * match - 2 * gap_open - 2 * gap_extend).

        Because mutational profiling causes many more substitutions than
        indels, the alignment with two substitutions should be preferred
        over the alignment with one insertion and one deletion (although
        both have two mutations).
        
        Thus, (2 * match - 2 * substitution) must be greater than
        (3 * match - 2 * gap_open - 2 * gap_extend), which simplifies to
        (2 * gap_open + 2 * gap_extend > match + 2 * substitution).
        """
        match, subst, ambig, gapop, gapex = self.parse_all_scores()
        self.assertGreater(2 * gapop + 2 * gapex, match + 2 * subst)


if __name__ == "__main__":
    ut.main()

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
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
