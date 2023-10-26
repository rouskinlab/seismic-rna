import unittest as ut

from ..relate import find_rels_line
from ...aux.iterread import iter_alignments
from ....align.sim import as_sam
from ....core.arg import opt_min_mapq
from ....core.ngs import OK_QUAL
from ....core.seq import DNA


class TestRelateRelateLineAmbrel(ut.TestCase):
    """ Test function `relate.relate_line`. """

    MAPQ = opt_min_mapq.default

    @staticmethod
    def relate(ref: str,
               refseq: DNA,
               read: DNA,
               qual: str,
               mapq: int,
               cigar: str,
               end5: int):
        """ Generate a SAM line from the given information, and use it
        to compute a relation vector. """
        sam_line = as_sam("read",
                          99,
                          ref,
                          end5,
                          mapq,
                          cigar,
                          "=",
                          1,
                          len(read),
                          read,
                          qual)
        return find_rels_line(sam_line,
                              "",
                              ref,
                              refseq,
                              opt_min_mapq.default,
                              OK_QUAL,
                              True)

    def iter_cases(self, refseq: DNA, max_ins: int = 2):
        """ Iterate through every test case. """
        for read, qual, cigar, end5, end3, rels in iter_alignments(refseq,
                                                                   max_ins,
                                                                   max_ins,
                                                                   max_ins):
            name, end5_, mid5, mid3, end3_, rels_ = self.relate("ref",
                                                                refseq,
                                                                read,
                                                                qual,
                                                                self.MAPQ,
                                                                cigar,
                                                                end5)
            with self.subTest(refseq=refseq,
                              read=read,
                              qual=qual,
                              end5=end5,
                              cigar=cigar,
                              rels=rels):
                self.assertEqual(end5_, end5)
                self.assertEqual(mid5, end5)
                self.assertEqual(mid3, end3)
                self.assertEqual(end3_, end3)
                self.assertEqual(rels_, rels)

    def test_aaaa_0ins(self):
        """ Test all possible reads with 0 insertions from AAAA. """
        self.iter_cases(DNA("AAAA"), 0)

    @ut.skip("Takes a long time to run; burdensome while debugging others")
    def test_aaaaaa_0ins(self):
        """ Test all possible reads with 0 insertions from AAAAAA. """
        self.iter_cases(DNA("AAAAAA"), 0)

    @ut.skip("Takes a long time to run; burdensome while debugging others")
    def test_aacc_1ins(self):
        """ Test all possible reads with ≤ 1 insertion from AACC. """
        self.iter_cases(DNA("AACC"), 1)

    def test_acgt_1ins(self):
        """ Test all possible reads with ≤ 1 insertion from ACGT. """
        self.iter_cases(DNA("ACGT"), 1)

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
