import unittest as ut
from sys import byteorder

import pandas as pd

from ..iterread import iter_alignments
from ..relate import relate_line
from ..seqpos import format_seq_pos
from ...align.sim import as_sam
from ...core.cli import opt_min_mapq
from ...core.qual import OK_QUAL
from ...core.relvect import NOCOV
from ...core.seq import DNA


class TestSeqposFormatSeqPos(ut.TestCase):
    """ Test function `seqpos.format_seq_pos`. """

    def test_acgt_index_1_acgt(self):
        """ Test with ACGT, 1-indexed. """
        index = format_seq_pos(DNA("ACGT"), [1, 2, 3, 4], 1)
        expect = pd.Index(["A1", "C2", "G3", "T4"])
        self.assertTrue(index.equals(expect))

    def test_acgt_index_1_cg(self):
        """ Test with ACGT, 1-indexed. """
        index = format_seq_pos(DNA("ACGT"), [2, 3], 1)
        expect = pd.Index(["C2", "G3"])
        self.assertTrue(index.equals(expect))

    def test_acgt_index_58_acgt(self):
        """ Test with ACGT, 58-indexed. """
        index = format_seq_pos(DNA("ACGT"), [58, 59, 60, 61], 58)
        expect = pd.Index(["A58", "C59", "G60", "T61"])
        self.assertTrue(index.equals(expect))

    def test_acgt_index_58_cg(self):
        """ Test with ACGT, 58-indexed. """
        index = format_seq_pos(DNA("ACGT"), [59, 60], 58)
        expect = pd.Index(["C59", "G60"])
        self.assertTrue(index.equals(expect))


class TestRelateRelateLineAmbrel(ut.TestCase):
    """ Test function `relate.relate_line`. """

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
        muts = bytearray(NOCOV.to_bytes(1, byteorder) * len(refseq))
        relate_line(muts,
                    sam_line,
                    ref,
                    refseq,
                    len(refseq),
                    opt_min_mapq.default,
                    OK_QUAL,
                    True)
        return muts

    def iter_cases(self, refseq: DNA, max_ins: int = 2):
        """ Iterate through every test case. """
        for read, qual, cigar, end5, end3, relvec in iter_alignments(refseq,
                                                                     max_ins,
                                                                     max_ins,
                                                                     max_ins):
            result = self.relate("ref",
                                 refseq,
                                 read,
                                 qual,
                                 opt_min_mapq.default,
                                 cigar,
                                 end5)
            with self.subTest(relvec=relvec, result=result):
                self.assertEqual(relvec.tobytes(), result)

    def test_aaaa_0ins(self):
        """ Test all possible reads with 0 insertions from AAAA. """
        self.iter_cases(DNA("AAAA"), 0)

    @ut.skip("Takes a long time to run")
    def test_aaaaaa_0ins(self):
        """ Test all possible reads with 0 insertions from AAAAAA. """
        self.iter_cases(DNA("AAAAAA"), 0)

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