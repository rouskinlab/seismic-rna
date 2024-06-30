import unittest as ut
from itertools import chain, product

from seismicrna.core.arg import opt_min_mapq
from seismicrna.core.ngs import LO_QUAL, OK_QUAL, MAX_FLAG, SAM_DELIM
from seismicrna.core.rel import DELET, IRREC, MATCH, NOCOV, SUB_G
from seismicrna.core.seq import DNA
from seismicrna.relate.py.cigar import CIG_ALIGN, CIG_DELET, CIG_SCLIP
from seismicrna.relate.py.encode import encode_relate
from seismicrna.relate.py.relate import _find_rels_read, _merge_mates, SamRead
from seismicrna.relate.aux.iterread import iter_alignments


def as_sam(name: str,
           flag: int,
           ref: str,
           end5: int,
           mapq: int,
           cigar: str,
           rnext: str,
           pnext: int,
           tlen: int,
           read: DNA,
           qual: str):
    """
    Return a line in SAM format from the given fields.

    Parameters
    ----------
    name: str
        Name of the read.
    flag: int
        SAM flag. Must be in [0, MAX_FLAG].
    ref: str
        Name of the reference.
    end5: int
        Most 5' position to which the read mapped (1-indexed).
    mapq: int
        Mapping quality score.
    cigar: str
        CIGAR string. Not checked for compatibility with the read.
    rnext: str
        Name of the mate's reference (if paired-end).
    pnext: int
        Most 5' position of the mate (if paired-end).
    tlen: int
        Length of the template.
    read: DNA
        Base calls in the read. Must be equal in length to `read`.
    qual: str
        Phred quality score string of the base calls. Must be equal in
        length to `read`.

    Returns
    -------
    str
        A line in SAM format containing the given fields.
    """
    if not name:
        raise ValueError("Read name is empty")
    if not 0 <= flag <= MAX_FLAG:
        raise ValueError(f"Invalid SAM flag: {flag}")
    if not ref:
        raise ValueError("Reference name is empty")
    if not end5 >= 1:
        raise ValueError(f"Invalid 5' mapping position: {end5}")
    if not cigar:
        raise ValueError("CIGAR string is empty")
    if not rnext:
        raise ValueError("Next reference name is empty")
    if not pnext >= 0:
        raise ValueError(f"Invalid next 5' mapping position: {pnext}")
    if not len(read) == len(qual):
        raise ValueError(
            f"Lengths of read ({len(read)}) and qual ({len(qual)}) disagree")
    return SAM_DELIM.join(map(str, (name, flag, ref, end5, mapq, cigar, rnext,
                                    pnext, tlen, read, f"{qual}\n")))


class TestFindRelsLine(ut.TestCase):
    """ Test function `relate.relate_line`. """

    @staticmethod
    def relate(ref: str,
               refseq: DNA,
               read: DNA,
               qual: str,
               cigar: str,
               end5: int,
               ambindel: bool,
               clip_end5: int,
               clip_end3: int):
        """ Generate a SAM line from the given information, and use it
        to compute a relation vector. """
        sam_read = SamRead(as_sam("read",
                                  99,
                                  ref,
                                  end5,
                                  opt_min_mapq.default,
                                  cigar,
                                  "=",
                                  1,
                                  len(read),
                                  read,
                                  qual))
        return _find_rels_read(sam_read,
                               refseq,
                               OK_QUAL,
                               ambindel,
                               clip_end5,
                               clip_end3)

    def iter_cases(self, refseq: DNA, max_ins: int = 2):
        """ Iterate through every test case. """
        for read, qual, cigar, end5, end3, rels in iter_alignments(refseq,
                                                                   max_ins,
                                                                   max_ins,
                                                                   max_ins):
            with self.subTest(refseq=refseq,
                              read=read,
                              qual=qual,
                              end5=end5,
                              cigar=cigar,
                              rels=rels):
                result = self.relate("ref",
                                     refseq,
                                     read,
                                     qual,
                                     cigar,
                                     end5,
                                     ambindel=True,
                                     clip_end5=0,
                                     clip_end3=0)
                expect = (end5, end3, rels)
                self.assertEqual(result, expect)

    def test_aaaa_0ins(self):
        """ Test all possible reads with 0 insertions from AAAA. """
        self.iter_cases(DNA("AAAA"), 0)

    def test_aaaaaa_0ins(self):
        """ Test all possible reads with 0 insertions from AAAAAA. """
        self.iter_cases(DNA("AAAAAA"), 0)

    def test_aacc_1ins(self):
        """ Test all possible reads with ≤ 1 insertion from AACC. """
        self.iter_cases(DNA("AACC"), 1)

    def test_acgt_1ins(self):
        """ Test all possible reads with ≤ 1 insertion from ACGT. """
        self.iter_cases(DNA("ACGT"), 1)

    def test_all_matches(self):
        for reflen in range(1, 10):
            refseq = DNA.random(reflen)
            for readlen in range(1, reflen + 1):
                for end5 in range(1, reflen - readlen + 1):
                    end3 = end5 + readlen - 1
                    read = refseq[end5 - 1: end3]
                    qual = OK_QUAL * readlen
                    cigar = f"{readlen}{CIG_ALIGN}"
                    for clip5 in range(10):
                        for clip3 in range(10):
                            with self.subTest(reflen=reflen,
                                              readlen=readlen,
                                              end5=end5,
                                              clip5=clip5,
                                              clip3=clip3):
                                result = self.relate("ref",
                                                     refseq,
                                                     read,
                                                     qual,
                                                     cigar,
                                                     end5,
                                                     True,
                                                     clip5,
                                                     clip3)
                                expect = (min(end5 + clip5, reflen + 1),
                                          max(end3 - clip3, 0),
                                          dict())
                                self.assertEqual(result, expect)

    def test_soft_clips(self):
        reflen = 10
        refseq = DNA.random(reflen)
        for readlen in range(1, reflen + 1):
            for soft5 in range(readlen + 1):
                cigar_s5 = f"{soft5}{CIG_SCLIP}" if soft5 else ""
                for soft3 in range(readlen - soft5 + 1):
                    cigar_s3 = f"{soft3}{CIG_SCLIP}" if soft3 else ""
                    soft = soft5 + soft3
                    for end5 in range(soft5 + 1, reflen - readlen + 2):
                        matches = readlen - soft
                        end3 = end5 + matches - 1
                        cigar_m = f"{matches}{CIG_ALIGN}" if matches else ""
                        read = sum([DNA("N") * soft5,
                                    refseq[end5 - 1: end3],
                                    DNA("N") * soft3],
                                   DNA(""))
                        qual = OK_QUAL * readlen
                        cigar = "".join([cigar_s5, cigar_m, cigar_s3])
                        for clip5 in range(3):
                            for clip3 in range(3):
                                with self.subTest(reflen=reflen,
                                                  readlen=readlen,
                                                  soft5=soft5,
                                                  soft3=soft3,
                                                  end5=end5,
                                                  clip5=clip5,
                                                  clip3=clip3):
                                    result = self.relate("ref",
                                                         refseq,
                                                         read,
                                                         qual,
                                                         cigar,
                                                         end5,
                                                         True,
                                                         clip5,
                                                         clip3)
                                    expect = (min(end5 + clip5, reflen + 1),
                                              max(end3 - clip3, 0),
                                              dict())
                                    self.assertEqual(result, expect)

    def test_ambig_delet_low_qual(self):
        """ Test ambiguous deletions with all low-quality positions. """
        reflen = 10
        refseq = DNA.random(reflen)
        for readlen in range(2, reflen):
            for soft5 in range(readlen - 1):
                cigar_s5 = f"{soft5}{CIG_SCLIP}" if soft5 else ""
                for soft3 in range(readlen - soft5 - 1):
                    cigar_s3 = f"{soft3}{CIG_SCLIP}" if soft3 else ""
                    soft = soft5 + soft3
                    for end5 in range(soft5 + 1, reflen - readlen + 1):
                        end3 = end5 + readlen - soft
                        for delpos in range(end5 + 1, end3):
                            cigar_md = "".join([f"{delpos - end5}{CIG_ALIGN}",
                                                f"{1}{CIG_DELET}",
                                                f"{end3 - delpos}{CIG_ALIGN}"])
                            read = sum([DNA("N") * soft5,
                                        refseq[end5 - 1: delpos - 1],
                                        refseq[delpos: end3],
                                        DNA("N") * soft3],
                                       DNA(""))
                            qual = LO_QUAL * readlen
                            cigar = "".join([cigar_s5, cigar_md, cigar_s3])
                            for clip5 in range(3):
                                for clip3 in range(3):
                                    with self.subTest(reflen=reflen,
                                                      readlen=readlen,
                                                      soft5=soft5,
                                                      soft3=soft3,
                                                      end5=end5,
                                                      clip5=clip5,
                                                      clip3=clip3):
                                        result = self.relate("ref",
                                                             refseq,
                                                             read,
                                                             qual,
                                                             cigar,
                                                             end5,
                                                             True,
                                                             clip5,
                                                             clip3)
                                        read5 = min(end5 + clip5, reflen + 1)
                                        read3 = max(end3 - clip3, 0)
                                        positions = list(range(read5,
                                                               read3 + 1))
                                        rels = {pos: 0 for pos in positions}
                                        for pos in positions:
                                            if end5 < pos < end3:
                                                rels[pos] |= DELET
                                            if not end5 + 1 == pos == end3 - 1:
                                                rels[pos] |= encode_relate(
                                                    refseq[pos - 1],
                                                    "N",
                                                    LO_QUAL,
                                                    OK_QUAL
                                                )
                                        expect = read5, read3, rels
                                        self.assertEqual(result, expect)


class TestMergeMates(ut.TestCase):

    def test_empty(self):
        result = _merge_mates(1, 10, {}, 1, 10, {}, True)
        expect = ([1, 1], [10, 10]), {}
        self.assertEqual(result, expect)

    def test_read1(self):
        end51 = 1
        end31 = 20
        end52 = 11
        end32 = 30
        for pos in range(end51, end31 + 1):
            for rel in range(MATCH + 1, NOCOV):
                result = _merge_mates(end51, end31, {pos: rel},
                                      end52, end32, {},
                                      True)
                if end52 <= pos <= end32:
                    # The relationship can be compensated by read 2.
                    if rel & MATCH:
                        # The match in read 2 compensated.
                        expect = ([1, 11], [20, 30]), {}
                    else:
                        # The match in read 2 is irreconcilable.
                        expect = ([1, 11], [20, 30]), {pos: IRREC}
                else:
                    # Read 2 cannot compensate.
                    expect = ([1, 11], [20, 30]), {pos: rel}
                self.assertEqual(result, expect)

    def test_read2(self):
        end51 = 1
        end31 = 20
        end52 = 11
        end32 = 30
        for pos in range(end52, end32 + 1):
            for rel in range(MATCH + 1, NOCOV):
                result = _merge_mates(end51, end31, {},
                                      end52, end32, {pos: rel},
                                      True)
                if end51 <= pos <= end31:
                    # The relationship can be compensated by read 1.
                    if rel & MATCH:
                        # The match in read 1 compensated.
                        expect = ([1, 11], [20, 30]), {}
                    else:
                        # The match in read 1 is irreconcilable.
                        expect = ([1, 11], [20, 30]), {pos: IRREC}
                else:
                    # Read 1 cannot compensate.
                    expect = ([1, 11], [20, 30]), {pos: rel}
                self.assertEqual(result, expect)

    def test_both_reads(self):
        end51 = 1
        end31 = 2
        end52 = 2
        end32 = 3
        for pos1 in range(end51, end31 + 1):
            for rel1 in range(MATCH + 1, NOCOV):
                rels1 = {pos1: rel1}
                for pos2 in range(end52, end32 + 1):
                    for rel2 in range(MATCH + 1, NOCOV):
                        rels2 = {pos2: rel2}
                        with self.subTest(pos1=pos1, rel1=rel1,
                                          pos2=pos2, rel2=rel2):
                            result = _merge_mates(end51, end31, rels1,
                                                  end52, end32, rels2,
                                                  True)
                            if pos1 == pos2:
                                merged = rel1 & rel2
                                if merged == MATCH:
                                    expect = ([1, 2], [2, 3]), {}
                                else:
                                    expect = ([1, 2], [2, 3]), {pos1: merged}
                            else:
                                expect = ([1, 2], [2, 3]), {}
                                merged1 = (rel1 & MATCH
                                           if end52 <= pos1 <= end32
                                           else rel1)
                                if merged1 != MATCH:
                                    expect[1][pos1] = merged1
                                merged2 = (rel2 & MATCH
                                           if end51 <= pos2 <= end31
                                           else rel2)
                                if merged2 != MATCH:
                                    expect[1][pos2] = merged2
                            self.assertEqual(result, expect)

    def test_both_blank(self):
        end51 = 1
        end31 = 2
        end52 = 2
        end32 = 3
        for pos1 in range(end51, end31 + 1):
            rels1 = {pos1: NOCOV}
            for pos2 in range(end52, end32 + 1):
                rels2 = {pos2: NOCOV}
                with self.subTest(pos1=pos1, pos2=pos2):
                    if end52 <= pos1 <= end32:
                        error = pos2
                    else:
                        error = pos1
                    self.assertRaisesRegex(
                        ValueError,
                        f"Cannot merge non-covered position {error}",
                        _merge_mates,
                        end51, end31, rels1,
                        end52, end32, rels2,
                        True
                    )

    def test_overhangs(self):
        for end5f, end5r, read_length in product(range(5), repeat=3):
            end3f = end5f + read_length
            end3r = end5r + read_length
            relsf = {pos: SUB_G for pos in range(end5f, end3f + 1)}
            relsr = {pos: SUB_G for pos in range(end5r, end3r + 1)}
            for overhangs in [True, False]:
                result = _merge_mates(end5f, end3f, relsf,
                                      end5r, end3r, relsr,
                                      overhangs)
                if overhangs:
                    ends = [end5f, end5r], [end3f, end3r]
                else:
                    ends = ([end5f, max(end5f, end5r)],
                            [min(end3f, end3r), end3r])
                rels = {pos: SUB_G
                        for pos in chain(range(ends[0][0], ends[1][0] + 1),
                                         range(ends[0][1], ends[1][1] + 1))}
                expect = ends, rels
                with self.subTest(overhangs=overhangs,
                                  end5f=end5f,
                                  end3f=end3f,
                                  end5r=end5r,
                                  end3r=end3r):
                    self.assertEqual(result, expect)


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
