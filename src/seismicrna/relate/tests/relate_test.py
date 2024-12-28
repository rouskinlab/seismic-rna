import unittest as ut
from itertools import chain, product
from typing import Any

from seismicrna.relate.cx.relate import (RelateError as RelateErrorCx,
                                         calc_rels_lines as calc_rels_lines_cx)

from seismicrna.core.ngs import LO_QUAL, OK_QUAL, MAX_FLAG, SAM_DELIM
from seismicrna.core.rel import (DELET,
                                 IRREC,
                                 MATCH,
                                 NOCOV,
                                 SUB_A,
                                 SUB_C,
                                 SUB_G,
                                 SUB_T,
                                 ANY_N)
from seismicrna.core.seq import DNA
from seismicrna.relate.aux.iterread import iter_alignments
from seismicrna.relate.py.cigar import CIG_ALIGN, CIG_DELET, CIG_SCLIP
from seismicrna.relate.py.encode import encode_relate
from seismicrna.relate.py.relate import (RelateError as RelateErrorPy,
                                         calc_rels_lines as calc_rels_lines_py,
                                         merge_mates)


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
           qual: str,
           validate: bool = True):
    """ Return a line in SAM format from the given fields.

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
    validate: bool
        Check that the fields are valid before assembling the line.

    Returns
    -------
    str
        A line in SAM format containing the given fields.
    """
    if validate:
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
                f"Lengths of read ({len(read)}) and qual ({len(qual)}) differ"
            )
    return SAM_DELIM.join(map(str, (name, flag, ref, end5, mapq, cigar, rnext,
                                    pnext, tlen, read, f"{qual}\n")))


class TestCalcRelsLinesSingle(ut.TestCase):

    def relate(self,
               ref: str,
               refseq: DNA,
               read: DNA,
               qual: str,
               cigar: str,
               end5: int,
               ambindel: bool,
               insert3: bool,
               clip_end5: int,
               clip_end3: int,
               paired: bool = False):
        """ Generate a SAM line from the given information, and use it
        to compute the relationships. """
        line1 = as_sam("read",
                       int(paired),
                       ref,
                       end5,
                       ord(OK_QUAL),
                       cigar,
                       "=",
                       1,
                       len(read),
                       read,
                       qual)
        line2 = line1 if paired else ""
        result_cx = calc_rels_lines_cx(line1,
                                       line2,
                                       ref,
                                       str(refseq),
                                       0,
                                       ord(OK_QUAL),
                                       insert3,
                                       ambindel,
                                       False,
                                       clip_end5,
                                       clip_end3)
        # Test calc_rels_lines_py second to ensure calc_rels_lines_cx
        # didn't corrupt any memory.
        result_py = calc_rels_lines_py(line1,
                                       line2,
                                       ref,
                                       str(refseq),
                                       0,
                                       ord(OK_QUAL),
                                       insert3,
                                       ambindel,
                                       False,
                                       clip_end5,
                                       clip_end3)
        self.assertEqual(result_cx, result_py)
        return result_cx

    def relate_error(self,
                     error_msg: str,
                     error_msg_py: str = "",
                     ref: str = "ref",
                     refseq: DNA = DNA("ACGT"),
                     read: DNA = DNA("ACGT"),
                     qual: str = "FFFF",
                     cigar: str = "4M",
                     end5: Any = 1,
                     sam_ref: str = "",
                     mapq: Any = None,
                     flag: Any = None,
                     ambindel: bool = True,
                     insert3: bool = True,
                     clip_end5: int = 0,
                     clip_end3: int = 0,
                     paired: bool = False):
        line1 = as_sam("read",
                       flag if flag is not None else int(paired),
                       ref,
                       end5,
                       mapq if mapq is not None else ord(OK_QUAL),
                       cigar,
                       "=",
                       1,
                       len(read),
                       read,
                       qual,
                       validate=False)
        line2 = line1 if paired else ""
        self.assertRaisesRegex(RelateErrorCx,
                               error_msg,
                               calc_rels_lines_cx,
                               line1,
                               line2,
                               sam_ref if sam_ref else ref,
                               str(refseq),
                               ord(OK_QUAL),
                               ord(OK_QUAL),
                               insert3,
                               ambindel,
                               False,
                               clip_end5,
                               clip_end3)
        # Test calc_rels_lines_py second to ensure calc_rels_lines_cx
        # didn't corrupt any memory.
        self.assertRaisesRegex(RelateErrorPy,
                               error_msg_py if error_msg_py else error_msg,
                               calc_rels_lines_py,
                               line1,
                               line2,
                               sam_ref if sam_ref else ref,
                               str(refseq),
                               ord(OK_QUAL),
                               ord(OK_QUAL),
                               insert3,
                               ambindel,
                               False,
                               clip_end5,
                               clip_end3)

    def relate_truncated(self,
                         num_fields: int,
                         error_msg: str,
                         ref: str = "ref",
                         refseq: DNA = DNA("ACGT"),
                         read: DNA = DNA("ACGT"),
                         qual: str = "FFFF",
                         cigar: str = "4M",
                         end5: Any = 1,
                         ambindel: bool = True,
                         insert3: bool = True,
                         clip_end5: int = 0,
                         clip_end3: int = 0,
                         paired: bool = False):
        """ Test errors caused by lines that are too short. """
        line1 = as_sam("read",
                       int(paired),
                       ref,
                       end5,
                       ord(OK_QUAL),
                       cigar,
                       "=",
                       1,
                       len(read),
                       read,
                       qual)
        line1 = SAM_DELIM.join(line1.split(SAM_DELIM)[:num_fields])
        line2 = line1 if paired else ""
        self.assertRaisesRegex(RelateErrorCx,
                               error_msg,
                               calc_rels_lines_cx,
                               line1,
                               line2,
                               ref,
                               str(refseq),
                               0,
                               ord(OK_QUAL),
                               insert3,
                               ambindel,
                               False,
                               clip_end5,
                               clip_end3)
        # Test calc_rels_lines_py second to ensure calc_rels_lines_cx
        # didn't corrupt any memory.
        self.assertRaisesRegex(RelateErrorPy,
                               error_msg,
                               calc_rels_lines_py,
                               line1,
                               line2,
                               ref,
                               str(refseq),
                               0,
                               ord(OK_QUAL),
                               insert3,
                               ambindel,
                               False,
                               clip_end5,
                               clip_end3)

    def iter_cases_insert3(self,
                           refseq: DNA,
                           max_ins: int,
                           insert3: bool,
                           paired: bool):
        """ Iterate through every test case. """
        for read, qual, cigar, end5, end3, rels in iter_alignments(
                refseq,
                insert3=insert3,
                max_ins=max_ins,
                max_ins_len=max_ins,
                max_ins_bases=max_ins
        ):
            with self.subTest(refseq=refseq,
                              insert3=insert3,
                              read=read,
                              qual=qual,
                              cigar=cigar,
                              end5=end5,
                              end3=end3,
                              rels=rels,
                              paired=paired):
                result = self.relate("ref",
                                     refseq,
                                     read,
                                     qual,
                                     cigar,
                                     end5,
                                     ambindel=True,
                                     insert3=insert3,
                                     clip_end5=0,
                                     clip_end3=0,
                                     paired=paired)
                if paired:
                    expect = ([end5, end5], [end3, end3]), rels
                else:
                    expect = ([end5], [end3]), rels
                self.assertEqual(result, expect)

    def iter_cases(self, refseq: DNA, max_ins: int, paired: bool = False):
        self.iter_cases_insert3(refseq, max_ins, False, paired)
        if max_ins > 0:
            self.iter_cases_insert3(refseq, max_ins, True, paired)

    def test_4nt_2ins(self):
        self.iter_cases(DNA("AGCT"), 2)

    def test_4nt_2ins_paired(self):
        self.iter_cases(DNA("CTAG"), 2, paired=True)

    def test_5nt_2ins(self):
        self.iter_cases(DNA("CAAAT"), 2)

    def test_6nt_2ins(self):
        self.iter_cases(DNA("GTATAC"), 2)

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
                                end5_expect = min(end5 + clip5, reflen + 1)
                                end3_expect = max(end3 - clip3, 0)
                                result = self.relate("ref",
                                                     refseq,
                                                     read,
                                                     qual,
                                                     cigar,
                                                     end5,
                                                     ambindel=True,
                                                     insert3=True,
                                                     clip_end5=clip5,
                                                     clip_end3=clip3)
                                expect = (([end5_expect], [end3_expect]),
                                          dict())
                                self.assertEqual(result, expect)

    def test_soft_clips(self):
        ref = "ref"
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
                        if matches > 0:
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
                                        end5_expect = min(end5 + clip5,
                                                          reflen + 1)
                                        end3_expect = max(end3 - clip3, 0)
                                        result = self.relate(ref,
                                                             refseq,
                                                             read,
                                                             qual,
                                                             cigar,
                                                             end5,
                                                             ambindel=True,
                                                             insert3=True,
                                                             clip_end5=clip5,
                                                             clip_end3=clip3)
                                        expect = (([end5_expect],
                                                   [end3_expect]),
                                                  dict())
                                        self.assertEqual(result, expect)
                        else:
                            cigar = f"{soft}{CIG_SCLIP}"
                            with self.subTest(reflen=reflen,
                                              readlen=readlen,
                                              soft5=soft5,
                                              soft3=soft3,
                                              end5=end5):
                                self.relate_error(
                                    ("CIGAR operations consumed 0 bases "
                                     "in the reference"),
                                    ref=ref,
                                    refseq=refseq,
                                    read=read,
                                    qual=qual,
                                    cigar=cigar,
                                    end5=end5,
                                )

    def test_ambig_delet_low_qual(self):
        """ Test ambiguous deletions with all low-quality positions. """
        reflen = 10
        refseq = DNA.random(reflen)
        for readlen in range(2, reflen):
            for soft5 in range(readlen - 1):
                cigar_s5 = f"{soft5}{CIG_SCLIP}" if soft5 else ""
                for soft3 in range(readlen - soft5 - 1):
                    cigar_s3 = f"{soft3}{CIG_SCLIP}" if soft3 else ""
                    for end5 in range(soft5 + 1, reflen - readlen + 1):
                        # The read has exactly 1 deletion, so readlen is
                        # already 1 less than the number of bases the
                        # read takes up of the reference sequence, so do
                        # not substract 1 from end3.
                        end3 = end5 + readlen - (soft5 + soft3)
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
                                                             ambindel=True,
                                                             insert3=True,
                                                             clip_end5=clip5,
                                                             clip_end3=clip3)
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
                                        expect = ([read5], [read3]), rels
                                        self.assertEqual(result, expect)

    def test_n_ref(self):
        """ Reference contains a non-ACGT base. """
        for n in "ACGTN":
            result = self.relate(ref="ref",
                                 refseq=DNA("N"),
                                 read=DNA(n),
                                 qual="F",
                                 cigar="1M",
                                 end5=1,
                                 ambindel=True,
                                 insert3=True,
                                 clip_end5=0,
                                 clip_end3=0)
            expect = (([1], [1]), {1: ANY_N})
            self.assertEqual(result, expect)

    def test_n_read(self):
        """ Read contains a non-ACGT base and reference is ACGT. """
        for n, sub in {"A": SUB_A, "C": SUB_C, "G": SUB_G, "T": SUB_T}.items():
            result = self.relate(ref="ref",
                                 refseq=DNA(n),
                                 read=DNA("N"),
                                 qual="F",
                                 cigar="1M",
                                 end5=1,
                                 ambindel=True,
                                 insert3=True,
                                 clip_end5=0,
                                 clip_end3=0)
            expect = (([1], [1]), {1: ANY_N - sub})
            self.assertEqual(result, expect)
            # Test it is equivalent to the read having low quality.
            self.assertEqual(result, self.relate(ref="ref",
                                                 refseq=DNA(n),
                                                 read=DNA(n),
                                                 qual="!",
                                                 cigar="1M",
                                                 end5=1,
                                                 ambindel=True,
                                                 insert3=True,
                                                 clip_end5=0,
                                                 clip_end3=0))

    def test_example_1(self):
        """ Soft clips in CIGAR plus clip_end5 and clip_end3.

        Seq  GGTATAG
        Qul  FFFFFFF
        CGR  SS====S
        Ref CAATATATC
        Pos 123456789
        """
        result = self.relate(ref="ref",
                             refseq=DNA("CAATATATC"),
                             read=DNA("GGTATAG"),
                             qual="FFFFFFF",
                             cigar="2S4=1S",
                             end5=4,
                             ambindel=True,
                             insert3=True,
                             clip_end5=1,
                             clip_end3=1)
        expect = (([5], [6]), {})
        self.assertEqual(result, expect)

    def test_example_2(self):
        """ Deletions cannot move out of soft-clipped regions.

        Seq  TATA--TAT
        Qul  FFFF--FFF
        CGR  SS==DD==S
        Ref ATATATATATA
        Pos 123456789ab
        """
        # No soft clips.
        result = self.relate(ref="ref",
                             refseq=DNA("ATATATATATA"),
                             read=DNA("TATATAT"),
                             qual="FFFFFFF",
                             cigar="4=2D3=",
                             end5=2,
                             ambindel=True,
                             insert3=True,
                             clip_end5=0,
                             clip_end3=0)
        expect = (([2], [10]), {3: 3, 4: 3, 5: 3, 6: 3, 7: 3, 8: 3, 9: 3})
        self.assertEqual(result, expect)
        # Soft clips.
        result = self.relate(ref="ref",
                             refseq=DNA("ATATATATATA"),
                             read=DNA("TATATAT"),
                             qual="FFFFFFF",
                             cigar="2S2=2D2=1S",
                             end5=4,
                             ambindel=True,
                             insert3=True,
                             clip_end5=0,
                             clip_end3=0)
        expect = (([4], [9]), {5: 3, 6: 3, 7: 3, 8: 3})
        self.assertEqual(result, expect)
        # Soft clips and clip_end5/clip_end3.
        result = self.relate(ref="ref",
                             refseq=DNA("ATATATATATA"),
                             read=DNA("TATATAT"),
                             qual="FFFFFFF",
                             cigar="2S2=2D2=1S",
                             end5=4,
                             ambindel=True,
                             insert3=True,
                             clip_end5=1,
                             clip_end3=1)
        expect = (([5], [8]), {5: 3, 6: 3, 7: 3, 8: 3})
        self.assertEqual(result, expect)

    def test_example_3(self):
        """ Insertions cannot move out of soft-clipped regions.

        Seq  TATATATAT
        Qul  FFFFFFFFF
        CGR  SS==II==S
        Ref ATATA--TATA
        Pos 12345--6789
        """
        # No soft clips.
        result = self.relate(ref="ref",
                             refseq=DNA("ATATATATA"),
                             read=DNA("TATATATAT"),
                             qual="FFFFFFFFF",
                             cigar="4=2I3=",
                             end5=2,
                             ambindel=True,
                             insert3=True,
                             clip_end5=0,
                             clip_end3=0)
        expect = (([2], [8]), {3: 9, 4: 9, 5: 9, 6: 9, 7: 9, 8: 9})
        self.assertEqual(result, expect)
        # Soft clips.
        result = self.relate(ref="ref",
                             refseq=DNA("ATATATATA"),
                             read=DNA("TATATATAT"),
                             qual="FFFFFFFFF",
                             cigar="2S2=2I2=1S",
                             end5=4,
                             ambindel=True,
                             insert3=True,
                             clip_end5=0,
                             clip_end3=0)
        expect = (([4], [7]), {5: 9, 6: 9, 7: 9})
        self.assertEqual(result, expect)
        # Soft clips and clip_end5/clip_end3.
        result = self.relate(ref="ref",
                             refseq=DNA("ATATATATA"),
                             read=DNA("TATATATAT"),
                             qual="FFFFFFFFF",
                             cigar="2S2=2I2=1S",
                             end5=4,
                             ambindel=True,
                             insert3=True,
                             clip_end5=1,
                             clip_end3=1)
        expect = (([5], [6]), {5: 9, 6: 9})
        self.assertEqual(result, expect)

    def test_error_name_missing(self):
        self.relate_truncated(0, "Failed to parse read name")

    def test_error_flag_missing(self):
        self.relate_truncated(1, "Failed to parse SAM flag")

    def test_error_flag_parse(self):
        self.relate_error("Failed to parse SAM flag",
                          flag="1X")

    def test_error_flag_large(self):
        self.relate_error("SAM flag is too large",
                          flag=MAX_FLAG + 1)

    def test_error_ref_missing(self):
        self.relate_truncated(2, "Failed to parse reference name")

    def test_error_pos_missing(self):
        self.relate_truncated(3, "Failed to parse mapping position")

    def test_error_pos_parse(self):
        self.relate_error("Failed to parse mapping position",
                          end5="2Y")

    def test_error_pos_zero(self):
        self.relate_error("Mapping position is 0",
                          end5=0)

    def test_error_pos_large(self):
        self.relate_error("Mapping position exceeds length of reference",
                          end5=5)

    def test_error_mapq_missing(self):
        self.relate_truncated(4, "Failed to parse mapping quality")

    def test_error_mapq(self):
        self.relate_error("Failed to parse mapping quality",
                          mapq="3Z")

    def test_error_cigar_missing(self):
        self.relate_truncated(5, "Failed to parse CIGAR string")

    def test_error_cigar_empty(self):
        self.relate_error("Failed to parse CIGAR string",
                          cigar="")

    def test_error_read_missing(self):
        self.relate_truncated(9, "Failed to parse read sequence")

    def test_error_qual_missing(self):
        self.relate_truncated(10, "Failed to parse read quality")

    def test_error_read_qual_diff(self):
        self.relate_error("Read sequence and quality strings differ in length",
                          qual="FFFFF")

    def test_error_ref_mismatch(self):
        self.relate_error("Reference name does not match name of SAM file",
                          sam_ref="other")

    def test_error_mapq_insufficient(self):
        self.relate_error("Mapping quality is insufficient",
                          mapq=ord(OK_QUAL) - 1)

    def test_error_line_paired_flag_unpaired(self):
        self.relate_error("Lines indicate read should be paired-end, "
                          "but it is marked as single-end",
                          paired=True,
                          flag=0)

    def test_error_line_unpaired_flag_paired(self):
        self.relate_error("Lines indicate read should be single-end, "
                          "but it is marked as paired-end",
                          paired=False,
                          flag=1)

    def test_error_line_improper_flag_proper(self):
        self.relate_error("Lines indicate read should be improperly paired, "
                          "but it is marked as properly paired",
                          paired=False,
                          flag=2)
        self.relate_error("Lines indicate read should be improperly paired, "
                          "but it is marked as properly paired",
                          paired=True,
                          flag=3)

    def test_error_cigar_parse(self):
        self.relate_error("Invalid CIGAR operation",
                          error_msg_py="Invalid CIGAR string",
                          cigar="M4M")
        self.relate_error("Invalid CIGAR operation",
                          error_msg_py="CIGAR operation has length 0",
                          cigar="0M")
        self.relate_error("Unsupported CIGAR operation",
                          error_msg_py="Invalid CIGAR string",
                          cigar="4A")

    def test_error_cigar_consecutive(self):
        for op in "M=XDIS":
            self.relate_error("Identical consecutive CIGAR operations",
                              cigar=f"1{op}2{op}")

    def test_error_cigar_adj_ins_del(self):
        self.relate_error("Adjacent insertion and deletion",
                          cigar="1D1I3M")
        self.relate_error("Adjacent insertion and deletion",
                          cigar="1I1D3M")
        self.relate_error("Adjacent insertion and deletion",
                          cigar="1M1D1I2M")
        self.relate_error("Adjacent insertion and deletion",
                          cigar="2M1I1D1M")
        self.relate_error("Adjacent insertion and deletion",
                          cigar="3M1D1I")
        self.relate_error("Adjacent insertion and deletion",
                          cigar="3M1I1D")

    def test_error_cigar_op_ref_zero(self):
        self.relate_error("CIGAR operations consumed 0 bases in the reference",
                          cigar="4S")

    def test_error_cigar_op_ref_long(self):
        for cigar in ["5M", "5=", "5X", "5D", "2M1D2M"]:
            self.relate_error("CIGAR operations extended out of the reference",
                              cigar=cigar)

    def test_error_cigar_op_read_diff(self):
        for cigar in ["3M", "5M", "5=", "5X", "2M1I2M", "2M1D1M"]:
            self.relate_error("CIGAR operations consumed a number of read "
                              "bases different from the read length",
                              refseq=DNA("ACGTA"),
                              cigar=cigar)

    def test_error_cigar_del_first_rel(self):
        self.relate_error("A deletion was the first relationship",
                          refseq=DNA("TACGT"),
                          end5=1,
                          cigar="1D4M")
        self.relate_error("A deletion was the first relationship",
                          refseq=DNA("GTACGT"),
                          end5=2,
                          cigar="1D4M")
        self.relate_error("A deletion was the first relationship",
                          refseq=DNA("TACGT"),
                          end5=1,
                          cigar="1D1S3M")
        self.relate_error("A deletion was the first relationship",
                          refseq=DNA("TACGT"),
                          end5=1,
                          cigar="1S1D3M")
        self.relate_error("A deletion was the first relationship",
                          refseq=DNA("GTACG"),
                          end5=2,
                          cigar="2S1D2M")

    def test_error_cigar_del_last_rel(self):
        self.relate_error("A deletion was the last relationship",
                          refseq=DNA("ACGTAC"),
                          end5=1,
                          cigar="4M1D")
        self.relate_error("A deletion was the last relationship",
                          refseq=DNA("TACGTA"),
                          end5=2,
                          cigar="4M1D")
        self.relate_error("A deletion was the last relationship",
                          refseq=DNA("ACGTA"),
                          end5=1,
                          cigar="3M1D1S")
        self.relate_error("A deletion was the last relationship",
                          refseq=DNA("TACGT"),
                          end5=2,
                          cigar="2M1D2S")

    def test_error_cigar_ins_first_rel(self):
        self.relate_error("An insertion was the first relationship",
                          refseq=DNA("TAC"),
                          end5=1,
                          cigar="1I3M")
        self.relate_error("An insertion was the first relationship",
                          refseq=DNA("GTAC"),
                          end5=2,
                          cigar="1I3M")
        self.relate_error("An insertion was the first relationship",
                          refseq=DNA("GT"),
                          end5=1,
                          cigar="1I1S2M")
        self.relate_error("An insertion was the first relationship",
                          refseq=DNA("GT"),
                          end5=1,
                          cigar="1S1I2M")
        self.relate_error("An insertion was the first relationship",
                          refseq=DNA("TA"),
                          end5=2,
                          cigar="2S1I1M")

    def test_error_cigar_ins_last_rel(self):
        self.relate_error("An insertion was the last relationship",
                          refseq=DNA("ACG"),
                          end5=1,
                          cigar="3M1I")
        self.relate_error("An insertion was the last relationship",
                          refseq=DNA("TACG"),
                          end5=2,
                          cigar="3M1I")
        self.relate_error("An insertion was the last relationship",
                          refseq=DNA("AC"),
                          end5=1,
                          cigar="2M1I1S")
        self.relate_error("An insertion was the last relationship",
                          refseq=DNA("TA"),
                          end5=2,
                          cigar="1M1I2S")

    def test_error_cigar_soft_clips(self):
        self.relate_error("A soft clip occurred in the middle",
                          refseq=DNA("AGT"),
                          end5=1,
                          cigar="1M1S2M")
        self.relate_error("A soft clip occurred in the middle",
                          refseq=DNA("TAGT"),
                          end5=2,
                          cigar="1M1S2M")
        self.relate_error("A soft clip occurred in the middle",
                          refseq=DNA("TACT"),
                          end5=2,
                          cigar="2M1S1M")
        self.relate_error("A soft clip occurred in the middle",
                          refseq=DNA("ACG"),
                          end5=1,
                          cigar="2M1S1I")
        self.relate_error("A soft clip occurred in the middle",
                          refseq=DNA("ACGT"),
                          end5=1,
                          cigar="3M1S1D")


class TestCalcRelsLinesPaired(ut.TestCase):

    def relate(self,
               ref: str,
               refseq: DNA,
               read1: DNA,
               qual1: str,
               cigar1: str,
               end51: int,
               read2: DNA,
               qual2: str,
               cigar2: str,
               end52: int,
               ambindel: bool = True,
               insert3: bool = True,
               clip_end5: int = 0,
               clip_end3: int = 0,
               read1rev: bool = False):
        """ Generate a SAM line from the given information, and use it
        to compute the relationships. """
        line1 = as_sam("read",
                       (83 if read1rev else 99),
                       ref,
                       end51,
                       ord(OK_QUAL),
                       cigar1,
                       "=",
                       1,
                       len(read1),
                       read1,
                       qual1)
        line2 = as_sam("read",
                       (163 if read1rev else 147),
                       ref,
                       end52,
                       ord(OK_QUAL),
                       cigar2,
                       "=",
                       1,
                       len(read2),
                       read2,
                       qual2)
        result_cx = calc_rels_lines_cx(line1,
                                       line2,
                                       ref,
                                       str(refseq),
                                       0,
                                       ord(OK_QUAL),
                                       insert3,
                                       ambindel,
                                       True,
                                       clip_end5,
                                       clip_end3)
        # Test calc_rels_lines_py second to ensure calc_rels_lines_cx
        # didn't corrupt any memory.
        result_py = calc_rels_lines_py(line1,
                                       line2,
                                       ref,
                                       str(refseq),
                                       0,
                                       ord(OK_QUAL),
                                       insert3,
                                       ambindel,
                                       True,
                                       clip_end5,
                                       clip_end3)
        self.assertEqual(result_py, result_cx)
        return result_cx

    def evaluate(self,
                 expect_ends1: tuple[int, int],
                 expect_ends2: tuple[int, int],
                 expect_rels: dict[int, int],
                 ref: str,
                 refseq: DNA,
                 read1: DNA,
                 qual1: str,
                 cigar1: str,
                 end51: int,
                 read2: DNA,
                 qual2: str,
                 cigar2: str,
                 end52: int):
        for swap_reads in [False, True]:
            if swap_reads:
                expect_ends1, expect_ends2 = expect_ends2, expect_ends1
                read1, read2 = read2, read1
                qual1, qual2 = qual2, qual1
                cigar1, cigar2 = cigar2, cigar1
                end51, end52 = end52, end51
            exp_end51, exp_end31 = expect_ends1
            exp_end52, exp_end32 = expect_ends2
            for read1_rev in [False, True]:
                result = self.relate(ref=ref,
                                     refseq=refseq,
                                     read1=read1,
                                     qual1=qual1,
                                     cigar1=cigar1,
                                     end51=end51,
                                     read2=read2,
                                     qual2=qual2,
                                     cigar2=cigar2,
                                     end52=end52,
                                     read1rev=read1_rev)
                if read1_rev:
                    expect = (([exp_end52, exp_end51],
                               [exp_end32, exp_end31]),
                              expect_rels)
                else:
                    expect = (([exp_end51, exp_end52],
                               [exp_end31, exp_end32]),
                              expect_rels)
                self.assertEqual(result, expect)

    def test_gap(self):
        """ Reads are separated by a gap.

        R1  AGTG
        R2       TCGT
        Ref AGTCAACGT
        Pos 123456789
        """
        self.evaluate((1, 4),
                      (6, 9),
                      {4: SUB_G, 6: SUB_T},
                      ref="ref",
                      refseq=DNA("AGTCAACGT"),
                      read1=DNA("AGTG"),
                      qual1="FFFF",
                      cigar1="4M",
                      end51=1,
                      read2=DNA("TCGT"),
                      qual2="FFFF",
                      cigar2="4M",
                      end52=6)

    def test_abut(self):
        """ Reads abut.

        R1   GTGC
        R2       TCGT
        Ref AGTCAACGT
        Pos 123456789
        """
        self.evaluate((2, 5),
                      (6, 9),
                      {4: SUB_G, 5: SUB_C, 6: SUB_T},
                      ref="ref",
                      refseq=DNA("AGTCAACGT"),
                      read1=DNA("GTGC"),
                      qual1="FFFF",
                      cigar1="4M",
                      end51=2,
                      read2=DNA("TCGT"),
                      qual2="FFFF",
                      cigar2="4M",
                      end52=6)

    def test_staggered(self):
        """ Reads overlap in a staggered manner.

        R1  aGTGgtA
        R2    AGATGc
        Ref AGTCAACGT
        Pos 123456789
        """
        self.evaluate((1, 7),
                      (3, 8),
                      {1: MATCH + SUB_C + SUB_G + SUB_T,
                       3: IRREC,
                       4: SUB_G,
                       6: SUB_T,
                       7: IRREC,
                       8: MATCH + SUB_A + SUB_C + SUB_T},
                      ref="ref",
                      refseq=DNA("AGTCAACGT"),
                      read1=DNA("aGTGgtA"),
                      qual1="!FFF!!F",
                      cigar1="7M",
                      end51=1,
                      read2=DNA("AGATGc"),
                      qual2="FFFFF!",
                      cigar2="6M",
                      end52=3)

    def test_contain_flush5(self):
        """ One read contains the other, with 5' ends flush.

        R1   ATcAggG
        R2   gTcaT
        Ref AGTCAACGT
        Pos 123456789
        """
        self.evaluate((2, 8),
                      (2, 6),
                      {2: SUB_A,
                       4: MATCH + SUB_A + SUB_G + SUB_T,
                       6: SUB_T,
                       7: MATCH + SUB_A + SUB_G + SUB_T},
                      ref="ref",
                      refseq=DNA("AGTCAACGT"),
                      read1=DNA("ATcAggG"),
                      qual1="FF!F!!F",
                      cigar1="7M",
                      end51=2,
                      read2=DNA("gTcaT"),
                      qual2="!F!!F",
                      cigar2="5M",
                      end52=2)

    def test_contain_flush53(self):
        """ Both reads start and end at the same positions.

        R1    tc-ACG
        R2    gGATCG
        Ref AGTCAACGT
        Pos 123456789
        """
        self.evaluate((3, 8),
                      (3, 8),
                      {3: MATCH + SUB_A + SUB_C + SUB_G,
                       4: SUB_G,
                       6: IRREC},
                      ref="ref",
                      refseq=DNA("AGTCAACGT"),
                      read1=DNA("tcACG"),
                      qual1="!!FFF",
                      cigar1="2M1D3M",
                      end51=3,
                      read2=DNA("gGATCG"),
                      qual2="!FFFFF",
                      cigar2="6M",
                      end52=3)

    def test_contain_flush3(self):
        """ One read contains the other, with 3' ends flush.

        R1     TAtCa
        R2   CggAtCc
        Ref AGTCAACGT
        Pos 123456789
        """
        self.evaluate((4, 8),
                      (2, 8),
                      {2: SUB_C,
                       3: MATCH + SUB_A + SUB_C + SUB_G,
                       4: SUB_T,
                       6: MATCH + SUB_C + SUB_G + SUB_T,
                       8: MATCH + SUB_A + SUB_C + SUB_T},
                      ref="ref",
                      refseq=DNA("AGTCAACGT"),
                      read1=DNA("TAtCa"),
                      qual1="FF!F!",
                      cigar1="5M",
                      end51=4,
                      read2=DNA("CggAtCc"),
                      qual2="F!!F!F!",
                      cigar2="7M",
                      end52=2)

    def test_contain(self):
        """ One read contains the other, with neither end flush.

        R1    TgcAT
        R2   gtaAACG
        Ref AGTCAACGT
        Pos 123456789
        """
        self.evaluate((3, 7),
                      (2, 8),
                      {2: MATCH + SUB_A + SUB_C + SUB_T,
                       4: MATCH + SUB_A + SUB_G + SUB_T,
                       7: IRREC},
                      ref="ref",
                      refseq=DNA("AGTCAACGT"),
                      read1=DNA("TgcAT"),
                      qual1="F!!FF",
                      cigar1="5M",
                      end51=3,
                      read2=DNA("gtaAACG"),
                      qual2="!!!FFFF",
                      cigar2="7M",
                      end52=2)

    def relate_error(self,
                     error_msg: str,
                     ref: str = "ref",
                     refseq: DNA = DNA("ACGT"),
                     name1: str = "read",
                     ref1: str = "ref",
                     cigar1: str = "4M",
                     flag1: int = 83,
                     end51: int = 1,
                     read1: DNA = DNA("ACGT"),
                     qual1: str = "FFFF",
                     name2: str = "read",
                     ref2: str = "ref",
                     flag2: int = 163,
                     end52: int = 1,
                     cigar2: str = "4M",
                     read2: DNA = DNA("ACGT"),
                     qual2: str = "FFFF",
                     ambindel: bool = True,
                     insert3: bool = True,
                     clip_end5: int = 0,
                     clip_end3: int = 0):
        line1 = as_sam(name1,
                       flag1,
                       ref1,
                       end51,
                       ord(OK_QUAL),
                       cigar1,
                       ref2,
                       1,
                       len(read1),
                       read1,
                       qual1)
        line2 = as_sam(name2,
                       flag2,
                       ref2,
                       end52,
                       ord(OK_QUAL),
                       cigar2,
                       ref1,
                       1,
                       len(read2),
                       read2,
                       qual2)
        self.assertRaisesRegex(RelateErrorCx,
                               error_msg,
                               calc_rels_lines_cx,
                               line1,
                               line2,
                               ref,
                               str(refseq),
                               0,
                               ord(OK_QUAL),
                               insert3,
                               ambindel,
                               True,
                               clip_end5,
                               clip_end3)
        # Test calc_rels_lines_py second to ensure calc_rels_lines_cx
        # didn't corrupt any memory.
        self.assertRaisesRegex(RelateErrorPy,
                               error_msg,
                               calc_rels_lines_py,
                               line1,
                               line2,
                               ref,
                               str(refseq),
                               0,
                               ord(OK_QUAL),
                               insert3,
                               ambindel,
                               True,
                               clip_end5,
                               clip_end3)

    def test_diff_names(self):
        self.relate_error("Mates 1 and 2 have different names",
                          name1="mate1",
                          name2="mate2")

    def test_unpaired(self):
        self.relate_error("Lines indicate read should be paired-end, "
                          "but it is marked as single-end",
                          flag1=83 ^ 1)
        self.relate_error("Lines indicate read should be paired-end, "
                          "but it is marked as single-end",
                          flag2=163 ^ 1)

    def test_improper(self):
        self.relate_error("Lines indicate read should be properly paired, "
                          "but it is marked as improperly paired",
                          flag1=83 ^ 2)
        self.relate_error("Lines indicate read should be properly paired, "
                          "but it is marked as improperly paired",
                          flag1=163 ^ 2)

    def test_read_marks(self):
        self.relate_error("Mate 1 is not marked as READ1",
                          flag1=83 ^ 64)
        self.relate_error("Mate 1 is not marked as READ1",
                          flag1=83 ^ 128)
        self.relate_error("Mate 2 is not marked as READ2",
                          flag2=163 ^ 128)
        self.relate_error("Mate 2 is not marked as READ2",
                          flag2=163 ^ 64)

    def test_read_orientation(self):
        self.relate_error("Mates 1 and 2 aligned in the same orientation",
                          flag1=83 ^ 16,
                          flag2=163)
        self.relate_error("Mates 1 and 2 aligned in the same orientation",
                          flag1=83,
                          flag2=163 ^ 16)


class TestMergeMates(ut.TestCase):

    def test_empty(self):
        result = merge_mates(1, 10, {}, 1, 10, {}, True)
        expect = ([1, 1], [10, 10]), {}
        self.assertEqual(result, expect)

    def test_read1(self):
        end51 = 1
        end31 = 20
        end52 = 11
        end32 = 30
        for pos in range(end51, end31 + 1):
            for rel in range(MATCH + 1, NOCOV):
                result = merge_mates(end51, end31, {pos: rel},
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
                result = merge_mates(end51, end31, {},
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
                            result = merge_mates(end51, end31, rels1,
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
                        RelateErrorPy,
                        f"Cannot merge non-covered position {error}",
                        merge_mates,
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
                result = merge_mates(end5f, end3f, relsf,
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
