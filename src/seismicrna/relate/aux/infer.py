from typing import Sequence

from .cigarop import CigarOp
from ..py.cigar import CIG_ALIGN, CIG_MATCH, CIG_SUBST, CIG_DELET, CIG_INSRT
from ..py.encode import SUBS_DECODINGS
from ...core.ngs import LO_QUAL, HI_QUAL
from ...core.rel import (MATCH,
                         DELET,
                         INS_5,
                         INS_3,
                         INSRT,
                         SUB_A,
                         SUB_C,
                         SUB_G,
                         SUB_T,
                         ANY_N,
                         IRREC,
                         NOCOV)
from ...core.seq import BASEN, DNA


def infer_read(refseq: DNA,
               end5: int,
               end3: int,
               muts: dict[int, int],
               hi_qual: str = HI_QUAL,
               lo_qual: str = LO_QUAL,
               ins_len: int | Sequence[int] = 1):
    """ Infer the sequence and quality string of a read from a reference
    sequence and relationships.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    end5: int
        5' end of the read with respect to the reference.
    end3: int
        3' end of the read with respect to the reference.
    muts: dict[int, int]
        Mutations in the read, keyed by their positions.
    hi_qual: str = MAX_QUAL
        Character to put in the quality string at every position that is
        high-quality.
    lo_qual: str = MIN_QUAL
        Character to put in the quality string at every position that is
        low-quality.
    ins_len: int | Sequence[int] = 1
        Number of bases to insert into the read and quality strings upon
        finding an insertion in the relationships. If an integer, then
        insert that number of bases for every insertion. If a sequence
        of integers, then the ith insertion gets a number of bases equal
        to the ith element of ins_len.

    Returns
    -------
    """
    # Validate the relationships and reference sequence.
    if end5 < 1:
        raise ValueError(f"end5 must be ≥ 1, but got {end5}")
    if end3 > len(refseq):
        raise ValueError(f"end3 must be ≤ {len(refseq)}, but got {end3}")
    if muts:
        if min(muts) < end5:
            raise ValueError(f"All positions must be ≥ {end5}, "
                             f"but got {[pos for pos in muts if pos < end5]}")
        if max(muts) > end3:
            raise ValueError(f"All positions must be ≤ {end3}, "
                             f"but got {[pos for pos in muts if pos > end3]}")
    # Validate the quality codes.
    if len(hi_qual) != 1:
        raise ValueError(f"hi_qual must be 1 character, but got {len(hi_qual)}")
    if len(lo_qual) != 1:
        raise ValueError(f"lo_qual must be 1 character, but got {len(hi_qual)}")
    if hi_qual < lo_qual:
        raise ValueError(f"The high-quality score ({hi_qual}) is less than "
                         f"than the low-quality score ({lo_qual})")
    # Build the read sequence, quality scores, and CIGAR string one base
    # at a time.
    read: list[str] = list()
    qual: list[str] = list()
    cigars: list[CigarOp] = list()
    ins_count = 0
    need_to_add_ins3 = False

    def add_to_cigar(op: str):
        """ Add one base of the relationships to the CIGAR string. """
        if cigars and cigars[-1].op == op:
            # The current operation matches that of the last CigarOp:
            # just increment its length.
            cigars[-1].lengthen()
        else:
            # Either no CigarOp objects exist or the current operation
            # does not match that of the last one: start a new CigarOp.
            cigars.append(CigarOp(op))

    for pos in range(end5, end3 + 1):
        ref_base = refseq[pos - 1]
        rel = muts.get(pos, MATCH)
        if rel == NOCOV:
            raise ValueError(f"Position {pos} in {end5}-{end3} is not covered")
        if need_to_add_ins3 or rel & INSRT:
            # Specially handle insertions because they may overlap any
            # other relation except for a deletion.
            if rel & DELET:
                # Being both a deletion and an insertion is forbidden.
                raise ValueError(f"Position {pos} in {muts} is del and ins")
            if need_to_add_ins3 or rel & INS_3:
                if pos <= end5:
                    # Insertions cannot occur before the beginning of
                    # the read.
                    raise ValueError(f"Position {pos} in {end5}-{end3} cannot "
                                     f"be 3' of an insertion")
                # The current position is 5' of an insertion, so the
                # inserted base must be added before the base at the
                # current position is added. Insert a number of bases
                # given by ins_len.
                if isinstance(ins_len, int):
                    # Add the same number of bases for every insertion.
                    n_ins = ins_len
                else:
                    # Add a number of bases for this specific insertion.
                    n_ins = ins_len[ins_count]
                read.append(BASEN * n_ins)
                qual.append(hi_qual * n_ins)
                for _ in range(n_ins):
                    add_to_cigar(CIG_INSRT)
                ins_count += 1
                # Reset need_to_add_ins3 to False once the insertion has
                # been added.
                need_to_add_ins3 = False
            if rel & INS_5:
                if pos >= end3:
                    raise ValueError(f"Position {pos} in {end5}-{end3} cannot "
                                     f"be 5' of an insertion")
                # The current position is 5' of an insertion, so the
                # inserted base must be added after the base at the
                # current position is added. Defer the responsibility of
                # adding the inserted base to the base that lies 3' of
                # the insertion; indicate that the next base must be 3'
                # of an insertion by setting need_to_add_ins3 to True.
                need_to_add_ins3 = True
            # Switch off the insertion flag(s) in the relation so that
            # the base at this position can be added below.
            rel = rel & ~INSRT
            if not rel:
                # If no other relationship was under the insertion, then
                # this position was otherwise a match.
                rel = MATCH
        if rel == MATCH:
            # Match: Add the reference base to the read.
            read.append(ref_base)
            qual.append(hi_qual)
            add_to_cigar(CIG_MATCH)
        elif rel == DELET:
            # Deletion from the read.
            if not end5 < pos < end3:
                raise ValueError(
                    f"Deletion cannot be at position {pos} in {end5}-{end3}"
                )
            add_to_cigar(CIG_DELET)
        elif rel ^ ANY_N in (SUB_A, SUB_C, SUB_G, SUB_T, IRREC):
            # Ambiguous substitution: Add any nucleotide as low quality.
            read.append(BASEN)
            qual.append(lo_qual)
            add_to_cigar(CIG_ALIGN)
        else:
            # Unambiguous substitution: Add the substituted nucleotide
            # as high quality.
            try:
                read_base = SUBS_DECODINGS[rel]
            except KeyError:
                raise ValueError(f"Invalid relation {rel} in {muts}[{pos}]")
            if ref_base == read_base:
                raise ValueError(
                    f"Cannot substitute {ref_base} to itself in {muts}[{pos}]"
                )
            read.append(read_base)
            qual.append(hi_qual)
            add_to_cigar(CIG_SUBST)
    if len(read) != len(qual):
        raise ValueError(
            f"Lengths of read ({len(read)}) and qual ({len(qual)}) differed"
        )
    if not read:
        raise ValueError("Read contained no bases")
    if need_to_add_ins3:
        raise ValueError("Failed to mark the base 3' of an insertion")
    # Assemble and return the read, quality, and CIGAR strings.
    return DNA("".join(read)), "".join(qual), "".join(map(str, cigars))
