"""

Relate Auxiliary Core Module

========================================================================

------------------------------------------------------------------------

"""

from collections import defaultdict
from itertools import (product, combinations,
                       combinations_with_replacement as cwr)
from typing import Sequence

import numpy as np
import pandas as pd

from .rel import (BASE_DECODINGS, BASE_ENCODINGS, MIN_INDEL_GAP, NP_TYPE,
                  CIG_ALIGN, CIG_MATCH, CIG_SUBST,
                  CIG_DELET, CIG_INSRT, CIG_SCLIP,
                  NOCOV, MATCH, DELET, INS_5, INS_3, INS_8,
                  ANY_N, SUB_A, SUB_C, SUB_G, SUB_T, IRREC,
                  MIN_QUAL, MAX_QUAL, blank_relvec, encode_relate, parse_cigar)
from .sect import Section
from .seq import BASEA, BASEC, BASEG, BASET, BASEN, DNA, expand_degenerate_seq
from .sim import rng


class CigarOp(object):
    """ Represent one operation in a CIGAR string. """

    def __init__(self, op: str):
        if op not in (CIG_ALIGN, CIG_MATCH, CIG_SUBST,
                      CIG_DELET, CIG_INSRT, CIG_SCLIP):
            raise ValueError(f"Invalid CIGAR operation: '{op}'")
        self._op = op
        self._len = 1

    @property
    def op(self):
        """ CIGAR operation as a character. """
        return self._op

    def lengthen(self):
        """ Lengthen the operation by 1 base call. """
        self._len += 1

    def __str__(self):
        """ Text that goes into the CIGAR string. """
        return f"{self._len}{self._op}"


def count_cigar_muts(cigar_string: str):
    """ Return the total number of mutations in a CIGAR string. """
    mutation_types = CIG_SUBST, CIG_DELET, CIG_INSRT
    return sum(olen for op, olen in parse_cigar(cigar_string)
               if op in mutation_types)


def find_cigar_op_pos(cigar_string: str, find_op: str):
    """ Yield the position in the read of every base with a particular
    type of operation specified by a CIGAR string. """
    consume_read = CIG_ALIGN, CIG_MATCH, CIG_SUBST, CIG_INSRT
    if find_op in consume_read:
        # Only these operations correspond to positions in the read.
        pos = 1
        # Check each CIGAR operation, starting at position 1.
        for op, olen in parse_cigar(cigar_string):
            if op == find_op:
                # If the operation matches the code, then yield the position
                # of every base consumed by that operation.
                yield from range(pos, pos := (pos + olen))
            elif op in consume_read:
                # Advance the position by the length of the operation.
                pos += olen


def validate_relvec(relvec: np.ndarray):
    """ Confirm that a relation vector is valid, then return it. """
    if not isinstance(relvec, np.ndarray):
        raise TypeError(f"Expected {np.ndarray}, but got {type(relvec)}")
    if relvec.dtype.type is not NP_TYPE:
        raise TypeError(f"Expected an array of type {NP_TYPE}, but got "
                        f"type {relvec.dtype.type}")
    if relvec.ndim != 1:
        raise ValueError("Expected an array with 1 dimension, but got "
                         f"{relvec.ndim} dimensions")
    called = np.flatnonzero(relvec != NOCOV)
    if called.size == 0:
        raise ValueError("Relation vector is blank")
    end5 = int(called[0]) + 1
    end3 = int(called[-1]) + 1
    if (nexp := end3 - end5 + 1) != called.size:
        raise ValueError(f"Expected {nexp} base calls "
                         f"({np.arange(end5, end3 + 1)}), but got "
                         f"{called.size} ({called})")
    return end5, end3


def random_relvecs(refseq: DNA, reads: pd.Index | np.ndarray | list | int,
                   ploq: pd.Series, pmut: pd.Series,
                   ptri: float = 0.50, ptrc: float = 0.25, ptrn: float = 0.25):
    """
    Return random relation vectors.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    reads: pd.Index | np.ndarray | list | int
        Names of the reads or number of reads (if int).
    ploq: Series
        Probability that each position is low-quality.
    pmut: Series
        Mutation rate of each position.
    ptri: float = 0.25
        Probability that a mutation is a transition.
    ptrc: float = 0.25
        Probability that a mutation is a complementary transversion.
    ptrn: float = 0.25
        Probability that a mutation is a non-complementary transversion.

    Returns
    -------
    """
    # Validate ploq and pmut.
    if not ploq.index.equals(pmut.index):
        raise ValueError(f"Indexes differ between ploq ({ploq.index}) "
                         f"and pmut ({pmut.index})")
    if np.any(ploq < 0.) or np.any(ploq >= 1.):
        raise ValueError(f"All ploq must be in [0, 1), but got {ploq}")
    if np.any(pmut < 0.) or np.any(pmut >= 1.):
        raise ValueError(f"All pmut must be in [0, 1), but got {pmut}")
    # Validate ptri, ptrc, and ptrn.
    if not 0. <= ptri <= 1.:
        raise ValueError(f"ptri must be in [0, 1], but got {ptri}")
    if not 0. <= ptrc <= 1.:
        raise ValueError(f"ptrc must be in [0, 1], but got {ptrc}")
    if not 0. <= ptrn <= 1.:
        raise ValueError(f"ptrn must be in [0, 1], but got {ptrn}")
    # The sum of the probabilities cannot exceed 1.
    if (qdel := ptri + ptrc + ptrn) > 1.:
        raise ValueError(f"(ptri + ptrc + ptrn) must be ≤ 1, but got {qdel}")
    # Compute the probability that a mutation is a deletion.
    pdel = 1. - qdel
    # Initialize a DataFrame of blank relation vectors.
    relvecs = blank_relvec(refseq, reads)
    n_reads, n_pos = relvecs.shape
    # Initially set all relationships within the section to matches.
    relvecs.loc[:, ploq.index] = MATCH
    # Simulate whether each position is low-quality.
    is_loq = rng.random((n_reads, ploq.size)) < ploq
    # Find the indexes of the low-quality positions.
    loq_rows, loq_cols = np.nonzero(is_loq)
    n_loqs = loq_rows.size
    if n_loqs != loq_cols.size:
        raise ValueError(f"Counts differ for low-quality rows ({n_loqs}) "
                         f"and columns ({loq_cols.size})")
    # Mark every low-quality base.
    for base, code in BASE_ENCODINGS.items():
        pass
    # Simulate whether each high-quality position is mutated.
    is_mut = np.logical_and(rng.random((n_reads, pmut.size)) < pmut,
                            np.logical_not(is_loq))
    # Find the indexes of the mutated positions.
    mut_rows, mut_cols = np.nonzero(is_mut)
    n_muts = mut_rows.size
    if n_muts != mut_cols.size:
        raise ValueError(f"Counts differ for mutated rows ({n_muts}) "
                         f"and columns ({mut_cols.size})")
    # Simulate the type of each mutation.
    pmut_types = [pdel, ptrc, ptri, ptrn]
    is_mut_types = rng.choice(np.arange(len(pmut_types)), n_muts, p=pmut_types)
    #


def iter_relvecs_q53(refseq: DNA, low_qual: Sequence[int] = (),
                     end5: int | None = None, end3: int | None = None,
                     max_ins: int = 2):
    """
    For a given reference sequence, yield every possible unambiguous
    relation vector between positions end5 and end3 that follows the
    given low-quality positions and has at most two insertions.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    low_qual: Sequence[int]
        List of positions in the read that are low-quality.
    end5: int | None = None
        5' end of the read; 1-indexed with respect to `refseq`.
    end3: int | None = None
        3' end of the read; 1-indexed with respect to `refseq`.
    max_ins: int = 2
        Maximum number of insertions in the read. Must be in [0, 2].
    """
    low_qual = set(low_qual)
    # Determine the section of the reference sequence that is occupied
    # by the read.
    sect = Section("", refseq, end5=end5, end3=end3)
    if low_qual - set(sect.range_int):
        raise ValueError(f"Invalid positions in low_qual: "
                         f"{sorted(low_qual - set(sect.range))}")
    if max_ins not in range(3):
        raise ValueError(f"Invalid value for max_ins: {max_ins}")
    # Find the possible relationships at each position in the section,
    # not including insertions.
    rel_opts: list[tuple[int, ...]] = [(NOCOV,)] * len(refseq)
    for pos in sect.range_int:
        # Find the base in the reference sequence (pos is 1-indexed).
        ref_base = refseq[pos - 1]
        if low_qual:
            # The only option is low-quality.
            opts = encode_relate(ref_base, ref_base, '0', '1'),
        else:
            # The options are a match and three substitutions.
            opts = tuple(encode_relate(ref_base, read_base, '1', '1')
                         for read_base in (BASEA, BASEC, BASEG, BASET))
        # A deletion is an option at every position except the ends of
        # the covered region.
        if sect.end5 < pos < sect.end3:
            opts = opts + (DELET,)
        rel_opts[pos - 1] = opts
    # Iterate through all possible relationships at each position.
    margin = MIN_INDEL_GAP + 1
    for rel in product(*rel_opts):
        # Generate a relation vector from the relationships.
        relvec = np.array(rel, dtype=NP_TYPE)
        yield relvec
        if max_ins == 0:
            # Do not introduce any insertions into the read.
            continue
        # Allow up to two insertions in one read, at every position
        # except the ends of the covered region and near deletions.
        for ins1_pos5, ins2_pos5 in cwr(range(sect.end5, sect.end3), 2):
            if max_ins == 1 and ins1_pos5 != ins2_pos5:
                # If at most one insertion can be in the read, then skip
                # cases where ins1_pos5 and ins2_pos5 are not equal.
                continue
            # Skip insertions nearby any deletion.
            mrg15 = max(ins1_pos5 - margin, 0)
            mrg13 = min(ins1_pos5 + margin, relvec.size)
            if np.any(relvec[mrg15: mrg13] == DELET):
                continue
            mrg25 = max(ins2_pos5 - margin, 0)
            mrg23 = min(ins2_pos5 + margin, relvec.size)
            if np.any(relvec[mrg25: mrg23] == DELET):
                continue
            # Yield the relation vector with these insertions.
            relvec_ins = relvec.copy()
            if ins1_pos5 >= 1:
                relvec_ins[ins1_pos5 - 1] |= INS_5
            if ins1_pos5 < relvec.size:
                relvec_ins[ins1_pos5] |= INS_3
            if ins2_pos5 >= 1:
                relvec_ins[ins2_pos5 - 1] |= INS_5
            if ins2_pos5 < relvec.size:
                relvec_ins[ins2_pos5] |= INS_3
            yield relvec_ins


def iter_relvecs_all(refseq: DNA, max_ins: int = 2):
    """
    For a given reference sequence, yield every possible unambiguous
    relation vector that has at most two insertions.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    max_ins: int = 2
        Maximum number of insertions in a read.
    """
    # Use every possible pair of 5' and 3' end positions.
    for end5, end3 in cwr(range(1, len(refseq) + 1), 2):
        # Use every possible number of low-quality base calls.
        for n_low_qual in range((end3 - end5 + 1) + 1):
            # Use every possible set of low-quality positions.
            for low_qual in combinations(range(end5, end3 + 1), n_low_qual):
                # Yield every possible relation vector.
                yield from iter_relvecs_q53(refseq, low_qual, end5, end3,
                                            max_ins)


def relvec_to_read(refseq: DNA, relvec: np.ndarray, hi_qual: str, lo_qual: str,
                   ins_len: int | Sequence[int] = 1):
    """
    Infer a read sequence and quality string from a reference sequence
    and relation vector.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    relvec: ndarray
        Relation vector.
    hi_qual: str
        Character to put in the quality string at every position that is
        high-quality according to the relation vector.
    lo_qual: str
        Character to put in the quality string at every position that is
        low-quality according to the relation vector.
    ins_len: int | Sequence[int] = 1
        Number of bases to insert into the read and quality strings upon
        finding an insertion in the relation vector. If an integer, then
        insert that number of bases for every insertion. If a sequence
        of integers, then the ith insertion gets a number of bases equal
        to the ith element of ins_len.

    Returns
    -------

    """
    # Validate the relation vector and reference sequence.
    end5, end3 = validate_relvec(relvec)
    if len(refseq) != len(relvec):
        raise ValueError(f"Reference sequence has {len(refseq)} nt, "
                         f"but relation vector has {len(relvec)} nt")
    if hi_qual < lo_qual:
        raise ValueError(f"The high-quality score ({hi_qual}) is less than "
                         f"than the low-quality score ({lo_qual})")
    # Validate the quality codes.
    if len(hi_qual) != 1:
        raise ValueError(f"hi_qual must be 1 character, but got {len(hi_qual)}")
    if len(lo_qual) != 1:
        raise ValueError(f"lo_qual must be 1 character, but got {len(hi_qual)}")
    # Build the read sequence, quality scores, and CIGAR string one base
    # at a time.
    read: list[str] = list()
    qual: list[str] = list()
    cigars: list[CigarOp] = list()
    ins3_next = False
    ins_count = 0

    def add_to_cigar(op: str):
        """ Add one base of the relation vector to the CIGAR string. """
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
        rel = relvec[pos - 1]
        if rel == NOCOV:
            raise ValueError(f"Position {pos} in {end5}-{end3} is not covered")
        if rel & INS_8:
            # Specially handle insertions because they may overlap any
            # other relation except for a deletion.
            if rel & DELET:
                # Being both a deletion and an insertion is forbidden.
                raise ValueError(f"Relation {rel} in {relvec} is del and ins")
            if rel & INS_3:
                if not pos > end5:
                    # Insertions cannot occur before the beginning of
                    # the read.
                    raise ValueError(f"Position {pos} in {end5}-{end3} cannot "
                                     f"be 3' of an insertion")
                if not ins3_next:
                    raise ValueError(f"Unexpected 3' ins at {pos} in {relvec}")
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
                # Being 3' of an insertion is not allowed until the next
                # position 5' of an insertion is reached.
                ins3_next = False
            elif ins3_next:
                # Check if this position should be 3' of an insertion.
                raise ValueError(f"Missing 3' ins at {pos} in {relvec}")
            if rel & INS_5:
                if not pos < end3:
                    raise ValueError(f"Position {pos} in {end5}-{end3} cannot "
                                     f"be 5' of an insertion")
                # The current position is 5' of an insertion, so the
                # inserted base must be added after the base at the
                # current position is added. Defer the responsibility of
                # adding the inserted base to the base that lies 3' of
                # the insertion and indicate that the next base must be
                # 3' of an insertion by setting ins3_next to True.
                ins3_next = True
            # Switch off the insertion flag(s) in the relation so that
            # the base at this position can be added below.
            rel = rel & ~INS_8
        elif ins3_next:
            # Check if this position should be 3' of an insertion.
            raise ValueError(f"Missing 3' ins at {pos} in {relvec}")
        if rel == MATCH:
            # Match: Add the reference base to the read.
            read.append(ref_base)
            qual.append(hi_qual)
            add_to_cigar(CIG_MATCH)
        elif rel == DELET:
            # Deletion from the read.
            if not end5 < pos < end3:
                raise ValueError(
                    f"Deletion cannot be at position {pos} in {end5}-{end3}")
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
                read_base = BASE_DECODINGS[rel]
            except KeyError:
                raise ValueError(f"Invalid relation {rel} in {relvec}")
            if ref_base == read_base:
                raise ValueError(
                    f"Cannot substitute {ref_base} to itself in {relvec}")
            read.append(read_base)
            qual.append(hi_qual)
            add_to_cigar(CIG_SUBST)
    # Check that the 5' end was found.
    if end5 == 0:
        raise ValueError(f"Relation vector had no 5' end: {relvec}")
    if len(read) != len(qual):
        raise ValueError(
            f"Lengths of read ({len(read)}) and qual ({len(qual)}) differed")
    if not read:
        raise ValueError("Read contained no bases")
    # Assemble and return the read, quality, and CIGAR strings.
    read_seq = DNA("".join(read))
    qual_string = "".join(qual)
    cigar_string = "".join(map(str, cigars))
    return read_seq, qual_string, cigar_string, end5, end3


def ref_to_alignments(refseq: DNA,
                      max_ins: int = 2,
                      max_ins_len: int = 1,
                      max_ins_bases: int | None = None):
    """
    For a given reference sequence, return maps from every possible read
    to the CIGAR string(s) and (possibly ambiguous) relation vector for
    the read.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    max_ins: int = 2
        Maximum number of insertions in the read. Must be ≥ 0.
    max_ins_len: int = 1
        Maximum length of (i.e. number of bases in) one insertion.
        Must be ≥ 1.
    max_ins_bases: int | None = None
        Maximum total number of bases inserted. Must be ≥ `max_ins`.
        If `None`, there is no limit.
    """
    # Initialize maps of reads to CIGAR strings and relation vectors.
    quals = dict()
    cigars = defaultdict(lambda: defaultdict(list))
    relvecs = defaultdict(lambda: defaultdict(lambda: np.zeros(len(refseq),
                                                               dtype=NP_TYPE)))
    if max_ins < 0:
        raise ValueError(f"max_ins must be ≥ 0, but got {max_ins}")
    if max_ins > 0:
        if max_ins_len < 1:
            raise ValueError(f"max_ins_len must be ≥ 1, but got {max_ins_len}")
        if max_ins_bases is not None and max_ins_bases < max_ins:
            raise ValueError(f"max_ins_bases ({max_ins_bases}) "
                             f"must be ≥ max_ins ({max_ins})")
    # Iterate through all possible relation vectors.
    for relvec in iter_relvecs_all(refseq, max_ins):
        # Check if there are insertions in the relation vector.
        n_ins = int(np.count_nonzero(np.logical_and(relvec & INS_5,
                                                    relvec != NOCOV)))
        n_ins3 = int(np.count_nonzero(np.logical_and(relvec & INS_3,
                                                     relvec != NOCOV)))
        if n_ins != n_ins3:
            raise ValueError(f"Got {n_ins} 5' and {n_ins3} 3' insertions")
        if n_ins > max_ins:
            # This should not happen. Checking just in case.
            raise ValueError(f"Got {n_ins} insertions, but max_ins = {max_ins}")
        # Iterate over all possible combinations of insertion lengths.
        for ins_len in product(range(1, max_ins_len + 1), repeat=n_ins):
            if max_ins_bases is not None and sum(ins_len) > max_ins_bases:
                # Skip insertion lengths whose sum exceeds the limit.
                continue
            # Determine the read(s) corresponding to this relation vector.
            degen, qual, cigar, end5, end3 = relvec_to_read(refseq, relvec,
                                                            MAX_QUAL, MIN_QUAL,
                                                            ins_len)
            if n_ins > 0:
                # If there are insertions, find their positions.
                ins_pos = list(find_cigar_op_pos(cigar, CIG_INSRT))
                # Remove quality codes of inserted bases because -- for
                # the purpose of aggregating reads based on sequence,
                # quality score, and position -- the "fake" quality
                # scores of inserted bases should not be considered.
                qual_no_ins = "".join(q for i, q in enumerate(qual, start=1)
                                      if i not in ins_pos)
            else:
                qual_no_ins = qual
            # Count the mutations in the CIGAR string.
            nmuts = count_cigar_muts(cigar)
            for read in expand_degenerate_seq(degen):
                key = read, qual_no_ins, end5, end3
                # Record the original quality score.
                quals[key] = qual
                # Gather every CIGAR string for the read.
                cigars[key][nmuts].append(cigar)
                # Accumulate the bitwise OR of all relation vectors.
                relvecs[key][nmuts] |= relvec
    # For every read-quality-end5-end3 key, keep only the CIGAR strings
    # and relation vector with the fewest mutations.
    cigars_best: dict[tuple[DNA, str, int, int], list[str]] = dict()
    relvec_best: dict[tuple[DNA, str, int, int], np.ndarray] = dict()
    for key, cig_key in cigars.items():
        # Find the minimum number of mutations for this key.
        min_nmuts = min(cig_key)
        # Export the results with the fewest mutations.
        cigars_best[key] = cigars[key][min_nmuts]
        relvec_best[key] = relvecs[key][min_nmuts]
    return quals, cigars_best, relvec_best


def iter_alignments(*args, **kwargs):
    """ For a given reference sequence, find every read that could come
    from the reference (with up to 2 bases inserted). For each read,
    yield the (possibly ambiguous) relation vector and every possible
    CIGAR string. """
    quals, cigars, relvecs = ref_to_alignments(*args, **kwargs)
    for key, relvec in relvecs.items():
        read, _, end5, end3 = key
        qual = quals[key]
        for cigar in cigars[key]:
            yield read, qual, cigar, end5, end3, relvec
