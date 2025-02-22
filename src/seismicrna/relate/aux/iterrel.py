from itertools import (product,
                       combinations,
                       combinations_with_replacement as cwr)
from typing import Sequence

from ..py.encode import encode_relate
from ...core.ngs import LO_QUAL, OK_QUAL, HI_QUAL
from ...core.rel import DELET, INS_5, INS_3, MATCH
from ...core.seq import DNA, Region


def iter_relvecs_q53(refseq: DNA,
                     low_qual: Sequence[int] = (),
                     end5: int | None = None,
                     end3: int | None = None,
                     insert3: bool = True,
                     max_ins: int | None = None):
    """
    For a given reference sequence, yield every possible unambiguous
    relationships between positions end5 and end3 that follows the
    given low-quality positions and has at most two insertions.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    low_qual: Sequence[int]
        List of positions in the read that are low-quality.
    end5: int | None
        5' end of the read; 1-indexed with respect to `refseq`.
    end3: int | None
        3' end of the read; 1-indexed with respect to `refseq`.
    insert3: bool
        Whether to mark the base 5' or 3' of an insertion.
    max_ins: int | None
        Maximum number of insertions in the read.
    """
    if max_ins is not None and max_ins < 0:
        raise ValueError(f"max_ins must be ≥ 0, but got {max_ins}")
    low_qual = set(low_qual)
    # Determine the region of the reference sequence that is occupied
    # by the read.
    region = Region("", refseq, end5=end5, end3=end3)
    all_positions = set(region.range_int)
    all5_positions = (set(region.range_int[:-1])
                      if region.length > 0
                      else set())
    if low_qual - all_positions:
        raise ValueError(f"Invalid positions in low_qual: "
                         f"{sorted(low_qual - set(region.range))}")
    # Find the possible relationships at each position in the region,
    # not including insertions.
    rel_opts = list()
    for pos in region.range_int:
        # Find the base in the reference sequence (pos is 1-indexed).
        ref_base = refseq[pos - 1]
        if pos in low_qual:
            # The only option is low-quality.
            opts = [encode_relate(ref_base, ref_base, LO_QUAL, OK_QUAL)]
        else:
            # The options are three substitutions.
            opts = [encode_relate(ref_base, read_base, HI_QUAL, OK_QUAL)
                    for read_base in DNA.four()]
        # A deletion is an option at every position except the ends of
        # the covered region.
        if region.end5 < pos < region.end3:
            opts.append(DELET)
        rel_opts.append([(pos, rel) for rel in opts])
    # Iterate through all possible relationships at each position.
    for rels in product(*rel_opts):
        # Generate one possible relationship.
        relvec = dict((pos, rel) for pos, rel in rels if rel != MATCH)
        yield region.end5, region.end3, relvec
        if (max_ins is None or max_ins > 0) and not low_qual:
            no_ins5_pos = {pos for del_pos, rel in relvec.items()
                           if rel == DELET
                           for pos in range(del_pos - 1,
                                            del_pos + 1)}
            ins5_pos = sorted(all5_positions - no_ins5_pos)
            for num_ins in range(1, 1 + min(max_ins
                                            if max_ins is not None
                                            else len(ins5_pos),
                                            len(ins5_pos))):
                for ins5s in combinations(ins5_pos, num_ins):
                    rv_ins = relvec.copy()
                    for i5 in ins5s:
                        if insert3:
                            i3 = i5 + 1
                            rv_ins[i3] = rv_ins.get(i3, 0) | INS_3
                        else:
                            rv_ins[i5] = rv_ins.get(i5, 0) | INS_5
                    yield region.end5, region.end3, rv_ins


def iter_relvecs_all(refseq: DNA,
                     insert3: bool = True,
                     max_ins: int | None = None):
    """
    For a given reference sequence, yield every possible unambiguous
    set of relationships that has at most two insertions.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    insert3: bool
        Whether to mark the base 5' or 3' of an insertion.
    max_ins: int | None
        Maximum number of insertions in a read.
    """
    # Use every possible pair of 5' and 3' end positions.
    for end5, end3 in cwr(range(1, len(refseq) + 1), 2):
        # Use every possible number of low-quality base calls.
        for n_low_qual in range((end3 - (end5 - 1)) + 1):
            # Use every possible set of low-quality positions.
            for low_qual in combinations(range(end5, end3 + 1), n_low_qual):
                # Yield every possible set of relationships.
                yield from iter_relvecs_q53(refseq,
                                            low_qual,
                                            end5,
                                            end3,
                                            insert3,
                                            max_ins)
