from collections import defaultdict
from functools import partial
from itertools import product
from typing import Sequence

import networkx as nx

from .cigarop import (count_cigar_muts,
                      find_cigar_op_pos_read,
                      find_cigar_op_pos_ref)
from .infer import infer_read
from .iterrel import iter_relvecs_all
from ..py.cigar import CIG_DELET, CIG_INSRT, parse_cigar
from ...core.rel import INS_5, INS_3, MATCH
from ...core.seq import DNA, expand_degenerate_seq


def _connected_indels(indels1: Sequence[int], indels2: Sequence[int]):
    if len(indels1) != len(indels2):
        return False
    for i1, i2 in zip(indels1, indels2, strict=True):
        if abs(i1 - i2) > 1:
            return False
    return True


def _indel_depths(cigar: str):
    depths = [0]
    for op, length in parse_cigar(cigar):
        if op == CIG_DELET:
            for _ in range(length):
                depths.append(depths[-1] - 1)
        elif op == CIG_INSRT:
            for _ in range(length):
                depths.append(depths[-1] + 1)
    return tuple(depths[1:])


def ref_to_alignments(refseq: DNA, *,
                      insert3: bool,
                      max_ins: int = 0,
                      max_ins_len: int = 1,
                      max_ins_bases: int | None = None):
    """ For a given reference sequence, map every possible read to its
    CIGAR string(s) and (possibly ambiguous) relationships.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    insert3: bool:
        Whether to mark the base 5' or 3' of an insertion.
    max_ins: int
        Maximum number of insertions in the read. Must be ≥ 0.
    max_ins_len: int
        Maximum length of (i.e. number of bases in) one insertion.
        Must be ≥ 1.
    max_ins_bases: int | None
        Maximum total number of bases inserted. Must be ≥ `max_ins`.
        If `None`, there is no limit.
    """
    # Initialize maps of reads to CIGAR strings and mutations.
    quals = dict()
    reads = defaultdict(partial(defaultdict, list))
    if max_ins < 0:
        raise ValueError(f"max_ins must be ≥ 0, but got {max_ins}")
    if max_ins > 0:
        if max_ins_len < 1:
            raise ValueError(f"max_ins_len must be ≥ 1, but got {max_ins_len}")
        if max_ins_bases is not None and max_ins_bases < max_ins:
            raise ValueError(f"max_ins_bases ({max_ins_bases}) "
                             f"must be ≥ max_ins ({max_ins})")
    # Iterate through all possible relationships.
    for end5, end3, muts in iter_relvecs_all(refseq, insert3, max_ins):
        # Check if there are insertions.
        if insert3:
            ins = INS_3
            anti_ins = INS_5
        else:
            ins = INS_5
            anti_ins = INS_3
        if any(rel & anti_ins for rel in muts.values()):
            # This should be impossible. Checking just in case.
            raise ValueError(f"muts cannot contain {anti_ins}, but got {muts}")
        n_ins = sum(bool(rel & ins) for rel in muts.values())
        if n_ins > max_ins:
            # This should be impossible. Checking just in case.
            raise ValueError(f"Got {n_ins} insertions, but max_ins = {max_ins}")
        # Iterate over all possible combinations of insertion lengths.
        for ins_len in product(range(1, max_ins_len + 1), repeat=n_ins):
            if max_ins_bases is not None and sum(ins_len) > max_ins_bases:
                # Skip insertion lengths whose sum exceeds the limit.
                continue
            # Determine the read(s) corresponding to these relationships.
            degen, qual, cigar = infer_read(refseq,
                                            end5,
                                            end3,
                                            muts,
                                            ins_len=ins_len)
            if n_ins > 0:
                # Remove quality codes of inserted bases because -- for
                # the purpose of aggregating reads based on sequence,
                # quality score, and position -- the "fake" quality
                # scores of inserted bases should not be considered.
                ins_pos = set(find_cigar_op_pos_read(cigar, CIG_INSRT))
                qual_no_ins = "".join(q for i, q in enumerate(qual, start=1)
                                      if i not in ins_pos)
            else:
                qual_no_ins = qual
            # Count the mutations in the CIGAR string.
            num_muts = count_cigar_muts(cigar)
            for read in expand_degenerate_seq(degen):
                key = read, qual_no_ins, end5, end3
                # Record the original quality score.
                quals[key] = qual
                # Gather every CIGAR string for the read.
                reads[key][num_muts].append((cigar, tuple(muts.items())))
    # Accumulate the bitwise OR of all relationships for each read,
    # quality, ends, and number of mutations.
    cigars_best = dict()
    muts_best = dict()
    for key, cigar_muts in reads.items():
        _, _, end5, end3 = key
        cigars_best[key] = list()
        muts_best[key] = list()
        # Find the minimum number of mutations for this key.
        min_muts = min(cigar_muts)
        # For all reads with that minimum number of mutations, group the
        # CIGAR strings that can be transformed into each other through
        # a series of moving each indel one step.
        transforms = nx.Graph()
        for (cigar, muts) in cigar_muts[min_muts]:
            dels = tuple(find_cigar_op_pos_ref(cigar, CIG_DELET, end5))
            inns = tuple(find_cigar_op_pos_read(cigar, CIG_INSRT))
            depths = _indel_depths(cigar)
            node = cigar, muts, dels, inns, depths
            add_edges = list()
            for other_node in transforms.nodes:
                # Consider two possible CIGAR strings "connected" if:
                # - the indels are in the same order, i.e. no deletions
                #   and insertions have swapped (depths == p)
                # - and either of the following is true:
                #   - their deletions are in the same positions and no
                #     insertion has moved more than one position
                #     (dels == d and _connected_indels(inns, i))
                #   - their insertions are in the same positions and no
                #     deletion has moved more than one position
                #     (inns == i and _connected_indels(dels, d))
                c, m, d, i, p = other_node
                if (depths == p
                        and ((dels == d and _connected_indels(inns, i))
                             or
                             (inns == i and _connected_indels(dels, d)))):
                    add_edges.append((node, other_node))
            if add_edges:
                for node1, node2 in add_edges:
                    transforms.add_edge(node1, node2)
            else:
                transforms.add_node(node)
        # For each group of CIGAR strings that can be transformed into
        # each other, take the bitwise union of the mutations.
        for nodes in nx.connected_components(transforms):
            cigars = list()
            muts_union = defaultdict(int)
            for cigar, muts, _, _, _ in nodes:
                muts = dict(muts)
                for pos in range(end5, end3 + 1):
                    muts_union[pos] |= muts.get(pos, MATCH)
                cigars.append(cigar)
            cigars_best[key].append(cigars)
            muts_best[key].append(
                {pos: rel for pos, rel in muts_union.items() if rel != MATCH}
            )
    return quals, cigars_best, muts_best


def iter_alignments(*args, **kwargs):
    """ For a given reference sequence, find every read that could come
    from the reference (with up to 2 bases inserted). For each read,
    yield the (possibly ambiguous) relationships and every possible
    CIGAR string. """
    quals, all_cigars, all_muts = ref_to_alignments(*args, **kwargs)
    for key in all_muts:
        read, _, end5, end3 = key
        qual = quals[key]
        for cigars, muts in zip(all_cigars[key], all_muts[key], strict=True):
            for cigar in cigars:
                yield read, qual, cigar, end5, end3, muts
