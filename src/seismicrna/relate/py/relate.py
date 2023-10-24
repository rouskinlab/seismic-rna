from __future__ import annotations

from collections import defaultdict

from .ambrel import Deletion, Insertion, find_ambrels
from .cigar import (CIG_ALIGN,
                    CIG_MATCH,
                    CIG_SUBST,
                    CIG_DELET,
                    CIG_INSRT,
                    CIG_SCLIP,
                    parse_cigar)
from .encode import encode_match, encode_relate
from .error import RelateValueError
from ...core.rel import MATCH, DELET, NOCOV
from ...core.seq import DNA
from ...core.ngs import MAX_FLAG


class SamFlag(object):
    """ Represents the set of 12 boolean flags for a SAM record. """

    # Define __slots__ to improve speed and memory performance.
    __slots__ = "flag", "paired", "rev", "first", "second"

    def __init__(self, flag: int):
        """
        Validate the integer value of the SAM flag, then set the flags
        that are needed.

        Parameters
        ----------
        flag: int
            The integer value of the SAM flag. For documentation, see
            https://samtools.github.io/hts-specs/

        """
        if not 0 <= flag <= MAX_FLAG:
            raise RelateValueError(f"Invalid flag: {repr(flag)}")
        self.flag = flag
        self.paired = bool(flag & 1)
        self.rev = bool(flag & 16)
        self.first = bool(flag & 64)
        self.second = bool(flag & 128)

    def __repr__(self):
        return f"{type(self).__name__}({self.flag})"


class SamRead(object):
    # Define __slots__ to improve speed and memory performance.
    __slots__ = "qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual"

    # Minimum number of fields in a valid SAM record
    MIN_FIELDS = 11

    def __init__(self, line: str):
        fields = line.rstrip().split('\t')
        if len(fields) < self.MIN_FIELDS:
            raise RelateValueError(f"Invalid SAM line:\n{line}")
        self.qname = fields[0]
        self.flag = SamFlag(int(fields[1]))
        self.rname = fields[2]
        self.pos = int(fields[3])
        if self.pos < 1:
            raise ValueError(f"Position must be ≥ 1, but got {self.pos}")
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.seq = DNA(fields[9])
        self.qual = fields[10]
        if len(self.seq) != len(self.qual):
            raise RelateValueError(f"Lengths of seq ({len(self.seq)}) and qual "
                                   f"string {len(self.qual)} did not match.")

    def __str__(self):
        attrs = {attr: getattr(self, attr) for attr in self.__slots__[1:]}
        return f"Read '{self.qname}' {attrs}"


def find_rels_read(read: SamRead, refseq: DNA, min_qual: str, ambrel: bool):
    """
    Find the relationships between a read and a reference.

    Parameters
    ----------
    read: SamRead
        Read from SAM file to be related
    refseq: str
        Reference sequence; refseq and muts must have the same length.
    min_qual: int
        ASCII encoding of the minimum Phred score to accept a base call
    ambrel: bool
        Whether to find and label all ambiguous insertions and deletions

    Returns
    -------
    tuple[int, int, dict[int, int]]
        - the 5' coordinate of the read
        - the 3' coordinate of the read
        - the relationship code for each mutation in the read
    """
    reflen = len(refseq)
    if read.pos > reflen:
        raise ValueError(f"Position must be ≤ reference length ({reflen}), "
                         f"but got {read.pos}")
    # Position in the reference: can be 0- or 1-indexed (initially 0).
    ref_pos = read.pos - 1
    # Position in the read: can be 0- or 1-indexed (initially 0).
    read_pos = 0
    # Record the relationship for each mutated position.
    rels: dict[int, int] = defaultdict(lambda: MATCH)
    # Record all deletions and insertions.
    dels: list[Deletion] = list()
    inns: list[Insertion] = list()
    # Read the CIGAR string one operation at a time.
    for cigar_op, op_length in parse_cigar(read.cigar):
        # Act based on the CIGAR operation and its length.
        if cigar_op == CIG_MATCH:
            # The read and reference sequences match over the entire
            # CIGAR operation.
            if ref_pos + op_length > reflen:
                raise ValueError("CIGAR operation overshot the reference")
            if read_pos + op_length > len(read.seq):
                raise ValueError("CIGAR operation overshot the read")
            for _ in range(op_length):
                rel = encode_match(read.seq[read_pos],
                                   read.qual[read_pos],
                                   min_qual)
                ref_pos += 1  # 1-indexed now until this iteration ends
                read_pos += 1  # 1-indexed now until this iteration ends
                if rel != MATCH:
                    rels[ref_pos] = rel
        elif cigar_op == CIG_ALIGN or cigar_op == CIG_SUBST:
            # There are only matches or substitutions over the entire
            # CIGAR operation.
            if ref_pos + op_length > reflen:
                raise ValueError("CIGAR operation overshot the reference")
            if read_pos + op_length > len(read.seq):
                raise ValueError("CIGAR operation overshot the read")
            for _ in range(op_length):
                rel = encode_relate(refseq[ref_pos],
                                    read.seq[read_pos],
                                    read.qual[read_pos],
                                    min_qual)
                ref_pos += 1  # 1-indexed now until this iteration ends
                read_pos += 1  # 1-indexed now until this iteration ends
                if rel != MATCH:
                    rels[ref_pos] = rel
        elif cigar_op == CIG_DELET:
            # The portion of the reference sequence corresponding to the
            # CIGAR operation is deleted from the read. Make a Deletion
            # object for each base in the reference sequence that has
            # been deleted from the read.
            if ref_pos + op_length > reflen:
                raise ValueError("CIGAR operation overshot the reference")
            read_pos += 1  # 1-indexed now until explicitly reset
            if not 1 < read_pos <= len(read.seq):
                raise ValueError(f"Deletion in {read}, pos {read_pos}")
            for _ in range(op_length):
                ref_pos += 1  # 1-indexed now until this iteration ends
                if not 1 < ref_pos < reflen:
                    raise ValueError(f"Deletion in {read}, ref {ref_pos}")
                dels.append(Deletion(ref_pos, read_pos))
                rels[ref_pos] = DELET
            read_pos -= 1  # 0-indexed now
        elif cigar_op == CIG_INSRT:
            # The read contains an insertion of one or more bases that
            # are not present in the reference sequence. Create one
            # Insertion object for each such base. Every mutation needs
            # a coordinate in the reference sequence. But each inserted
            # base, being absent from the reference, does not correspond
            # to a single reference coordinate but rather lies between
            # two coordinates. Either of these coordinates could be used
            # for an insertion; this algorithm uses the 3' coordinate.
            # For example, if two bases are inserted between coordinates
            # 45 and 46 in the reference, then both inserted bases will
            # be assigned coordinate 46. The reason for this convention
            # is that the math is simpler than it would be if using the
            # 5' coordinate. Because all inserted bases lie between the
            # previous and subsequent CIGAR operations, and ref_pos is
            # the position immediately 3' of the previous operation,
            # ref_pos is naturally the coordinate 3' of the insertion.
            if read_pos + op_length > len(read.seq):
                raise ValueError("CIGAR operation overshot the read")
            ref_pos += 1  # 1-indexed now until explicitly reset
            if not 1 < ref_pos <= reflen:
                raise ValueError(f"Insertion in {read}, ref {ref_pos}")
            for _ in range(op_length):
                read_pos += 1  # 1-indexed now until this iteration ends
                if not 1 < read_pos < len(read.seq):
                    raise ValueError(f"Insertion in {read}, pos {read_pos}")
                inns.append(Insertion(read_pos, ref_pos))
            # Insertions do not consume the reference, so do not add
            # any information to mut_vectors yet; it will be added later
            # via the method Insertion.stamp().
            ref_pos -= 1  # 0-indexed now
        elif cigar_op == CIG_SCLIP:
            # Bases were soft-clipped from the 5' or 3' end of the
            # read during alignment. Like insertions, they consume
            # the read but not the reference. Unlike insertions,
            # they are not mutations, so they do not require any
            # processing or boundary checking.
            if read_pos + op_length > len(read.seq):
                raise ValueError("CIGAR operation overshot the read")
            read_pos += op_length
        else:
            raise RelateValueError(
                f"Invalid CIGAR operation: '{cigar_op.decode()}'")
    # Verify that the sum of all CIGAR operations that consumed the read
    # equals the length of the read. The former equals read_pos because
    # for each CIGAR operation that consumed the read, the length of the
    # operation was added to read_pos.
    if read_pos != len(read.seq):
        raise RelateValueError(
            f"CIGAR string '{read.cigar}' consumed {read_pos} bases "
            f"from read, but read is {len(read.seq)} bases long.")
    # Add insertions to muts.
    for ins in inns:
        ins.stamp(rels, len(refseq))
    # Find and label all relationships that are ambiguous due to indels.
    if ambrel and (dels or inns):
        find_ambrels(rels, refseq, read.seq, read.qual, min_qual, dels, inns)
    return read.pos, ref_pos, dict(rels)


def validate_read(read: SamRead, ref: str, min_mapq: int):
    if read.rname != ref:
        raise ValueError(f"Read {repr(read.qname)} mapped to a reference named "
                         f"{repr(read.rname)} but is in an alignment map file "
                         f"for a reference named {repr(ref)}")
    if read.mapq < min_mapq:
        raise ValueError(f"Read {repr(read.qname)} mapped with quality score "
                         f"{read.mapq}, less than the minimum of {min_mapq}")


def validate_pair(read1: SamRead, read2: SamRead):
    """ Ensure that reads 1 and 2 are compatible mates. """
    if not read1.flag.paired:
        raise RelateValueError(f"Read 1 ({read1.qname}) was not paired, "
                               f"but read 2 ('{read2.qname}') was given")
    if not read2.flag.paired:
        raise RelateValueError(f"Read 2 ({read2.qname}) was not paired, "
                               f"but read 1 ({read1.qname}) was given")
    if read1.qname != read2.qname:
        raise RelateValueError(f"Reads 1 ({read1.qname}) and 2 "
                               f"({read2.qname}) had different names")
    if not (read1.flag.first and read2.flag.second):
        raise RelateValueError(f"Read '{read1.qname}' had mate 1 "
                               f"labeled {2 - read1.flag.first} and mate 2 "
                               f"labeled {1 + read2.flag.second}")
    if read1.flag.rev == read2.flag.rev:
        raise RelateValueError(f"Read '{read1.qname}' had "
                               "mates 1 and 2 facing the same way")


def merge_ends(end51: int, end31: int, end52: int, end32: int):
    return (min(end51, end52),
            max(end51, end52),
            min(end31, end32),
            max(end31, end32))


def merge_rels(rels1: dict[int, int], rels2: dict[int, int]):
    return {ref_idx: rels1.get(ref_idx, NOCOV) & rels2.get(ref_idx, NOCOV)
            for ref_idx in rels1 | rels2}


def find_rels_line(line1: str,
                   line2: str,
                   ref: str,
                   refseq: DNA,
                   min_mapq: int,
                   qmin: str,
                   ambrel: bool):
    read1 = SamRead(line1)
    validate_read(read1, ref, min_mapq)
    end5, end3, rels = find_rels_read(read1, refseq, qmin, ambrel)
    if line2:
        read2 = SamRead(line2)
        validate_read(read2, ref, min_mapq)
        validate_pair(read1, read2)
        end52, end32, rels2 = find_rels_read(read2, refseq, qmin, ambrel)
        end5, mid5, mid3, end3 = merge_ends(end5, end3, end52, end32)
        rels = merge_rels(rels, rels2)
    else:
        mid5, mid3 = end5, end3
    return read1.qname, end5, mid5, mid3, end3, rels


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
