from __future__ import annotations

from .ambrel import Deletion, Insertion, find_ambrels
from .cigar import (CIG_ALIGN, CIG_MATCH, CIG_SUBST,
                    CIG_DELET, CIG_INSRT, CIG_SCLIP,
                    parse_cigar)
from .encode import encode_match, encode_relate
from .error import RelateValueError
from ..core.rel import DELET
from ..core.seq import DNA
from ..core.xam import MAX_FLAG


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
            raise RelateValueError(f"Invalid flag: '{flag}'")
        self.flag = flag
        self.paired = bool(flag & 1)
        self.rev = bool(flag & 16)
        self.first = bool(flag & 64)
        self.second = bool(flag & 128)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.flag})"


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
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.seq = DNA(fields[9])
        self.qual = fields[10]
        if len(self.seq) != len(self.qual):
            raise RelateValueError(f"Lengths of seq ({len(self.seq)}) and qual "
                                   f"string {len(self.qual)} did not match.")

    def __str__(self):
        attrs = {attr: self.__getattribute__(attr) for attr in self.__slots__[1:]}
        return f"Read '{self.qname}' {attrs}"


def relate_read(relvec: bytearray,
                read: SamRead,
                refseq: DNA,
                length: int,
                min_qual: str,
                ambrel: bool):
    """
    Generate a relation vector of a read aligned to a reference.

    Parameters
    ----------
    read: SamRead
        Read from SAM file to be related
    relvec: bytearray
        Relation vector (initially blank) into which to write bytes;
        muts and refseq must have the same length.
    refseq: str
        Reference sequence; refseq and muts must have the same length.
    length: int (≥ 1)
        Length of the reference; must equal len(refseq) and len(muts)
    min_qual: int
        ASCII encoding of the minimum Phred score to accept a base call
    ambrel: bool
        Whether to find and label all ambiguous insertions and deletions

    """
    if len(relvec) != length:
        raise ValueError(
            f"Expected relvec to have length {length}, but got {len(relvec)}")
    if len(refseq) != length:
        raise ValueError(
            f"Expected refseq to have length {length}, but got {len(refseq)}")
    if length == 0:
        raise ValueError(f"Length of reference cannot be 0")
    # Current position in the reference (0-indexed)
    ref_idx = read.pos - 1
    if ref_idx < 0:
        raise ValueError(
            f"Read {read} mapped to a coordinate 5' of the reference")
    if ref_idx > length:
        raise ValueError(
            f"Read {read} mapped to a coordinate 3' of the reference")
    # Current position in the read (0-indexed)
    read_idx = 0
    # Record all deletions and insertions.
    dels: list[Deletion] = list()
    inns: list[Insertion] = list()
    # Read the CIGAR string one operation at a time.
    for cigar_op, op_length in parse_cigar(read.cigar):
        # Act based on the CIGAR operation and its length.
        if cigar_op == CIG_MATCH:
            # The read and reference sequences match over the entire
            # CIGAR operation.
            if ref_idx + op_length > length:
                raise ValueError("CIGAR operation overshot the reference")
            if read_idx + op_length > len(read.seq):
                raise ValueError("CIGAR operation overshot the read")
            for _ in range(op_length):
                relvec[ref_idx] &= encode_match(read.seq[read_idx],
                                                read.qual[read_idx],
                                                min_qual)
                ref_idx += 1
                read_idx += 1
        elif cigar_op == CIG_ALIGN or cigar_op == CIG_SUBST:
            # The read contains only matches or substitutions (no
            # indels) relative to the reference over the entire
            # CIGAR operation.
            if ref_idx + op_length > length:
                raise ValueError("CIGAR operation overshot the reference")
            if read_idx + op_length > len(read.seq):
                raise ValueError("CIGAR operation overshot the read")
            for _ in range(op_length):
                relvec[ref_idx] &= encode_relate(refseq[ref_idx],
                                                 read.seq[read_idx],
                                                 read.qual[read_idx],
                                                 min_qual)
                ref_idx += 1
                read_idx += 1
        elif cigar_op == CIG_DELET:
            # The portion of the reference sequence corresponding
            # to the CIGAR operation is deleted from the read.
            # Create one Deletion object for each base in the
            # reference sequence that is missing from the read.
            if ref_idx + op_length > length:
                raise ValueError("CIGAR operation overshot the reference")
            if not 1 <= read_idx < len(read.seq):
                raise ValueError(f"Deletion in {read}, pos {read_idx + 1}")
            for _ in range(op_length):
                if not 1 <= ref_idx < length - 1:
                    raise ValueError(f"Deletion in {read}, ref {ref_idx + 1}")
                dels.append(Deletion(ref_idx, read_idx))
                relvec[ref_idx] &= DELET
                ref_idx += 1
        elif cigar_op == CIG_INSRT:
            # The read contains an insertion of one or more bases
            # that are not present in the reference sequence.
            # Create one Insertion object for each base in the read
            # sequence that is not present in the reference. Every
            # mutation needs to be assigned a coordinate in the
            # region in order to appear at that coordinate in the
            # relation vector. But each inserted base, being absent
            # from the reference, does not correspond to a single
            # coordinate in the region; instead, each inserted base
            # lies between two coordinates in the region. Either of
            # these coordinates could be chosen; this code assigns
            # the 3' coordinate to the insertion. For example, if
            # two bases are inserted between coordinates 45 and 46
            # of the region, then both will be given coordinate 46.
            # The reason for this convention is that the math is
            # simpler than it would be if using the 5' coordinate.
            # Because region_idx5 is, by definition, immediately 3'
            # of the previous CIGAR operation; and the bases that
            # are inserted lie between the previous and subsequent
            # CIGAR operations; region_idx5 is the coordinate
            # immediately 3' of the inserted bases. In this special
            # case, region_idx5 also equals region_idx3 (because
            # the insertion does not consume the reference, so
            # region_idx3 += op_length was not run at the beginning
            # of this loop iteration), as well as the length of mut_vectors
            # (because Python is 0-indexed, so the length of a range
            # of indexes such as [0, 1, ... , 45] equals the value
            # of the next index in the range, 46). Thus, there are
            # three variables that already equal the 3' coordinate
            # and none that equal the 5' coordinate.
            if read_idx + op_length > len(read.seq):
                raise ValueError("CIGAR operation overshot the read")
            if not 1 <= ref_idx < length:
                raise ValueError(f"Insertion in {read}, ref {ref_idx}")
            for _ in range(op_length):
                if not 1 <= read_idx < len(read.seq) - 1:
                    raise ValueError(f"Insertion in {read}, pos {read_idx + 1}")
                inns.append(Insertion(read_idx, ref_idx))
                read_idx += 1
            # Insertions do not consume the reference, so do not add
            # any information to mut_vectors yet; it will be added later
            # via the method Insertion.stamp().
        elif cigar_op == CIG_SCLIP:
            # Bases were soft-clipped from the 5' or 3' end of the
            # read during alignment. Like insertions, they consume
            # the read but not the reference. Unlike insertions,
            # they are not mutations, so they do not require any
            # processing or boundary checking.
            if read_idx + op_length > len(read.seq):
                raise ValueError("CIGAR operation overshot the read")
            read_idx += op_length
        else:
            raise RelateValueError(
                f"Invalid CIGAR operation: '{cigar_op.decode()}'")
    # Verify that the sum of all CIGAR operations that consumed the read
    # equals the length of the read. The former equals read_idx because
    # for each CIGAR operation that consumed the read, the length of the
    # operation was added to read_idx.
    if read_idx != len(read.seq):
        raise RelateValueError(
            f"CIGAR string '{read.cigar}' consumed {read_idx} bases "
            f"from read, but read is {len(read.seq)} bases long.")
    # Add insertions to muts.
    for ins in inns:
        ins.stamp(relvec)
    # Find and label all relationships that are ambiguous due to indels.
    if ambrel and (dels or inns):
        find_ambrels(relvec, refseq, read.seq, read.qual, min_qual, dels, inns)


def validate_refs(read: SamRead, ref: str):
    if read.rname != ref:
        raise ValueError(f"Reference mismatch: read '{read.qname}' aligned to "
                         f"'{read.rname}' but is in file for reference '{ref}'")


def relate_line(relvec: bytearray, line: str,
                ref: str, refseq: DNA, length: int, qmin: str, ambrel: bool):
    read = SamRead(line)
    validate_refs(read, ref)
    relate_read(relvec, read, refseq, length, qmin, ambrel)


def relate_pair(relvec: bytearray, line1: str, line2: str,
                ref: str, refseq: DNA, length: int, qmin: str, ambrel: bool):
    # Parse lines 1 and 2 into SAM reads.
    read1 = SamRead(line1)
    validate_refs(read1, ref)
    read2 = SamRead(line2)
    validate_refs(read2, ref)
    # Ensure that reads 1 and 2 are compatible mates.
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
    # Vectorize read 1.
    relate_read(relvec, read1, refseq, length, qmin, ambrel)
    # Vectorize read 2.
    relate_read(relvec, read2, refseq, length, qmin, ambrel)

########################################################################
#                                                                      #
# ©2023, the Rouskin Lab.                                              #
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
