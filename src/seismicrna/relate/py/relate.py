from .ambindel import (IndelPod,
                       DeletionPod,
                       InsertionPod,
                       find_ambindels,
                       get_ins_rel)
from .cigar import (CIG_ALIGN,
                    CIG_MATCH,
                    CIG_SUBST,
                    CIG_DELET,
                    CIG_INSRT,
                    CIG_SCLIP,
                    parse_cigar)
from .encode import encode_match, encode_relate
from .error import RelateError
from ...core.ngs import MAX_FLAG
from ...core.rel import MATCH, DELET, NOCOV, IRREC


class SamFlag(object):
    """ Represents the set of 12 boolean flags for a SAM record. """

    # Define __slots__ to improve speed and memory performance.
    __slots__ = ["flag", "paired", "rev", "first", "second"]

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
            raise RelateError(f"Invalid flag: {repr(flag)}")
        self.flag = flag
        self.paired = bool(flag & 1)
        self.rev = bool(flag & 16)
        self.first = bool(flag & 64)
        self.second = bool(flag & 128)

    def __repr__(self):
        return f"{type(self).__name__}({self.flag})"


class SamRead(object):
    """ One read in a SAM file. """

    # Define __slots__ to improve speed and memory performance.
    __slots__ = ["name", "flag", "ref", "pos", "mapq", "cigar", "seq", "qual"]

    # Minimum number of fields in a valid SAM record
    MIN_FIELDS = 11

    def __init__(self, line: str):
        fields = line.rstrip().split('\t')
        if len(fields) < self.MIN_FIELDS:
            raise RelateError(
                f"Each line in the SAM file must have ≥ {self.MIN_FIELDS} "
                f"fields, but got {len(fields)} in\n{line}"
            )
        self.name = fields[0]
        self.flag = SamFlag(int(fields[1]))
        self.ref = fields[2]
        self.pos = int(fields[3])
        if self.pos < 1:
            raise RelateError(f"Position must be ≥ 1, but got {self.pos}")
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.seq = fields[9]
        self.qual = fields[10]
        if len(self.seq) != len(self.qual):
            raise RelateError(f"Lengths of seq ({len(self.seq)}) and qual "
                              f"string {len(self.qual)} did not match")

    def __str__(self):
        attrs = {attr: getattr(self, attr) for attr in self.__slots__[1:]}
        return f"Read {repr(self.name)} {attrs}"


def _add_indel(pods: list[IndelPod],
               pod_type: type[IndelPod],
               opposite: int,
               lateral3: int):
    if not pods or not isinstance(pods[-1], pod_type):
        pods.append(pod_type(len(pods)))
    pod = pods[-1]
    indel_type = pod.indel_type()
    return indel_type(opposite, lateral3, pod)


def _calc_rels_read(read: SamRead,
                    ref_seq: str,
                    min_qual: str,
                    insert3: bool,
                    ambindel: bool,
                    clip_end5: int,
                    clip_end3: int):
    """
    Find the relationships between a read and a reference.

    Parameters
    ----------
    read: SamRead
        Read from SAM file to be related
    ref_seq: str
        Reference sequence; refseq and muts must have the same length.
    min_qual: bytes
        ASCII encoding of the minimum Phred score to accept a base call
    ambindel: bool
        Whether to find and label all ambiguous insertions and deletions
    insert3: bool
        Whether to mark an insertion on the base immediately 3' (True)
        or 5' (False) of the insertion.
    clip_end5: int
        Number of bases to clip from the 5' end of the read.
    clip_end3: int
        Number of bases to clip from the 3' end of the read.

    Returns
    -------
    tuple[int, int, dict[int, int]]
        - the 5' coordinate of the read
        - the 3' coordinate of the read
        - the relationship code for each mutation in the read
    """
    read_length = len(read.seq)
    ref_length = len(ref_seq)
    if not 1 <= read.pos <= ref_length:
        raise RelateError(f"Mapping position must be ≥ 1 and ≤ reference "
                          f"length ({ref_length}), but got {read.pos}")
    # Position in the reference: can be 0- or 1-indexed (initially 0).
    ref_pos = read.pos - 1
    # Position in the read: can be 0- or 1-indexed (initially 0).
    read_pos = 0
    # Ends of the read, in the read coordinates (1-indexed).
    read_end5 = 1
    read_end3 = read_length
    # Record the relationship for each mutated position.
    rels: dict[int, int] = dict()
    # Record the positions of deletions and insertions.
    pods: list[IndelPod] = list()
    # Read the CIGAR string one operation at a time.
    for cigar_op, op_length in parse_cigar(read.cigar):
        # Act based on the CIGAR operation and its length.
        if cigar_op == CIG_MATCH:
            # The read and reference sequences match over the entire
            # CIGAR operation.
            if ref_pos + op_length > ref_length:
                raise RelateError("CIGAR operation overshot the reference")
            if read_pos + op_length > len(read.seq):
                raise RelateError("CIGAR operation overshot the read")
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
            if ref_pos + op_length > ref_length:
                raise RelateError("CIGAR operation overshot the reference")
            if read_pos + op_length > len(read.seq):
                raise RelateError("CIGAR operation overshot the read")
            for _ in range(op_length):
                rel = encode_relate(ref_seq[ref_pos],
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
            if ref_pos + op_length > ref_length:
                raise RelateError("CIGAR operation overshot the reference")
            read_pos += 1  # 1-indexed now until explicitly reset
            if not 1 < read_pos <= len(read.seq):
                raise RelateError(f"Deletion in {read}, pos {read_pos}")
            for _ in range(op_length):
                ref_pos += 1  # 1-indexed now until this iteration ends
                if not 1 < ref_pos < ref_length:
                    raise RelateError(f"Deletion in {read}, ref {ref_pos}")
                rels[ref_pos] = DELET
                _add_indel(pods, DeletionPod, ref_pos, read_pos)
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
                raise RelateError("CIGAR operation overshot the read")
            ref_pos += 1  # 1-indexed now until explicitly reset
            if not 1 < ref_pos <= ref_length:
                raise RelateError(f"Insertion in {read}, ref {ref_pos}")
            for _ in range(op_length):
                read_pos += 1  # 1-indexed now until this iteration ends
                if not 1 < read_pos < len(read.seq):
                    raise RelateError(f"Insertion in {read}, pos {read_pos}")
                _add_indel(pods, InsertionPod, read_pos, ref_pos)
            # Insertions can lie on top of other relationships, so do
            # not the information to rels yet; it will be added later.
            ref_pos -= 1  # 0-indexed now
        elif cigar_op == CIG_SCLIP:
            # Bases were soft-clipped from the 5' or 3' end of the
            # read during alignment. Like insertions, they consume
            # the read but not the reference. Unlike insertions,
            # they are not mutations, so they do not require any
            # processing or boundary checking.
            next_read_pos = read_pos + op_length
            if next_read_pos > len(read.seq):
                raise RelateError("CIGAR operation overshot the read")
            if read_pos == 0:
                # This is the soft clip from the 5' end of the read.
                if read_end5 != 1:
                    raise RelateError(f"{repr(read.cigar)} has >1 5' soft clip")
                read_end5 += op_length
            else:
                # This is the soft clip from the 3' end of the read.
                if next_read_pos != read_length:
                    raise RelateError(
                        f"{repr(read.cigar)} has a soft clip in the middle"
                    )
                if read_end3 != read_length:
                    raise RelateError(f"{repr(read.cigar)} has >1 3' soft clip")
                assert read_pos == read_end3 - op_length
                read_end3 = read_pos
            read_pos = next_read_pos
        else:
            raise RelateError(
                f"Unsupported CIGAR operation: {repr(cigar_op.decode())}"
            )
    # For consistency with the C version, the read must consume at least
    # one position in the reference.
    if ref_pos == read.pos - 1:
        raise RelateError(f"CIGAR string {repr(read.cigar)} consumed "
                          f"0 bases in the reference")
    # Verify that the sum of all CIGAR operations that consumed the read
    # equals the length of the read.
    if read_pos != len(read.seq):
        raise RelateError(
            f"CIGAR string {repr(read.cigar)} consumed {read_pos} bases "
            f"in the read, but the read is {len(read.seq)} bases long"
        )
    # Add insertions to rels.
    ins_rel = get_ins_rel(insert3)
    for pod in pods:
        if isinstance(pod, InsertionPod):
            for ins in pod.indels:
                ins_pos = ins.get_lateral(insert3)
                if 1 <= ins_pos <= ref_length:
                    # The position at which the insertion is located may
                    # have already been assigned another mutation;
                    # if so, then add the insertion on top.
                    rels[ins_pos] = rels.get(ins_pos, IRREC) | ins_rel
    # Find and label all relationships that are ambiguous due to indels.
    if ambindel and pods:
        find_ambindels(rels=rels,
                       pods=pods,
                       insert3=insert3,
                       ref_seq=ref_seq,
                       read_seq=read.seq,
                       read_qual=read.qual,
                       min_qual=min_qual,
                       ref_end5=read.pos,
                       ref_end3=ref_pos,
                       read_end5=read_end5,
                       read_end3=read_end3)
    # Clip bases from the 5' and 3' ends of the read.
    if clip_end5 < 0:
        raise RelateError(f"clip_end5 must be ≥ 0, but got {clip_end5}")
    if clip_end3 < 0:
        raise RelateError(f"clip_end3 must be ≥ 0, but got {clip_end3}")
    ref_end5 = min(read.pos + clip_end5, ref_length + 1)
    ref_end3 = max(ref_pos - clip_end3, 0)
    # Convert rels to a non-default dict and select only the positions
    # between ref_end5 and ref_end3.
    rels = {pos: rel for pos, rel in rels.items()
            if ref_end5 <= pos <= ref_end3}
    return ref_end5, ref_end3, rels


def _validate_read(read: SamRead, ref: str, min_mapq: int):
    if read.ref != ref:
        raise RelateError(f"Read {repr(read.name)} mapped to a reference named "
                          f"{repr(read.ref)} but is in an alignment map file "
                          f"for a reference named {repr(ref)}")
    if read.mapq < min_mapq:
        raise RelateError(f"Read {repr(read.name)} mapped with quality score "
                          f"{read.mapq}, less than the minimum of {min_mapq}")


def _validate_pair(read1: SamRead, read2: SamRead):
    """ Ensure that reads 1 and 2 are compatible mates. """
    if not read1.flag.paired:
        raise RelateError(f"Read 1 ({read1.name}) was not paired, "
                          f"but read 2 ({read2.name}) was given")
    if not read2.flag.paired:
        raise RelateError(f"Read 2 ({read2.name}) was not paired, "
                          f"but read 1 ({read1.name}) was given")
    if read1.name != read2.name:
        raise RelateError(f"Got different names for reads "
                          f"1 ({read1.name}) and 2 ({read2.name})")
    if not (read1.flag.first and read2.flag.second):
        raise RelateError(f"Read {repr(read1.name)} had mate 1 "
                          f"labeled {2 - read1.flag.first} and mate 2 "
                          f"labeled {1 + read2.flag.second}")
    if read1.flag.rev == read2.flag.rev:
        raise RelateError(f"Read {repr(read1.name)} had "
                          "mates 1 and 2 facing the same way")


def _merge_rels(end5f: int,
                end3f: int,
                relsf: dict[int, int],
                end5r: int,
                end3r: int,
                relsr: dict[int, int]):
    merged_rels = dict()
    for pos in relsf | relsr:
        relf = relsf.get(pos, MATCH if end5f <= pos <= end3f else NOCOV)
        relr = relsr.get(pos, MATCH if end5r <= pos <= end3r else NOCOV)
        rel = relf & relr
        if rel != MATCH:
            if rel == NOCOV:
                raise RelateError(f"Cannot merge non-covered position {pos}")
            merged_rels[pos] = rel
    return merged_rels


def _merge_mates(end5f: int,
                 end3f: int,
                 relsf: dict[int, int],
                 end5r: int,
                 end3r: int,
                 relsr: dict[int, int],
                 overhangs: bool):
    if not overhangs:
        # The 5' end of the reverse mate cannot extend past the 5' end
        # of the forward mate.
        if end5r < end5f:
            end5r = end5f
            relsr = {pos: rel for pos, rel in relsr.items() if pos >= end5r}
        # The 3' end of the forward mate cannot extend past the 3' end
        # of the reverse mate.
        if end3f > end3r:
            end3f = end3r
            relsf = {pos: rel for pos, rel in relsf.items() if pos <= end3f}
    rels = _merge_rels(end5f, end3f, relsf, end5r, end3r, relsr)
    return ([end5f, end5r], [end3f, end3r]), rels


def calc_rels_lines(line1: str,
                    line2: str,
                    ref: str,
                    refseq: str,
                    min_mapq: int,
                    min_qual: int,
                    insert3: bool,
                    ambindel: bool,
                    overhangs: bool,
                    clip_end5: int = 0,
                    clip_end3: int = 0):
    # Generate the relationships for read 1.
    read1 = SamRead(line1)
    _validate_read(read1, ref, min_mapq)
    end51, end31, rels1 = _calc_rels_read(read1,
                                          refseq,
                                          chr(min_qual),
                                          insert3,
                                          ambindel,
                                          clip_end5,
                                          clip_end3)
    if line2:
        if line2 == line1:
            # This read is paired-end but comprises only one mate, which
            # can occur only if Bowtie2 is run in mixed mode.
            ends = [end51, end51], [end31, end31]
            rels = rels1
        else:
            # Generate the relationships for read 2.
            read2 = SamRead(line2)
            _validate_read(read2, ref, min_mapq)
            _validate_pair(read1, read2)
            end52, end32, rels2 = _calc_rels_read(read2,
                                                  refseq,
                                                  chr(min_qual),
                                                  insert3,
                                                  ambindel,
                                                  clip_end5,
                                                  clip_end3)
            # Determine which read (1 or 2) faces forward and reverse.
            if read2.flag.rev:
                end5f, end3f, relsf = end51, end31, rels1
                end5r, end3r, relsr = end52, end32, rels2
            else:
                end5f, end3f, relsf = end52, end32, rels2
                end5r, end3r, relsr = end51, end31, rels1
            # Merge the relationships and end coordinates.
            ends, rels = _merge_mates(
                end5f, end3f, relsf, end5r, end3r, relsr, overhangs,
            )
    else:
        # The read is single-ended.
        ends = [end51], [end31]
        rels = rels1
    return ends, rels

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
