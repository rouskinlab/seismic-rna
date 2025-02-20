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
                    CIG_INTRN,
                    parse_cigar)
from .encode import encode_relate
from .error import RelateError
from ...core.ngs import MAX_FLAG, SAM_DELIM
from ...core.rel import MATCH, DELET, NOCOV, IRREC


class SamFlag(object):
    """ Represents the set of 12 boolean flags for a SAM record. """

    # Define __slots__ to improve speed and memory performance.
    __slots__ = ["flag", "paired", "proper", "rev", "read1", "read2"]

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
            raise RelateError("SAM flag is too large")
        self.flag = flag
        self.paired = bool(flag & 1)
        self.proper = bool(flag & 2)
        self.rev = bool(flag & 16)
        self.read1 = bool(flag & 64)
        self.read2 = bool(flag & 128)

    def __repr__(self):
        return f"{type(self).__name__}({self.flag})"


class SamRead(object):
    """ One read in a SAM file. """

    # Define __slots__ to improve speed and memory performance.
    __slots__ = ["name", "flag", "ref", "pos", "mapq", "cigar", "seq", "quals"]

    def __init__(self, line: str):
        fields = line.rstrip().split(SAM_DELIM)
        # Read name
        try:
            self.name = fields[0]
            if not self.name:
                raise IndexError
        except IndexError:
            raise RelateError("Failed to parse read name") from None
        # Bitwise flag
        try:
            self.flag = SamFlag(int(fields[1]))
        except (IndexError, ValueError):
            raise RelateError("Failed to parse SAM flag") from None
        # Reference name
        try:
            self.ref = fields[2]
            if not self.ref:
                raise IndexError
        except IndexError:
            raise RelateError("Failed to parse reference name") from None
        # Mapping position
        try:
            self.pos = int(fields[3])
        except (IndexError, ValueError):
            raise RelateError("Failed to parse mapping position") from None
        if self.pos < 1:
            raise RelateError("Mapping position is 0")
        # Mapping quality
        try:
            self.mapq = int(fields[4])
        except (IndexError, ValueError):
            raise RelateError("Failed to parse mapping quality") from None
        # CIGAR string
        try:
            self.cigar = fields[5]
            if not self.cigar:
                raise IndexError
        except IndexError:
            raise RelateError("Failed to parse CIGAR string") from None
        # Read sequence
        try:
            self.seq = fields[9]
            if not self.seq:
                raise IndexError
        except IndexError:
            raise RelateError("Failed to parse read sequence") from None
        # Read quality
        try:
            self.quals = fields[10]
            if not self.quals:
                raise IndexError
        except IndexError:
            raise RelateError("Failed to parse read quality") from None
        if len(self.seq) != len(self.quals):
            raise RelateError(
                "Read sequence and quality strings differ in length"
            )

    def __str__(self):
        attrs = {attr: getattr(self, attr) for attr in self.__slots__[1:]}
        return f"Read {repr(self.name)} {attrs}"


def _add_indel(pods: list[IndelPod],
               pod_type: type[IndelPod],
               opposite: int,
               lateral3: int):
    if not pods or not isinstance(pods[-1], pod_type):
        pods.append(pod_type())
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
    min_qual: str
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
    assert read.pos >= 1
    if read.pos > ref_length:
        raise RelateError("Mapping position exceeds length of reference")
    # Position in the reference: can be 0- or 1-indexed (initially 0).
    init_ref_pos = read.pos - 1
    ref_pos = init_ref_pos
    # Position in the read: can be 0- or 1-indexed (initially 0).
    read_pos = 0
    # Ends of the read, in the read coordinates (1-indexed).
    read_end5 = 1
    read_end3 = read_length
    # Lists of segment ends, in reference coordinates (1-indexed)
    seg_ends5 = list()
    seg_ends3 = list()
    # Record the relationship for each mutated position.
    rels: dict[int, int] = dict()
    # Record the positions of deletions and insertions.
    pods: list[IndelPod] = list()

    # Count reference and read bases consumed by the CIGAR string.
    ref_bases = 0
    read_bases = 0
    for cigar_op, op_length in parse_cigar(read.cigar):
        if (cigar_op == CIG_ALIGN
                or cigar_op == CIG_MATCH
                or cigar_op == CIG_SUBST):
            # These operations consume the reference and read.
            ref_bases += op_length
            read_bases += op_length
        elif cigar_op == CIG_DELET or cigar_op == CIG_INTRN:
            # These operations consume only the reference.
            ref_bases += op_length
        elif cigar_op == CIG_INSRT or cigar_op == CIG_SCLIP:
            # These operations consume only the read.
            read_bases += op_length
        else:
            # It should be impossible for parse_cigar() to yield any
            # operation besides those already handled.
            assert False, "Unsupported CIGAR operation"
    # The number of reference bases consumed must be > 0 but not extend
    # out of the reference sequence.
    if ref_bases == 0:
        raise RelateError("CIGAR operations consumed 0 bases in the reference")
    if ref_bases + init_ref_pos > ref_length:
        raise RelateError("CIGAR operations extended out of the reference")
    # The number of read bases consumed must equal the read length.
    if read_bases != len(read.seq):
        raise RelateError("CIGAR operations consumed a number of read bases "
                          "different from the read length")

    # Read the CIGAR string one operation at a time.
    for cigar_op, op_length in parse_cigar(read.cigar):
        assert ref_pos - init_ref_pos <= ref_bases
        assert read_pos <= len(read.seq)
        # Act based on the CIGAR operation and its length.
        if (cigar_op == CIG_ALIGN
                or cigar_op == CIG_MATCH
                or cigar_op == CIG_SUBST):
            # There are only matches or substitutions over the entire
            # CIGAR operation.
            assert ref_pos + op_length <= ref_length, ("CIGAR operations "
                                                       "extended out of the "
                                                       "reference")
            assert read_pos + op_length <= len(read.seq), ("CIGAR operations "
                                                           "extended out of "
                                                           "the read")
            for _ in range(op_length):
                rel = encode_relate(ref_seq[ref_pos],
                                    read.seq[read_pos],
                                    read.quals[read_pos],
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
            assert ref_pos + op_length <= ref_length, ("CIGAR operations "
                                                       "extended out of the "
                                                       "reference")
            read_pos += 1  # 1-indexed now until explicitly reset
            for _ in range(op_length):
                if ref_pos == init_ref_pos:
                    raise RelateError("A deletion was the first relationship")
                ref_pos += 1  # 1-indexed now until this iteration ends
                if ref_pos - init_ref_pos >= ref_bases:
                    raise RelateError("A deletion was the last relationship")
                rels[ref_pos] = DELET
                _add_indel(pods, DeletionPod, ref_pos, read_pos)
            read_pos -= 1  # 0-indexed now
            assert read_pos > 0
            assert read_pos < len(read.seq)
        elif cigar_op == CIG_INTRN:
            # The portion of the reference sequence corresponding to the
            # CIGAR operation is spliced out of the read. Treat each side of
            # the splice junction as its own segment.
            assert ref_pos + op_length <= ref_length, ("CIGAR operations "
                                                       "extended out of the "
                                                       "reference")
            if ref_pos == init_ref_pos:
                raise RelateError("An intron was the first relationship")
            # The start of an intron delimits the end of segment.
            # The segment ends at the position before the intron starts.
            seg_ends3.append(ref_pos)
            ref_pos += op_length
            if ref_pos - init_ref_pos >= ref_bases:
                raise RelateError("An intron was the last relationship")
            # The end of an intron delimits the start of a segment.
            # The segment starts at the position after the intron ends.
            seg_ends5.append(ref_pos + 1)
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
            assert read_pos + op_length <= len(read.seq), ("CIGAR operations "
                                                           "extended out of "
                                                           "the read")
            if ref_pos == init_ref_pos:
                raise RelateError("An insertion was the first relationship")
            assert read_pos > 0
            if ref_pos - init_ref_pos >= ref_bases:
                raise RelateError("An insertion was the last relationship")
            ref_pos += 1  # 1-indexed now until explicitly reset
            for _ in range(op_length):
                assert read_pos < len(read.seq)
                read_pos += 1  # 1-indexed now until this iteration ends
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
            assert next_read_pos <= len(read.seq), ("CIGAR operations "
                                                    "extended out of "
                                                    "the read")
            if read_pos == 0:
                # This is the soft clip from the 5' end of the read.
                assert ref_pos == init_ref_pos
                assert read_end5 == 1
                read_end5 += op_length
                assert read_end5 <= len(read.seq) + 1
            else:
                # This is the soft clip from the 3' end of the read.
                if (ref_pos - init_ref_pos < ref_bases
                        or next_read_pos < read_length):
                    raise RelateError("A soft clip occurred in the middle")
                assert read_end3 == read_length == read_pos + op_length
                read_end3 = read_pos
            read_pos = next_read_pos
        else:
            # It should be impossible for parse_cigar() to yield any
            # operation besides those already handled.
            assert False, "Unsupported CIGAR operation"
    # The second pass through the CIGAR string should have resulted in
    # the same counts as did the first pass.
    assert ref_pos - init_ref_pos == ref_bases, ("ref_pos - init_ref_pos "
                                                 "≠ ref_bases")
    assert read_pos == read_bases, "read_pos ≠ read_bases"

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
    if ambindel:
        find_ambindels(rels=rels,
                       pods=pods,
                       insert3=insert3,
                       ref_seq=ref_seq,
                       read_seq=read.seq,
                       read_qual=read.quals,
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

    seg_ends5.insert(0, ref_end5)
    seg_ends3.append(ref_end3)

    assert len(seg_ends5) == len(seg_ends3), ("The number of 5' segment ends "
                                              "must match the number of "
                                              "3' segment ends")

    # Convert rels to a non-default dict and select only the positions
    # contained within a segment.
    rels = {pos: rel for pos, rel in rels.items()
            if any(end5 <= pos <= end3 for end5, end3 in zip(seg_ends5,
                                                             seg_ends3))}
    return seg_ends5, seg_ends3, rels


def _validate_read(read: SamRead,
                   ref: str,
                   min_mapq: int,
                   paired: bool,
                   proper: bool):
    if read.ref != ref:
        raise RelateError("Reference name does not match name of SAM file")
    if read.mapq < min_mapq:
        raise RelateError("Mapping quality is insufficient")
    if read.flag.paired != paired:
        expect = "paired" if paired else "single"
        marked = "paired" if read.flag.paired else "single"
        raise RelateError(f"Lines indicate read should be {expect}-end, "
                          f"but it is marked as {marked}-end")
    if read.flag.proper != proper:
        expect = "properly" if proper else "improperly"
        marked = "properly" if read.flag.proper else "improperly"
        raise RelateError(f"Lines indicate read should be {expect} paired, "
                          f"but it is marked as {marked} paired")


def _validate_pair(read1: SamRead, read2: SamRead):
    """ Ensure that reads 1 and 2 are compatible mates. """
    if read1.name != read2.name:
        raise RelateError("Mates 1 and 2 have different names")
    if not read1.flag.read1 or read1.flag.read2:
        raise RelateError("Mate 1 is not marked as READ1")
    if not read2.flag.read2 or read2.flag.read1:
        raise RelateError("Mate 2 is not marked as READ2")
    if read1.flag.rev == read2.flag.rev:
        raise RelateError("Mates 1 and 2 aligned in the same orientation")


def trim_segs_start(seg5s: list[int],
                    seg3s: list[int],
                    min_start: int) -> tuple[list[int], list[int]]:
    """
    For each segment (start, end), the segment's start is replaced by
    max(start, min_start).
    """
    trimmed = [(max(start, min_start), end)
               for start, end in zip(seg5s, seg3s)]
    new_seg5s, new_seg3s = zip(*trimmed)
    return list(new_seg5s), list(new_seg3s)


def trim_segs_end(seg5s: list[int],
                  seg3s: list[int],
                  max_end: int) -> tuple[list[int], list[int]]:
    """
    For each segment (start, end), the segment's end is replaced by
    min(start, max_end).
    """
    trimmed = [(start, min(end, max_end))
               for start, end in zip(seg5s, seg3s)]
    new_seg5s, new_seg3s = zip(*trimmed)
    return list(new_seg5s), list(new_seg3s)


def _merge_rels(end5sf: list[int],
                end3sf: list[int],
                relsf: dict[int, int],
                end5sr: list[int],
                end3sr: list[int],
                relsr: dict[int, int]):
    merged_rels = dict()
    for pos in relsf | relsr:
        relf = relsf.get(pos, MATCH if any(end5f <= pos <= end3f
                                           for end5f, end3f
                                           in zip(end5sf, end3sf, strict=True)) else NOCOV)
        relr = relsr.get(pos, MATCH if any(end5r <= pos <= end3r
                                           for end5r, end3r
                                           in zip(end5sr, end3sr, strict=True)) else NOCOV)
        rel = relf & relr
        if rel != MATCH:
            if rel == NOCOV:
                raise RelateError(f"Cannot merge non-covered position {pos}")
            merged_rels[pos] = rel
    return merged_rels


def merge_mates(end5sf: list[int],
                end3sf: list[int],
                relsf: dict[int, int],
                end5sr: list[int],
                end3sr: list[int],
                relsr: dict[int, int],
                overhangs: bool):
    if not overhangs:
        # The 5' end of the reverse mate cannot extend past the 5' end
        # of the forward mate.
        min_end5r = min(end5sr)
        min_end5f = min(end5sf)
        if min_end5r < min_end5f:
            end5sr, end3sr = trim_segs_start(end5sr, end3sr, min_end5f)
            relsr = {pos: rel for pos, rel in relsr.items() if pos >= min_end5f}
        # The 3' end of the forward mate cannot extend past the 3' end
        # of the reverse mate.
        max_end3r = max(end3sr)
        max_end3f = max(end3sf)
        if max_end3f > max_end3r:
            end5sf, end3sf = trim_segs_end(end5sf, end3sf, max_end3r)
            relsf = {pos: rel for pos, rel in relsf.items() if pos <= max_end3r}
    rels = _merge_rels(end5sf, end3sf, relsf, end5sr, end3sr, relsr)
    return ([*end5sf, *end5sr], [*end3sf, *end3sr]), rels


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
    # Determine if the read is paired-end and properly paired.
    paired = bool(line2)
    proper = paired and line1 != line2
    # Generate the relationships for read 1.
    read1 = SamRead(line1)
    _validate_read(read1, ref, min_mapq, paired, proper)
    end5s1, end3s1, rels1 = _calc_rels_read(read1,
                                            refseq,
                                            chr(min_qual),
                                            insert3,
                                            ambindel,
                                            clip_end5,
                                            clip_end3)
    if paired:
        # The read is paired-end.
        if proper:
            # The read is properly paired.
            # Generate the relationships for read 2.
            read2 = SamRead(line2)
            _validate_read(read2, ref, min_mapq, paired, proper)
            _validate_pair(read1, read2)
            end5s2, end3s2, rels2 = _calc_rels_read(read2,
                                                    refseq,
                                                    chr(min_qual),
                                                    insert3,
                                                    ambindel,
                                                    clip_end5,
                                                    clip_end3)
            # Determine which read (1 or 2) faces forward and reverse.
            if read2.flag.rev:
                end5sf, end3sf, relsf = end5s1, end3s1, rels1
                end5sr, end3sr, relsr = end5s2, end3s2, rels2
            else:
                end5sf, end3sf, relsf = end5s2, end3s2, rels2
                end5sr, end3sr, relsr = end5s1, end3s1, rels1

            # Merge the relationships and end coordinates.
            ends, rels = merge_mates(
                end5sf, end3sf, relsf, end5sr, end3sr, relsr, overhangs,
            )
        else:
            # The read is improperly paired (paired-end but comprises
            # only one mate), which can occur only if Bowtie2 is run in
            # mixed mode.
            ends = [*end5s1, *end5s1], [*end3s1, *end3s1]
            rels = rels1
    else:
        # The read is single-end.
        ends = end5s1, end3s1
        rels = rels1
    return ends, rels
