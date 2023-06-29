from __future__ import annotations

from ..core.rel import (DELET, INS_5, INS_3, SUB_N,
                        CIG_ALIGN, CIG_MATCH, CIG_SUBST,
                        CIG_DELET, CIG_INSRT, CIG_SCLIP,
                        parse_cigar, encode_match, encode_relate)


class RelateError(Exception):
    """ Any error that occurs during relating. """


class RelateValueError(RelateError, ValueError):
    """ Any ValueError that occurs during relating. """


class RelateNotImplementedError(RelateError, NotImplementedError):
    """ Any NotImplementedError that occurs during relating. """


class Indel(object):
    """
    Base class for an Insertion or Deletion (collectively, "indel")
    It is used to find alternative positions for indels by keeping track
    of an indel's current coordinates (as it is moved) and determining
    whether a specific move is valid.
    Arguments
    ---------
    rel_ins_idx: int
        The 0-indexed position of the indel with respect to the sequence
        (ref or read) with the relative insertion. This position points
        to one specific base. If the mutation is labeled an insertion,
        then the read is the sequence with the relative insertion (since
        it has a base that is not in the reference), and rel_ins_idx is
        the 0-based index of the inserted base in the coordinates of the
        read sequence. If the mutation is labeled a deletion, then the
        reference is the sequence with the relative insertion (since it
        has a base that is not in the read), and rel_ins_idx is the
        0-based index of the deleted base in the coordinates of the
        reference sequence.
    rel_del_idx (int):  The opposite of rel_ins_idx: the 0-indexed position
                        of the indel with respect to the sequence with a
                        relative deletion (that is, the read if the mutation
                        is denoted a deletion, and the ref if an insertion).
                        Because the deleted base does not actually exist
                        within the sequence whose coordinates it is based on,
                        rel_del_idx does not refer to a specific position in
                        the sequence, but instead to the two extant positions
                        in the sequence that flank the deleted position.
                        It is most convenient for the algorithm to have this
                        argument refer to the position 3' of the deleted base
                        and define the 5' position as a property.
    """

    # Define __slots__ to improve speed and memory performance.
    __slots__ = ["_ins_idx", "_ins_init", "_del_idx", "_del_init", "_tunneled"]

    # Minimum distance between an insertion and a deletion
    MIN_INDEL_DIST = 2

    def __init__(self, rel_ins_idx: int, rel_del_idx: int) -> None:
        self._ins_idx = rel_ins_idx
        self._ins_init = rel_ins_idx
        self._del_idx = rel_del_idx
        self._del_init = rel_del_idx
        self._tunneled = False

    @property
    def ins_idx(self):
        return self._ins_idx

    @property
    def del_idx5(self):
        return self._del_idx - 1

    @property
    def del_idx3(self):
        return self._del_idx

    @property
    def tunneled(self):
        return self._tunneled

    @property
    def rank(self) -> int:
        raise RelateNotImplementedError

    def reset(self):
        """ Reset the indel to its initial position, and erase its
        history of tunneling. """
        self._ins_idx = self._ins_init
        self._del_idx = self._del_init
        self._tunneled = False

    @staticmethod
    def _get_indel_by_idx(indels: list[Indel], idx: int):
        for indel in indels:
            if indel.ins_idx == idx:
                return indel

    def _peek_out_of_indel(self, indels: list[Indel], from3to5: bool):
        inc = -1 if from3to5 else 1
        idx = self.ins_idx + inc
        tunneled_indels: list[Indel] = list()
        while indel := (self._get_indel_by_idx(indels, idx)):
            idx += inc
            tunneled_indels.append(indel)
        self._tunneled = bool(tunneled_indels)
        return idx, tunneled_indels

    def _collision(self, other: Indel, swap_idx: int):
        return self.MIN_INDEL_DIST > (min(abs(swap_idx - other.del_idx5),
                                          abs(swap_idx - other.del_idx3)))

    def _collisions(self, indels: list[Indel], swap_idx: int):
        return any(self._collision(indel, swap_idx) for indel in indels)

    def step_del_idx(self, swap_idx: int):
        # Move the indel's position (self._ins_idx) to swap_idx.
        # Move self._del_idx one step in the same direction.
        if swap_idx == self.ins_idx:
            raise RelateValueError(f"swap ({swap_idx}) = ins ({self.ins_idx})")
        self._del_idx += 1 if swap_idx > self.ins_idx else -1

    def _step(self, swap_idx: int):
        self.step_del_idx(swap_idx)
        self._ins_idx = swap_idx

    @staticmethod
    def _consistent_rels(curr_rel: int, swap_rel: int):
        if curr_rel & swap_rel or (curr_rel & SUB_N
                                   and swap_rel & SUB_N):
            # Relationship between reference and read base (read_code) and
            # relationship between reference and swap base (swap_code)
            # are consistent, meaning either
            # - both match the reference
            # - one matches and the other potentially matches (i.e. low qual)
            # - one is a substitution and the other could be a substitution
            # - both are substitutions (for each code, code & SUB_N == code)
            return curr_rel
        # Otherwise, i.e.g if one base matches and the other is a substitution,
        # then the relationships are not consistent.
        return 0

    def _encode_swap(self, *args, **kwargs) -> bool:
        raise RelateNotImplementedError

    def _try_swap(self, *args, **kwargs) -> bool:
        raise RelateNotImplementedError

    def sweep(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
              min_qual: int, dels: list[Deletion], inns: list[Insertion],
              from3to5: bool, tunnel: bool):
        # Move the indel as far as possible in either the 5' or 3' direction.
        while self._try_swap(muts, ref, read, qual, min_qual, dels, inns,
                             from3to5, tunnel):
            # All actions happen in _try_swap, so loop body is empty.
            pass


class Deletion(Indel):
    @property
    def rank(self):
        return self._ins_idx

    @classmethod
    def _encode_swap(cls, ref_base: int, swap_base: int, read_base: int,
                     read_qual: int, min_qual: int):
        curr_rel = encode_relate(ref_base, read_base, read_qual, min_qual)
        swap_rel = encode_relate(swap_base, read_base, read_qual, min_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self, muts: bytearray, swap_idx: int, relation: int):
        """
        Arguments
        ---------
        muts: bytearray
            Mutation vector
        swap_idx: int
            Index in the reference to which the deletion moves during
            this swap
        relation: int
            Relationship (match, sub, etc.) between the base located at
            swap_idx and the base in the read
        """
        # The base at swap_idx moves to self.ref_idx, so after the swap, the
        # relationship between self.ref_idx and the read base will be swap_code.
        muts[self.ins_idx] |= relation
        # The base at self.ref_idx is marked as a deletion (by definition), so
        # mark the position it moves to (swap_idx) as a deletion too.
        muts[swap_idx] |= DELET
        self._step(swap_idx)

    def _try_swap(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                  min_qual: int, dels: list[Deletion], inns: list[Insertion],
                  from3to5: bool, tunnel: bool) -> bool:
        swap_idx, tunneled_indels = self._peek_out_of_indel(dels, from3to5)
        read_idx = self.del_idx5 if from3to5 else self.del_idx3
        if (1 <= swap_idx < len(ref) - 1 and 1 <= read_idx < len(read) - 1
                and (tunnel or not self.tunneled)
                and not self._collisions(inns, swap_idx)):
            relation = self._encode_swap(ref[self.ins_idx], ref[swap_idx],
                                         read[read_idx], qual[read_idx],
                                         min_qual)
            if relation:
                self._swap(muts, swap_idx, relation)
                for indel in tunneled_indels:
                    indel.step_del_idx(swap_idx)
                return True
        return False


class Insertion(Indel):
    @property
    def rank(self):
        return self._del_idx

    def stamp(self, muts: bytearray):
        """ Stamp the relation vector with a 5' and a 3' insertion. """
        if 0 <= self.del_idx5 < len(muts):
            muts[self.del_idx5] |= INS_5
        if 0 <= self.del_idx3 < len(muts):
            muts[self.del_idx3] |= INS_3

    @classmethod
    def _encode_swap(cls, ref_base: int, read_base: int, read_qual: int,
                     swap_base: int, swap_qual: int, min_qual: int):
        curr_rel = encode_relate(ref_base, read_base, read_qual, min_qual)
        swap_rel = encode_relate(ref_base, swap_base, swap_qual, min_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self, muts: bytearray, ref_idx: int,
              swap_idx: int, relation: int):
        """
        Arguments
        mut_vectors (bytearray): mutation vector
        swap_idx (int): the index in the read to which the deletion moves
                        during this swap
        swap_code (int): the relationship (match, sub, etc.) between the
                         base located at swap_idx and the base in the ref
        """
        # The base at ins_idx moves to swap_idx, so after the swap, the
        # relationship between ref_idx and the read base will be relation.
        muts[ref_idx] |= relation
        self._step(swap_idx)
        # Mark the new positions of the insertion.
        self.stamp(muts)

    def _try_swap(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                  min_qual: int, dels: list[Deletion], inns: list[Insertion],
                  from3to5: bool, tunnel: bool) -> bool:
        swap_idx, tunneled_indels = self._peek_out_of_indel(inns, from3to5)
        ref_idx = self.del_idx5 if from3to5 else self.del_idx3
        if (1 <= swap_idx < len(read) - 1 and 1 <= ref_idx < len(ref) - 1
                and (tunnel or not self.tunneled)
                and not self._collisions(dels, swap_idx)):
            relation = self._encode_swap(ref[ref_idx], read[self.ins_idx],
                                         qual[self.ins_idx], read[swap_idx],
                                         qual[swap_idx], min_qual)
            if relation:
                self._swap(muts, ref_idx, swap_idx, relation)
                for indel in tunneled_indels:
                    indel.step_del_idx(swap_idx)
                return True
        return False


def sweep_indels(muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                 min_qual: int, dels: list[Deletion], inns: list[Insertion],
                 from3to5: bool, tunnel: bool):
    """
    For every insertion and deletion,
    Parameters
    ----------
    muts: bytearray
        Mutation vector
    ref: bytes
        Reference sequence
    read: bytes
        Sequence of the read
    qual: bytes
        Phred quality scores of the read, encoded as ASCII characters
    min_qual: int
        The minimum Phred quality score needed to consider a base call
        informative: integer value of the ASCII character
    dels: list[Deletion]
        List of deletions identified by `vectorize_read`
    inns: list[Insertion]
        List of insertions identified by `vectorize_read`
    from3to5: bool
        Whether to move indels in the 3' -> 5' direction (True) or the
        5' -> 3' direction (False)
    tunnel: bool
        Whether to allow tunneling
    """
    # Collect all indels into one list.
    indels: list[Indel] = list()
    indels.extend(dels)
    indels.extend(inns)
    # Reset each indel to its initial state. This operation does nothing
    # the first time sweep_indels is called because all indels start in
    # their initial state (by definition). But the indels may move when
    # this function runs, so resetting is necessary at the beginning of
    # the second and subsequent calls to sweep_indels to ensure that the
    # algorithm starts from the initial state every time.
    for indel in indels:
        indel.reset()
    # Sort the indels by their rank, which is
    sort_rev = from3to5 != tunnel
    indels.sort(key=lambda idl: idl.rank, reverse=sort_rev)
    while indels:
        indel = indels.pop()
        indel.sweep(muts, ref, read, qual, min_qual,
                    dels, inns, from3to5, tunnel)
        i = len(indels)
        if sort_rev:
            while i > 0 and indel.rank > indels[i - 1].rank:
                i -= 1
        else:
            while i > 0 and indel.rank < indels[i - 1].rank:
                i -= 1
        if i < len(indels):
            indels.insert(i, indel)


def find_ambrels(muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                 min_qual: int, dels: list[Deletion], inns: list[Insertion]):
    """
    Find and label all positions in the vector that are ambiguous due to
    insertions and deletions.
    Parameters
    ----------
    muts: bytearray
        Mutation vector
    ref: bytes
        Reference sequence
    read: bytes
        Sequence of the read
    qual: bytes
        Phred quality scores of the read, encoded as ASCII characters
    min_qual: int
        The minimum Phred quality score needed to consider a base call
        informative: integer value of the ASCII character
    dels: list[Deletion]
        List of deletions identified by `vectorize_read`
    inns: list[Insertion]
        List of insertions identified by `vectorize_read`
    """
    # Each indel might be able to be moved in the 5' -> 3' direction
    # (from3to5 is False) or 3' -> 5' direction (from3to5 is True).
    # Test both directions.
    for from3to5 in (False, True):
        # For each indel, try to move it as far as it can go in the
        # direction indicated by from3to5. Allow tunneling so that any
        # runs of consecutive insertions or consecutive deletions can
        # effectively move together.
        sweep_indels(muts, ref, read, qual, min_qual,
                     dels, inns, from3to5, tunnel=True)
        if any(d.tunneled for d in dels) or any(i.tunneled for i in inns):
            # If any indel tunneled,
            sweep_indels(muts, ref, read, qual, min_qual,
                         dels, inns, from3to5, tunnel=False)


def op_consumes_ref(op: bytes) -> bool:
    """ Return whether the CIGAR operation consumes the reference. """
    return op != CIG_INSRT and op != CIG_SCLIP


def op_consumes_read(op: bytes) -> bool:
    """ Return whether the CIGAR operation consumes the read. """
    return op != CIG_DELET


class SamFlag(object):
    """ Represents the set of 12 boolean flags for a SAM record. """

    # Define __slots__ to improve speed and memory performance.
    __slots__ = "flag", "paired", "rev", "first", "second"

    # Maximum value of a valid SAM flag representation, corresponding
    # to all 12 flags set to 1: 111111111111 (binary) = 4095 (decimal)
    MAX_FLAG = 4095

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
        if not 0 <= flag <= self.MAX_FLAG:
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

    def __init__(self, line: bytes):
        fields = line.rstrip().split(b"\t")
        if len(fields) < self.MIN_FIELDS:
            raise RelateValueError(f"Invalid SAM line:\n{line}")
        self.qname = fields[0]
        self.flag = SamFlag(int(fields[1]))
        self.rname = fields[2].decode()
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.seq = fields[9]
        self.qual = fields[10]
        if len(self.seq) != len(self.qual):
            raise RelateValueError(f"Lengths of seq ({len(self.seq)}) and qual "
                                   f"string {len(self.qual)} did not match.")

    def __str__(self):
        attrs = {attr: self.__getattribute__(attr) for attr in self.__slots__[1:]}
        return f"Read '{self.qname.decode()}' {attrs}"


def relate_read(read: SamRead,
                muts: bytearray,
                refseq: bytes,
                length: int,
                min_qual: int,
                ambrel: bool):
    """
    Generate a relation vector of a read aligned to a reference.

    Parameters
    ----------
    read: SamRead
        Read from SAM file to be related
    muts: bytearray
        Relation vector (initially blank) into which to write bytes;
        muts and refseq must have the same length.
    refseq: bytes
        Reference sequence; refseq and muts must have the same length.
    length: int (≥ 1)
        Length of the reference; must equal len(refseq) and len(muts)
    min_qual: int
        ASCII encoding of the minimum Phred score to accept a base call
    ambrel: bool
        Whether to find and label all ambiguous insertions and deletions
    """
    if len(muts) != length:
        raise ValueError(
            f"Expected muts to have length {length}, but got {len(muts)}")
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
                muts[ref_idx] &= encode_match(read.seq[read_idx],
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
                muts[ref_idx] &= encode_relate(refseq[ref_idx],
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
                muts[ref_idx] &= DELET
                ref_idx += 1
        elif cigar_op == CIG_INSRT:
            # The read contains an insertion of one or more bases
            # that are not present in the reference sequence.
            # Create one Insertion object for each base in the read
            # sequence that is not present in the reference. Every
            # mutation needs to be assigned a coordinate in the
            # region in order to appear at that coordinate in the
            # mutation vector. But each inserted base, being absent
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
            f"CIGAR string '{read.cigar.decode()}' consumed {read_idx} "
            f"bases from read, but read is {len(read.seq)} bases long.")
    # Add insertions to muts.
    for ins in inns:
        ins.stamp(muts)
    # Find and label all relationships that are ambiguous due to indels.
    if ambrel and (dels or inns):
        find_ambrels(muts, refseq, read.seq, read.qual, min_qual, dels, inns)


def relate_line(line: bytes, muts: bytearray, refseq: bytes,
                length: int, ref: str, qmin: int, ambrel: bool):
    read = SamRead(line)
    if read.rname != ref:
        raise RelateValueError(f"Read '{read.qname.decode()}' had reference "
                               f"'{read.rname}' (≠ '{ref}')")
    relate_read(read, muts, refseq, length, qmin, ambrel)


def relate_pair(line1: bytes, line2: bytes, muts: bytearray, refseq: bytes,
                length: int, ref: str, qmin: int, ambrel: bool):
    # Parse lines 1 and 2 into SAM reads.
    read1 = SamRead(line1)
    read2 = SamRead(line2)
    # Ensure that reads 1 and 2 are compatible mates.
    if not read1.flag.paired:
        raise RelateValueError(f"Read 1 ({read1.qname.decode()}) was not "
                               f"paired, but read 2 ('{read2.qname.decode()}') "
                               f"was given")
    if not read2.flag.paired:
        raise RelateValueError(f"Read 2 ({read2.qname.decode()}) was not "
                               f"paired, but read 1 ({read1.qname.decode()}) "
                               f"was given")
    if read1.qname != read2.qname:
        raise RelateValueError(f"Reads 1 ({read1.qname.decode()}) and 2 "
                               f"({read2.qname.decode()}) had different names")
    if read1.rname != read2.rname:
        raise RelateValueError(f"Read '{read1.qname.decode()}' had "
                               "different references for mates 1 "
                               f"('{read1.rname}') and 2 "
                               f"('{read2.rname}')")
    if read1.rname != ref:
        raise RelateValueError(f"Read '{read1.qname.decode()}' had reference "
                               f"'{read1.rname}' (≠ '{ref}')")
    if not (read1.flag.first and read2.flag.second):
        raise RelateValueError(f"Read '{read1.qname.decode()}' had mate 1 "
                               f"labeled {2 - read1.flag.first} and mate 2 "
                               f"labeled {1 + read2.flag.second}")
    if read1.flag.rev == read2.flag.rev:
        raise RelateValueError(f"Read '{read1.qname.decode()}' had "
                               "mates 1 and 2 facing the same way")
    # Vectorize read 1.
    relate_read(read1, muts, refseq, length, qmin, ambrel)
    # Vectorize read 2.
    relate_read(read2, muts, refseq, length, qmin, ambrel)
