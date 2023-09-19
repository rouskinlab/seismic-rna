from __future__ import annotations

from .encode import encode_relate
from .error import RelateNotImplementedError, RelateValueError
from ..core.rel import DELET, INS_5, INS_3, SUB_N
from ..core.seq import DNA


class Indel(object):
    """
    Base class for an Insertion or Deletion (collectively, "indel")
    It is used to find alternative positions for indels by keeping track
    of an indel's current coordinates (as it is moved) and determining
    whether a specific move is valid.

    Parameters
    ----------

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
    rel_del_idx (int):
        The opposite of rel_ins_idx: the 0-indexed position of the indel
        with respect to the sequence with a relative deletion (that is,
        the read if the mutation is denoted a deletion, and the ref if
        an insertion). Because the deleted base does not actually exist
        in the sequence whose coordinates it is based on, rel_del_idx
        does not refer to a specific position in the sequence, rather to
        the two extant positions in the sequence that flank the deleted
        position. It is most convenient for the algorithm to have this
        argument refer to the position 3' of the deleted base and define
        the 5' position as a property.

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
        if curr_rel & swap_rel or (curr_rel & SUB_N and swap_rel & SUB_N):
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

    def sweep(self, muts: bytearray, ref: DNA, read: DNA, qual: str,
              min_qual: str, dels: list[Deletion], inns: list[Insertion],
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
    def _encode_swap(cls, ref_base: str, swap_base: str, read_base: str,
                     read_qual: str, min_qual: str):
        curr_rel = encode_relate(ref_base, read_base, read_qual, min_qual)
        swap_rel = encode_relate(swap_base, read_base, read_qual, min_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self, muts: bytearray, swap_idx: int, relation: int):
        """
        Parameters
        ----------
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

    def _try_swap(self, relvec: bytearray, refseq: DNA, read: DNA, qual: str,
                  min_qual: str, dels: list[Deletion], inns: list[Insertion],
                  from3to5: bool, tunnel: bool) -> bool:
        swap_idx, tunneled_indels = self._peek_out_of_indel(dels, from3to5)
        read_idx = self.del_idx5 if from3to5 else self.del_idx3
        if (1 <= swap_idx < len(refseq) - 1 and 1 <= read_idx < len(read) - 1
                and (tunnel or not self.tunneled)
                and not self._collisions(inns, swap_idx)):
            relation = self._encode_swap(refseq[self.ins_idx],
                                         refseq[swap_idx],
                                         read[read_idx],
                                         qual[read_idx],
                                         min_qual)
            if relation:
                self._swap(relvec, swap_idx, relation)
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
    def _encode_swap(cls, ref_base: str, read_base: str, read_qual: str,
                     swap_base: str, swap_qual: str, min_qual: str):
        curr_rel = encode_relate(ref_base, read_base, read_qual, min_qual)
        swap_rel = encode_relate(ref_base, swap_base, swap_qual, min_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self, muts: bytearray, ref_idx: int,
              swap_idx: int, relation: int):
        """
        Parameters
        ----------
        muts: bytearray:
            Relation vector
        swap_idx: int:
            Index in the read to which the deletion is swapped
        relation: int
            Relationship (match, sub, etc.) between the base located at
            swap_idx and the base in the ref

        """
        # The base at ref_idx moves to swap_idx, so after the swap, the
        # relationship between ref_idx and the read base is relation.
        muts[ref_idx] |= relation
        self._step(swap_idx)
        # Mark the new positions of the insertion.
        self.stamp(muts)

    def _try_swap(self, muts: bytearray, refseq: DNA, read: DNA, qual: str,
                  min_qual: str, dels: list[Deletion], inns: list[Insertion],
                  from3to5: bool, tunnel: bool) -> bool:
        swap_idx, tunneled_indels = self._peek_out_of_indel(inns, from3to5)
        ref_idx = self.del_idx5 if from3to5 else self.del_idx3
        if (1 <= swap_idx < len(read) - 1 and 1 <= ref_idx < len(refseq) - 1
                and (tunnel or not self.tunneled)
                and not self._collisions(dels, swap_idx)):
            relation = self._encode_swap(refseq[ref_idx], read[self.ins_idx],
                                         qual[self.ins_idx], read[swap_idx],
                                         qual[swap_idx], min_qual)
            if relation:
                self._swap(muts, ref_idx, swap_idx, relation)
                for indel in tunneled_indels:
                    indel.step_del_idx(swap_idx)
                return True
        return False


def sweep_indels(muts: bytearray, refseq: DNA, read: DNA, qual: str,
                 min_qual: str, dels: list[Deletion], inns: list[Insertion],
                 from3to5: bool, tunnel: bool):
    """
    For every insertion and deletion,

    Parameters
    ----------
    muts: bytearray
        Mutation vector
    refseq: bytes
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
        indel.sweep(muts, refseq, read, qual, min_qual,
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


def find_ambrels(relvec: bytearray, refseq: DNA, read: DNA, qual: str,
                 min_qual: str, dels: list[Deletion], inns: list[Insertion]):
    """
    Find and label all positions in the vector that are ambiguous due to
    insertions and deletions.

    Parameters
    ----------
    relvec: bytearray
        Mutation vector
    refseq: DNA
        Reference sequence
    read: DNA
        Sequence of the read
    qual: str
        Phred quality scores of the read, encoded as ASCII characters
    min_qual: str
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
        sweep_indels(relvec, refseq, read, qual, min_qual,
                     dels, inns, from3to5, tunnel=True)
        if any(d.tunneled for d in dels) or any(i.tunneled for i in inns):
            # If any indel tunneled,
            sweep_indels(relvec, refseq, read, qual, min_qual,
                         dels, inns, from3to5, tunnel=False)

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
