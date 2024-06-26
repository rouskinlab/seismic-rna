from __future__ import annotations

from abc import ABC, abstractmethod

from .encode import encode_relate
from .error import RelateValueError
from ...core.rel import DELET, INS_5, INS_3, SUB_N
from ...core.seq import DNA

# Minimum distance between an insertion and a deletion
MIN_INDEL_DIST = 2


class Indel(ABC):
    """
    Base class for an Insertion or Deletion (collectively, "indel")
    It is used to find alternative positions for indels by keeping track
    of an indel's current coordinates (as it is moved) and determining
    whether a specific move is valid.

    Parameters
    ----------

    rel_ins_pos: int
        The 0-indexed position of the indel with respect to the sequence
        (ref or read) with the relative insertion. This position points
        to one specific base. If the mutation is labeled an insertion,
        then the read is the sequence with the relative insertion (since
        it has a base that is not in the reference), and rel_ins_pos is
        the 0-based index of the inserted base in the coordinates of the
        read sequence. If the mutation is labeled a deletion, then the
        reference is the sequence with the relative insertion (since it
        has a base that is not in the read), and rel_ins_pos is the
        0-based index of the deleted base in the coordinates of the
        reference sequence.
    rel_del_pos (int):
        The opposite of rel_ins_pos: the 0-indexed position of the indel
        with respect to the sequence with a relative deletion (that is,
        the read if the mutation is denoted a deletion, and the ref if
        an insertion). Because the deleted base does not actually exist
        in the sequence whose coordinates it is based on, rel_del_pos
        does not refer to a specific position in the sequence, rather to
        the two extant positions in the sequence that flank the deleted
        position. It is most convenient for the algorithm to have this
        argument refer to the position 3' of the deleted base and define
        the 5' position as a property.

    """

    __slots__ = "_ins_pos", "_ins_init", "_del_pos", "_del_init", "_tunneled"

    def __init__(self, rel_ins_pos: int, rel_del_pos: int):
        self._ins_pos = rel_ins_pos
        self._ins_init = rel_ins_pos
        self._del_pos = rel_del_pos
        self._del_init = rel_del_pos
        self._tunneled = False

    @property
    def ins_pos(self):
        return self._ins_pos

    @property
    def del_pos5(self):
        return self._del_pos - 1

    @property
    def del_pos3(self):
        return self._del_pos

    @property
    def tunneled(self):
        return self._tunneled

    @property
    @abstractmethod
    def rank(self) -> int:
        """ Rank of the indel. """

    def reset(self):
        """ Reset the indel to its initial position, and erase its
        history of tunneling. """
        self._ins_pos = self._ins_init
        self._del_pos = self._del_init
        self._tunneled = False

    @staticmethod
    def _get_indel_by_pos(indels: list[Indel], ins_pos: int):
        for indel in indels:
            if indel.ins_pos == ins_pos:
                return indel

    def _peek_out_of_indel(self, indels: list[Indel], from3to5: bool):
        increment = -1 if from3to5 else 1
        ins_pos = self.ins_pos + increment
        tunneled_indels: list[Indel] = list()
        while indel := (self._get_indel_by_pos(indels, ins_pos)):
            ins_pos += increment
            tunneled_indels.append(indel)
        self._tunneled = bool(tunneled_indels)
        return ins_pos, tunneled_indels

    @staticmethod
    def _collision(other: Indel, swap_pos: int):
        return MIN_INDEL_DIST > (min(abs(swap_pos - other.del_pos5),
                                     abs(swap_pos - other.del_pos3)))

    @classmethod
    def _collisions(cls, indels: list[Indel], swap_pos: int):
        return any(cls._collision(indel, swap_pos) for indel in indels)

    def step_del_pos(self, swap_pos: int):
        # Move the indel's position (self._ins_pos) to swap_pos.
        # Move self._del_pos one step in the same direction.
        if swap_pos == self.ins_pos:
            raise RelateValueError(f"swap ({swap_pos}) = ins ({self.ins_pos})")
        self._del_pos += 1 if swap_pos > self.ins_pos else -1

    def _step(self, swap_pos: int):
        self.step_del_pos(swap_pos)
        self._ins_pos = swap_pos

    @staticmethod
    def _consistent_rels(curr_rel: int, swap_rel: int):
        if curr_rel & swap_rel or (curr_rel & SUB_N and swap_rel & SUB_N):
            # Relationship of the reference and read base (curr_rel) and
            # relationship of the reference and swap base (swap_rel) are
            # consistent, meaning one of the following is true:
            # - both match the reference
            # - one matches and the other maybe matches (i.e. low qual)
            # - one is a substitution and the other could also be
            # - both are substitutions
            return curr_rel
        # Otherwise, e.g. if one base matches and the other is a
        # substitution, then the relationships are not consistent.
        return 0

    @classmethod
    @abstractmethod
    def _encode_swap(cls, *args, **kwargs) -> bool:
        """ Encode a swap. """

    @abstractmethod
    def _try_swap(self, *args, **kwargs) -> bool:
        """ Perform a swap if possible. """

    @abstractmethod
    def sweep(self,
              muts: dict[int, int],
              end5_ref: int,
              end3_ref: int,
              end5_read: int,
              end3_read: int,
              ref: DNA,
              read: DNA,
              qual: str,
              min_qual: str,
              dels: list[Deletion],
              inns: list[Insertion],
              from3to5: bool,
              tunnel: bool):
        """ Move the indel as far as possible in one direction. """


class Deletion(Indel):
    @property
    def rank(self):
        return self._ins_pos

    @classmethod
    def _encode_swap(cls,
                     ref_base: str,
                     swap_base: str,
                     read_base: str,
                     read_qual: str,
                     min_qual: str):
        curr_rel = encode_relate(ref_base, read_base, read_qual, min_qual)
        swap_rel = encode_relate(swap_base, read_base, read_qual, min_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self, muts: dict[int, int], swap_pos: int, relation: int):
        """
        Parameters
        ----------
        muts: dict
            Mutation vector
        swap_pos: int
            Position in the reference to which the deletion moves during
            this swap
        relation: int
            Relationship (match, sub, etc.) between the base located at
            swap_pos and the base in the read

        """
        # The base at swap_pos moves to self.ins_pos, so after the swap,
        # the relationship between self.ins_pos and the read base will
        # be relation.
        muts[self.ins_pos] |= relation
        # The base at self.ref_pos is a deletion (by definition), so
        # mark the position it moves to (swap_pos) as a deletion, too.
        muts[swap_pos] |= DELET
        self._step(swap_pos)

    def _try_swap(self,
                  muts: dict[int, int],
                  end5_ref: int,
                  end3_ref: int,
                  refseq: DNA,
                  read: DNA,
                  qual: str,
                  min_qual: str,
                  dels: list[Deletion],
                  inns: list[Insertion],
                  from3to5: bool,
                  tunnel: bool):
        swap_pos, tunneled_indels = self._peek_out_of_indel(dels, from3to5)
        read_pos = self.del_pos5 if from3to5 else self.del_pos3
        if (end5_ref < swap_pos < end3_ref
                and 1 < read_pos < len(read)
                and (tunnel or not self.tunneled)
                and not self._collisions(inns, swap_pos)):
            relation = self._encode_swap(refseq[self.ins_pos - 1],
                                         refseq[swap_pos - 1],
                                         read[read_pos - 1],
                                         qual[read_pos - 1],
                                         min_qual)
            if relation:
                self._swap(muts, swap_pos, relation)
                for indel in tunneled_indels:
                    indel.step_del_pos(swap_pos)
                return True
        return False

    def sweep(self,
              muts: dict[int, int],
              end5_ref: int,
              end3_ref: int,
              end5_read: int,
              end3_read: int,
              ref: DNA,
              read: DNA,
              qual: str,
              min_qual: str,
              dels: list[Deletion],
              inns: list[Insertion],
              from3to5: bool,
              tunnel: bool):
        """ Move the indel as far as possible in one direction. """
        while self._try_swap(muts,
                             end5_ref,
                             end3_ref,
                             ref,
                             read,
                             qual,
                             min_qual,
                             dels,
                             inns,
                             from3to5,
                             tunnel):
            # All actions happen in _try_swap, so loop body is empty.
            pass


class Insertion(Indel):
    @property
    def rank(self):
        return self._del_pos

    def stamp(self, muts: dict[int, int], reflen: int):
        """ Stamp the relation vector with a 5' and a 3' insertion. """
        if 1 <= self.del_pos5 <= reflen:
            muts[self.del_pos5] |= INS_5
        if 1 <= self.del_pos3 <= reflen:
            muts[self.del_pos3] |= INS_3

    @classmethod
    def _encode_swap(cls,
                     ref_base: str,
                     read_base: str,
                     read_qual: str,
                     swap_base: str,
                     swap_qual: str,
                     min_qual: str):
        curr_rel = encode_relate(ref_base, read_base, read_qual, min_qual)
        swap_rel = encode_relate(ref_base, swap_base, swap_qual, min_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self,
              muts: dict[int, int],
              ref_pos: int,
              swap_pos: int,
              relation: int,
              reflen: int):
        """
        Parameters
        ----------
        muts: dict
            Mutations
        swap_pos: int
            Position in the read to which the deletion is swapped
        relation: int
            Relationship (match, sub, etc.) between the base located at
            swap_pos and the base in the ref
        reflen: int
            Length of the reference sequence.

        """
        # The base at ref_pos moves to swap_pos, so after the swap, the
        # relationship between ref_pos and the read base is relation.
        muts[ref_pos] |= relation
        self._step(swap_pos)
        # Mark the new positions of the insertion.
        self.stamp(muts, reflen)

    def _try_swap(self,
                  muts: dict[int, int],
                  end5_read: int,
                  end3_read: int,
                  refseq: DNA,
                  read: DNA,
                  qual: str,
                  min_qual: str,
                  dels: list[Deletion],
                  inns: list[Insertion],
                  from3to5: bool,
                  tunnel: bool):
        swap_pos, tunneled_indels = self._peek_out_of_indel(inns, from3to5)
        ref_pos = self.del_pos5 if from3to5 else self.del_pos3
        if (end5_read < swap_pos < end3_read
                and 1 < ref_pos < len(refseq)
                and (tunnel or not self.tunneled)
                and not self._collisions(dels, swap_pos)):
            relation = self._encode_swap(refseq[ref_pos - 1],
                                         read[self.ins_pos - 1],
                                         qual[self.ins_pos - 1],
                                         read[swap_pos - 1],
                                         qual[swap_pos - 1],
                                         min_qual)
            if relation:
                self._swap(muts, ref_pos, swap_pos, relation, len(refseq))
                for indel in tunneled_indels:
                    indel.step_del_pos(swap_pos)
                return True
        return False

    def sweep(self,
              muts: dict[int, int],
              end5_ref: int,
              end3_ref: int,
              end5_read: int,
              end3_read: int,
              ref: DNA,
              read: DNA,
              qual: str,
              min_qual: str,
              dels: list[Deletion],
              inns: list[Insertion],
              from3to5: bool,
              tunnel: bool):
        """ Move the indel as far as possible in one direction. """
        while self._try_swap(muts,
                             end5_read,
                             end3_read,
                             ref,
                             read,
                             qual,
                             min_qual,
                             dels,
                             inns,
                             from3to5,
                             tunnel):
            # All actions happen in _try_swap, so loop body is empty.
            pass


def sweep_indels(muts: dict[int, int],
                 end5_ref: int,
                 end3_ref: int,
                 end5_read: int,
                 end3_read: int,
                 refseq: DNA,
                 read: DNA,
                 qual: str,
                 min_qual: str,
                 dels: list[Deletion],
                 inns: list[Insertion],
                 from3to5: bool,
                 tunnel: bool):
    """
    For every insertion and deletion,

    Parameters
    ----------
    muts: dict
        Mutations
    end5_ref: int
        5' most position of the read that is not soft-clipped, using
        reference coordinates.
    end3_ref: int
        3' most position of the read that is not soft-clipped, using
        reference coordinates.
    end5_read: int
        5' most position of the read that is not soft-clipped, using
        read coordinates.
    end3_read: int
        3' most position of the read that is not soft-clipped, using
        read coordinates.
    refseq: DNA
        Reference sequence
    read: DNA
        Sequence of the read
    qual: str
        Phred quality scores of the read, encoded as ASCII characters
    min_qual: int
        The minimum Phred quality score needed to consider a base call
        informative: integer value of the ASCII character
    dels: list[Deletion]
        List of deletions.
    inns: list[Insertion]
        List of insertions.
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
    # Sort the indels by their rank.
    sort_rev = from3to5 != tunnel
    indels.sort(key=lambda idl: idl.rank, reverse=sort_rev)
    while indels:
        indel = indels.pop()
        indel.sweep(muts,
                    end5_ref,
                    end3_ref,
                    end5_read,
                    end3_read,
                    refseq,
                    read,
                    qual,
                    min_qual,
                    dels,
                    inns,
                    from3to5,
                    tunnel)
        idx = len(indels)
        if sort_rev:
            while idx > 0 and indel.rank > indels[idx - 1].rank:
                idx -= 1
        else:
            while idx > 0 and indel.rank < indels[idx - 1].rank:
                idx -= 1
        if idx < len(indels):
            indels.insert(idx, indel)


def find_ambindels(muts: dict[int, int],
                   end5_ref: int,
                   end3_ref: int,
                   end5_read: int,
                   end3_read: int,
                   refseq: DNA,
                   read: DNA,
                   qual: str,
                   min_qual: str,
                   dels: list[Deletion],
                   inns: list[Insertion]):
    """
    Find and label all positions in the vector that are ambiguous due to
    insertions and deletions.

    Parameters
    ----------
    muts: dict
        Mutations
    end5_ref: int
        5' most position of the read that is not soft-clipped, using
        reference coordinates.
    end3_ref: int
        3' most position of the read that is not soft-clipped, using
        reference coordinates.
    end5_read: int
        5' most position of the read that is not soft-clipped, using
        read coordinates.
    end3_read: int
        3' most position of the read that is not soft-clipped, using
        read coordinates.
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
        List of deletions.
    inns: list[Insertion]
        List of insertions.

    """
    # Each indel might be able to be moved in the 5' -> 3' direction
    # (from3to5 is False) or 3' -> 5' direction (from3to5 is True).
    # Test both directions.
    for from3to5 in (False, True):
        # For each indel, try to move it as far as it can go in the
        # direction indicated by from3to5. Allow tunneling so that any
        # runs of consecutive insertions or consecutive deletions can
        # effectively move together.
        sweep_indels(muts,
                     end5_ref,
                     end3_ref,
                     end5_read,
                     end3_read,
                     refseq,
                     read,
                     qual,
                     min_qual,
                     dels,
                     inns,
                     from3to5,
                     tunnel=True)
        if any(d.tunneled for d in dels) or any(i.tunneled for i in inns):
            # If any indel tunneled,
            sweep_indels(muts,
                         end5_ref,
                         end3_ref,
                         end5_read,
                         end3_read,
                         refseq,
                         read,
                         qual,
                         min_qual,
                         dels,
                         inns,
                         from3to5,
                         tunnel=False)

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
