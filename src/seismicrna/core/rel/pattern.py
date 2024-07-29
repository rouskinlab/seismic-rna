"""

Core -- Pattern Module

========================================================================

"""

from __future__ import annotations

import re
from functools import cache
from itertools import product
from logging import getLogger
from typing import Iterable

from .code import (MATCH,
                   DELET,
                   INS_5,
                   INS_3,
                   SUB_A,
                   SUB_C,
                   SUB_G,
                   SUB_T,
                   REL_TYPE)
from ..seq import BASEA, BASEC, BASEG, BASET, DNA

logger = getLogger(__name__)

READ_DEL = "D"
READ_INS = "I"


class HalfRelPattern(object):
    """ """

    __slots__ = "a", "c", "g", "t"

    ref_bases = "".join(DNA.four())
    read_bases = "".join((ref_bases, READ_DEL, READ_INS))
    mut_bits = bytes([SUB_A, SUB_C, SUB_G, SUB_T, DELET, INS_5 | INS_3])
    fmt_plain = "{}{}"
    fmt_fancy = "{} -> {}"
    ptrn_plain = re.compile(f"([{ref_bases.lower()}])([{read_bases.lower()}])")
    ptrn_fancy = re.compile(f"([{ref_bases}]) -> ([{read_bases}])")

    @classmethod
    def as_match(cls, code: str) -> re.Match[str]:
        """
        Return a re.Match object if the code matches either the plain or
        the fancy format. Otherwise, raise ValueError.
        """
        # If code matches ptrn_plain, cls.ptrn_plain.match(code.lower())
        # is truthy, so short-circuit the OR and return the plain match.
        # If code matches ptrn_fancy, cls.ptrn_fancy.match(code.upper())
        # is truthy, so match becomes truthy and is returned.
        match = (cls.ptrn_plain.match(code.lower())
                 or cls.ptrn_fancy.match(code.upper()))
        if match:
            return match
        raise ValueError(f"Failed to match code: {repr(code)}")

    @classmethod
    def as_plain(cls, code: str):
        """
        Convert a ref-read code into plain format, as follows:

        - 2-character lowercase string
        - 1st character is reference base
        - 2nd character is read base, 'd' (deletion), or 'i' (insertion)

        Examples:

        - 'ag' means an A in the reference is a G in the read
        - 'cd' means a C in the reference is deleted in the read
        """
        return cls.fmt_plain.format(*cls.as_match(code).groups()).lower()

    @classmethod
    def as_fancy(cls, code: str):
        """
        Convert a ref-read code into fancy format, as follows:

        - 6-character uppercase string
        - 1st character is reference base
        - 6th character is read base, 'D' (deletion), or 'I' (insertion)
        - 2nd to 5th characters form an arrow: ' -> '

        Examples:

        - 'A -> G' means an A in the reference is a G in the read
        - 'C -> D' means a C in the reference is deleted in the read
        """
        return cls.fmt_fancy.format(*cls.as_match(code).groups()).upper()

    @classmethod
    def compile(cls, codes: Iterable[str]):
        """
        Given one or more codes in plain or fancy format, return a dict
        that maps each reference base to a pattern that will match all
        and only the codes given for that reference base.

        This function is the inverse of `cls.decompile`.
        """
        # Create a dict that maps each reference base to a query byte,
        # which is an integer in the range [0, 256). Initialize to 0.
        queries: dict[str, int] = {ref: 0 for ref in cls.ref_bases}
        # For each code given, get the ref and read bases by converting
        # the code to plain format, then to uppercase, then to a tuple.
        for ref, read in map(str.upper, map(cls.as_plain, codes)):
            # Update the query byte for the reference base. If the read
            # and reference bases are equal, then this code represents
            # a match, so update using the match byte, MATCH.
            # Otherwise, update using the mutation bit that corresponds
            # to the read base (it is at the same index in cls.mut_bytes
            # as the read base is in cls.read_bases). Update by taking
            # the bitwise OR so that all query bytes are accumulated.
            queries[ref] |= (MATCH if read == ref
                             else cls.mut_bits[cls.read_bases.index(read)])
        return queries

    @classmethod
    def decompile(cls, patterns: dict[str, int]):
        """
        For each reference base and its one-byte pattern, yield all
        codes that the pattern will count.

        This function is the inverse of `cls.compile`.
        """
        # Check each pattern.
        for ref, pattern in patterns.items():
            if ref not in cls.ref_bases:
                raise ValueError(f"Invalid reference base: {repr(ref)}")
            if pattern & MATCH:
                # The pattern has its match bit set to 1, so the code
                # in which this ref base matches the read base counts.
                yield cls.as_fancy(f"{ref}{ref}")
            # For each mutation bit, check whether the pattern has the bit
            # set to 1.
            for mut_bit, read in zip(cls.mut_bits, cls.read_bases, strict=True):
                if pattern & mut_bit:
                    # If the mutation bit is set to 1 in the pattern,
                    # then the code where the ref base becomes the read
                    # base (or deletion/insertion) counts.
                    yield cls.as_fancy(f"{ref}{read}")

    @classmethod
    def from_report_format(cls, mut_codes: Iterable[str]):
        return cls(*list(mut_codes))

    @classmethod
    def from_counts(cls, *,
                    count_ref: bool = False,
                    count_sub: bool = False,
                    count_del: bool = False,
                    count_ins: bool = False,
                    discount: Iterable[str] = ()):
        """
        Return a new SemiBitCaller by specifying which general types of
        relationships are to be counted.

        Parameters
        ----------
        count_ref: bool = False
            Whether to call True all matches between the read and ref.
        count_sub: bool = False
            Whether to call True all substitutions in the read.
        count_del: bool = False
            Whether to call True all deletions in the read.
        count_ins: bool = False
            Whether to call True all insertions in the read.
        discount: Iterable[str] = ()
            Do not count any of these relationships between the read and
            the reference, even if they would be counted according to
            any of the other parameters. Should be an iterable of str in
            either plain or fancy format (except case-insensitive).

        Returns
        -------
        HalfRelPattern
            New HalfRefPattern instance that counts the specified bytes.
        """
        codes: set[str] = set()
        if count_ref:
            # Count all matches between the read and reference.
            codes.update(cls.as_plain(2 * base) for base in cls.ref_bases)
        if count_sub:
            # Count all substitutions in the read.
            codes.update(cls.as_plain(f"{base1}{base2}")
                         for base1, base2 in product(cls.ref_bases, repeat=2)
                         if base1 != base2)
        if count_del:
            # Count all deletions in the read.
            codes.update(cls.as_plain(f"{base}D") for base in cls.ref_bases)
        if count_ins:
            # Count all insertions in the read.
            codes.update(cls.as_plain(f"{base}I") for base in cls.ref_bases)
        # Remove all the codes to be discounted.
        codes -= set(map(cls.as_plain, discount))
        return cls(*codes)

    @classmethod
    @cache
    def allc(cls):
        """ Fits every relationship except for no coverage. """
        return cls.from_counts(count_ref=True,
                               count_sub=True,
                               count_del=True,
                               count_ins=True)

    @classmethod
    @cache
    def muts(cls):
        """ Fits every mutation. """
        return cls.from_counts(count_sub=True,
                               count_del=True,
                               count_ins=True)

    @classmethod
    @cache
    def refs(cls):
        """ Fits every match. """
        return cls.from_counts(count_ref=True)

    @classmethod
    @cache
    def none(cls):
        """ Fits nothing. """
        return cls()

    def __init__(self, *codes: str):
        # Compile the codes into patterns.
        patterns = self.compile(codes)
        self.a = patterns.pop(BASEA, 0)
        self.c = patterns.pop(BASEC, 0)
        self.g = patterns.pop(BASEG, 0)
        self.t = patterns.pop(BASET, 0)
        if patterns:
            raise ValueError(f"Unexpected reference base(s): {patterns}")

    @property
    def patterns(self):
        return {BASEA: self.a, BASEC: self.c, BASEG: self.g, BASET: self.t}

    @property
    def codes(self):
        """ Return the codes of the relationships counted. """
        return sorted(self.decompile(self.patterns))

    def fits(self, base: str, rel: int):
        """ Test whether a relationship code fits the pattern. """
        return ((pattern := self.patterns.get(base)) is not None
                and (REL_TYPE(rel) | pattern) == pattern)

    def to_report_format(self):
        """ Return the types of counted relationships as a list. """
        return self.codes

    def intersect(self, other: HalfRelPattern):
        """ Intersect the HalfRelPattern with another. """
        return self.__class__(*(set(self.codes) & set(other.codes)))

    def __str__(self):
        return f"{type(self).__name__} {self.patterns}"

    def __hash__(self):
        return hash(tuple(getattr(self, x) for x in self.__slots__))

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.patterns == other.patterns


class RelPattern(object):
    __slots__ = "yes", "nos"

    @classmethod
    def from_counts(cls,
                    count_del: bool = False,
                    count_ins: bool = False,
                    discount: Iterable[str] = ()):
        """ Return a new RelPattern by specifying which general types of
        mutations are to be counted, with optional ones to discount. """
        discount = list(discount)
        return cls(HalfRelPattern.from_counts(count_sub=True,
                                              count_del=count_del,
                                              count_ins=count_ins,
                                              discount=discount),
                   HalfRelPattern.from_counts(count_ref=True,
                                              discount=discount))

    @classmethod
    @cache
    def allc(cls):
        """ Fits every relationship except for no coverage. """
        return cls(HalfRelPattern.allc(), HalfRelPattern.none())

    @classmethod
    @cache
    def muts(cls):
        """ Fits every mutation. """
        return cls(HalfRelPattern.muts(), HalfRelPattern.refs())

    def __init__(self, yes: HalfRelPattern, nos: HalfRelPattern):
        self.yes = yes
        self.nos = nos

    def fits(self, base: str, rel: int):
        """ Check whether the base and relationship give a definitive
        result and whether they fit the pattern. """
        is_yes = self.yes.fits(base, rel)
        is_nos = self.nos.fits(base, rel)
        return is_yes != is_nos, is_yes

    def intersect(self, other: RelPattern | None, invert: bool = False):
        """ Intersect the pattern with another. """
        if other is not None:
            yes = self.yes.intersect(other.yes)
            nos = self.nos.intersect(other.nos)
        else:
            yes = self.yes
            nos = self.nos
        return self.__class__(nos, yes) if invert else self.__class__(yes, nos)

    def invert(self):
        """ Swap the `yes` and `nos` patterns. """
        return self.__class__(self.nos, self.yes)

    def __str__(self):
        return f"{type(self).__name__}  +[{self.yes}]  -[{self.nos}]"

    def __hash__(self):
        return hash(tuple(getattr(self, x) for x in self.__slots__))

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.yes == other.yes and self.nos == other.nos

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
