"""

Sequence Core Module.

========================================================================

Define alphabets and classes for nucleic acid sequences, and functions
for reading them from and writing them to FASTA files.

"""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from functools import cache, cached_property
from itertools import chain, product
from string import printable
from typing import Any

import numpy as np

from ..types import BITS_PER_BYTE, get_uint_type

# Nucleic acid sequence alphabets.
BASEA = "A"
BASEC = "C"
BASEG = "G"
BASET = "T"
BASEU = "U"
BASEN = "N"

# IUPAC extended nucleic acid alphabet.
IUPAC_AC = "M"
IUPAC_AG = "R"
IUPAC_AT = "W"
IUPAC_CG = "S"
IUPAC_CT = "Y"
IUPAC_GT = "K"
IUPAC_ACG = "V"
IUPAC_ACT = "H"
IUPAC_AGT = "D"
IUPAC_CGT = "B"
IUPAC_CODES = frozenset([BASEA,
                         BASEC,
                         BASEG,
                         BASET,
                         BASEU,
                         IUPAC_AC,
                         IUPAC_AG,
                         IUPAC_AT,
                         IUPAC_CG,
                         IUPAC_CT,
                         IUPAC_GT,
                         IUPAC_ACG,
                         IUPAC_ACT,
                         IUPAC_AGT,
                         IUPAC_CGT,
                         BASEN])

# Nucleic acid compression symbols.
COMPRESS_TYPE = get_uint_type(1)
BITS_PER_BASE = 2
NUM_BASES = 2 ** BITS_PER_BASE
BLOCK_SIZE = BITS_PER_BYTE // BITS_PER_BASE
BLOCK_FORMAT = ":".join(["{}", f"{BASEN}<{BLOCK_SIZE}"])

# Nucleic acid pictogram characters.
PICTA = "▲"
PICTC = "⌠"
PICTN = "○"
PICTG = "⌡"
PICTT = PICTU = "▼"


class XNA(ABC):
    __slots__ = "_seq",

    @classmethod
    @abstractmethod
    def alph(cls) -> tuple[str, str, str, str, str]:
        """ Sequence alphabet. """

    @classmethod
    @abstractmethod
    def pict(cls) -> tuple[str, str, str, str, str]:
        """ Sequence pictograms. """

    @classmethod
    @cache
    def four(cls):
        """ Get the four standard bases. """
        four = tuple(n for n in cls.alph() if n != BASEN)
        if len(four) != NUM_BASES:
            raise ValueError(f"Expected {NUM_BASES} bases, but got {four}")
        return four

    @classmethod
    @cache
    def get_alphaset(cls):
        """ Get the alphabet as a set. """
        return frozenset(cls.alph())

    @classmethod
    @cache
    def get_nonalphaset(cls):
        """ Get the printable characters not in the alphabet. """
        return frozenset(printable) - cls.get_alphaset()

    @classmethod
    @cache
    def get_other_iupac(cls):
        """ Get the IUPAC extended characters not in the alphabet. """
        return frozenset(IUPAC_CODES - cls.get_alphaset())

    @classmethod
    @cache
    def get_comp(cls):
        """ Get the complementary alphabet as a tuple. """
        return tuple(reversed(cls.alph()))

    @classmethod
    @cache
    def get_comptrans(cls):
        """ Get the translation table for complementary bases. """
        return str.maketrans(dict(zip(cls.alph(), cls.get_comp(), strict=True)))

    @classmethod
    @cache
    def get_pictrans(cls):
        """ Get the translation table for pictogram characters. """
        return str.maketrans(dict(zip(cls.alph(), cls.pict(), strict=True)))

    @classmethod
    def t_or_u(cls):
        """ Get the base that is complementary to A. """
        return max(cls.four())

    @classmethod
    def random(cls,
               nt: int,
               a: float = 0.25,
               c: float = 0.25,
               g: float = 0.25,
               t: float = 0.25):
        """
        Return a random sequence of the given length.

        Parameters
        ----------
        nt: int
            Number of nucleotides to simulate. Must be ≥ 0.
        a: float = 0.25
            Expected proportion of A.
        c: float = 0.25
            Expected proportion of C.
        g: float = 0.25
            Expected proportion of G.
        t: float = 0.25
            Expected proportion of T (if DNA) or U (if RNA).

        Returns
        -------
        DNA | RNA
            A random sequence.
        """
        # Calculate expected proportion of N.
        n = 1. - (a + c + g + t)
        if not 0. <= n <= 1.:
            raise ValueError(f"Sum of A, C, G, and {cls.t_or_u()} proportions "
                             f"must be in [0, 1], but got {1. - n}")
        return cls("".join(np.random.default_rng().choice(cls.alph(),
                                                          size=nt,
                                                          p=(a, c, n, g, t))))

    def __init__(self, seq: Any):
        self._seq = str(seq)
        # Check for invalid characters, including lowercase letters.
        if inv := set(self._seq) - self.get_alphaset():
            # If there are any invalid characters, check whether they
            # are valid if converted to uppercase.
            if inv := {i for i in inv if i.upper() not in self.get_alphaset()}:
                # If there are invalid uppercase characters, then raise.
                raise ValueError(f"Invalid {type(self).__name__} bases: {inv}")
            self._seq = self._seq.upper()

    @cached_property
    def rc(self):
        """ Reverse complement. """
        return self.__class__(str(self)[::-1].translate(self.get_comptrans()))

    @cached_property
    def picto(self):
        """ Pictogram string. """
        return str(self).translate(self.get_pictrans())

    @cached_property
    def array(self):
        """ NumPy array of Unicode characters for the sequence. """
        return np.array(list(self))

    def compress(self):
        """ Compress the sequence. """
        return CompressedSeq(self)

    def kmers(self, k: int):
        """ Every subsequence of length k (k-mer). """
        if k < 0:
            raise ValueError(f"k must be ≥ 0, but got {k}")
        for i in range(len(self) - k + 1):
            yield self[i: i + k]

    def __str__(self):
        return self._seq

    def __repr__(self):
        """ Encapsulate the sequence string with the class name. """
        return f"{type(self).__name__}({str(self).__repr__()})"

    def __hash__(self):
        """ Define __hash__ so that Seq subclasses can be used as keys
        for dict-like mappings. Use the hash of the plain string. """
        return str(self).__hash__()

    def __getitem__(self, item):
        """ If item is a slice, then return an instance of the class.
        Otherwise, return an instance of str. """
        value = str(self).__getitem__(item)
        return self.__class__(value) if isinstance(item, slice) else value

    def __contains__(self, item):
        """ Check if a sequence is contained in this sequence. """
        return isinstance(item, type(self)) and str(item) in str(self)

    def __iter__(self):
        return str(self).__iter__()

    def __len__(self):
        return len(str(self))

    def __bool__(self):
        """ Empty sequences return False; all else, True. """
        return bool(str(self))

    def __add__(self, other):
        """ Allow addition (concatenation) of two sequences only if the
        sequences have the same class. """
        if type(self) is type(other):
            return self.__class__(str(self).__add__(str(other)))
        return NotImplemented

    def __mul__(self, other):
        """ Multiply a sequence by an int like a str times an int. """
        return self.__class__(str(self).__mul__(other))

    def __eq__(self, other):
        """ Return True if both the type of the sequence and the bases
        in the sequence match, otherwise False. """
        return type(self) is type(other) and str(self) == str(other)

    def __ne__(self, other):
        return not self.__eq__(other)


class DNA(XNA):

    @classmethod
    def alph(cls):
        return BASEA, BASEC, BASEN, BASEG, BASET

    @classmethod
    def pict(cls):
        return PICTA, PICTC, PICTN, PICTG, PICTT

    def tr(self):
        """ Transcribe DNA to RNA. """
        return RNA(str(self).replace(BASET, BASEU))


class RNA(XNA):

    @classmethod
    def alph(cls):
        return BASEA, BASEC, BASEN, BASEG, BASEU

    @classmethod
    def pict(cls):
        return PICTA, PICTC, PICTN, PICTG, PICTU

    def rt(self):
        """ Reverse transcribe RNA to DNA. """
        return DNA(str(self).replace(BASEU, BASET))


def expand_degenerate_seq(seq: DNA):
    """ Given a (possibly degenerate) sequence, yield every definite
    sequence that could derive from it. Only the degenerate base N is
    supported by this function; other IUPAC codes (e.g. R) are not. """
    # Split the sequence into every segment that does not have an N.
    segs = str(seq).split(BASEN)
    # The number of N bases is one less than the number of segments.
    if ns := len(segs) - 1:
        # If the sequence contains at least one N, then yield every
        # possible sequence by replacing each N with each base.
        for bases in product(DNA.four(), repeat=ns):
            yield DNA("".join(chain((segs[0],), *zip(bases, segs[1:],
                                                     strict=True))))
    else:
        # If the sequence contains no N bases, then yield it as DNA.
        yield DNA(segs[0])


class CompressedSeq(object):
    """ Compress a sequence into two bits per base. """

    def __init__(self, seq: XNA):
        self.r = isinstance(seq, RNA)
        self.s = len(seq)
        self.b = _compress_seq(seq)
        self.n = _find_ns(seq)

    @property
    def type(self):
        return RNA if self.r else DNA

    def decompress(self):
        """ Restore the original sequence. """
        return decompress(self)

    def __eq__(self, other):
        if not isinstance(other, CompressedSeq):
            return NotImplemented
        return (self.r == other.r
                and self.s == other.s
                and self.b == other.b
                and self.n == other.n)


@cache
def _base_to_index(base: str, alph: tuple[str, str, str, str]):
    try:
        return alph.index(base)
    except ValueError:
        return 0


@cache
def _compress_block(block: str, alph: tuple[str, str, str, str]):
    """ Compress one block of a sequence. """
    if len(alph) != NUM_BASES:
        raise ValueError(f"Expected {NUM_BASES} bases, but got {alph}")
    if len(block) > BLOCK_SIZE:
        raise ValueError(f"Expected block of no more than {BLOCK_SIZE} nt, "
                         f"but got {repr(block)}")
    if len(block) < BLOCK_SIZE:
        # Pad the end of the block with Ns until it is long enough.
        block = BLOCK_FORMAT.format(block)
    return sum(_base_to_index(base, alph) << (i * BITS_PER_BASE)
               for i, base in enumerate(block))


def _get_blocks(seq: str):
    return (seq[i: i + BLOCK_SIZE] for i in range(0, len(seq), BLOCK_SIZE))


def _compress_seq(seq: XNA):
    return bytes(_compress_block(block, seq.four())
                 for block in _get_blocks(str(seq)))


def _find_ns(seq: XNA):
    return tuple(match.start() for match in re.finditer(BASEN, str(seq)))


@cache
def _decompress_block(byte: int, alph: tuple[str, str, str, str]):
    """ Decompress one block of a sequence. """
    if len(alph) != NUM_BASES:
        raise ValueError(f"Expected {NUM_BASES} bases, but got {alph}")
    byte = COMPRESS_TYPE(byte)
    return "".join(alph[(byte >> (i * BITS_PER_BASE)) % len(alph)]
                   for i in range(BLOCK_SIZE))


def decompress(seq: CompressedSeq):
    """ Restore the original sequence from a CompressedSeq object. """
    bases = [base for byte in seq.b
             for base in _decompress_block(byte, seq.type.four())][:seq.s]
    for n in seq.n:
        bases[n] = BASEN
    return seq.type("".join(bases))

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
