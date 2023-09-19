"""

Sequence Core Module.

========================================================================

Define alphabets and classes for nucleic acid sequences, and functions
for reading them from and writing them to FASTA files.

"""

from functools import cache, cached_property
from itertools import chain, product
from string import printable
from typing import Any

import numpy as np

from .rand import rng

# Nucleic acid sequence alphabets.
BASEA = 'A'
BASEC = 'C'
BASEG = 'G'
BASET = 'T'
BASEU = 'U'
BASEN = 'N'

# Nucleic acid pictogram characters.
PICTA = '▲'
PICTC = '⌠'
PICTN = '○'
PICTG = '⌡'
PICTT = PICTU = '▼'


class Seq(object):
    __slots__ = "_seq",

    alph: tuple[str, str, str, str, str]
    pict: tuple[str, str, str, str, str]

    def __init__(self, seq: Any):
        self._seq = str(seq)
        if invalid := set(self._seq) - self.get_alphaset():
            raise ValueError(
                f"Invalid {self.__class__.__name__} bases: {sorted(invalid)}")

    @cached_property
    def rc(self):
        """ Reverse complement. """
        return self.__class__(str(self)[::-1].translate(self.get_comptrans()))

    @cached_property
    def picto(self):
        """ Pictogram string. """
        return str(self).translate(self.get_pictrans())

    @cache
    def to_array(self):
        """ NumPy array of Unicode characters for the sequence. """
        return np.array(list(self))

    @classmethod
    def t_or_u(cls):
        """ Get the base that is complementary to A. """
        return cls.alph[-1]

    @classmethod
    @cache
    def get_unambig(cls):
        """ Get the unambiguous bases. """
        return tuple(n for n in cls.alph if n != BASEN)

    @classmethod
    @cache
    def get_alphaset(cls):
        """ Get the alphabet as a set. """
        return frozenset(cls.alph)

    @classmethod
    @cache
    def get_nonalphaset(cls):
        """ Get the printable characters not in the alphabet. """
        return frozenset(printable) - cls.get_alphaset()

    @classmethod
    @cache
    def get_comp(cls):
        """ Get the complementary alphabet as a tuple. """
        return tuple(reversed(cls.alph))

    @classmethod
    @cache
    def get_comptrans(cls):
        """ Get the translation table for complementary bases. """
        return str.maketrans(dict(zip(cls.alph, cls.get_comp(), strict=True)))

    @classmethod
    @cache
    def get_pictrans(cls):
        """ Get the translation table for pictogram characters. """
        return str.maketrans(dict(zip(cls.alph, cls.pict, strict=True)))

    @classmethod
    def random(cls, nt: int,
               a: float = 0.25, c: float = 0.25,
               g: float = 0.25, t: float = 0.25):
        """
        Return a random sequence of the given length.

        Parameters
        ----------
        nt: int
            Number of nucleotides to simulate. Must be ≥ 1.
        a: float = 0.25
            Expected proportion of A.
        c: float = 0.25
            Expected proportion of C.
        g: float = 0.25
            Expected proportion of G.
        t: float = 0.25
            Expected proportion of T (DNA) or U (RNA).

        Returns
        -------
        Seq
            A random sequence.
        """
        # Calculate expected proportion of N.
        n = 1. - (a + c + g + t)
        if not 0. <= n <= 1.:
            raise ValueError(f"Sum of A, C, G, and {cls.t_or_u()} proportions "
                             f"must be in [0, 1], but got {1. - n}")
        return cls("".join(rng.choice(cls.alph, size=nt, p=(a, c, n, g, t))))

    def __str__(self):
        return self._seq

    def __repr__(self):
        """ Encapsulate the sequence string with the class name. """
        return f"{self.__class__.__name__}({str(self).__repr__()})"

    def __hash__(self):
        """ Define __hash__ so that Seq subclasses can be used as keys
        for dict-like mappings. Use the hash of the plain string. """
        return str(self).__hash__()

    def __getitem__(self, item):
        """ If item is a slice, then return an instance of the class.
        Otherwise, return an instance of str. """
        value = str(self).__getitem__(item)
        return self.__class__(value) if isinstance(item, slice) else value

    def __iter__(self):
        return self._seq.__iter__()

    def __len__(self):
        return len(str(self))

    def __bool__(self):
        """ Empty sequences return False; all else, True. """
        return bool(str(self))

    def __add__(self, other):
        """ Allow addition (concatenation) of two sequences only if the
        sequences have the same class. """
        if self.__class__ is other.__class__:
            return self.__class__(str(self).__add__(str(other)))
        return NotImplemented

    def __mul__(self, other):
        """ Multiply a sequence by an int like a str times an int. """
        return self.__class__(str(self).__mul__(other))

    def __eq__(self, other):
        """ Return True if both the type of the sequence and the bases
        in the sequence match, otherwise False. """
        return self.__class__ is other.__class__ and str(self) == str(other)

    def __ne__(self, other):
        return not self.__eq__(other)


class DNA(Seq):
    alph = BASEA, BASEC, BASEN, BASEG, BASET
    pict = PICTA, PICTC, PICTN, PICTG, PICTT

    @cache
    def tr(self):
        """ Transcribe DNA to RNA. """
        return RNA(str(self).replace(BASET, BASEU))


class RNA(Seq):
    alph = BASEA, BASEC, BASEN, BASEG, BASEU
    pict = PICTA, PICTC, PICTN, PICTG, PICTU

    @cache
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
        for bases in product(DNA.get_unambig(), repeat=ns):
            yield DNA("".join(chain((segs[0],), *zip(bases, segs[1:],
                                                     strict=True))))
    else:
        # If the sequence contains no N bases, then yield it as DNA.
        yield DNA(segs[0])

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
