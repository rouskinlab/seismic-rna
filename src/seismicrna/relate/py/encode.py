"""

Relate Code Module

========================================================================

Convert the relationships between reads and a reference from SAM format
(which encodes relationships implicitly as CIGAR strings) to vectorized
format (which encodes relationships explicitly as elements of arrays).

------------------------------------------------------------------------

"""

from ...core.rel import MATCH, SUB_A, SUB_C, SUB_G, SUB_T, ANY_N, IRREC
from ...core.seq import DNA, BASEA, BASEC, BASEG, BASET, BASEN


# Map bases to integer encodings and vice versa.
BASE_ENCODINGS = {BASEA: SUB_A, BASEC: SUB_C, BASEG: SUB_G, BASET: SUB_T,
                  BASEN: IRREC}
BASE_DECODINGS = {code: base for base, code in BASE_ENCODINGS.items()
                  if base != BASEN}


def encode_relate(ref_base: str, read_base: str, read_qual: str, min_qual: str):
    """
    Encode the relation between a base in the read and a base in the
    reference sequence. If the read quality is sufficient, then return
    the match encoding if the read and reference bases match, otherwise
    the encoding of the substitution for the base in the read.
    If the read quality is insufficient, then return the fully ambiguous
    base encoding, that is a match or substitution to any base except
    the reference base, since a "substitution to the reference base"
    would be a match, not a substitution.

    Parameters
    ----------
    ref_base: DNA
        Base in the reference sequence.
    read_base: DNA
        Base in the read sequence.
    read_qual: str
        ASCII encoding for the Phred quality score of the read base.
    min_qual: str
        Minimum value of `read_qual` to not call the relation ambiguous.
    """
    return ((MATCH if ref_base == read_base
             else BASE_ENCODINGS[read_base])
            if (read_qual >= min_qual
                and read_base != BASEN
                and ref_base != BASEN)
            else ANY_N ^ BASE_ENCODINGS[ref_base])


def encode_match(read_base: str, read_qual: str, min_qual: str):
    """
    A more efficient version of `encode_compare` given prior knowledge
    from the CIGAR string that the read and reference match at this
    position. Note that there is no analagous version when there is a
    known substitution because substitutions are relatively infrequent,
    so optimizing their processing would speed the program only slightly
    while making the source code more complex and harder to maintain.
    """
    return (MATCH if (read_qual >= min_qual and read_base != BASEN)
            else ANY_N ^ BASE_ENCODINGS[read_base])

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
