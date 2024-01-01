from ..arg import opt_min_phred, opt_phred_enc


def encode_phred(phred_score: int, phred_encoding: int):
    """
    Encode a numeric Phred quality score as an ASCII character.

    Parameters
    ----------
    phred_score : int
        The Phred score as an integer.
    phred_encoding : int
        The encoding offset for Phred scores. A Phred score is encoded
        as the character whose ASCII value is the sum of the phred score
        and the encoding offset.

    Returns
    -------
    str
        The character whose ASCII code, in the encoding scheme of the
        FASTQ file, represents valid quality.
    """
    return chr(phred_score + phred_encoding)


def decode_phred(quality_code: str, phred_encoding: int):
    """
    Decode the ASCII character for a Phred quality score to an integer.

    Parameters
    ----------
    quality_code : str
        The Phred score encoded as an ASCII character.
    phred_encoding : int
        The encoding offset for Phred scores. A Phred score is encoded
        as the character whose ASCII value is the sum of the phred score
        and the encoding offset.

    Returns
    -------
    str
        The character whose ASCII code, in the encoding scheme of the
        FASTQ file, represents valid quality.
    """
    return ord(quality_code) - phred_encoding


# Default high, medium, and low quality codes.
LO_PHRED = 0
OK_PHRED = opt_min_phred.default
HI_PHRED = 40
LO_QUAL = encode_phred(LO_PHRED, opt_phred_enc.default)
OK_QUAL = encode_phred(OK_PHRED, opt_phred_enc.default)
HI_QUAL = encode_phred(HI_PHRED, opt_phred_enc.default)

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
