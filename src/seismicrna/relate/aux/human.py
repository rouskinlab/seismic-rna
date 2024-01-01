from sys import byteorder

import numpy as np

from ...core.rel import (MATCH,
                         DELET,
                         INS_5,
                         INS_3,
                         SUB_A,
                         SUB_C,
                         SUB_G,
                         SUB_T,
                         NOCOV,
                         IRREC)

# Number of unique bytes
N_BYTES = 2 ** 8
# Map each common byte in the vector encoding (e.g. MATCH) to a byte of
# human-readable text.
BYTE_SYMBOLS = {NOCOV: b'_',
                MATCH: b'=',
                DELET: b'.',
                INS_5 | MATCH: b'{',
                INS_5 | MATCH | INS_3: b'+',
                INS_3 | MATCH: b'}',
                SUB_A: b'A',
                SUB_C: b'C',
                SUB_G: b'G',
                SUB_T: b'T',
                SUB_A | SUB_C | SUB_G | MATCH: b'N',
                SUB_T | SUB_A | SUB_C | MATCH: b'N',
                SUB_G | SUB_T | SUB_A | MATCH: b'N',
                SUB_C | SUB_G | SUB_T | MATCH: b'N',
                IRREC: b'!'}
# Default for uncommon bytes in the vector encoding.
OTHER_SYMBOL = b'?'
# Create an array that maps each vector byte to its readable character.
map_array = bytearray(OTHER_SYMBOL) * N_BYTES
for byte, symbol in BYTE_SYMBOLS.items():
    map_array[byte] = int.from_bytes(symbol, byteorder)
# Create a translation table from vector to human-readable encodings.
map_table = bytes.maketrans(bytes(range(N_BYTES)), map_array)


def humanize_relvec(relvec: np.ndarray):
    """ Translate a binary relation vector to human-readable text. """
    return relvec.tobytes().translate(map_table)

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
