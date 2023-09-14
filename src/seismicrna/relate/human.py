from sys import byteorder

import numpy as np

from .encode import (MATCH, DELET, INS_5, INS_3,
                     SUB_A, SUB_C, SUB_G, SUB_T,
                     NOCOV, IRREC)

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
