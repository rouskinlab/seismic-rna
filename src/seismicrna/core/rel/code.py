from ..types import get_uint_type

REL_SIZE = 1
REL_TYPE = get_uint_type(REL_SIZE)

# Integer encodings for relationships
IRREC = int("00000000", 2)  # 000: irreconcilable mates
MATCH = int("00000001", 2)  # 001: match
DELET = int("00000010", 2)  # 002: deletion
INS_5 = int("00000100", 2)  # 004: 5' of insertion
INS_3 = int("00001000", 2)  # 008: 3' of insertion
SUB_A = int("00010000", 2)  # 016: substitution to A
SUB_C = int("00100000", 2)  # 032: substitution to C
SUB_G = int("01000000", 2)  # 064: substitution to G
SUB_T = int("10000000", 2)  # 128: substitution to T
NOCOV = int("11111111", 2)  # 255: not covered by read
SUB_N = SUB_A | SUB_C | SUB_G | SUB_T
SUB_B = SUB_N ^ SUB_A
SUB_D = SUB_N ^ SUB_C
SUB_H = SUB_N ^ SUB_G
SUB_V = SUB_N ^ SUB_T
ANY_N = SUB_N | MATCH
ANY_B = SUB_B | MATCH
ANY_D = SUB_D | MATCH
ANY_H = SUB_H | MATCH
ANY_V = SUB_V | MATCH
INSRT = INS_5 | INS_3
INDEL = DELET | INSRT
