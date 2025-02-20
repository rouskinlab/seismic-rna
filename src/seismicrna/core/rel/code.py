from ..types import get_uint_type

REL_SIZE = 1
REL_TYPE = get_uint_type(REL_SIZE)

# Integer encodings for relation vectors
IRREC = REL_TYPE(int("00000000", 2))  # 000: irreconcilable mates
MATCH = REL_TYPE(int("00000001", 2))  # 001: match
DELET = REL_TYPE(int("00000010", 2))  # 002: deletion
INS_5 = REL_TYPE(int("00000100", 2))  # 004: 5' of insertion
INS_3 = REL_TYPE(int("00001000", 2))  # 008: 3' of insertion
SUB_A = REL_TYPE(int("00010000", 2))  # 016: substitution to A
SUB_C = REL_TYPE(int("00100000", 2))  # 032: substitution to C
SUB_G = REL_TYPE(int("01000000", 2))  # 064: substitution to G
SUB_T = REL_TYPE(int("10000000", 2))  # 128: substitution to T
NOCOV = REL_TYPE(int("11111111", 2))  # 255: not covered by read
SUB_N = REL_TYPE(SUB_A | SUB_C | SUB_G | SUB_T)
SUB_B = REL_TYPE(SUB_N ^ SUB_A)
SUB_D = REL_TYPE(SUB_N ^ SUB_C)
SUB_H = REL_TYPE(SUB_N ^ SUB_G)
SUB_V = REL_TYPE(SUB_N ^ SUB_T)
ANY_N = REL_TYPE(SUB_N | MATCH)
ANY_B = REL_TYPE(SUB_B | MATCH)
ANY_D = REL_TYPE(SUB_D | MATCH)
ANY_H = REL_TYPE(SUB_H | MATCH)
ANY_V = REL_TYPE(SUB_V | MATCH)
INSRT = REL_TYPE(INS_5 | INS_3)
INDEL = REL_TYPE(DELET | INSRT)
