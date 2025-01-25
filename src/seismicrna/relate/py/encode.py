from ...core.rel import MATCH, SUB_A, SUB_C, SUB_G, SUB_T, ANY_N
from ...core.seq import DNA, BASEA, BASEC, BASEG, BASET

# Encode bases as substitutions and vice versa.
SUBS_ENCODINGS = {BASEA: SUB_A,
                  BASEC: SUB_C,
                  BASEG: SUB_G,
                  BASET: SUB_T}
SUBS_DECODINGS = {code: base for base, code in SUBS_ENCODINGS.items()}


def is_acgt(base: str):
    """ Check whether a character is a standard DNA base. """
    return base == BASEA or base == BASEC or base == BASEG or base == BASET


def encode_relate(ref_base: str, read_base: str, read_qual: str, min_qual: str):
    """ Encode the relationship between a base in the read and a base in
    the reference sequence.

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
    if not is_acgt(ref_base):
        # If the reference base is unknown, then the read base could be
        # a match or a substitution to any base.
        return ANY_N
    if read_qual < min_qual or not is_acgt(read_base):
        # If the read base is unknown or low quality, then it could be
        # a match or substitution to any base but the reference base.
        return ANY_N ^ SUBS_ENCODINGS[ref_base]
    # Both reference and read bases are known and high-quality.
    return MATCH if ref_base == read_base else SUBS_ENCODINGS[read_base]
