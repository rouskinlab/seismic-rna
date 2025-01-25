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
    int
        The Phred quality score represented by the ASCII character.
    """
    return ord(quality_code) - phred_encoding


# Default high, medium, and low quality codes.
LO_PHRED = 0
OK_PHRED = opt_min_phred.default
HI_PHRED = 40
LO_QUAL = encode_phred(LO_PHRED, opt_phred_enc.default)
OK_QUAL = encode_phred(OK_PHRED, opt_phred_enc.default)
HI_QUAL = encode_phred(HI_PHRED, opt_phred_enc.default)
