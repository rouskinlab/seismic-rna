from functools import cache
from itertools import chain, product
from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np

from . import path

logger = getLogger(__name__)

# FASTA format
FASTA_NAMESYM = b">"

# Byte encodings for nucleic acid alphabets
BASES = b"ACGT"
COMPS = b"TGCA"
RBASE = b"ACGU"
RCOMP = b"UGCA"
BASEN = b"N"
BASES_LIST = [chr(base).encode() for base in BASES]
BASES_ARR = np.frombuffer(BASES, dtype=np.uint8)

# Integer encodings for nucleic acid alphabets
A_INT = BASES[0]
C_INT = BASES[1]
G_INT = BASES[2]
T_INT = BASES[3]
N_INT = BASEN[0]


class Seq(bytes):
    __slots__ = []

    alph = b""
    comp = b""
    alphaset = set(alph)
    comptrans = alph.maketrans(alph, comp)

    def __init__(self, seq: bytes):
        self.validate_seq(seq)
        super().__init__()

    @classmethod
    def validate_seq(cls, seq: bytes | bytearray):
        if not isinstance(seq, (bytes, bytearray)):
            raise TypeError(
                f"Expected bytes or bytearray, but got '{type(seq).__name__}'")
        if not seq:
            raise ValueError("seq is empty")
        if set(seq) - cls.alphaset:
            raise ValueError(f"Invalid characters in seq: '{seq.decode()}'")

    @classmethod
    def is_valid(cls, seq: bytes | bytearray):
        """ Whether the given sequence is valid for the class. """
        # Check if validate_seq raises ValueError (but not TypeError).
        try:
            cls.validate_seq(seq)
        except ValueError:
            # Instantiation raised ValueError: sequence is invalid.
            return False
        # Instantiation succeeded: sequence is valid.
        return True

    @property
    def rc(self):
        return self.__class__(self[::-1].translate(self.comptrans))

    @cache
    def to_int_array(self):
        """ NumPy array of ASCII integers for the sequence. """
        return np.frombuffer(self, dtype=np.uint8)

    @cache
    def to_str_array(self):
        """ NumPy array of Unicode characters for the sequence. """
        return np.array(list(self.decode()))

    def __getitem__(self, item):
        value = super().__getitem__(item)
        # If a single index is selected, then value will be an int.
        # If a slice is selected, then value will be bytes.
        return value if isinstance(value, int) else self.__class__(value)

    @classmethod
    def parse(cls, seq: str):
        return cls(seq.encode())

    def __str__(self):
        return self.decode()


class DNA(Seq):
    alph = BASES
    comp = COMPS
    alphaset = set(alph)
    comptrans = alph.maketrans(alph, comp)

    def tr(self):
        """ Transcribe DNA to RNA. """
        return RNA(self.replace(b"T", b"U"))


class RNA(Seq):
    alph = RBASE
    comp = RCOMP
    alphaset = set(alph)
    comptrans = alph.maketrans(alph, comp)

    def rt(self):
        """ Reverse transcribe RNA to DNA. """
        return DNA(self.replace(b"U", b"T"))


def expand_degenerate_seq(seq: bytes):
    """ Given a (possibly degenerate) sequence, yield every definite
    sequence that could derive from it. Only the degenerate base N is
    supported by this function; other IUPAC codes (e.g. R) are not. """
    # Split the sequence into every segment that does not have an N.
    segs = seq.split(BASEN)
    seg0 = segs[0]
    # The number of N bases is one less than the number of segments.
    ns = len(segs) - 1
    if ns:
        # If the sequence contains at least one N, then yield every
        # possible sequence by replacing each N with each base.
        for bases in product(BASES_LIST, repeat=ns):
            yield DNA(b"".join(chain((seg0,), *zip(bases, segs[1:],
                                                   strict=True))))
    else:
        # If the sequence contains no N bases, then yield it as DNA.
        yield DNA(seg0)


def parse_fasta(fasta: Path, rna: bool = False):
    """ Parse a FASTA file and iterate through the reference names and
    sequences. """
    if not fasta:
        raise TypeError("No FASTA file given")
    seq_type = "RNA" if rna else "DNA"
    logger.info(f"Began parsing FASTA of {seq_type}: {fasta}")
    # Get the name of the set of references.
    refset = path.parse(fasta, path.FastaSeg)[path.REF]
    has_ref_named_refset = False
    # Record the names of all the references.
    names = set()
    with open(fasta, "rb") as f:
        line = f.readline()
        while line:
            # Read the name from the current line.
            if not line.startswith(FASTA_NAMESYM):
                logger.error(f"Name line '{line.strip()}' in {fasta} does not "
                             f"start with name symbol '{FASTA_NAMESYM}'")
                continue
            # Get the name of the reference up to the first whitespace.
            name = line.split(maxsplit=1)[0][len(FASTA_NAMESYM):].decode()
            # Read the sequence of the reference up until the next
            # reference or the end of the file, whichever comes first.
            seqarray = bytearray()
            while (line := f.readline()) and not line.startswith(FASTA_NAMESYM):
                seqarray.extend(line.rstrip().upper())
            # Confirm that the sequence is valid.
            try:
                seq = RNA(seqarray) if rna else DNA(seqarray)
            except Exception as error:
                logger.error(
                    f"Failed to parse {seq_type} sequence in {fasta}: {error}")
                continue
            # Confirm that the name is not blank.
            if not name:
                logger.error(f"Blank name line '{line.strip()}' in {fasta}")
                continue
            # If there are two or more references with the same name,
            # then the sequence of only the first is used.
            if name in names:
                logger.warning(f"Duplicate reference '{name}' in {fasta}")
                continue
            # If any reference has the same name as the file, then the
            # file is not allowed to have any additional references
            # because, if it did, then the files of all references and
            # of only the self-named reference would have the same names
            # and thus be indistinguishable by their paths.
            if name == refset:
                has_ref_named_refset = True
                logger.debug(f"Reference '{name}' had same name as {fasta}")
            if has_ref_named_refset and names:
                raise ValueError(f"Because {fasta} had a reference with the "
                                 f"same name as the file ('{name}'), it was "
                                 f"not allowed to have any other references, "
                                 f"but it also had {', '.join(names)}")
            # Yield the validated name and sequence.
            names.add(name)
            logger.debug(f"Read {seq_type} reference '{name}' of length "
                         f"{len(seq)} from {fasta}")
            yield name, seq
    logger.info(f"Ended parsing {len(names)} {seq_type} sequences from {fasta}")


def get_ref_seq(fasta: Path, ref: str):
    """ Return the sequence of the reference named `ref` in the FASTA
    file; raise ValueError if not found. """
    for curr_ref, curr_seq in parse_fasta(fasta):
        if curr_ref == ref:
            # The reference was found.
            return curr_seq
    # The reference sequence was not found.
    raise ValueError(f"Reference '{ref}' is not in {fasta}")


def write_fasta(fasta: Path, refs: Iterable[tuple[str, Seq]]):
    """ Write an iterable of reference names and DNA sequences to a
    FASTA file. """
    if not fasta:
        raise TypeError("No FASTA file given")
    logger.info(f"Began writing FASTA file: {fasta}")
    # Get the name of the set of references.
    refset = path.parse(fasta, path.FastaSeg)[path.REF]
    has_ref_named_refset = False
    # Record the names of all the references.
    names = set()
    with open(fasta, "xb") as f:
        for name, seq in refs:
            # Confirm that the name is not blank.
            if not name:
                logger.error(f"Blank reference name")
                continue
            # If there are two or more references with the same name,
            # then the sequence of only the first is used.
            if name in names:
                logger.warning(f"Duplicate reference '{name}'")
                continue
            # If any reference has the same name as the file, then the
            # file is not allowed to have any additional references
            # because, if it did, then the files of all references and
            # of only the self-named reference would have the same names
            # and thus be indistinguishable by their paths.
            if name == refset:
                has_ref_named_refset = True
                logger.debug(f"Reference '{name}' had same name as {fasta}")
            if has_ref_named_refset and names:
                raise ValueError(f"Because {fasta} got a reference with the "
                                 f"same name as the file ('{name}'), it was "
                                 f"not allowed to get any other references, "
                                 f"but it also got {', '.join(names)}")
            try:
                f.write(b"".join([FASTA_NAMESYM, name.encode(), b"\n",
                                  seq, b"\n"]))
            except Exception as error:
                logger.error(
                    f"Error writing reference '{name}' to {fasta}: {error}")
            else:
                logger.debug(f"Wrote reference '{name}' of length {len(seq)}"
                             f"to {fasta}")
                names.add(name)
    logger.info(f"Wrote {len(names)} references to {fasta}")
