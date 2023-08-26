"""

Sequence Core Module.

========================================================================

Define alphabets and classes for nucleic acid sequences, and functions
for reading them from and writing them to FASTA files.

"""

from functools import cache, cached_property
from itertools import chain, product
from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np

from . import path
from .sim import rng

logger = getLogger(__name__)

# FASTA record mark.
FASTA_RECORD = '>'

# Nucleic acid sequence alphabets.
BASEA = 'A'
BASEC = 'C'
BASEG = 'G'
BASET = 'T'
BASEU = 'U'
BASEN = 'N'


class Seq(str):
    __slots__ = []

    alph: tuple[str, str, str, str]

    def __init__(self, seq: str):
        self.validate_seq(seq)
        super().__init__()

    @cached_property
    def rc(self):
        """ Reverse complement. """
        return self.__class__(self[::-1].translate(self.get_comptrans()))

    @cache
    def to_array(self):
        """ NumPy array of Unicode characters for the sequence. """
        return np.array(list(self))

    @classmethod
    @cache
    def get_alphaset(cls):
        """ Get the alphabet as a set. """
        return set(cls.alph)

    @classmethod
    @cache
    def validate_seq(cls, seq: str):
        if not isinstance(seq, str):
            raise TypeError(f"Expected str, but got {type(seq).__name__}")
        if invalid := set(seq) - cls.get_alphaset():
            raise ValueError(f"Invalid {cls.__name__} bases: {sorted(invalid)}")

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
    def random(cls, nt: int,
               a: float = 0.25, c: float = 0.25,
               g: float = 0.25, tu: float = 0.25):
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
        tu: float = 0.25
            Expected proportion of T (for DNA) or U (for RNA).

        Returns
        -------
        Seq
            A random sequence.
        """
        return cls("".join(rng.choice(cls.alph, size=nt, p=(a, c, g, tu))))

    def __getitem__(self, item):
        if isinstance(item, slice):
            # For slices, return an instance of the class, not of str.
            return self.__class__(super().__getitem__(item))
        # Otherwise, return an instance of str.
        return super().__getitem__(item)

    def __repr__(self):
        """ Encapsulate the sequence string with DNA(). """
        return f"{self.__class__.__name__}({super().__repr__()})"

    def __eq__(self, other):
        """ Return True if both the type of the sequence and the bases
        in the sequence match, otherwise False. """
        return other.__class__ is self.__class__ and super().__eq__(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        """ Define __hash__ so that Seq subclasses can be used as keys
        for dict-like mappings. Use the hash of the plain string. """
        return super().__hash__()

    def __bool__(self):
        """ Empty sequences return False; all else, True. """
        return bool(len(self))


class DNA(Seq):
    alph = BASEA, BASEC, BASEG, BASET

    @cache
    def tr(self):
        """ Transcribe DNA to RNA. """
        return RNA(self.replace(BASET, BASEU))


class RNA(Seq):
    alph = BASEA, BASEC, BASEG, BASEU

    @cache
    def rt(self):
        """ Reverse transcribe RNA to DNA. """
        return DNA(self.replace(BASEU, BASET))


class DNAmbig(DNA):
    alph = BASEA, BASEC, BASEN, BASEG, BASET

    @cache
    def tr(self):
        """ Transcribe DNAmbig to RNAmbig. """
        return RNAmbig(self.replace(BASET, BASEU))


class RNAmbig(RNA):
    alph = BASEA, BASEC, BASEN, BASEG, BASEU

    @cache
    def rt(self):
        """ Reverse transcribe RNAmbig to DNAmbig. """
        return DNAmbig(self.replace(BASEU, BASET))


def expand_degenerate_seq(seq: DNAmbig):
    """ Given a (possibly degenerate) sequence, yield every definite
    sequence that could derive from it. Only the degenerate base N is
    supported by this function; other IUPAC codes (e.g. R) are not. """
    # Split the sequence into every segment that does not have an N.
    segs = seq.split(BASEN)
    # The number of N bases is one less than the number of segments.
    if ns := len(segs) - 1:
        # If the sequence contains at least one N, then yield every
        # possible sequence by replacing each N with each base.
        for bases in product(DNA.alph, repeat=ns):
            yield DNA("".join(chain((segs[0],), *zip(bases, segs[1:],
                                                     strict=True))))
    else:
        # If the sequence contains no N bases, then yield it as DNA.
        yield DNA(segs[0])


def parse_fasta(fasta: Path, rna: bool = False):
    """ Parse a FASTA file and iterate through the reference names and
    sequences. """
    if not fasta:
        raise TypeError("No FASTA file given")
    seq_type = RNA if rna else DNA
    logger.info(f"Began parsing FASTA of {seq_type}: {fasta}")
    # Get the name of the set of references.
    refset = path.parse(fasta, path.FastaSeg)[path.REF]
    has_ref_named_refset = False
    # Record the names of all the references.
    names = set()
    with open(fasta) as f:
        line = f.readline()
        while line:
            # Read the name from the current line.
            if not line.startswith(FASTA_RECORD):
                logger.error(f"Name line '{line.strip()}' in {fasta} does not "
                             f"start with name symbol '{FASTA_RECORD}'")
                continue
            # Get the name of the reference up to the first whitespace.
            name = line.split(maxsplit=1)[0][len(FASTA_RECORD):]
            # Read the sequence of the reference up until the next
            # reference or the end of the file, whichever comes first.
            segments = list()
            while (line := f.readline()) and not line.startswith(FASTA_RECORD):
                segments.append(line.rstrip().upper())
            # Confirm that the sequence is valid.
            try:
                seq = seq_type("".join(segments))
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
            logger.debug(f"Read {seq_type.__name__} reference '{name}' "
                         f"of length {len(seq)} from {fasta}")
            yield name, seq
    logger.info(f"Ended parsing {len(names)} {seq_type.__name__} sequence(s) "
                f"from {fasta}")


def get_ref_seq(fasta: Path, ref: str):
    """ Return the sequence of the reference named `ref` in the FASTA
    file; raise ValueError if not found. """
    for curr_ref, curr_seq in parse_fasta(fasta):
        if curr_ref == ref:
            # The reference was found.
            return curr_seq
    # The reference sequence was not found.
    raise ValueError(f"Reference '{ref}' is not in {fasta}")


def write_fasta(fasta: Path, refs: Iterable[tuple[str, Seq]],
                overwrite: bool = False):
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
    with open(fasta, 'w' if overwrite else 'x') as f:
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
                f.write("".join([FASTA_RECORD, name, '\n', seq, '\n']))
            except Exception as error:
                logger.error(
                    f"Error writing reference '{name}' to {fasta}: {error}")
            else:
                logger.debug(f"Wrote reference '{name}' of length {len(seq)}"
                             f"to {fasta}")
                names.add(name)
    logger.info(f"Wrote {len(names)} sequences(s) to {fasta}")
