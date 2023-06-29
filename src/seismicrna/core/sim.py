"""
Core -- Simulation Module
========================================================================
Auth: Matty
Date: 2023-06-06
"""

from itertools import product
from typing import Sequence

import numpy as np

from .rel import DELET, INS_5, INS_3, NOCOV, encode_relate
from .sect import Section
from .seq import DNA, BASES, BASES_ARR

READ_SIM_COUNT = 0



def rand_dna(length: int):
    """ Return a random DNA sequence of the given length ≥ 1. """
    return DNA(np.random.choice(BASES_ARR, size=length).tobytes())


def rand_rel(length: int, n_miss: int = 0, n_lowq: int = 0,
             n_sub: int = 0, n_del: int = 0, n_ins: int = 0):
    """ Return a random relation vector. """


def rand_read(refseq: DNA):
    """ Return a random read and relation vector given a reference. """


def rand_sam_line(refseq: DNA, ref: str):
    """
    Generate the SAM line of a random read given a reference sequence.

    Parameters
    ----------
    refseq: DNA
        Sequence of the reference.
    ref: str
        Name of the reference.

    """
