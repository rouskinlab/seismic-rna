import numpy as np
from numba import jit

from ..core.ngs import HI_QUAL, LO_QUAL
from ..core.rel import (MATCH,
                        NOCOV,
                        SUB_A,
                        SUB_C,
                        SUB_G,
                        SUB_T,
                        DELET)
from ..core.seq import DNA, BASEA, BASEC, BASEG, BASET, BASEN
from ..relate.data import QnamesDataset, RelateDataset


@jit()
def _generate_fastq_read_qual(seq: str,
                              rels: np.ndarray,
                              hi_qual: str,
                              lo_qual: str):
    """ Generate a FASTQ line for a read. """
    length = rels.size
    read = np.empty(length, dtype="U1")
    qual = np.empty(length, dtype="U1")
    empty = chr(0)
    for i in range(length):
        rel = rels[i]
        if rel == MATCH:
            read[i] = seq[i]
            qual[i] = hi_qual
        elif rel == SUB_T:
            read[i] = BASET
            qual[i] = hi_qual
        elif rel == SUB_G:
            read[i] = BASEG
            qual[i] = hi_qual
        elif rel == SUB_C:
            read[i] = BASEC
            qual[i] = hi_qual
        elif rel == SUB_A:
            read[i] = BASEA
            qual[i] = hi_qual
        elif rel == DELET:
            read[i] = empty
            qual[i] = empty
        else:
            read[i] = BASEN
            qual[i] = lo_qual
    nonzero = np.flatnonzero(read)
    if nonzero.size != length:
        read = read[nonzero]
        qual = qual[nonzero]
    return "".join(read), "".join(qual)


def generate_fastq_record(seq: str,
                          name: str,
                          rels: np.ndarray,
                          reverse: bool = False,
                          hi_qual: str = HI_QUAL,
                          lo_qual: str = LO_QUAL):
    """ Generate a FASTQ line for a read. """
    # Find the first and last positions with coverage.
    covered, = np.where(rels != NOCOV)
    try:
        first = covered[0]
        last = covered[-1]
    except IndexError:
        return DNA("")
    read, qual = _generate_fastq_read_qual(seq[first: last + 1],
                                           rels[first: last + 1],
                                           hi_qual,
                                           lo_qual)
    if reverse:
        read = DNA(read).rc
        qual = qual[::-1]
    return f"@{name}\n{read}\n+\n{qual}\n"


def generate_fastq(rdata: RelateDataset, ndata: QnamesDataset):
    """ Generate FASTQ file(s) from a dataset. """
    seq = str(rdata.section.seq)
    for rbatch, nbatch in zip(rdata.iter_batches(),
                              ndata.iter_batches(),
                              strict=True):
        for rels, name in zip(rbatch.matrix.values, nbatch.names, strict=True):
            yield generate_fastq_record(seq, name, rels)
