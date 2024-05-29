import gzip
import os
from pathlib import Path
from typing import Any

import numpy as np
from click import command
from numba import jit

from ..core import path
from ..core.arg import (docdef,
                        arg_input_path,
                        opt_max_procs,
                        opt_parallel,
                        opt_force)
from ..core.ngs import HI_QUAL, LO_QUAL
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.rel import (MATCH,
                        NOCOV,
                        SUB_A,
                        SUB_C,
                        SUB_G,
                        SUB_T,
                        DELET)
from ..core.seq import DNA, BASEA, BASEC, BASEG, BASET, BASEN
from ..core.write import need_write, write_mode
from ..pool.data import load_relate_dataset
from ..relate.data import QnamesDataset, RelateDataset

rng = np.random.default_rng()

COMMAND = __name__.split(os.path.extsep)[-1]


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
        read = ""
        qual = ""
    else:
        read, qual = _generate_fastq_read_qual(seq[first: last + 1],
                                               rels[first: last + 1],
                                               hi_qual,
                                               lo_qual)
        if reverse:
            read = str(DNA(read).rc)
            qual = qual[::-1]
    return f"@{name}\n{read}\n+\n{qual}\n"


def _get_common_attr(a: Any, b: Any, attr: str):
    rval = getattr(a, attr)
    nval = getattr(b, attr)
    if rval != nval:
        raise ValueError(f"Got different values of {repr(attr)} for {a} "
                         f"({repr(rval)}) and {b} ({repr(nval)})")
    return rval


def generate_fastq(rdata: RelateDataset,
                   ndata: QnamesDataset,
                   force: bool = False):
    """ Generate FASTQ file(s) from a dataset. """
    top = _get_common_attr(rdata, ndata, "top")
    sample = _get_common_attr(rdata, ndata, "sample")
    ref = _get_common_attr(rdata, ndata, "ref")
    fastq = path.buildpar(*path.DMFASTQ_SEGS,
                          top=top,
                          sample=sample,
                          ref=ref,
                          ext=path.FQ_EXTS[0])
    if need_write(fastq, force):
        seq = str(rdata.section.seq)
        if fastq.suffix.endswith(".gz"):
            open_func = gzip.open
            binary = True
        else:
            open_func = open
            binary = False
        with open_func(fastq, write_mode(force, binary=binary)) as fq:
            for rbatch, nbatch in zip(rdata.iter_batches(),
                                      ndata.iter_batches(),
                                      strict=True):
                num_reads = _get_common_attr(rbatch, nbatch, "num_reads")
                reverse = rng.integers(2, size=num_reads).astype(bool,
                                                                 copy=False)
                for rels, name, rev in zip(rbatch.matrix.values,
                                           nbatch.names,
                                           reverse,
                                           strict=True):
                    record = generate_fastq_record(seq, name, rels, rev)
                    if binary:
                        record = record.encode()
                    fq.write(record)
    return fastq


def run_report(report_file: Path, force: bool):
    rdata = RelateDataset.load(report_file)
    ndata = QnamesDataset.load(report_file)
    return generate_fastq(rdata, ndata, force)


@docdef.auto()
def run(input_path: tuple[str, ...],
        max_procs: int,
        parallel: bool,
        force: bool):
    report_files = path.find_files_chain(
        input_path,
        load_relate_dataset.report_path_seg_types
    )
    return list(map(Path, dispatch(run_report,
                                   max_procs=max_procs,
                                   parallel=parallel,
                                   pass_n_procs=False,
                                   args=as_list_of_tuples(report_files),
                                   kwargs=dict(force=force))))


params = [arg_input_path,
          opt_max_procs,
          opt_parallel,
          opt_force]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate a FASTQ file. """
    run(*args, **kwargs)
