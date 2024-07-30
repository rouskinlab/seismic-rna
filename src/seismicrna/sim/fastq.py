import gzip
import os
from itertools import chain
from logging import getLogger
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from click import command
from numba import jit

from .relate import get_param_dir_fields, load_param_dir
from ..core import path
from ..core.arg import (ADAPTER_SEQ_ILLUMINA_3P,
                        arg_input_path,
                        opt_param_dir,
                        opt_profile_name,
                        opt_sample,
                        opt_read_length,
                        opt_paired_end,
                        opt_reverse_fraction,
                        opt_min_mut_gap,
                        opt_fq_gzip,
                        opt_num_reads,
                        opt_batch_size,
                        opt_max_procs,
                        opt_parallel,
                        opt_force)
from ..core.array import get_length
from ..core.ngs import HI_QUAL, LO_QUAL
from ..core.rel import (MATCH,
                        NOCOV,
                        SUB_A,
                        SUB_C,
                        SUB_G,
                        SUB_T,
                        DELET)
from ..core.report import SampleF
from ..core.run import run_func
from ..core.seq import DNA, BASEA, BASEC, BASEG, BASET, BASEN
from ..core.task import as_list_of_tuples, dispatch
from ..core.write import need_write, write_mode
from ..pool.data import load_relate_dataset
from ..relate.batch import QnamesBatch, RelateBatch
from ..relate.data import QnamesDataset, RelateDataset
from ..relate.report import RelateReport
from ..relate.sim import simulate_batches

logger = getLogger(__name__)

rng = np.random.default_rng()

COMMAND = __name__.split(os.path.extsep)[-1]


@jit()
def _complement(base: str):
    if base == BASEA:
        return BASET
    if base == BASEC:
        return BASEG
    if base == BASEG:
        return BASEC
    if base == BASET:
        return BASEA
    return BASEN


@jit()
def _generate_fastq_read_qual(rels: np.ndarray,
                              refseq: str,
                              read_length: int,
                              revcomp: bool,
                              hi_qual: str,
                              lo_qual: str):
    """ Generate a FASTQ line for a read. """
    ref_length, = rels.shape
    # Map each type of relationship to a type of base in the read.
    if revcomp:
        ref_pos = ref_length - 1
        ref_inc = -1
        sub_a = BASET
        sub_c = BASEG
        sub_g = BASEC
        sub_t = BASEA
    else:
        ref_pos = 0
        ref_inc = 1
        sub_a = BASEA
        sub_c = BASEC
        sub_g = BASEG
        sub_t = BASET
    # Fill the read with high-quality Gs (the default base in Illumina).
    read = np.full(read_length, BASEG)
    qual = np.full(read_length, hi_qual)
    # Add bases from the read.
    read_pos = 0
    while read_pos < read_length and 0 <= ref_pos < ref_length:
        rel = rels[ref_pos]
        if rel == MATCH:
            ref_base = refseq[ref_pos]
            read[read_pos] = _complement(ref_base) if revcomp else ref_base
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel == SUB_T:
            read[read_pos] = sub_t
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel == SUB_G:
            read[read_pos] = sub_g
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel == SUB_C:
            read[read_pos] = sub_c
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel == SUB_A:
            read[read_pos] = sub_a
            qual[read_pos] = hi_qual
            read_pos += 1
        elif rel != DELET:
            read[read_pos] = BASEN
            qual[read_pos] = lo_qual
            read_pos += 1
        ref_pos += ref_inc
    # Add the adapter to the end of the read.
    adapter_pos = 0
    adapter_length = len(ADAPTER_SEQ_ILLUMINA_3P)
    while read_pos < read_length and adapter_pos < adapter_length:
        read[read_pos] = ADAPTER_SEQ_ILLUMINA_3P[adapter_pos]
        read_pos += 1
        adapter_pos += 1
    return "".join(read), "".join(qual)


def generate_fastq_record(name: str,
                          rels: np.ndarray,
                          refseq: str,
                          read_length: int,
                          reverse: bool = False,
                          hi_qual: str = HI_QUAL,
                          lo_qual: str = LO_QUAL):
    """ Generate a FASTQ line for a read. """
    if len(refseq) < get_length(rels, "rels"):
        raise ValueError(f"Length of the reference sequence ({len(refseq)}) "
                         f"is less than that of rels ({get_length(rels)})")
    if np.any(rels == NOCOV):
        raise ValueError(f"rels contains {NOCOV}: {rels}")
    read, qual = _generate_fastq_read_qual(rels,
                                           refseq,
                                           read_length,
                                           reverse,
                                           hi_qual,
                                           lo_qual)
    return f"@{name}\n{read}\n+\n{qual}\n"


def _get_common_attr(a: Any, b: Any, attr: str):
    rval = getattr(a, attr)
    nval = getattr(b, attr)
    if rval != nval:
        raise ValueError(f"Got different values of {repr(attr)} for {a} "
                         f"({repr(rval)}) and {b} ({repr(nval)})")
    return rval


def generate_fastq(top: Path,
                   sample: str,
                   ref: str,
                   refseq: DNA,
                   paired: bool,
                   read_length: int,
                   batches: Iterable[tuple[RelateBatch, QnamesBatch]],
                   p_rev: float = 0.5,
                   fq_gzip: bool = True,
                   force: bool = False):
    """ Generate FASTQ file(s) from a dataset. """
    seq_str = str(refseq)
    if paired:
        segs = [path.DMFASTQ1_SEGS, path.DMFASTQ2_SEGS]
        exts = [path.FQ1_EXTS[0], path.FQ2_EXTS[0]]
    else:
        segs = [path.DMFASTQ_SEGS]
        exts = [path.FQ_EXTS[0]]
    if fq_gzip:
        exts = [(ext if ext.endswith(path.GZIP_EXT)
                 else f"{ext}{path.GZIP_EXT}")
                for ext in exts]
        open_func = gzip.open
    else:
        exts = [(ext if not ext.endswith(path.GZIP_EXT)
                 else ext[:-len(path.GZIP_EXT)])
                for ext in exts]
        open_func = open
    fastq_paths = [path.buildpar(*seg, top=top, sample=sample, ref=ref, ext=ext)
                   for seg, ext in zip(segs, exts, strict=True)]
    if any(need_write(fastq, force, warn=False) for fastq in fastq_paths):
        fastq_files = list()
        try:
            for fastq in fastq_paths:
                fastq_files.append(open_func(fastq, write_mode(force, fq_gzip)))
            for rbatch, nbatch in batches:
                num_reads = _get_common_attr(rbatch, nbatch, "num_reads")
                seg_end5s = np.asarray(rbatch.seg_end5s)
                seg_end3s = np.asarray(rbatch.seg_end3s)
                if rbatch.num_segments == 1:
                    reverse = rng.random(num_reads) < p_rev
                elif rbatch.num_segments == 2:
                    reverse = seg_end5s[:, 0] > seg_end5s[:, 1]
                else:
                    raise ValueError(f"Each batch must have 1 or 2 segments, "
                                     f"but got {rbatch.num_segments}")
                for rels, end5s, end3s, name, rev in zip(rbatch.matrix.values,
                                                         seg_end5s,
                                                         seg_end3s,
                                                         nbatch.names,
                                                         reverse,
                                                         strict=True):
                    for i, (fq, end5, end3) in enumerate(zip(fastq_files,
                                                             end5s,
                                                             end3s,
                                                             strict=True)):
                        record = generate_fastq_record(name,
                                                       rels[end5 - 1: end3],
                                                       seq_str[end5 - 1: end3],
                                                       read_length,
                                                       bool((rev + i) % 2))
                        fq.write(record.encode() if fq_gzip else record)
        finally:
            # Close the FASTQ files.
            for fq in fastq_files:
                try:
                    fq.close()
                except Exception as error:
                    logger.warning(f"Failed to close file {fq}: {error}")
    else:
        # Warn that the FASTQ file(s) already exist(s).
        for fastq in fastq_paths:
            need_write(fastq, force, warn=True)
    return fastq_paths


def from_report(report_file: Path, *,
                read_length: int,
                p_rev: float,
                fq_gzip: bool,
                force: bool):
    """ Simulate a FASTQ file from a Relate report. """
    report = RelateReport.load(report_file)
    sample = report.get_field(SampleF)
    rdata = RelateDataset.load(report_file)
    ndata = QnamesDataset.load(report_file)
    sim_dir = _get_common_attr(rdata, ndata, "top")
    section = rdata.section
    batches = zip(rdata.iter_batches(), ndata.iter_batches())
    return generate_fastq(sim_dir,
                          sample,
                          section.ref,
                          section.seq,
                          rdata.paired,
                          read_length,
                          batches,
                          p_rev=p_rev,
                          fq_gzip=fq_gzip,
                          force=force)


def from_param_dir(param_dir: Path, *,
                   sample: str,
                   profile: str,
                   read_length: int,
                   paired: bool,
                   p_rev: float,
                   fq_gzip: bool,
                   force: bool,
                   **kwargs):
    """ Simulate a FASTQ file from parameter files. """
    sim_dir, _, _ = get_param_dir_fields(param_dir)
    section, pmut, u5s, u3s, pends, pclust = load_param_dir(param_dir, profile)
    batches = simulate_batches(sample=sample,
                               ref=section.ref,
                               pmut=pmut,
                               uniq_end5s=u5s,
                               uniq_end3s=u3s,
                               pends=pends,
                               pclust=pclust,
                               paired=paired,
                               read_length=read_length,
                               p_rev=p_rev,
                               batch_size=opt_batch_size.default,
                               **kwargs)
    return generate_fastq(sim_dir.joinpath(path.SIM_SAMPLES_DIR),
                          sample,
                          section.ref,
                          section.seq,
                          paired,
                          read_length,
                          batches,
                          p_rev=p_rev,
                          fq_gzip=fq_gzip,
                          force=force)


@run_func(logger.critical)
def run(*,
        input_path: tuple[str, ...],
        param_dir: tuple[str, ...],
        profile_name: str,
        sample: str,
        paired_end: bool,
        read_length: int,
        reverse_fraction: float,
        min_mut_gap: int,
        fq_gzip: bool,
        num_reads: int,
        max_procs: int,
        parallel: bool,
        force: bool):
    report_files = as_list_of_tuples(path.find_files_chain(
        input_path,
        load_relate_dataset.report_path_seg_types
    ))
    param_dirs = as_list_of_tuples(map(Path, param_dir))
    fastqs = list()
    if report_files:
        fastqs.extend(chain(*dispatch(from_report,
                                      max_procs=max_procs,
                                      parallel=parallel,
                                      pass_n_procs=False,
                                      args=report_files,
                                      kwargs=dict(read_length=read_length,
                                                  p_rev=reverse_fraction,
                                                  fq_gzip=fq_gzip,
                                                  force=force))))
    if param_dirs:
        fastqs.extend(chain(*dispatch(from_param_dir,
                                      max_procs=max_procs,
                                      parallel=parallel,
                                      pass_n_procs=False,
                                      args=param_dirs,
                                      kwargs=dict(sample=sample,
                                                  profile=profile_name,
                                                  paired=paired_end,
                                                  read_length=read_length,
                                                  p_rev=reverse_fraction,
                                                  min_mut_gap=min_mut_gap,
                                                  fq_gzip=fq_gzip,
                                                  num_reads=num_reads,
                                                  force=force))))
    if not fastqs:
        logger.warning("No FASTQ files were generated")
    return fastqs


params = [arg_input_path,
          opt_param_dir,
          opt_profile_name,
          opt_sample,
          opt_paired_end,
          opt_read_length,
          opt_reverse_fraction,
          opt_min_mut_gap,
          opt_fq_gzip,
          opt_num_reads,
          opt_max_procs,
          opt_parallel,
          opt_force]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate a FASTQ file. """
    run(*args, **kwargs)
