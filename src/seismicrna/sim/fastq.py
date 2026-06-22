from __future__ import annotations
import gzip
import os
from itertools import chain
from pathlib import Path
from typing import Any, Iterable

from click import command

from .idmut import (
    _get_param_dir_fields,
    _load_param_dir,
    parse_injected_mut_probs,
    parse_min_mut_gap_weights,
    set_sim_mut_params,
)
from ..core import path
from ..core.arg.cli import (
    ILLUMINA_TRUSEQ_ADAPTER_R1,
    ILLUMINA_TRUSEQ_ADAPTER_R2,
    arg_input_path,
    opt_param_dir,
    opt_profile_name,
    opt_sample_sim,
    opt_read_length,
    opt_paired_end,
    opt_reverse_fraction,
    opt_probe,
    opt_min_mut_gap_weights,
    opt_mut_collisions,
    opt_injected_mut_probs,
    opt_fq_gzip,
    opt_num_reads,
    opt_batch_size,
    opt_num_cpus,
    opt_seed,
    opt_force,
)
from ..core.array import get_length
from ..core.logs import logger
from ..core.random import get_random_integer_generator
from ..core.ngs.phred import HI_QUAL, LO_QUAL
from ..core.rel.code import NOCOV
from ..core.report import SampleF
from ..core.run import run_func
from ..core.seq.xna import DNA
from ..core.task import as_list_of_tuples, dispatch
from ..core.write import need_write, write_mode
from ..idmut.batch import ReadNamesBatch, IDmutRegionMutsBatch
from ..idmut.dataset import ReadNamesDataset, IDmutMutsDataset, load_idmut_dataset
from ..idmut.report import IDmutReport
from ..idmut.sim import simulate_batches

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

COMMAND = __name__.split(os.path.extsep)[-1]


def generate_fastq_record(
    name: str,
    rels: np.ndarray,
    refseq: str,
    adapter: str,
    read_length: int,
    reverse: bool = False,
    hi_qual: str = HI_QUAL,
    lo_qual: str = LO_QUAL,
):
    """Generate a FASTQ line for a read."""
    import numpy as np

    if len(refseq) < get_length(rels, "rels"):
        raise ValueError(
            f"Length of the reference sequence ({len(refseq)}) "
            f"is less than that of rels ({get_length(rels)})"
        )
    if np.any(rels == NOCOV):
        raise ValueError(f"rels contains {NOCOV}: {rels}")
    # Imported here (not at module level) so importing this module does not
    # import numba via fastq_jit.
    from .fastq_jit import generate_fastq_read_qual

    read, qual = generate_fastq_read_qual(
        rels, refseq, adapter, read_length, reverse, hi_qual, lo_qual
    )
    return f"@{name}\n{read}\n+\n{qual}\n"


def _get_common_attr(a: Any, b: Any, attr: str):
    """
    Return the value of an attribute that must be equal on two objects.

    Parameters
    ----------
    a: Any
        First object.
    b: Any
        Second object.
    attr: str
        Name of the attribute to compare.

    Returns
    -------
    Any
        The shared value of the attribute.

    Raises
    ------
    ValueError
        If the attribute values differ between `a` and `b`.
    """
    rval = getattr(a, attr)
    nval = getattr(b, attr)
    if rval != nval:
        raise ValueError(
            f"Got different values of {repr(attr)} for {a} "
            f"({repr(rval)}) and {b} ({repr(nval)})"
        )
    return rval


def generate_fastq(
    top: Path,
    sample: str,
    ref: str,
    refseq: DNA,
    paired: bool,
    read_length: int,
    batches: Iterable[tuple[IDmutRegionMutsBatch, ReadNamesBatch]],
    p_rev: float = 0.5,
    fq_gzip: bool = True,
    force: bool = False,
    seed: int | None = None,
):
    """Generate FASTQ file(s) from a dataset."""
    import numpy as np

    rng = np.random.default_rng(seed)
    if paired:
        segs_list = [path.DMFASTQ1_SEGS, path.DMFASTQ2_SEGS]
        exts = [path.FQ1_EXTS[0], path.FQ2_EXTS[0]]
        adapters = [ILLUMINA_TRUSEQ_ADAPTER_R1, ILLUMINA_TRUSEQ_ADAPTER_R2]
    else:
        segs_list = [path.DMFASTQ_SEGS]
        exts = [path.FQ_EXTS[0]]
        adapters = [ILLUMINA_TRUSEQ_ADAPTER_R1]
    if fq_gzip:
        exts = [
            (ext if ext.endswith(path.GZIP_EXT) else f"{ext}{path.GZIP_EXT}")
            for ext in exts
        ]
        open_func = gzip.open
    else:
        exts = [
            (ext if not ext.endswith(path.GZIP_EXT) else ext[: -len(path.GZIP_EXT)])
            for ext in exts
        ]
        open_func = open
    fastq_paths = [
        path.buildpar(
            segs,
            {
                path.TOP: top,
                path.SAMPLE: sample,
                path.STEP: path.DEMULT_STEP,
                path.BRANCHES: dict(),
                path.REF: ref,
                path.EXT: ext,
            },
        )
        for segs, ext in zip(segs_list, exts, strict=True)
    ]
    if any(need_write(fastq, force, warn=False) for fastq in fastq_paths):
        seq_str = str(refseq)
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
                    raise ValueError(
                        "Each batch must have 1 or 2 segments, "
                        f"but got {rbatch.num_segments}"
                    )
                for rels, end5s, end3s, name, rev in zip(
                    rbatch.matrix.values,
                    seg_end5s,
                    seg_end3s,
                    nbatch.names,
                    reverse,
                    strict=True,
                ):
                    for i, (fq, end5, end3, adapter) in enumerate(
                        zip(fastq_files, end5s, end3s, adapters, strict=True)
                    ):
                        record = generate_fastq_record(
                            name,
                            rels[end5 - 1 : end3],
                            seq_str[end5 - 1 : end3],
                            adapter,
                            read_length,
                            bool((rev + i) % 2),
                        )
                        fq.write(record.encode() if fq_gzip else record)
        finally:
            # Close the FASTQ files.
            for fq in fastq_files:
                try:
                    fq.close()
                except Exception as error:
                    logger.warning(error)
    else:
        # Warn that the FASTQ file(s) already exist(s).
        for fastq in fastq_paths:
            need_write(fastq, force, warn=True)
    return fastq_paths


def from_report(
    report_file: Path,
    *,
    read_length: int,
    p_rev: float,
    fq_gzip: bool,
    force: bool,
    seed: int | None,
):
    """Simulate a FASTQ file from an IDmut report."""
    report = IDmutReport.load(report_file)
    sample = report.get_field(SampleF)
    rdata = IDmutMutsDataset(report_file)
    ndata = ReadNamesDataset(report_file)
    sim_dir = _get_common_attr(rdata, ndata, "top")
    region = rdata.region
    batches = zip(rdata.iter_batches(), ndata.iter_batches())
    return generate_fastq(
        sim_dir,
        sample,
        region.ref,
        region.seq,
        rdata.paired,
        read_length,
        batches,
        p_rev=p_rev,
        fq_gzip=fq_gzip,
        force=force,
        seed=seed,
    )


def from_param_dir(
    param_dir: Path,
    *,
    sample: str,
    profile: str,
    read_length: int,
    paired: bool,
    p_rev: float,
    fq_gzip: bool,
    force: bool,
    seed: int | None,
    **kwargs,
):
    """Simulate a FASTQ file from parameter files."""
    seeds = get_random_integer_generator(seed)
    sim_dir, _, _ = _get_param_dir_fields(param_dir)
    region, pmut, u5s, u3s, pends, pclust = _load_param_dir(param_dir, profile)
    batches = simulate_batches(
        sample=sample,
        branches=dict(),
        ref=region.ref,
        pmut=pmut,
        uniq_end5s=u5s,
        uniq_end3s=u3s,
        pends=pends,
        pclust=pclust,
        paired=paired,
        read_length=read_length,
        p_rev=p_rev,
        batch_size=opt_batch_size.default,
        write_read_names=True,
        seed=next(seeds),
        **kwargs,
    )
    # Convert each IDmutBatchIO into an IDmutRegionMutsBatch, which is
    # required by generate_fastq().
    batches = [(rbatch.to_region_batch(region), nbatch) for rbatch, nbatch in batches]
    return generate_fastq(
        sim_dir.joinpath(path.SIM_SAMPLES_DIR),
        sample,
        region.ref,
        region.seq,
        paired,
        read_length,
        batches,
        p_rev=p_rev,
        fq_gzip=fq_gzip,
        force=force,
        seed=next(seeds),
    )


@run_func(COMMAND)
def run(
    *,
    input_path: Iterable[str | Path],
    param_dir: Iterable[str | Path],
    profile_name: str,
    sample: str,
    paired_end: bool,
    read_length: int,
    reverse_fraction: float,
    probe: str,
    min_mut_gap_weights: str | None,
    mut_collisions: str,
    injected_mut_probs: str | None,
    fq_gzip: bool,
    num_reads: int,
    num_cpus: int,
    force: bool,
    seed: int | None,
):
    """
    Simulate FASTQ file(s) from idmut reports or parameter directories.

    Parameters
    ----------
    input_path: Iterable[str | Path]
        Paths to idmut report files or directories containing them;
        used to generate FASTQ files from existing IDmut data.
    param_dir: Iterable[str | Path]
        Paths to simulation parameter directories; used to generate
        FASTQ files from CT/parameter files.
    profile_name: str
        Name of the mutation profile to use from the parameter directory.
    sample: str
        Sample name to embed in the output FASTQ paths.
    paired_end: bool
        Whether to simulate paired-end reads.
    read_length: int
        Length of each simulated read.
    reverse_fraction: float
        Fraction of reads where mate 1 is reverse-complemented.
    probe: str
        Probe type (e.g. DMS); used to set probe-specific defaults.
    min_mut_gap_weights: str | None
        Comma-separated gap:weight pairs for a bias mixture; None to use
        the probe default, empty string to use the probe-default single gap.
    mut_collisions: str
        How to handle reads with close mutations: "drop", "merge", or
        "auto" to select based on the probe.
    injected_mut_probs: str | None
        Comma-separated offset:prob pairs (offset ≥ 1) for injecting a
        mutation at each offset 5' of an existing mutation; None to use
        the probe default, empty string to disable injection.
    fq_gzip: bool
        Whether to gzip-compress the output FASTQ files.
    num_reads: int
        Total number of reads to simulate per param_dir run.
    num_cpus: int
        Number of CPU cores to use.
    force: bool
        Whether to overwrite existing output files.
    seed: int | None
        Random seed for reproducibility; None for no fixed seed.

    Returns
    -------
    list[Path]
        Paths of all generated FASTQ files.
    """
    seeds = get_random_integer_generator(seed)
    mut_collisions, min_mut_gap_weights, injected_mut_probs = set_sim_mut_params(
        probe, mut_collisions, min_mut_gap_weights, injected_mut_probs
    )
    min_mut_gap_weights_dict = parse_min_mut_gap_weights(min_mut_gap_weights)
    injected_mut_probs_dict = parse_injected_mut_probs(injected_mut_probs)
    report_files = as_list_of_tuples(
        path.find_files_chain(input_path, load_idmut_dataset.report_path_seg_types)
    )
    param_dirs = as_list_of_tuples(map(Path, param_dir))
    fastqs = list()
    if report_files:
        fastqs.extend(
            chain(
                *dispatch(
                    from_report,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=False,
                    ordered=False,
                    raise_on_error=False,
                    args=report_files,
                    kwargs=dict(
                        read_length=read_length,
                        p_rev=reverse_fraction,
                        fq_gzip=fq_gzip,
                        force=force,
                        seed=next(seeds),
                    ),
                )
            )
        )
    if param_dirs:
        fastqs.extend(
            chain(
                *dispatch(
                    from_param_dir,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=False,
                    ordered=False,
                    raise_on_error=False,
                    args=param_dirs,
                    kwargs=dict(
                        sample=sample,
                        profile=profile_name,
                        paired=paired_end,
                        read_length=read_length,
                        p_rev=reverse_fraction,
                        min_mut_gap_weights=min_mut_gap_weights_dict,
                        injected_mut_probs=injected_mut_probs_dict,
                        mut_collisions=mut_collisions,
                        fq_gzip=fq_gzip,
                        num_reads=num_reads,
                        force=force,
                        seed=next(seeds),
                    ),
                )
            )
        )
    if not fastqs:
        logger.warning("No FASTQ files were generated")
    return fastqs


params = [
    arg_input_path,
    opt_param_dir,
    opt_profile_name,
    opt_sample_sim,
    opt_paired_end,
    opt_read_length,
    opt_reverse_fraction,
    opt_probe,
    opt_min_mut_gap_weights,
    opt_mut_collisions,
    opt_injected_mut_probs,
    opt_fq_gzip,
    opt_num_reads,
    opt_num_cpus,
    opt_force,
    opt_seed,
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """Simulate a FASTQ file."""
    run(*args, **kwargs)
