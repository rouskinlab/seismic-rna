import os
from pathlib import Path
from typing import Iterable

from click import command

from . import (fastq as fastq_mod,
               fold as fold_mod,
               params as params_mod,
               ref as ref_mod,
               relate as relate_mod)
from ..core import path
from ..core.arg import (arg_fasta,
                        arg_input_path,
                        opt_branch,
                        opt_batch_size,
                        opt_brotli_level,
                        opt_ct_file,
                        opt_param_dir,
                        opt_write_read_names,
                        merge_params,
                        extra_defaults)
from ..core.run import run_func
from ..core.seq import DNA

COMMAND = __name__.split(os.path.extsep)[-1]


@run_func(COMMAND,
          default=None,
          extra_defaults=extra_defaults)
def run(*,
        sim_dir: str | Path,
        tmp_pfx: str | Path,
        sample: str,
        refs: str,
        ref: str,
        reflen: int,
        profile_name: str,
        fold_coords: Iterable[tuple[str, int, int]],
        fold_primers: Iterable[tuple[str, DNA, DNA]],
        fold_regions_file: str | None,
        fold_constraint: str | None,
        fold_temp: float,
        fold_md: int,
        fold_mfe: bool,
        fold_max: int,
        fold_percent: float,
        pmut_paired: Iterable[tuple[str, float]],
        pmut_unpaired: Iterable[tuple[str, float]],
        vmut_paired: float,
        vmut_unpaired: float,
        center_fmean: float,
        center_fvar: float,
        length_fmean: float,
        length_fvar: float,
        clust_conc: float,
        paired_end: bool,
        read_length: int,
        reverse_fraction: float,
        min_mut_gap: int,
        fq_gzip: bool,
        num_reads: int,
        keep_tmp: bool,
        force: bool,
        num_cpus: int):
    """ Simulate FASTQ files from scratch. """
    fasta = str(ref_mod.run(
        sim_dir=sim_dir,
        refs=refs,
        ref=ref,
        reflen=reflen,
        force=force)
    )
    ct_file = fold_mod.run(
        fasta=fasta,
        sim_dir=sim_dir,
        tmp_pfx=tmp_pfx,
        profile_name=profile_name,
        fold_coords=fold_coords,
        fold_primers=fold_primers,
        fold_regions_file=fold_regions_file,
        fold_constraint=fold_constraint,
        fold_temp=fold_temp,
        fold_md=fold_md,
        fold_mfe=fold_mfe,
        fold_max=fold_max,
        fold_percent=fold_percent,
        keep_tmp=keep_tmp,
        force=force,
        num_cpus=num_cpus
    )
    params_mod.run(
        ct_file=ct_file,
        pmut_paired=pmut_paired,
        pmut_unpaired=pmut_unpaired,
        vmut_paired=vmut_paired,
        vmut_unpaired=vmut_unpaired,
        center_fmean=center_fmean,
        center_fvar=center_fvar,
        length_fmean=length_fmean,
        length_fvar=length_fvar,
        clust_conc=clust_conc,
        force=force,
        num_cpus=num_cpus
    )
    return fastq_mod.run(
        input_path=(),
        param_dir=path.deduplicate(f.parent for f in ct_file),
        profile_name=profile_name,
        sample=sample,
        paired_end=paired_end,
        read_length=read_length,
        reverse_fraction=reverse_fraction,
        min_mut_gap=min_mut_gap,
        fq_gzip=fq_gzip,
        num_reads=num_reads,
        force=force,
        num_cpus=num_cpus
    )


params = merge_params(fastq_mod.params,
                      fold_mod.params,
                      params_mod.params,
                      ref_mod.params,
                      relate_mod.params,
                      exclude=[arg_fasta,
                               arg_input_path,
                               opt_branch,
                               opt_batch_size,
                               opt_brotli_level,
                               opt_ct_file,
                               opt_param_dir,
                               opt_write_read_names])


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate FASTQ files from scratch. """
    return run(*args, **kwargs)
