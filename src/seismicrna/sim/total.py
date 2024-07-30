import os
from logging import getLogger
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
                        opt_batch_size,
                        opt_brotli_level,
                        opt_ct_file,
                        opt_param_dir,
                        merge_params,
                        extra_defaults)
from ..core.run import run_func
from ..core.seq import DNA

logger = getLogger(__name__)

COMMAND = __name__.split(os.path.extsep)[-1]


def as_tuple_str(items: Iterable):
    return tuple(map(str, items))


@run_func(logger.critical,
          default=None,
          extra_defaults=extra_defaults)
def run(*,
        sim_dir: str,
        tmp_pfx: str,
        sample: str,
        refs: str,
        ref: str,
        reflen: int,
        profile_name: str,
        fold_coords: tuple[tuple[str, int, int], ...],
        fold_primers: tuple[tuple[str, DNA, DNA], ...],
        fold_sections_file: str | None,
        fold_constraint: str | None,
        fold_temp: float,
        fold_md: int,
        fold_mfe: bool,
        fold_max: int,
        fold_percent: float,
        pmut_paired: tuple[tuple[str, float], ...],
        pmut_unpaired: tuple[tuple[str, float], ...],
        vmut_paired: float,
        vmut_unpaired: float,
        end3_fmean: float,
        insert_fmean: float,
        ends_var: float,
        clust_conc: float,
        paired_end: bool,
        read_length: int,
        reverse_fraction: float,
        min_mut_gap: int,
        fq_gzip: bool,
        num_reads: int,
        keep_tmp: bool,
        force: bool,
        max_procs: int,
        parallel: bool):
    """ Simulate FASTQ files from scratch. """
    fasta = str(ref_mod.run(
        sim_dir=sim_dir,
        refs=refs,
        ref=ref,
        reflen=reflen,
        force=force)
    )
    ct_file = as_tuple_str(fold_mod.run(
        fasta=fasta,
        sim_dir=sim_dir,
        tmp_pfx=tmp_pfx,
        profile_name=profile_name,
        fold_coords=fold_coords,
        fold_primers=fold_primers,
        fold_sections_file=fold_sections_file,
        fold_constraint=fold_constraint,
        fold_temp=fold_temp,
        fold_md=fold_md,
        fold_mfe=fold_mfe,
        fold_max=fold_max,
        fold_percent=fold_percent,
        keep_tmp=keep_tmp,
        force=force,
        max_procs=max_procs,
        parallel=parallel)
    )
    params_mod.run(
        ct_file=ct_file,
        pmut_paired=pmut_paired,
        pmut_unpaired=pmut_unpaired,
        vmut_paired=vmut_paired,
        vmut_unpaired=vmut_unpaired,
        end3_fmean=end3_fmean,
        insert_fmean=insert_fmean,
        ends_var=ends_var,
        clust_conc=clust_conc,
        force=force,
        parallel=parallel,
        max_procs=max_procs
    )
    param_dir = as_tuple_str(path.deduplicate(map(os.path.dirname, ct_file)))
    return fastq_mod.run(
        input_path=(),
        param_dir=param_dir,
        profile_name=profile_name,
        sample=sample,
        paired_end=paired_end,
        read_length=read_length,
        reverse_fraction=reverse_fraction,
        min_mut_gap=min_mut_gap,
        fq_gzip=fq_gzip,
        num_reads=num_reads,
        force=force,
        parallel=parallel,
        max_procs=max_procs
    )


params = merge_params(fastq_mod.params,
                      fold_mod.params,
                      params_mod.params,
                      ref_mod.params,
                      relate_mod.params,
                      exclude=[arg_fasta,
                               arg_input_path,
                               opt_batch_size,
                               opt_brotli_level,
                               opt_ct_file,
                               opt_param_dir])


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """ Simulate FASTQ files from scratch. """
    return run(*args, **kwargs)

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
