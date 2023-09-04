"""

SEISMIC-RNA Main Module

========================================================================

This module is the entry point for the command line interface. Running

$ seismic [OPTIONS] command [OPTIONS] [ARGS]

calls the function main_cli() defined in this module.

"""

import cProfile
import os

from click import Context, command, group, pass_context, version_option

from . import (__version__,
               demult as demultiplex_mod,
               align as align_mod,
               relate as relate_mod,
               mask as mask_mod,
               cluster as cluster_mod,
               table as table_mod,
               fold as fold_mod,
               graph as graph_mod,
               test as test_mod)
from .core import docdef, logs
from .core.cli import (merge_params, opt_demultiplex,
                       opt_verbose, opt_quiet, opt_log, opt_profile, opt_fold)
from .core.seq import DNA

all_params = merge_params([opt_demultiplex],
                          demultiplex_mod.params,
                          align_mod.params,
                          relate_mod.params,
                          mask_mod.params,
                          cluster_mod.params,
                          table_mod.params,
                          [opt_fold],
                          fold_mod.params)


@command("all", params=all_params)
def all_cli(*args, **kwargs):
    """ Run 'align', 'relate', 'mask', (optionally) 'cluster', 'table',
    (optionally) 'fold', and (optionally) 'graph', in that order. """
    return run(*args, **kwargs)


@docdef.auto()
def run(*,
        # Arguments
        input_path: tuple[str, ...],
        # General options
        out_dir: str,
        temp_dir: str,
        save_temp: bool,
        rerun: bool,
        max_procs: int,
        parallel: bool,
        # FASTQ options
        fasta: str,
        fastqs: tuple[str, ...],
        fastqi: tuple[str, ...],
        fastqp: tuple[str, ...],
        phred_enc: int,
        # Demultiplexing options
        demulti_overwrite: bool,
        demult_on: bool,
        parallel_demultiplexing: bool,
        clipped: int,
        mismatch_tolerence: int,
        index_tolerance: int,
        barcode_start: int,
        barcode_length: int,
        # Align options
        dmfastqs: tuple[str, ...],
        dmfastqi: tuple[str, ...],
        dmfastqp: tuple[str, ...],
        fastqc: bool,
        qc_extract: bool,
        cut: bool,
        cut_q1: int,
        cut_q2: int,
        cut_g1: tuple[str, ...],
        cut_a1: tuple[str, ...],
        cut_g2: tuple[str, ...],
        cut_a2: tuple[str, ...],
        cut_o: int,
        cut_e: float,
        cut_indels: bool,
        cut_nextseq: bool,
        cut_discard_trimmed: bool,
        cut_discard_untrimmed: bool,
        cut_m: int,
        bt2_local: bool,
        bt2_discordant: bool,
        bt2_mixed: bool,
        bt2_dovetail: bool,
        bt2_contain: bool,
        bt2_unal: bool,
        bt2_score_min_e2e: str,
        bt2_score_min_loc: str,
        bt2_i: int,
        bt2_x: int,
        bt2_gbar: int,
        bt2_l: int,
        bt2_s: str,
        bt2_d: int,
        bt2_r: int,
        bt2_dpad: int,
        bt2_orient: str,
        min_mapq: int,
        # Relate options
        min_phred: int,
        min_reads: int,
        ambrel: bool,
        batch_size: float,
        # Mask options
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, DNA, DNA], ...],
        primer_gap: int,
        sections_file: str,
        count_del: bool,
        count_ins: bool,
        discount_mut: tuple[str, ...],
        exclude_polya: int,
        exclude_gu: bool,
        exclude_pos: tuple[tuple[str, int], ...],
        min_finfo_read: float,
        max_fmut_read: int,
        min_mut_gap: int,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        # Cluster options
        max_clusters: int,
        min_nmut_read: int,
        em_runs: int,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        # Table options
        rels: tuple[str, ...],
        # Fold options
        fold: bool,
        quantile: float):
    """ Run entire pipeline. """
    # Demultiplexing
    if demult_on:
        for dms, dmi, dmm in demultiplex_mod.run(
                fasta=fasta,
                sections_file=sections_file,
                out_dir=out_dir,
                temp_dir=temp_dir,
                demulti_overwrite=demulti_overwrite,
                fastqp=fastqp,
                clipped=clipped,
                index_tolerance=index_tolerance,
                mismatch_tolerence=mismatch_tolerence,
                parallel_demultiplexing=parallel_demultiplexing,
                barcode_start=barcode_start,
                barcode_length=barcode_length,
                phred_enc=phred_enc):
            dmfastqs = dmfastqs + dms
            dmfastqi = dmfastqi + dmi
            dmfastqp = dmfastqp + dmm
    # Alignment
    input_path += tuple(map(str, align_mod.run(
        out_dir=out_dir,
        temp_dir=temp_dir,
        save_temp=save_temp,
        rerun=rerun,
        max_procs=max_procs,
        parallel=parallel,
        fasta=fasta,
        fastqs=fastqs,
        fastqi=fastqi,
        fastqp=fastqp,
        dmfastqs=dmfastqs,
        dmfastqi=dmfastqi,
        dmfastqp=dmfastqp,
        phred_enc=phred_enc,
        fastqc=fastqc,
        qc_extract=qc_extract,
        cut=cut,
        cut_q1=cut_q1,
        cut_q2=cut_q2,
        cut_g1=cut_g1,
        cut_a1=cut_a1,
        cut_g2=cut_g2,
        cut_a2=cut_a2,
        cut_o=cut_o,
        cut_e=cut_e,
        cut_indels=cut_indels,
        cut_nextseq=cut_nextseq,
        cut_discard_trimmed=cut_discard_trimmed,
        cut_discard_untrimmed=cut_discard_untrimmed,
        cut_m=cut_m,
        bt2_local=bt2_local,
        bt2_discordant=bt2_discordant,
        bt2_mixed=bt2_mixed,
        bt2_dovetail=bt2_dovetail,
        bt2_contain=bt2_contain,
        bt2_unal=bt2_unal,
        bt2_score_min_e2e=bt2_score_min_e2e,
        bt2_score_min_loc=bt2_score_min_loc,
        bt2_i=bt2_i,
        bt2_x=bt2_x,
        bt2_gbar=bt2_gbar,
        bt2_l=bt2_l,
        bt2_s=bt2_s,
        bt2_d=bt2_d,
        bt2_r=bt2_r,
        bt2_dpad=bt2_dpad,
        bt2_orient=bt2_orient,
        min_mapq=min_mapq,
        min_reads=min_reads,
    )))
    # Relating
    input_path += tuple(map(str, relate_mod.run(
        fasta=fasta,
        input_path=input_path,
        out_dir=out_dir,
        temp_dir=temp_dir,
        min_mapq=min_mapq,
        min_reads=min_reads,
        batch_size=batch_size,
        phred_enc=phred_enc,
        min_phred=min_phred,
        ambrel=ambrel,
        max_procs=max_procs,
        parallel=parallel,
        rerun=rerun,
        save_temp=save_temp,
    )))
    # Masking
    input_path += tuple(map(str, mask_mod.run(
        input_path=input_path,
        coords=coords,
        primers=primers,
        primer_gap=primer_gap,
        sections_file=sections_file,
        count_del=count_del,
        count_ins=count_ins,
        discount_mut=discount_mut,
        exclude_polya=exclude_polya,
        exclude_gu=exclude_gu,
        exclude_pos=exclude_pos,
        min_finfo_read=min_finfo_read,
        max_fmut_read=max_fmut_read,
        min_mut_gap=min_mut_gap,
        min_ninfo_pos=min_ninfo_pos,
        max_fmut_pos=max_fmut_pos,
        max_procs=max_procs,
        parallel=parallel,
        rerun=rerun,
    )))
    # Clustering
    input_path += tuple(map(str, cluster_mod.run(
        input_path=input_path,
        max_clusters=max_clusters,
        min_nmut_read=min_nmut_read,
        em_runs=em_runs,
        min_em_iter=min_em_iter,
        max_em_iter=max_em_iter,
        em_thresh=em_thresh,
        max_procs=max_procs,
        parallel=parallel,
        rerun=rerun,
    )))
    # Table
    input_path += tuple(map(str, table_mod.run(
        input_path=input_path,
        rels=rels,
        max_procs=max_procs,
        parallel=parallel,
        rerun=rerun,
    )))
    # Fold
    if fold:
        fold_mod.run(
            input_path=input_path,
            fasta=fasta,
            sections_file=sections_file,
            coords=coords,
            primers=primers,
            primer_gap=primer_gap,
            quantile=quantile,
            temp_dir=temp_dir,
            save_temp=save_temp,
            max_procs=max_procs,
            parallel=parallel,
            rerun=rerun,
        )
    # Graph


main_params = [
    opt_verbose,
    opt_quiet,
    opt_log,
    opt_profile,
]


# Group for main commands
@group(params=main_params, context_settings={"show_default": True})
@version_option(__version__)
@pass_context
def main_cli(ctx: Context, verbose: int, quiet: int, log: str, profile: str,
             **kwargs):
    """ SEISMIC-RNA main command line interface """
    # Configure logging.
    if log:
        log_file = os.path.abspath(log)
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
    else:
        log_file = None
    logs.config(verbose, quiet, log_file=log_file)
    # If no subcommand was given, then run the entire pipeline.
    if ctx.invoked_subcommand is None:
        if profile:
            profile_path = os.path.abspath(profile)
            # Profile the program as it runs and write results to the
            # file given in the parameter profile.
            os.makedirs(os.path.dirname(profile_path), exist_ok=True)
            cProfile.runctx("run(**kwargs)",
                            globals=globals(),
                            locals=locals(),
                            filename=profile_path,
                            sort="time")
        else:
            # Run without profiling.
            run(**kwargs)


# Add all commands to the main CLI command group.
main_cli.add_command(all_cli)
main_cli.add_command(demultiplex_mod.cli)
main_cli.add_command(align_mod.cli)
main_cli.add_command(relate_mod.cli)
main_cli.add_command(mask_mod.cli)
main_cli.add_command(cluster_mod.cli)
main_cli.add_command(table_mod.cli)
main_cli.add_command(fold_mod.cli)
main_cli.add_command(graph_mod.cli)
main_cli.add_command(test_mod.cli)


########################################################################
#                                                                      #
# ©2023, the Rouskin Lab.                                              #
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
