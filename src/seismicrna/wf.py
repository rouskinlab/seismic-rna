from typing import Iterable

from click import command

from . import (demult as demultiplex_mod,
               align as align_mod,
               relate as relate_mod,
               mask as mask_mod,
               cluster as cluster_mod,
               fold as fold_mod,
               draw as draw_mod,
               export as export_mod)
from .core.arg import (CMD_WORKFLOW,
                       merge_params,
                       opt_demultiplex,
                       opt_cluster,
                       opt_min_clusters,
                       opt_max_clusters,
                       opt_fold,
                       opt_export,
                       opt_cgroup,
                       opt_hist_bins,
                       opt_hist_margin,
                       opt_struct_file,
                       opt_fold_full,
                       opt_draw,
                       opt_window,
                       opt_winmin,
                       opt_csv,
                       opt_html,
                       opt_svg,
                       opt_pdf,
                       opt_png,
                       opt_graph_mprof,
                       opt_graph_tmprof,
                       opt_graph_ncov,
                       opt_graph_mhist,
                       opt_graph_giniroll,
                       opt_graph_roc,
                       opt_graph_aucroll,
                       extra_defaults)
from .core.run import run_func
from .core.seq import DNA
from .core.table import (DELET_REL,
                         INSRT_REL,
                         MUTAT_REL,
                         REL_NAMES,
                         SUB_A_REL,
                         SUB_C_REL,
                         SUB_G_REL,
                         SUB_T_REL,
                         INFOR_REL)
from .graph.aucroll import RollingAUCRunner
from .graph.giniroll import RollingGiniRunner
from .graph.histread import ReadHistogramRunner
from .graph.profile import ProfileRunner
from .graph.roc import ROCRunner

MUTAT_RELS = "".join(REL_NAMES[code] for code in [SUB_A_REL,
                                                  SUB_C_REL,
                                                  SUB_G_REL,
                                                  SUB_T_REL,
                                                  DELET_REL,
                                                  INSRT_REL])


def as_tuple_str(items: Iterable):
    return tuple(map(str, items))


@run_func(CMD_WORKFLOW,
          default=None,
          extra_defaults=extra_defaults)
def run(fasta: str,
        input_path: tuple[str, ...], *,
        # General options
        out_dir: str,
        tmp_pfx: str,
        keep_tmp: bool,
        brotli_level: int,
        force: bool,
        max_procs: int,
        # FASTQ options
        fastqz: tuple[str, ...],
        fastqy: tuple[str, ...],
        fastqx: tuple[str, ...],
        phred_enc: int,
        # Demultiplexing options
        demulti_overwrite: bool,
        demult_on: bool,
        parallel_demultiplexing: bool,
        clipped: int,
        mismatch_tolerence: int,
        index_tolerance: int,
        barcode_start: int,
        barcode_end: int,
        # Align options
        dmfastqz: tuple[str, ...],
        dmfastqy: tuple[str, ...],
        dmfastqx: tuple[str, ...],
        fastp: bool,
        fastp_5: bool,
        fastp_3: bool,
        fastp_w: int,
        fastp_m: int,
        fastp_poly_g: str,
        fastp_poly_g_min_len: int,
        fastp_poly_x: bool,
        fastp_poly_x_min_len: int,
        fastp_adapter_trimming: bool,
        fastp_adapter_1: str,
        fastp_adapter_2: str,
        fastp_adapter_fasta: str | None,
        fastp_detect_adapter_for_pe: bool,
        fastp_min_length: int,
        bt2_local: bool,
        bt2_discordant: bool,
        bt2_mixed: bool,
        bt2_dovetail: bool,
        bt2_contain: bool,
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
        bt2_un: bool,
        min_mapq: int,
        sep_strands: bool,
        f1r2_fwd: bool,
        rev_label: str,
        # Relate options
        min_phred: int,
        min_reads: int,
        insert3: bool,
        ambindel: bool,
        overhangs: bool,
        clip_end5: int,
        clip_end3: int,
        batch_size: int,
        relate_pos_table: bool,
        relate_read_table: bool,
        relate_cx: bool,
        # Mask options
        mask_coords: tuple[tuple[str, int, int], ...],
        mask_primers: tuple[tuple[str, DNA, DNA], ...],
        primer_gap: int,
        mask_regions_file: str | None,
        mask_del: bool,
        mask_ins: bool,
        mask_mut: tuple[str, ...],
        mask_polya: int,
        mask_gu: bool,
        mask_pos: tuple[tuple[str, int], ...],
        mask_pos_file: str | None,
        mask_read: tuple[str, ...],
        mask_read_file: str | None,
        mask_discontig: bool,
        min_ncov_read: int,
        min_finfo_read: float,
        max_fmut_read: float,
        min_mut_gap: int,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        quick_unbias: bool,
        quick_unbias_thresh: float,
        max_mask_iter: int,
        mask_pos_table: bool,
        mask_read_table: bool,
        # Cluster options
        cluster: bool,
        min_clusters: int,
        max_clusters: int,
        em_runs: int,
        jackpot: bool,
        jackpot_conf_level: float,
        max_jackpot_quotient: float,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        min_nrmsd_run: float,
        max_pearson_run: float,
        max_loglike_vs_best: float,
        min_pearson_vs_best: float,
        max_nrmsd_vs_best: float,
        try_all_ks: bool,
        write_all_ks: bool,
        cluster_pos_table: bool,
        cluster_abundance_table: bool,
        verify_times: bool,
        # Fold options
        fold: bool,
        fold_coords: tuple[tuple[str, int, int], ...],
        fold_primers: tuple[tuple[str, DNA, DNA], ...],
        fold_regions_file: str | None,
        fold_full: bool,
        quantile: float,
        fold_temp: float,
        fold_constraint: str | None,
        fold_md: int,
        fold_mfe: bool,
        fold_max: int,
        fold_percent: float,
        # Draw options,
        draw: bool,
        struct_num: tuple[int, ...],
        color: bool,
        # Export options,
        export: bool,
        samples_meta: str,
        refs_meta: str,
        all_pos: bool,
        # Graph options
        cgroup: str,
        hist_bins: int,
        hist_margin: float,
        struct_file: tuple[str, ...],
        window: int,
        winmin: int,
        csv: bool,
        html: bool,
        svg: bool,
        pdf: bool,
        png: bool,
        graph_mprof: bool,
        graph_tmprof: bool,
        graph_ncov: bool,
        graph_mhist: bool,
        graph_giniroll: bool,
        graph_roc: bool,
        graph_aucroll: bool):
    """ Run the entire workflow. """
    # Demultiplex
    if demult_on:
        for dms, dmi, dmm in demultiplex_mod.run_dm(
                fasta=fasta,
                refs_meta=refs_meta,
                out_dir=out_dir,
                tmp_pfx=tmp_pfx,
                demulti_overwrite=demulti_overwrite,
                fastqx=fastqx,
                clipped=clipped,
                index_tolerance=index_tolerance,
                mismatch_tolerence=mismatch_tolerence,
                parallel_demultiplexing=parallel_demultiplexing,
                barcode_start=barcode_start,
                barcode_end=barcode_end,
                phred_enc=phred_enc,
                keep_tmp=keep_tmp):
            dmfastqz = dmfastqz + dms
            dmfastqy = dmfastqy + dmi
            dmfastqx = dmfastqx + dmm
        # Clear the input FASTQ files once the demultiplexed FASTQ files
        # have been generated.
        fastqx = tuple()
    # Align
    input_path += as_tuple_str(align_mod.run(
        out_dir=out_dir,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        force=force,
        max_procs=max_procs,
        fasta=fasta,
        fastqz=fastqz,
        fastqy=fastqy,
        fastqx=fastqx,
        dmfastqz=dmfastqz,
        dmfastqy=dmfastqy,
        dmfastqx=dmfastqx,
        phred_enc=phred_enc,
        fastp=fastp,
        fastp_5=fastp_5,
        fastp_3=fastp_3,
        fastp_w=fastp_w,
        fastp_m=fastp_m,
        fastp_poly_g=fastp_poly_g,
        fastp_poly_g_min_len=fastp_poly_g_min_len,
        fastp_poly_x=fastp_poly_x,
        fastp_poly_x_min_len=fastp_poly_x_min_len,
        fastp_adapter_trimming=fastp_adapter_trimming,
        fastp_adapter_1=fastp_adapter_1,
        fastp_adapter_2=fastp_adapter_2,
        fastp_adapter_fasta=fastp_adapter_fasta,
        fastp_detect_adapter_for_pe=fastp_detect_adapter_for_pe,
        fastp_min_length=fastp_min_length,
        bt2_local=bt2_local,
        bt2_discordant=bt2_discordant,
        bt2_mixed=bt2_mixed,
        bt2_dovetail=bt2_dovetail,
        bt2_contain=bt2_contain,
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
        bt2_un=bt2_un,
        min_mapq=min_mapq,
        min_reads=min_reads,
        sep_strands=sep_strands,
        f1r2_fwd=f1r2_fwd,
        rev_label=rev_label,
    ))
    # Relate
    input_path += as_tuple_str(relate_mod.run(
        fasta=fasta,
        input_path=input_path,
        out_dir=out_dir,
        tmp_pfx=tmp_pfx,
        min_mapq=min_mapq,
        min_reads=min_reads,
        batch_size=batch_size,
        phred_enc=phred_enc,
        min_phred=min_phred,
        insert3=insert3,
        ambindel=ambindel,
        overhangs=overhangs,
        clip_end5=clip_end5,
        clip_end3=clip_end3,
        sep_strands=sep_strands,
        rev_label=rev_label,
        relate_pos_table=relate_pos_table,
        relate_read_table=relate_read_table,
        relate_cx=relate_cx,
        max_procs=max_procs,
        brotli_level=brotli_level,
        force=force,
        keep_tmp=keep_tmp,
    ))
    # Mask
    input_path += as_tuple_str(mask_mod.run(
        input_path=input_path,
        tmp_pfx=tmp_pfx,
        mask_coords=mask_coords,
        mask_primers=mask_primers,
        primer_gap=primer_gap,
        mask_regions_file=mask_regions_file,
        mask_del=mask_del,
        mask_ins=mask_ins,
        mask_mut=mask_mut,
        mask_polya=mask_polya,
        mask_gu=mask_gu,
        mask_pos=mask_pos,
        mask_pos_file=mask_pos_file,
        mask_read=mask_read,
        mask_read_file=mask_read_file,
        mask_discontig=mask_discontig,
        min_ncov_read=min_ncov_read,
        min_finfo_read=min_finfo_read,
        max_fmut_read=max_fmut_read,
        min_mut_gap=min_mut_gap,
        min_ninfo_pos=min_ninfo_pos,
        max_fmut_pos=max_fmut_pos,
        quick_unbias=quick_unbias,
        quick_unbias_thresh=quick_unbias_thresh,
        max_mask_iter=max_mask_iter,
        mask_pos_table=mask_pos_table,
        mask_read_table=mask_read_table,
        brotli_level=brotli_level,
        max_procs=max_procs,
        force=force,
    ))
    # Cluster
    if (cluster
            or min_clusters != opt_min_clusters.default
            or max_clusters != opt_max_clusters.default):
        input_path += as_tuple_str(cluster_mod.run(
            input_path=input_path,
            tmp_pfx=tmp_pfx,
            min_clusters=min_clusters,
            max_clusters=max_clusters,
            em_runs=em_runs,
            jackpot=jackpot,
            jackpot_conf_level=jackpot_conf_level,
            max_jackpot_quotient=max_jackpot_quotient,
            min_em_iter=min_em_iter,
            max_em_iter=max_em_iter,
            em_thresh=em_thresh,
            min_nrmsd_run=min_nrmsd_run,
            max_pearson_run=max_pearson_run,
            max_loglike_vs_best=max_loglike_vs_best,
            min_pearson_vs_best=min_pearson_vs_best,
            max_nrmsd_vs_best=max_nrmsd_vs_best,
            try_all_ks=try_all_ks,
            write_all_ks=write_all_ks,
            cluster_pos_table=cluster_pos_table,
            cluster_abundance_table=cluster_abundance_table,
            verify_times=verify_times,
            brotli_level=brotli_level,
            max_procs=max_procs,
            force=force,
        ))
    # Fold
    if fold:
        input_path += as_tuple_str(fold_mod.run(
            input_path=input_path,
            fold_regions_file=fold_regions_file,
            fold_coords=fold_coords,
            fold_primers=fold_primers,
            fold_full=fold_full,
            quantile=quantile,
            fold_temp=fold_temp,
            fold_constraint=fold_constraint,
            fold_md=fold_md,
            fold_mfe=fold_mfe,
            fold_max=fold_max,
            fold_percent=fold_percent,
            tmp_pfx=tmp_pfx,
            keep_tmp=keep_tmp,
            max_procs=max_procs,
            force=force,
        ))
    # Draw
    if draw:
        draw_mod.run(
            input_path=input_path,
            struct_num=struct_num,
            color=color,
            tmp_pfx=tmp_pfx,
            keep_tmp=keep_tmp,
            max_procs=max_procs,
            force=force,
        )
    if graph_mprof or graph_tmprof:
        rels = ()
        if graph_mprof:
            rels += REL_NAMES[MUTAT_REL],
        if graph_tmprof:
            rels += MUTAT_RELS,
        # Graph mutational profiles.
        ProfileRunner.run(input_path=input_path,
                          rels=(REL_NAMES[MUTAT_REL], MUTAT_RELS),
                          use_ratio=True,
                          quantile=0.,
                          cgroup=cgroup,
                          csv=csv,
                          html=html,
                          svg=svg,
                          pdf=pdf,
                          png=png,
                          max_procs=max_procs,
                          force=force)
    if graph_ncov:
        # Graph information per position.
        ProfileRunner.run(input_path=input_path,
                          rels=(REL_NAMES[INFOR_REL],),
                          use_ratio=False,
                          quantile=0.,
                          cgroup=cgroup,
                          csv=csv,
                          html=html,
                          svg=svg,
                          pdf=pdf,
                          png=png,
                          max_procs=max_procs,
                          force=force)
    if graph_mhist:
        # Graph mutations per read.
        ReadHistogramRunner.run(input_path=input_path,
                                rels=(REL_NAMES[MUTAT_REL],),
                                use_ratio=False,
                                quantile=0.,
                                cgroup=cgroup,
                                hist_bins=hist_bins,
                                hist_margin=hist_margin,
                                csv=csv,
                                html=html,
                                svg=svg,
                                pdf=pdf,
                                png=png,
                                max_procs=max_procs,
                                force=force)
    if graph_giniroll:
        # Graph Gini coefficient.
        RollingGiniRunner.run(input_path=input_path,
                              rels=(REL_NAMES[MUTAT_REL],),
                              use_ratio=True,
                              quantile=0.,
                              window=window,
                              winmin=winmin,
                              cgroup=cgroup,
                              csv=csv,
                              html=html,
                              svg=svg,
                              pdf=pdf,
                              png=png,
                              max_procs=max_procs,
                              force=force)
    if fold:
        if graph_roc:
            # Graph ROC curves.
            ROCRunner.run(input_path=input_path,
                          rels=(REL_NAMES[MUTAT_REL],),
                          use_ratio=True,
                          quantile=0.,
                          struct_file=struct_file,
                          fold_regions_file=fold_regions_file,
                          fold_coords=fold_coords,
                          fold_primers=fold_primers,
                          fold_full=fold_full,
                          cgroup=cgroup,
                          csv=csv,
                          html=html,
                          svg=svg,
                          pdf=pdf,
                          png=png,
                          max_procs=max_procs,
                          force=force)
        if graph_aucroll:
            # Graph rolling AUC-ROC.
            RollingAUCRunner.run(input_path=input_path,
                                 rels=(REL_NAMES[MUTAT_REL],),
                                 use_ratio=True,
                                 quantile=0.,
                                 struct_file=struct_file,
                                 fold_regions_file=fold_regions_file,
                                 fold_coords=fold_coords,
                                 fold_primers=fold_primers,
                                 fold_full=fold_full,
                                 window=window,
                                 winmin=winmin,
                                 cgroup=cgroup,
                                 csv=csv,
                                 html=html,
                                 svg=svg,
                                 pdf=pdf,
                                 png=png,
                                 max_procs=max_procs,
                                 force=force)
    # Export
    if export:
        export_mod.run(
            input_path=input_path,
            samples_meta=samples_meta,
            refs_meta=refs_meta,
            all_pos=all_pos,
            max_procs=max_procs,
            force=force,
        )


graph_options = [opt_cgroup,
                 opt_hist_bins,
                 opt_hist_margin,
                 opt_struct_file,
                 opt_fold_full,
                 opt_window,
                 opt_winmin,
                 opt_csv,
                 opt_html,
                 opt_svg,
                 opt_pdf,
                 opt_png,
                 opt_graph_mprof,
                 opt_graph_tmprof,
                 opt_graph_ncov,
                 opt_graph_mhist,
                 opt_graph_giniroll,
                 opt_graph_roc,
                 opt_graph_aucroll]

params = merge_params([opt_demultiplex],
                      demultiplex_mod.params,
                      align_mod.params,
                      relate_mod.params,
                      mask_mod.params,
                      [opt_cluster],
                      cluster_mod.params,
                      [opt_fold],
                      fold_mod.params,
                      [opt_draw],
                      draw_mod.params,
                      [opt_export],
                      export_mod.params,
                      graph_options)


@command(CMD_WORKFLOW, params=params)
def cli(*args, **kwargs):
    """ Run the entire workflow. """
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
