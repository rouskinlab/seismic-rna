from logging import getLogger
from typing import Iterable

from click import command

from . import (demult as demultiplex_mod,
               align as align_mod,
               relate as relate_mod,
               mask as mask_mod,
               cluster as cluster_mod,
               table as table_mod,
               fold as fold_mod,
               export as export_mod)
from .core.arg import (CMD_WORKFLOW,
                       merge_params,
                       opt_demultiplex,
                       opt_cluster,
                       opt_fold,
                       opt_export,
                       opt_cgroup,
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
                       opt_graph_aucroll,
                       extra_defaults)
from .core.run import run_func
from .core.seq import DNA
from .graph.aucroll import RollingAUCRunner
from .graph.giniroll import RollingGiniRunner
from .graph.histread import ReadHistogramRunner
from .graph.profile import ProfileRunner
from .graph.roc import ROCRunner
from .table.base import (DELET_REL,
                         INSRT_REL,
                         MUTAT_REL,
                         REL_NAMES,
                         SUB_A_REL,
                         SUB_C_REL,
                         SUB_G_REL,
                         SUB_T_REL,
                         UNAMB_REL)

logger = getLogger(__name__)

MUTAT_RELS = "".join(REL_NAMES[code] for code in [SUB_A_REL,
                                                  SUB_C_REL,
                                                  SUB_G_REL,
                                                  SUB_T_REL,
                                                  DELET_REL,
                                                  INSRT_REL])


def as_tuple_str(items: Iterable):
    return tuple(map(str, items))


@run_func(logger.critical,
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
        parallel: bool,
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
        f1r2_plus: bool,
        minus_label: str,
        # Relate options
        min_phred: int,
        min_reads: int,
        ambindel: bool,
        overhangs: bool,
        clip_end5: int,
        clip_end3: int,
        batch_size: int,
        # Mask options
        mask_coords: tuple[tuple[str, int, int], ...],
        mask_primers: tuple[tuple[str, DNA, DNA], ...],
        primer_gap: int,
        mask_sections_file: str | None,
        mask_del: bool,
        mask_ins: bool,
        mask_mut: tuple[str, ...],
        mask_polya: int,
        mask_gu: bool,
        mask_pos_file: str | None,
        mask_pos: tuple[tuple[str, int], ...],
        mask_discontig: bool,
        min_ncov_read: int,
        min_finfo_read: float,
        max_fmut_read: int,
        min_mut_gap: int,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        quick_unbias: bool,
        quick_unbias_thresh: float,
        # Cluster options
        cluster: bool,
        min_clusters: int,
        max_clusters: int,
        em_runs: int,
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
        # Table options
        table_pos: bool,
        table_read: bool,
        table_clust: bool,
        # Fold options
        fold: bool,
        fold_coords: tuple[tuple[str, int, int], ...],
        fold_primers: tuple[tuple[str, DNA, DNA], ...],
        fold_sections_file: str | None,
        fold_full: bool,
        quantile: float,
        fold_temp: float,
        fold_constraint: str | None,
        fold_md: int,
        fold_mfe: bool,
        fold_max: int,
        fold_percent: float,
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
        parallel=parallel,
        fasta=fasta,
        fastqz=fastqz,
        fastqy=fastqy,
        fastqx=fastqx,
        dmfastqz=dmfastqz,
        dmfastqy=dmfastqy,
        dmfastqx=dmfastqx,
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
        f1r2_plus=f1r2_plus,
        minus_label=minus_label,
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
        ambindel=ambindel,
        overhangs=overhangs,
        clip_end5=clip_end5,
        clip_end3=clip_end3,
        sep_strands=sep_strands,
        minus_label=minus_label,
        max_procs=max_procs,
        parallel=parallel,
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
        mask_sections_file=mask_sections_file,
        mask_del=mask_del,
        mask_ins=mask_ins,
        mask_mut=mask_mut,
        mask_polya=mask_polya,
        mask_gu=mask_gu,
        mask_pos_file=mask_pos_file,
        mask_pos=mask_pos,
        mask_discontig=mask_discontig,
        min_ncov_read=min_ncov_read,
        min_finfo_read=min_finfo_read,
        max_fmut_read=max_fmut_read,
        min_mut_gap=min_mut_gap,
        min_ninfo_pos=min_ninfo_pos,
        max_fmut_pos=max_fmut_pos,
        quick_unbias=quick_unbias,
        quick_unbias_thresh=quick_unbias_thresh,
        brotli_level=brotli_level,
        max_procs=max_procs,
        parallel=parallel,
        force=force,
    ))
    # Cluster
    if cluster:
        input_path += as_tuple_str(cluster_mod.run(
            input_path=input_path,
            tmp_pfx=tmp_pfx,
            min_clusters=min_clusters,
            max_clusters=max_clusters,
            em_runs=em_runs,
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
            brotli_level=brotli_level,
            max_procs=max_procs,
            parallel=parallel,
            force=force,
        ))
    # Table
    input_path += as_tuple_str(table_mod.run(
        input_path=input_path,
        table_pos=table_pos,
        table_read=table_read,
        table_clust=table_clust,
        max_procs=max_procs,
        parallel=parallel,
        force=force,
    ))
    # Fold
    if fold:
        fold_mod.run(
            input_path=input_path,
            fold_sections_file=fold_sections_file,
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
            parallel=parallel,
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
                          parallel=parallel,
                          force=force)
    if graph_ncov:
        # Graph information per position.
        ProfileRunner.run(input_path=input_path,
                          rels=(REL_NAMES[UNAMB_REL],),
                          use_ratio=False,
                          quantile=0.,
                          cgroup=cgroup,
                          csv=csv,
                          html=html,
                          svg=svg,
                          pdf=pdf,
                          png=png,
                          max_procs=max_procs,
                          parallel=parallel,
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
                                parallel=parallel,
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
                              parallel=parallel,
                              force=force)
    if fold:
        if graph_roc:
            # Graph ROC curves.
            ROCRunner.run(input_path=input_path,
                          rels=(REL_NAMES[MUTAT_REL],),
                          use_ratio=True,
                          quantile=0.,
                          struct_file=struct_file,
                          fold_sections_file=fold_sections_file,
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
                          parallel=parallel,
                          force=force)
        if graph_aucroll:
            # Graph rolling AUC-ROC.
            RollingAUCRunner.run(input_path=input_path,
                                 rels=(REL_NAMES[MUTAT_REL],),
                                 use_ratio=True,
                                 quantile=0.,
                                 struct_file=struct_file,
                                 fold_sections_file=fold_sections_file,
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
                                 parallel=parallel,
                                 force=force)
    # Export
    if export:
        export_mod.run(
            input_path=input_path,
            samples_meta=samples_meta,
            refs_meta=refs_meta,
            all_pos=all_pos,
            max_procs=max_procs,
            parallel=parallel,
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
                      table_mod.params,
                      [opt_fold],
                      fold_mod.params,
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
