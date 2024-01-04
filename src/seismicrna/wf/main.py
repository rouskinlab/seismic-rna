"""

Whole Pipeline Main Module
========================================================================

"""

from typing import Iterable

from click import command

from .. import (demult as demultiplex_mod,
                align as align_mod,
                relate as relate_mod,
                pool as pool_mod,
                mask as mask_mod,
                cluster as cluster_mod,
                table as table_mod,
                fold as fold_mod,
                export as export_mod)
from ..core.arg import (CMD_WORKFLOW,
                        docdef,
                        merge_params,
                        opt_demultiplex,
                        opt_fold,
                        opt_export,
                        opt_cgroup,
                        opt_hist_bins,
                        opt_hist_margin,
                        opt_window,
                        opt_winmin,
                        opt_csv,
                        opt_html,
                        opt_pdf)
from ..core.seq import DNA
from ..graph.aucroll import RollingAUCRunner
from ..graph.histread import ReadHistogramRunner
from ..graph.profile import ProfileRunner
from ..graph.roc import ROCRunner

graph_options = [opt_cgroup,
                 opt_hist_bins,
                 opt_hist_margin,
                 opt_window,
                 opt_winmin,
                 opt_csv,
                 opt_html,
                 opt_pdf]

params = merge_params([opt_demultiplex],
                      demultiplex_mod.params,
                      align_mod.params,
                      relate_mod.params,
                      pool_mod.params,
                      mask_mod.params,
                      cluster_mod.params,
                      table_mod.params,
                      [opt_fold],
                      fold_mod.params,
                      [opt_export],
                      export_mod.params,
                      graph_options)


def as_tuple_str(items: Iterable):
    return tuple(map(str, items))


@command(CMD_WORKFLOW, params=params)
def cli(*args, **kwargs):
    """ Run the entire workflow. """
    return run(*args, **kwargs)


@docdef.auto()
def run(*,
        # Arguments
        input_path: tuple[str, ...],
        # General options
        out_dir: str,
        temp_dir: str,
        keep_temp: bool,
        brotli_level: int,
        force: bool,
        max_procs: int,
        parallel: bool,
        # FASTQ options
        fasta: str,
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
        cram: bool,
        # Relate options
        min_phred: int,
        min_reads: int,
        ambrel: bool,
        batch_size: float,
        # Pool options
        pool: str,
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
        em_runs: int,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        # Table options
        table_pos: bool,
        table_read: bool,
        table_clust: bool,
        # Fold options
        fold: bool,
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
        window: int,
        winmin: int,
        csv: bool,
        html: bool,
        pdf: bool):
    """ Run the entire workflow. """
    # Demultiplex
    if demult_on:
        for dms, dmi, dmm in demultiplex_mod.run(
                fasta=fasta,
                refs_meta=refs_meta,
                out_dir=out_dir,
                temp_dir=temp_dir,
                demulti_overwrite=demulti_overwrite,
                fastqx=fastqx,
                clipped=clipped,
                index_tolerance=index_tolerance,
                mismatch_tolerence=mismatch_tolerence,
                parallel_demultiplexing=parallel_demultiplexing,
                barcode_start=barcode_start,
                barcode_end=barcode_end,
                phred_enc=phred_enc,
                keep_temp=keep_temp):
            dmfastqz = dmfastqz + dms
            dmfastqy = dmfastqy + dmi
            dmfastqx = dmfastqx + dmm
        # Clear the input FASTQ files once the demultiplexed FASTQ files
        # have been generated.
        fastqx = tuple()
        fastqy = tuple()
        fastqz = tuple()
    # Align
    input_path += as_tuple_str(align_mod.run(
        out_dir=out_dir,
        temp_dir=temp_dir,
        keep_temp=keep_temp,
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
        cram=cram,
    ))
    # Relate
    input_path += as_tuple_str(relate_mod.run(
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
        brotli_level=brotli_level,
        force=force,
        keep_temp=keep_temp,
    ))
    # Pool
    input_path += as_tuple_str(pool_mod.run(
        input_path=input_path,
        pool=pool,
        max_procs=max_procs,
        parallel=parallel,
        force=force,
    ))
    # Mask
    input_path += as_tuple_str(mask_mod.run(
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
        brotli_level=brotli_level,
        max_procs=max_procs,
        parallel=parallel,
        force=force,
    ))
    # Cluster
    input_path += as_tuple_str(cluster_mod.run(
        input_path=input_path,
        max_clusters=max_clusters,
        em_runs=em_runs,
        min_em_iter=min_em_iter,
        max_em_iter=max_em_iter,
        em_thresh=em_thresh,
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
            sections_file=sections_file,
            coords=coords,
            primers=primers,
            primer_gap=primer_gap,
            quantile=quantile,
            fold_temp=fold_temp,
            fold_constraint=fold_constraint,
            fold_md=fold_md,
            fold_mfe=fold_mfe,
            fold_max=fold_max,
            fold_percent=fold_percent,
            temp_dir=temp_dir,
            keep_temp=keep_temp,
            max_procs=max_procs,
            parallel=parallel,
            force=force,
        )
    # Graph mutational profiles.
    ProfileRunner.run(input_path=input_path,
                      rels=("m", "acgtdi"),
                      use_ratio=True,
                      quantile=0.,
                      cgroup=cgroup,
                      out_dir=out_dir,
                      csv=csv,
                      html=html,
                      pdf=pdf,
                      max_procs=max_procs,
                      parallel=parallel,
                      force=force)
    # Graph information per position.
    ProfileRunner.run(input_path=input_path,
                      rels=("n",),
                      use_ratio=False,
                      quantile=0.,
                      cgroup=cgroup,
                      out_dir=out_dir,
                      csv=csv,
                      html=html,
                      pdf=pdf,
                      max_procs=max_procs,
                      parallel=parallel,
                      force=force)
    # Graph mutations per read.
    ReadHistogramRunner.run(input_path=input_path,
                            rels=("m",),
                            use_ratio=False,
                            quantile=0.,
                            cgroup=cgroup,
                            hist_bins=hist_bins,
                            hist_margin=hist_margin,
                            out_dir=out_dir,
                            csv=csv,
                            html=html,
                            pdf=pdf,
                            max_procs=max_procs,
                            parallel=parallel,
                            force=force)
    if fold:
        # Graph ROC curves.
        ROCRunner.run(input_path=input_path,
                      rels=("m",),
                      use_ratio=True,
                      quantile=0.,
                      cgroup=cgroup,
                      out_dir=out_dir,
                      csv=csv,
                      html=html,
                      pdf=pdf,
                      max_procs=max_procs,
                      parallel=parallel,
                      force=force)
        # Graph rolling AUC-ROC.
        RollingAUCRunner.run(input_path=input_path,
                             rels=("m",),
                             use_ratio=True,
                             quantile=0.,
                             window=window,
                             winmin=winmin,
                             cgroup=cgroup,
                             out_dir=out_dir,
                             csv=csv,
                             html=html,
                             pdf=pdf,
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
