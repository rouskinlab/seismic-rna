from pathlib import Path
from typing import Iterable
from click import command
from . import (
    demult as demultiplex_mod,
    align as align_mod,
    idmut as idmut_mod,
    filter as filter_mod,
    cluster as cluster_mod,
    fold as fold_mod,
    draw as draw_mod,
    collate as collate_mod,
    export as export_mod,
)
from .core.arg import (
    CMD_WORKFLOW,
    merge_params,
    opt_branch,
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
    opt_terminal_pairs,
    opt_draw,
    opt_window,
    opt_winmin,
    opt_csv,
    opt_html,
    opt_svg,
    opt_pdf,
    opt_png,
    opt_graph_quantile,
    opt_graph_mprof,
    opt_graph_tmprof,
    opt_graph_ncov,
    opt_graph_mhist,
    opt_graph_abundance,
    opt_graph_giniroll,
    opt_graph_roc,
    opt_graph_aucroll,
    opt_graph_poscorr,
    opt_graph_mutdist,
    opt_mutdist_null,
    opt_collate,
)
from .core.run import run_func
from .core.seq import DNA
from .core.table import (
    DELET_REL,
    INSRT_REL,
    MUTAT_REL,
    REL_NAMES,
    SUB_A_REL,
    SUB_C_REL,
    SUB_G_REL,
    SUB_T_REL,
    INFOR_REL,
)
from .graph.abundance import ClusterAbundanceRunner
from .graph.aucroll import RollingAUCRunner
from .graph.giniroll import RollingGiniRunner
from .graph.histread import ReadHistogramRunner
from .graph.mutdist import MutationDistanceRunner
from .graph.poscorr import PositionCorrelationRunner
from .graph.profile import ProfileRunner
from .graph.roc import ROCRunner

MUTAT_RELS = "".join(
    REL_NAMES[code]
    for code in [SUB_A_REL, SUB_C_REL, SUB_G_REL, SUB_T_REL, DELET_REL, INSRT_REL]
)


def flatten(nested):
    for item in nested:
        if isinstance(item, (list, tuple)):
            yield from flatten(item)
        else:
            yield item


@run_func(CMD_WORKFLOW, default=None)
def run(
    fasta: str | Path,
    input_path: Iterable[str | Path],
    *,
    # General options
    out_dir: str | Path,
    tmp_pfx: str | Path,
    keep_tmp: bool,
    brotli_level: int,
    force: bool,
    num_cpus: int,
    # FASTQ options
    fastqz: Iterable[str | Path],
    fastqy: Iterable[str | Path],
    fastqx: Iterable[str | Path],
    phred_enc: int,
    # Demultiplexing options
    demult: bool,
    barcode_start: int,
    barcode_end: int,
    read_pos: int,
    barcode: tuple[tuple[str, DNA, int]],
    mismatch_tolerance: int,
    index_tolerance: int,
    allow_n: bool,
    # Align options
    dmfastqz: Iterable[str | Path],
    dmfastqy: Iterable[str | Path],
    dmfastqx: Iterable[str | Path],
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
    # IDmut options
    min_phred: int,
    min_reads: int,
    insert3: bool,
    ambindel: bool,
    overhangs: bool,
    clip_end5: int,
    clip_end3: int,
    batch_size: int,
    write_read_names: bool,
    idmut_pos_table: bool,
    idmut_read_table: bool,
    idmut_cx: bool,
    # Filter options
    region_coords: Iterable[tuple[str, int, int]],
    region_primers: Iterable[tuple[str, DNA, DNA]],
    primer_gap: int,
    regions_file: str | None,
    count_del: bool,
    count_ins: bool,
    no_mut: Iterable[str],
    only_mut: Iterable[str],
    probe: str,
    mask_a: bool | None,
    mask_c: bool | None,
    mask_g: bool | None,
    mask_u: bool | None,
    mask_polya: int,
    mask_pos: Iterable[tuple[str, int]],
    mask_pos_file: Iterable[str | Path],
    drop_read: Iterable[str],
    drop_read_file: Iterable[str | Path],
    drop_discontig: bool,
    min_ncov_read: int,
    min_fcov_read: float,
    min_finfo_read: float,
    max_fmut_read: float,
    min_mut_gap: int,
    mut_collisions: str,
    min_ninfo_pos: int,
    max_fmut_pos: float,
    quick_unbias: bool,
    quick_unbias_thresh: float,
    max_filter_iter: int,
    filter_pos_table: bool,
    filter_read_table: bool,
    # Cluster options
    cluster: bool,
    min_clusters: int,
    max_clusters: int,
    min_em_runs: int,
    max_em_runs: int,
    jackpot: bool,
    jackpot_conf_level: float,
    max_jackpot_quotient: float,
    max_jackpot_sims: int,
    jackpot_max_data: int,
    min_em_iter: int,
    max_em_iter: int,
    em_thresh: float,
    min_marcd_run: float,
    max_pearson_run: float,
    max_arcd_vs_ens_avg: float,
    max_gini_run: float,
    max_loglike_vs_best: float,
    min_pearson_vs_best: float,
    max_marcd_vs_best: float,
    try_all_ks: bool,
    write_all_ks: bool,
    cluster_pos_table: bool,
    cluster_abundance_table: bool,
    verify_times: bool,
    # Self-contained batches
    self_contained: bool,
    # Fold options
    fold: bool,
    fold_coords: Iterable[tuple[str, int, int]],
    fold_primers: Iterable[tuple[str, DNA, DNA]],
    fold_regions_file: str | None,
    fold_full: bool,
    fold_dry_run: bool,
    fold_backend: str,
    pseudoknots: bool,
    fold_temp: float,
    fold_energy_method: str,
    deigan_slope: float,
    deigan_intercept: float,
    fold_quantile: float,
    fold_constraint: str | None,
    fold_commands: str | None,
    eddy_prior_paired_file: str | None,
    eddy_prior_unpaired_file: str | None,
    fold_md: int,
    fold_mfe: bool,
    fold_max: int,
    fold_percent: float,
    fold_edelta: float,
    fold_isolated: bool,
    # Draw options,
    draw: bool,
    struct_num: Iterable[int],
    color: bool,
    draw_svg: bool,
    draw_png: bool,
    update_rnartistcore: bool,
    # Export options,
    export: bool,
    samples_meta: str,
    refs_meta: str,
    all_pos: bool,
    # Graph options
    cgroup: str,
    graph_quantile: float,
    hist_bins: int,
    hist_margin: float,
    struct_file: Iterable[str | Path],
    terminal_pairs: bool,
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
    graph_abundance: bool,
    graph_giniroll: bool,
    graph_roc: bool,
    graph_aucroll: bool,
    graph_poscorr: bool,
    graph_mutdist: bool,
    mutdist_null: bool,
    # Collate options
    collate: bool,
    name: str,
    verbose_name: bool,
    include_svg: bool,
    include_graph: bool,
    group: str,
    portable: bool,
    collate_out_dir: str | Path | None = None,
    seed: int | None,
):
    """Run the entire workflow."""
    # Ensure that each iterable argument is a list rather than an
    # exhaustible generator.
    input_path = list(input_path)
    fastqx = list(fastqx)
    fastqy = list(fastqy)
    fastqz = list(fastqz)
    dmfastqx = list(dmfastqx)
    dmfastqy = list(dmfastqy)
    dmfastqz = list(dmfastqz)
    region_coords = list(region_coords)
    region_primers = list(region_primers)
    mask_pos = list(mask_pos)
    drop_read = list(drop_read)
    fold_coords = list(fold_coords)
    fold_primers = list(fold_primers)
    struct_num = list(struct_num)
    struct_file = list(struct_file)
    if demult:
        for dmz, dmy, dmx in demultiplex_mod.run(
            fasta=fasta,
            refs_meta=refs_meta,
            out_dir=out_dir,
            tmp_pfx=tmp_pfx,
            keep_tmp=keep_tmp,
            barcode_start=barcode_start,
            barcode_end=barcode_end,
            read_pos=read_pos,
            barcode=barcode,
            mismatch_tolerance=mismatch_tolerance,
            index_tolerance=index_tolerance,
            allow_n=allow_n,
            num_cpus=num_cpus,
            force=force,
            phred_enc=phred_enc,
        ):
            dmfastqz = dmfastqz + dmz
            dmfastqy = dmfastqy + dmy
            dmfastqx = dmfastqx + dmx
        # Clear the input FASTQ files once the demultiplexed FASTQ files
        # have been generated.
        fastqx = list()
    input_path.extend(
        flatten(
            align_mod.run(
                out_dir=out_dir,
                tmp_pfx=tmp_pfx,
                keep_tmp=keep_tmp,
                force=force,
                num_cpus=num_cpus,
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
            )
        )
    )
    input_path.extend(
        flatten(
            idmut_mod.run(
                fasta=fasta,
                input_path=input_path,
                out_dir=out_dir,
                tmp_pfx=tmp_pfx,
                keep_tmp=keep_tmp,
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
                write_read_names=write_read_names,
                idmut_pos_table=idmut_pos_table,
                idmut_read_table=idmut_read_table,
                idmut_cx=idmut_cx,
                num_cpus=num_cpus,
                brotli_level=brotli_level,
                force=force,
            )
        )
    )
    input_path.extend(
        flatten(
            filter_mod.run(
                input_path=input_path,
                tmp_pfx=tmp_pfx,
                keep_tmp=keep_tmp,
                region_coords=region_coords,
                region_primers=region_primers,
                primer_gap=primer_gap,
                regions_file=regions_file,
                count_del=count_del,
                count_ins=count_ins,
                no_mut=no_mut,
                only_mut=only_mut,
                probe=probe,
                mask_a=mask_a,
                mask_c=mask_c,
                mask_g=mask_g,
                mask_u=mask_u,
                mask_polya=mask_polya,
                mask_pos=mask_pos,
                mask_pos_file=mask_pos_file,
                drop_read=drop_read,
                drop_read_file=drop_read_file,
                drop_discontig=drop_discontig,
                min_ncov_read=min_ncov_read,
                min_fcov_read=min_fcov_read,
                min_finfo_read=min_finfo_read,
                max_fmut_read=max_fmut_read,
                min_mut_gap=min_mut_gap,
                mut_collisions=mut_collisions,
                min_ninfo_pos=min_ninfo_pos,
                max_fmut_pos=max_fmut_pos,
                quick_unbias=quick_unbias,
                quick_unbias_thresh=quick_unbias_thresh,
                max_filter_iter=max_filter_iter,
                filter_pos_table=filter_pos_table,
                filter_read_table=filter_read_table,
                brotli_level=brotli_level,
                self_contained=self_contained,
                num_cpus=num_cpus,
                force=force,
            )
        )
    )
    if (
        cluster
        or min_clusters != opt_min_clusters.default
        or max_clusters != opt_max_clusters.default
    ):
        input_path.extend(
            flatten(
                cluster_mod.run(
                    input_path=input_path,
                    tmp_pfx=tmp_pfx,
                    keep_tmp=keep_tmp,
                    min_clusters=min_clusters,
                    max_clusters=max_clusters,
                    min_em_runs=min_em_runs,
                    max_em_runs=max_em_runs,
                    jackpot=jackpot,
                    jackpot_conf_level=jackpot_conf_level,
                    max_jackpot_quotient=max_jackpot_quotient,
                    max_jackpot_sims=max_jackpot_sims,
                    jackpot_max_data=jackpot_max_data,
                    min_em_iter=min_em_iter,
                    max_em_iter=max_em_iter,
                    em_thresh=em_thresh,
                    min_marcd_run=min_marcd_run,
                    max_pearson_run=max_pearson_run,
                    max_arcd_vs_ens_avg=max_arcd_vs_ens_avg,
                    max_gini_run=max_gini_run,
                    max_loglike_vs_best=max_loglike_vs_best,
                    min_pearson_vs_best=min_pearson_vs_best,
                    max_marcd_vs_best=max_marcd_vs_best,
                    try_all_ks=try_all_ks,
                    write_all_ks=write_all_ks,
                    cluster_pos_table=cluster_pos_table,
                    cluster_abundance_table=cluster_abundance_table,
                    verify_times=verify_times,
                    brotli_level=brotli_level,
                    self_contained=self_contained,
                    num_cpus=num_cpus,
                    force=force,
                    seed=seed,
                )
            )
        )
    if fold:
        input_path.extend(
            flatten(
                fold_mod.run(
                    input_path=input_path,
                    branch="",
                    fold_regions_file=fold_regions_file,
                    fold_coords=fold_coords,
                    fold_primers=fold_primers,
                    fold_full=fold_full,
                    fold_dry_run=fold_dry_run,
                    fold_backend=fold_backend,
                    pseudoknots=pseudoknots,
                    fold_energy_method=fold_energy_method,
                    deigan_slope=deigan_slope,
                    deigan_intercept=deigan_intercept,
                    fold_temp=fold_temp,
                    fold_quantile=fold_quantile,
                    fold_constraint=fold_constraint,
                    fold_commands=fold_commands,
                    eddy_prior_paired_file=eddy_prior_paired_file,
                    eddy_prior_unpaired_file=eddy_prior_unpaired_file,
                    fold_md=fold_md,
                    fold_mfe=fold_mfe,
                    fold_max=fold_max,
                    fold_percent=fold_percent,
                    fold_edelta=fold_edelta,
                    fold_isolated=fold_isolated,
                    tmp_pfx=tmp_pfx,
                    keep_tmp=keep_tmp,
                    num_cpus=num_cpus,
                    force=force,
                )
            )
        )
        if graph_roc:
            input_path.extend(
                flatten(
                    ROCRunner.run(
                        input_path=input_path,
                        branch="",
                        rels=[REL_NAMES[MUTAT_REL]],
                        use_ratio=True,
                        struct_file=struct_file,
                        fold_regions_file=fold_regions_file,
                        fold_coords=fold_coords,
                        fold_primers=fold_primers,
                        fold_full=fold_full,
                        terminal_pairs=terminal_pairs,
                        cgroup=cgroup,
                        graph_quantile=graph_quantile,
                        csv=csv,
                        html=html,
                        svg=svg,
                        pdf=pdf,
                        png=png,
                        verify_times=verify_times,
                        num_cpus=num_cpus,
                        force=force,
                    )
                )
            )
        if graph_aucroll:
            input_path.extend(
                flatten(
                    RollingAUCRunner.run(
                        input_path=input_path,
                        rels=[REL_NAMES[MUTAT_REL]],
                        use_ratio=True,
                        struct_file=struct_file,
                        fold_regions_file=fold_regions_file,
                        fold_coords=fold_coords,
                        fold_primers=fold_primers,
                        fold_full=fold_full,
                        terminal_pairs=terminal_pairs,
                        window=window,
                        winmin=winmin,
                        cgroup=cgroup,
                        graph_quantile=graph_quantile,
                        csv=csv,
                        html=html,
                        svg=svg,
                        pdf=pdf,
                        png=png,
                        verify_times=verify_times,
                        num_cpus=num_cpus,
                        force=force,
                    )
                )
            )
    if draw:
        input_path.extend(
            flatten(
                draw_mod.run(
                    input_path=input_path,
                    struct_num=struct_num,
                    color=color,
                    draw_svg=draw_svg,
                    draw_png=draw_png,
                    update_rnartistcore=update_rnartistcore,
                    tmp_pfx=tmp_pfx,
                    keep_tmp=keep_tmp,
                    num_cpus=num_cpus,
                    force=force,
                )
            )
        )
    if graph_mprof or graph_tmprof:
        rels = list()
        if graph_mprof:
            rels.append(REL_NAMES[MUTAT_REL])
        if graph_tmprof:
            rels.append(MUTAT_RELS)
        input_path.extend(
            flatten(
                ProfileRunner.run(
                    input_path=input_path,
                    rels=rels,
                    use_ratio=True,
                    graph_quantile=graph_quantile,
                    cgroup=cgroup,
                    csv=csv,
                    html=html,
                    svg=svg,
                    pdf=pdf,
                    png=png,
                    verify_times=verify_times,
                    num_cpus=num_cpus,
                    force=force,
                )
            )
        )
    if graph_ncov:
        input_path.extend(
            flatten(
                ProfileRunner.run(
                    input_path=input_path,
                    rels=[REL_NAMES[INFOR_REL]],
                    use_ratio=False,
                    graph_quantile=graph_quantile,
                    cgroup=cgroup,
                    csv=csv,
                    html=html,
                    svg=svg,
                    pdf=pdf,
                    png=png,
                    verify_times=verify_times,
                    num_cpus=num_cpus,
                    force=force,
                )
            )
        )
    if graph_mhist:
        input_path.extend(
            flatten(
                ReadHistogramRunner.run(
                    input_path=input_path,
                    rels=[REL_NAMES[MUTAT_REL]],
                    use_ratio=False,
                    graph_quantile=graph_quantile,
                    cgroup=cgroup,
                    hist_bins=hist_bins,
                    hist_margin=hist_margin,
                    csv=csv,
                    html=html,
                    svg=svg,
                    pdf=pdf,
                    png=png,
                    verify_times=verify_times,
                    num_cpus=num_cpus,
                    force=force,
                )
            )
        )
    if graph_abundance:
        input_path.extend(
            flatten(
                ClusterAbundanceRunner.run(
                    input_path=input_path,
                    use_ratio=True,
                    csv=csv,
                    html=html,
                    svg=svg,
                    pdf=pdf,
                    png=png,
                    verify_times=verify_times,
                    num_cpus=num_cpus,
                    force=force,
                )
            )
        )
    if graph_giniroll:
        input_path.extend(
            flatten(
                RollingGiniRunner.run(
                    input_path=input_path,
                    rels=[REL_NAMES[MUTAT_REL]],
                    use_ratio=True,
                    graph_quantile=graph_quantile,
                    window=window,
                    winmin=winmin,
                    cgroup=cgroup,
                    csv=csv,
                    html=html,
                    svg=svg,
                    pdf=pdf,
                    png=png,
                    verify_times=verify_times,
                    num_cpus=num_cpus,
                    force=force,
                )
            )
        )
    if graph_poscorr:
        input_path.extend(
            flatten(
                PositionCorrelationRunner.run(
                    input_path=input_path,
                    rels=[REL_NAMES[MUTAT_REL]],
                    cgroup=cgroup,
                    verify_times=verify_times,
                    csv=csv,
                    html=html,
                    svg=svg,
                    pdf=pdf,
                    png=png,
                    num_cpus=num_cpus,
                    force=force,
                )
            )
        )
    if graph_mutdist:
        input_path.extend(
            flatten(
                MutationDistanceRunner.run(
                    input_path=input_path,
                    rels=[REL_NAMES[MUTAT_REL]],
                    cgroup=cgroup,
                    verify_times=verify_times,
                    mutdist_null=mutdist_null,
                    csv=csv,
                    html=html,
                    svg=svg,
                    pdf=pdf,
                    png=png,
                    num_cpus=num_cpus,
                    force=force,
                )
            )
        )
    if collate:
        collate_mod.run(
            input_path=input_path,
            name=name,
            verbose_name=verbose_name,
            include_svg=include_svg,
            include_graph=include_graph,
            group=group,
            portable=portable,
            collate_out_dir=collate_out_dir,
            force=force,
        )
    if export:
        export_mod.run(
            input_path=input_path,
            samples_meta=samples_meta,
            refs_meta=refs_meta,
            verify_times=verify_times,
            all_pos=all_pos,
            num_cpus=num_cpus,
            force=force,
        )


graph_options = [
    opt_cgroup,
    opt_hist_bins,
    opt_hist_margin,
    opt_struct_file,
    opt_fold_full,
    opt_terminal_pairs,
    opt_window,
    opt_winmin,
    opt_graph_quantile,
    opt_csv,
    opt_html,
    opt_svg,
    opt_pdf,
    opt_png,
    opt_graph_mprof,
    opt_graph_tmprof,
    opt_graph_ncov,
    opt_graph_mhist,
    opt_graph_abundance,
    opt_graph_giniroll,
    opt_graph_roc,
    opt_graph_aucroll,
    opt_graph_poscorr,
    opt_graph_mutdist,
    opt_mutdist_null,
]

params = merge_params(
    [opt_demultiplex],
    demultiplex_mod.params,
    align_mod.params,
    idmut_mod.params,
    filter_mod.params,
    [opt_cluster],
    cluster_mod.params,
    [opt_fold],
    fold_mod.params,
    [opt_draw],
    draw_mod.params,
    [opt_export],
    export_mod.params,
    graph_options,
    [opt_collate],
    collate_mod.params,
    exclude=[opt_branch],
)


@command(CMD_WORKFLOW, params=params)
def cli(*args, **kwargs):
    """Run the entire workflow."""
    return run(*args, **kwargs)
