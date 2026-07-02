from __future__ import annotations
from datetime import datetime
from pathlib import Path

from .. import cluster as cluster_mod
from ..core import path
from ..core.logs import logger
from ..core.report import BestKF, SampleF, RefF, RegF, BranchesF, DomainCoordsF
from ..core.seq.region import hyphenate_ends
from ..core.write import need_write
from ..cluster.report import ClusterReport
from ..filter.report import FilterReport
from ..filterscan.report import FilterScanReport
from .report import ClusterScanReport


def clusterscan(
    filterscan_report_file: Path,
    *,
    # General options
    branch: str,
    tmp_pfx: str | Path,
    keep_tmp: bool,
    brotli_level: int,
    force: bool,
    num_cpus: int,
    # Cluster options
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
    self_contained: bool,
    seed: int | None,
):
    """Cluster the domains detected by filterscan for one FilterScanReport.

    Locate the domain filter results that filterscan produced, cluster
    each domain, and write a ClusterScanReport recording the cluster
    directories (relative to the output directory) and the best number
    of clusters per domain.
    """
    began = datetime.now()
    # Load the FilterScanReport and the output (top) directory.
    top, _ = FilterScanReport.parse_path(filterscan_report_file)
    report = FilterScanReport.load(filterscan_report_file)
    sample = report.get_field(SampleF)
    ref = report.get_field(RefF)
    reg = report.get_field(RegF)
    domain_coords = report.get_field(DomainCoordsF)
    # Recover the filterscan branch and the branches of the domain filter
    # results that filterscan produced.
    branches = report.get_field(BranchesF)
    filterscan_branch = branches[path.FILTERSCAN_STEP]
    idmut_branches = {
        step: name for step, name in branches.items() if step != path.FILTERSCAN_STEP
    }
    filter_branches = path.add_branch(
        path.FILTER_STEP, filterscan_branch, idmut_branches
    )
    # Build the path of the ClusterScanReport.
    report_branches = path.add_branch(path.CLUSTERSCAN_STEP, branch, idmut_branches)
    report_file = ClusterScanReport.build_path(
        {
            path.TOP: top,
            path.SAMPLE: sample,
            path.BRANCHES: report_branches,
            path.REF: ref,
            path.REG: reg,
        }
    )
    if need_write(report_file, force):
        # Locate each domain's filter results (named by their coordinates).
        domain_dirs = [
            FilterReport.build_path(
                {
                    path.TOP: top,
                    path.SAMPLE: sample,
                    path.BRANCHES: filter_branches,
                    path.REF: ref,
                    path.REG: hyphenate_ends(end5, end3),
                }
            ).parent
            for end5, end3 in domain_coords
        ]
        if domain_dirs:
            cluster_dirs = cluster_mod.run(
                input_path=domain_dirs,
                branch=branch,
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
                self_contained=self_contained,
                brotli_level=brotli_level,
                num_cpus=num_cpus,
                force=force,
                seed=seed,
            )
            best_ks = [
                ClusterReport.load(f).get_field(BestKF)
                for f in path.find_files_chain(
                    cluster_dirs, ClusterReport.get_path_seg_types()
                )
            ]
        else:
            logger.warning(
                "No domains in {}; skipping clustering", filterscan_report_file
            )
            cluster_dirs = list()
            best_ks = list()
        # Store the cluster directories relative to the output directory.
        rel_cluster_dirs = [str(Path(d).relative_to(top)) for d in cluster_dirs]
        report_file.parent.mkdir(parents=True, exist_ok=True)
        ClusterScanReport(
            sample=sample,
            ref=ref,
            reg=reg,
            branches=report_branches,
            cluster_dirs=rel_cluster_dirs,
            best_ks=best_ks,
            began=began,
            ended=datetime.now(),
        ).save(top, force=force)
    return report_file.parent
