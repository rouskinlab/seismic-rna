from __future__ import annotations
from datetime import datetime
from pathlib import Path

from .. import cluster as cluster_mod
from .. import filter as filter_mod
from ..core import path
from ..core.logs import logger
from ..core.report import (
    BestKF,
    SampleF,
    RefF,
    RegF,
    BranchesF,
    DomainCoordsF,
    # Filter parameters recorded in each domain's FilterReport, reused to
    # re-filter merged regions consistently with how filterscan filtered.
    ProbeF,
    MaskAF,
    MaskCF,
    MaskGF,
    MaskUF,
    MaskPolyAF,
    MaskPosF,
    MinNInfoPosF,
    MaxFMutPosF,
    MinNCovReadF,
    MinFCovReadF,
    DiscontigF,
    MinFInfoReadF,
    MaxFMutReadF,
    MinMutGapF,
    MutCollisionsF,
    QuickUnbiasF,
    QuickUnbiasThreshF,
    MaxFilterIterF,
)
from ..core.rel.pattern import RelPattern
from ..core.seq.region import hyphenate_ends
from ..core.write import need_write
from ..cluster.report import ClusterReport
from ..cluster.uniq import UniqReads
from ..filter.dataset import FilterMutsDataset
from ..filter.report import FilterReport
from ..filterscan.report import FilterScanReport
from .report import ClusterScanReport
from .gap import evaluate_gap


def _load_domain_clustering(cluster_report_file: Path):
    """Load a clustered domain's best number of clusters and the best run's
    mutation rates (position-indexed DataFrame, one column per cluster).

    Returns (best_k, mus) where mus is None if the parameter file is missing
    (in which case the gap test on that domain is skipped)."""
    import pandas as pd

    report = ClusterReport.load(cluster_report_file)
    best_k = max(int(report.get_field(BestKF)), 1)
    mus_file = cluster_report_file.parent.joinpath(
        "parameters", f"k{best_k}-r0_mus.csv"
    )
    if not mus_file.is_file():
        logger.warning("Missing cluster mutation rates {}", mus_file)
        return best_k, None
    mus = pd.read_csv(mus_file, header=[0, 1], index_col=0)
    mus.columns = [int(cluster) for _, cluster in mus.columns]
    mus.index = mus.index.astype(int)
    mus.index.name = "pos"
    return best_k, mus


def _reconstruct_filter_kwargs(filter_report: FilterReport, ref: str) -> dict:
    """Rebuild the filter kwargs filterscan used, from a domain's FilterReport,
    so a merged region can be re-filtered identically."""
    return dict(
        probe=filter_report.get_field(ProbeF),
        mask_a=filter_report.get_field(MaskAF),
        mask_c=filter_report.get_field(MaskCF),
        mask_g=filter_report.get_field(MaskGF),
        mask_u=filter_report.get_field(MaskUF),
        mask_polya=filter_report.get_field(MaskPolyAF),
        mask_pos=[(ref, pos) for pos in filter_report.get_field(MaskPosF)],
        min_ninfo_pos=filter_report.get_field(MinNInfoPosF),
        max_fmut_pos=filter_report.get_field(MaxFMutPosF),
        min_ncov_read=filter_report.get_field(MinNCovReadF),
        min_fcov_read=filter_report.get_field(MinFCovReadF),
        drop_discontig=filter_report.get_field(DiscontigF),
        min_finfo_read=filter_report.get_field(MinFInfoReadF),
        max_fmut_read=filter_report.get_field(MaxFMutReadF),
        min_mut_gap=filter_report.get_field(MinMutGapF),
        mut_collisions=filter_report.get_field(MutCollisionsF),
        quick_unbias=filter_report.get_field(QuickUnbiasF),
        quick_unbias_thresh=filter_report.get_field(QuickUnbiasThreshF),
        max_filter_iter=filter_report.get_field(MaxFilterIterF),
        # filterscan/clusterscan use default mutation counting; verified below.
        count_del=True,
        count_ins=True,
        no_mut=(),
        only_mut=(),
        # Options with no persisted effect for a re-filter.
        region_primers=(),
        primer_gap=0,
        regions_file=None,
        mask_pos_file=(),
        drop_read=(),
        drop_read_file=(),
    )


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
    # Gap-validation options
    validate_gaps: bool,
    gap_min_assoc: float,
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
    html: bool,
    verify_times: bool,
    self_contained: bool,
    seed: int | None,
):
    """Cluster the domains detected by filterscan for one FilterScanReport.
    If ``validate_gaps`` is set, also validate every gap between adjacent
    domains: test on the reads spanning the gap whether the two domains'
    clusters are independent. If a gap fails (the domains are coupled),
    merge the two domains, re-cluster the merged region, and repeat until
    every remaining gap passes. Write a ClusterScanReport recording the
    final domains, their numbers of clusters, and the original filterscan
    domains comprising each final domain.
    """
    began = datetime.now()
    top, _ = FilterScanReport.parse_path(filterscan_report_file)
    report = FilterScanReport.load(filterscan_report_file)
    sample = report.get_field(SampleF)
    ref = report.get_field(RefF)
    reg = report.get_field(RegF)
    domain_coords = sorted(report.get_field(DomainCoordsF).keys())
    # Branches: filterscan's branches locate the per-domain filter results;
    # clusterscan adds its own branch to build the ClusterScanReport path.
    filterscan_branches = report.get_field(BranchesF)
    report_branches = path.add_branch(
        path.CLUSTERSCAN_STEP, branch, filterscan_branches
    )
    report_file = ClusterScanReport.build_path(
        {
            path.TOP: top,
            path.SAMPLE: sample,
            path.BRANCHES: report_branches,
            path.REF: ref,
            path.REG: reg,
        }
    )
    if not need_write(report_file, force):
        return report_file.parent

    def domain_filter_dir(end5: int, end3: int) -> Path:
        return FilterReport.build_path(
            {
                path.TOP: top,
                path.SAMPLE: sample,
                path.BRANCHES: filterscan_branches,
                path.REF: ref,
                path.REG: hyphenate_ends(end5, end3),
            }
        ).parent

    if not domain_coords:
        logger.warning("No domains in {}; skipping clustering", filterscan_report_file)
        report_file.parent.mkdir(parents=True, exist_ok=True)
        ClusterScanReport(
            sample=sample,
            ref=ref,
            reg=reg,
            branches=report_branches,
            best_ks={},
            merged_domains={},
            began=began,
            ended=datetime.now(),
        ).save(top, force=force)
        return report_file.parent

    cluster_kwargs = dict(
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
        html=html,
        verify_times=verify_times,
        self_contained=self_contained,
        brotli_level=brotli_level,
        num_cpus=num_cpus,
        force=force,
        seed=seed,
    )

    def cluster_regions(filter_dirs: list[Path]) -> dict[tuple[int, int], Path]:
        """Cluster the given filter dirs; return {(end5, end3): cluster_report}."""
        cluster_dirs = cluster_mod.run(input_path=filter_dirs, **cluster_kwargs)
        reports: dict[tuple[int, int], Path] = {}
        for f in path.find_files_chain(
            cluster_dirs, ClusterReport.get_path_seg_types()
        ):
            name = ClusterReport.load(f).get_field(RegF)
            end5, end3 = map(int, name.split("-"))
            reports[(end5, end3)] = f
        return reports

    # 1. Cluster every initial filterscan domain.
    domain_dirs = [domain_filter_dir(e5, e3) for e5, e3 in domain_coords]
    domain_reports = cluster_regions(domain_dirs)

    # 2. Build the working list of domains.
    domains: list[dict] = []
    for e5, e3 in domain_coords:
        k, mus = _load_domain_clustering(domain_reports[(e5, e3)])
        domains.append(dict(end5=e5, end3=e3, k=k, mus=mus, members=[(e5, e3)]))

    # 3. Agglomeratively test gaps and merge coupled pairs.
    if validate_gaps:
        # Reach the upstream idmut report and the filter parameters (from any
        # domain's FilterReport) so merged regions can be re-filtered
        # identically.
        a_filter_report_file = domain_filter_dir(*domain_coords[0]).joinpath(
            "filter-report.json"
        )
        idmut_report_file = FilterMutsDataset.get_dataset1_report_file(
            a_filter_report_file, verify_times
        )
        a_filter_report = FilterReport.load(a_filter_report_file)
        if FilterMutsDataset(a_filter_report_file).pattern != RelPattern.from_counts():
            raise ValueError(
                "clusterscan gap validation supports only default mutation "
                "counting (count_del, count_ins, no_mut, only_mut)"
            )
        filter_kwargs = _reconstruct_filter_kwargs(a_filter_report, ref)

        merged_filter_dirs: dict[tuple[int, int], Path] = {}

        def merged_filter_dir(end5: int, end3: int) -> Path:
            """Re-filter the merged (end5, end3) region from idmut, memoized."""
            key = (end5, end3)
            if key not in merged_filter_dirs:
                dirs = filter_mod.run(
                    input_path=[idmut_report_file],
                    branch=filterscan_branches[path.FILTERSCAN_STEP],
                    tmp_pfx=tmp_pfx,
                    keep_tmp=keep_tmp,
                    region_coords=[(ref, end5, end3)],
                    filter_pos_table=False,
                    filter_read_table=False,
                    self_contained=self_contained,
                    brotli_level=brotli_level,
                    num_cpus=num_cpus,
                    force=force,
                    **filter_kwargs,
                )
                merged_filter_dirs[key] = Path(dirs[0])
            return merged_filter_dirs[key]

        # Cache gap tests by the pair's coordinates so only pairs touching a
        # newly merged domain are recomputed.
        tested: dict[tuple, tuple[bool, float]] = {}

        def gap_test(left: dict, right: dict) -> tuple[bool, float]:
            key = (left["end5"], left["end3"], right["end5"], right["end3"])
            if key in tested:
                return tested[key]
            if (
                left["k"] < 2
                or right["k"] < 2
                or left["mus"] is None
                or right["mus"] is None
            ):
                result = (True, float("-inf"))
            else:
                ur = UniqReads.from_dataset_contig(
                    FilterMutsDataset(
                        merged_filter_dir(left["end5"], right["end3"]).joinpath(
                            "filter-report.json"
                        ),
                        verify_times=verify_times,
                    ),
                    branch,
                )
                passed, info = evaluate_gap(
                    ur,
                    left["mus"],
                    right["mus"],
                    left["end3"],
                    right["end5"],
                    gap_min_assoc,
                )
                result = (passed, info.get("coupling", float("-inf")))
            tested[key] = result
            return result

        while len(domains) > 1:
            worst_i = None
            worst_coupling = float("-inf")
            for i in range(len(domains) - 1):
                passed, coupling = gap_test(domains[i], domains[i + 1])
                if not passed and coupling > worst_coupling:
                    worst_i = i
                    worst_coupling = coupling
            if worst_i is None:
                break
            # Merge the most-coupled failing pair and re-cluster the merged
            # region.
            left = domains[worst_i]
            right = domains[worst_i + 1]
            merged5, merged3 = left["end5"], right["end3"]
            logger.info(
                "Merging coupled domains {}-{} and {}-{} (coupling={:.1f})",
                left["end5"],
                left["end3"],
                right["end5"],
                right["end3"],
                worst_coupling,
            )
            merged_dir = merged_filter_dir(merged5, merged3)
            merged_report = cluster_regions([merged_dir])[(merged5, merged3)]
            k, mus = _load_domain_clustering(merged_report)
            domains[worst_i : worst_i + 2] = [
                dict(
                    end5=merged5,
                    end3=merged3,
                    k=k,
                    mus=mus,
                    members=left["members"] + right["members"],
                )
            ]

    # 4. Write the ClusterScanReport with the final domains.
    report_file.parent.mkdir(parents=True, exist_ok=True)
    ClusterScanReport(
        sample=sample,
        ref=ref,
        reg=reg,
        branches=report_branches,
        best_ks={(d["end5"], d["end3"]): d["k"] for d in domains},
        merged_domains={
            (d["end5"], d["end3"]): d["members"]
            for d in domains
            if len(d["members"]) >= 2
        },
        began=began,
        ended=datetime.now(),
    ).save(top, force=force)
    return report_file.parent
