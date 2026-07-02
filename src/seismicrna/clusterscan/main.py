from pathlib import Path
from typing import Iterable

from click import command

from .write import clusterscan
from .. import cluster as cluster_mod
from ..core import path
from ..core.arg.cmd import CMD_CLUSTERSCAN
from ..core.arg.cli import merge_params
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..filterscan.report import FilterScanReport


@run_func(CMD_CLUSTERSCAN)
def run(
    input_path: Iterable[str | Path],
    *,
    branch: str,
    tmp_pfx: str | Path,
    keep_tmp: bool,
    # Cluster options
    min_clusters: int,
    max_clusters: int,
    min_em_runs: int,
    max_em_runs: int,
    min_em_iter: int,
    max_em_iter: int,
    em_thresh: float,
    min_marcd_run: float,
    max_pearson_run: float,
    max_arcd_vs_ens_avg: float,
    max_gini_run: float,
    jackpot: bool,
    jackpot_conf_level: float,
    max_jackpot_quotient: float,
    max_jackpot_sims: int,
    jackpot_max_data: int,
    max_loglike_vs_best: float,
    min_pearson_vs_best: float,
    max_marcd_vs_best: float,
    try_all_ks: bool,
    write_all_ks: bool,
    cluster_pos_table: bool,
    cluster_abundance_table: bool,
    verify_times: bool,
    brotli_level: int,
    self_contained: bool,
    num_cpus: int,
    force: bool,
    seed: int | None,
):
    """Cluster the domains detected by filterscan."""
    seg_types = FilterScanReport.get_path_seg_types()
    filterscan_report_files = list(path.find_files_chain(input_path, seg_types))
    kwargs = dict(
        branch=branch,
        tmp_pfx=tmp_pfx,
        keep_tmp=keep_tmp,
        brotli_level=brotli_level,
        force=force,
        min_clusters=min_clusters,
        max_clusters=max_clusters,
        min_em_runs=min_em_runs,
        max_em_runs=max_em_runs,
        min_em_iter=min_em_iter,
        max_em_iter=max_em_iter,
        em_thresh=em_thresh,
        min_marcd_run=min_marcd_run,
        max_pearson_run=max_pearson_run,
        max_arcd_vs_ens_avg=max_arcd_vs_ens_avg,
        max_gini_run=max_gini_run,
        jackpot=jackpot,
        jackpot_conf_level=jackpot_conf_level,
        max_jackpot_quotient=max_jackpot_quotient,
        max_jackpot_sims=max_jackpot_sims,
        jackpot_max_data=jackpot_max_data,
        max_loglike_vs_best=max_loglike_vs_best,
        min_pearson_vs_best=min_pearson_vs_best,
        max_marcd_vs_best=max_marcd_vs_best,
        try_all_ks=try_all_ks,
        write_all_ks=write_all_ks,
        cluster_pos_table=cluster_pos_table,
        cluster_abundance_table=cluster_abundance_table,
        verify_times=verify_times,
        self_contained=self_contained,
        seed=seed,
    )
    return dispatch(
        clusterscan,
        num_cpus=num_cpus,
        pass_num_cpus=True,
        as_list=True,
        ordered=False,
        raise_on_error=False,
        args=as_list_of_tuples(filterscan_report_files),
        kwargs=kwargs,
    )


params = merge_params(cluster_mod.params)


@command(CMD_CLUSTERSCAN, params=params)
def cli(*args, **kwargs):
    """Cluster the domains detected by filterscan."""
    return run(*args, **kwargs)
