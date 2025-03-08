from pathlib import Path
from typing import Iterable

from click import command

from .write import cluster
from ..core.arg import (CMD_CLUSTER,
                        arg_input_path,
                        opt_branch,
                        opt_tmp_pfx,
                        opt_keep_tmp,
                        opt_min_clusters,
                        opt_max_clusters,
                        opt_min_em_runs,
                        opt_max_em_runs,
                        opt_em_thresh,
                        opt_min_em_iter,
                        opt_max_em_iter,
                        opt_jackpot,
                        opt_jackpot_conf_level,
                        opt_max_jackpot_quotient,
                        opt_max_pearson_run,
                        opt_min_marcd_run,
                        opt_max_arcd_vs_ens_avg,
                        opt_max_gini_run,
                        opt_max_loglike_vs_best,
                        opt_min_pearson_vs_best,
                        opt_max_marcd_vs_best,
                        opt_try_all_ks,
                        opt_write_all_ks,
                        opt_cluster_pos_table,
                        opt_cluster_abundance_table,
                        opt_verify_times,
                        opt_brotli_level,
                        opt_num_cpus,
                        opt_force)
from ..core.run import run_func
from ..core.task import as_list_of_tuples, dispatch
from ..mask.dataset import load_mask_dataset


@run_func(CMD_CLUSTER)
def run(input_path: Iterable[str | Path], *,
        branch: str,
        tmp_pfx: str | Path,
        keep_tmp: bool,
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
        max_loglike_vs_best: float,
        min_pearson_vs_best: float,
        max_marcd_vs_best: float,
        try_all_ks: bool,
        write_all_ks: bool,
        cluster_pos_table: bool,
        cluster_abundance_table: bool,
        verify_times: bool,
        brotli_level: int,
        num_cpus: int,
        force: bool) -> list[Path]:
    """ Infer alternative structures by clustering reads' mutations. """
    datasets = load_mask_dataset.iterate(input_path, verify_times=verify_times)
    return dispatch(cluster,
                    num_cpus=num_cpus,
                    pass_num_cpus=True,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=as_list_of_tuples(datasets),
                    kwargs=dict(tmp_pfx=tmp_pfx,
                                branch=branch,
                                keep_tmp=keep_tmp,
                                min_clusters=min_clusters,
                                max_clusters=max_clusters,
                                min_em_runs=min_em_runs,
                                max_em_runs=max_em_runs,
                                min_iter=min_em_iter,
                                max_iter=max_em_iter,
                                em_thresh=em_thresh,
                                min_marcd_run=min_marcd_run,
                                max_pearson_run=max_pearson_run,
                                max_arcd_vs_ens_avg=max_arcd_vs_ens_avg,
                                max_gini_run=max_gini_run,
                                jackpot=jackpot,
                                jackpot_conf_level=jackpot_conf_level,
                                max_jackpot_quotient=max_jackpot_quotient,
                                max_loglike_vs_best=max_loglike_vs_best,
                                min_pearson_vs_best=min_pearson_vs_best,
                                max_marcd_vs_best=max_marcd_vs_best,
                                try_all_ks=try_all_ks,
                                write_all_ks=write_all_ks,
                                cluster_pos_table=cluster_pos_table,
                                cluster_abundance_table=cluster_abundance_table,
                                verify_times=verify_times,
                                brotli_level=brotli_level,
                                force=force))


params = [
    # Input and output files
    arg_input_path,
    opt_branch,
    # Clustering options
    opt_min_clusters,
    opt_max_clusters,
    opt_min_em_runs,
    opt_max_em_runs,
    opt_em_thresh,
    opt_min_em_iter,
    opt_max_em_iter,
    opt_max_pearson_run,
    opt_min_marcd_run,
    opt_max_arcd_vs_ens_avg,
    opt_max_gini_run,
    opt_jackpot,
    opt_jackpot_conf_level,
    opt_max_jackpot_quotient,
    opt_max_loglike_vs_best,
    opt_min_pearson_vs_best,
    opt_max_marcd_vs_best,
    opt_try_all_ks,
    opt_write_all_ks,
    # Table options
    opt_cluster_pos_table,
    opt_cluster_abundance_table,
    # Validation
    opt_verify_times,
    # Compression
    opt_brotli_level,
    # Parallelization
    opt_num_cpus,
    # Effort
    opt_force,
    opt_tmp_pfx,
    opt_keep_tmp,
]


@command(CMD_CLUSTER, params=params)
def cli(*args, **kwargs):
    """ Infer structure ensembles by clustering reads' mutations. """
    return run(*args, **kwargs)
