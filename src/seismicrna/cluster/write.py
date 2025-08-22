from datetime import datetime
from math import inf
from pathlib import Path
from typing import Iterable

import numpy as np

from .data import ClusterDatasetTabulator, ClusterMutsDataset
from .em import EMRun
from .emk import EMRunsK, find_best_k, sort_runs
from .io import ClusterBatchWriter
from .obsexp import write_obs_exp_counts
from .params import write_mus, write_pis
from .report import ClusterReport
from .summary import write_summaries
from .uniq import UniqReads
from ..core import path
from ..core.header import validate_ks
from ..core.logs import logger
from ..core.task import dispatch
from ..core.tmp import release_to_out, with_tmp_dir
from ..core.types import get_max_uint
from ..core.write import need_write
from ..mask.dataset import MaskMutsDataset, JoinMaskMutsDataset

SEED_DTYPE = np.uint32


def run_k(uniq_reads: UniqReads,
          k: int,
          num_runs: int, *,
          num_cpus: int,
          **kwargs) -> list[EMRun]:
    """ Run EM with a specific number of clusters. """
    if k < 1:
        raise ValueError(f"k must be ≥ 1, but got {k}")
    if num_runs < 1:
        logger.warning(
            f"Expected num_runs to be ≥ 1, but got {num_runs}: setting to 1"
        )
        num_runs = 1
    rng = np.random.default_rng()
    # On some but not all platforms, using this central source of seeds
    # is necessary because otherwise the runs would all take identical
    # trajectories, defeating the purpose of replicates.
    seeds = rng.integers(get_max_uint(SEED_DTYPE),
                         size=num_runs,
                         dtype=SEED_DTYPE)
    args = [(uniq_reads, k, seed) for seed in seeds]
    runs = dispatch(EMRun,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=args,
                    kwargs=kwargs)
    if len(runs) < num_runs:
        logger.warning(f"Obtained only {len(runs)} (of {num_runs}) "
                       f"run(s) of {uniq_reads} with {k} cluster(s)")
    return sort_runs(runs)


def run_ks(uniq_reads: UniqReads,
           ks: Iterable[int], *,
           min_em_runs: int,
           max_em_runs: int,
           try_all_ks: bool,
           min_marcd_run: float,
           max_pearson_run: float,
           max_arcd_vs_ens_avg: float,
           max_gini_run: float,
           max_jackpot_quotient: float,
           max_loglike_vs_best: float,
           min_pearson_vs_best: float,
           max_marcd_vs_best: float,
           min_iter: int,
           max_iter: int,
           em_thresh: float,
           num_cpus: int,
           top: Path,
           **kwargs):
    """ Run EM with multiple numbers of clusters. """
    if min_em_runs < 1:
        logger.warning(
            f"min_em_runs must be ≥ 1, but got {min_em_runs}; setting to 1"
        )
        min_em_runs = 1
    if max_em_runs < min_em_runs:
        logger.warning(
            "max_em_runs must be ≥ min_em_runs, "
            f"but got min_em_runs={min_em_runs} and max_em_runs={max_em_runs}; "
            f"setting to {min_em_runs}"
        )
        max_em_runs = min_em_runs
    if num_cpus < 1:
        logger.warning(
            f"num_cpus must be ≥ 1, but got {num_cpus}; setting to 1"
        )
        num_cpus = 1
    path_kwargs = {path.TOP: top,
                   path.BRANCHES: uniq_reads.branches,
                   path.SAMPLE: uniq_reads.sample,
                   path.REF: uniq_reads.ref,
                   path.REG: uniq_reads.region.name}
    runs_ks = dict()
    ks = validate_ks(ks)
    # Loop through every K if try_all_ks is True, otherwise go until the
    # current K is worse than the previous K or raises an error.
    for k in ks:
        assert k >= 1
        try:
            if k > 1:
                min_runs_k = min_em_runs
                max_runs_k = max_em_runs
            else:
                min_runs_k = 1
                max_runs_k = 1
            logger.routine(f"Began {min_runs_k} - {max_runs_k} run(s) of EM "
                           f"with {k} cluster(s)")
            # Accumulate EM runs for this K.
            runs_k = list()
            num_runs_k = 0
            passing = False
            while (num_runs_k < max_runs_k
                   and (len(runs_k) < min_runs_k or not passing)):
                if len(runs_k) < min_runs_k:
                    # The minimum number of runs has not been reached,
                    # so run as many times as needed to reach it.
                    num_runs = min_runs_k - len(runs_k)
                else:
                    # The minimum number of runs has been reached, but
                    # the runs do not pass all filters, so run up to
                    # max_runs_k - num_runs_k more. To be efficient,
                    # run up to num_cpus at a time.
                    assert not passing
                    num_runs = min(max_runs_k - num_runs_k, num_cpus)
                num_runs_k += num_runs
                # Run a batch of EM.
                runs_k.extend(run_k(uniq_reads,
                                    k,
                                    num_runs=num_runs,
                                    em_thresh=(em_thresh if k > 1 else inf),
                                    min_iter=(min_iter if k > 1 else 2),
                                    max_iter=(max_iter if k > 1 else 2),
                                    max_jackpot_quotient=max_jackpot_quotient,
                                    num_cpus=num_cpus,
                                    **kwargs))
                # Collect all runs so far for this K.
                runs_ks[k] = EMRunsK(runs_k,
                                     max_pearson_run=max_pearson_run,
                                     min_marcd_run=min_marcd_run,
                                     max_arcd_vs_ens_avg=max_arcd_vs_ens_avg,
                                     max_gini_run=max_gini_run,
                                     max_jackpot_quotient=max_jackpot_quotient,
                                     max_loglike_vs_best=max_loglike_vs_best,
                                     min_pearson_vs_best=min_pearson_vs_best,
                                     max_marcd_vs_best=max_marcd_vs_best)
                # Check if the runs pass. Use allow_underclustered=False
                # because this is the criteria used to choose the best
                # number of clusters at the end, so make the best effort
                # now to find clusters that meet those criteria.
                passing = runs_ks[k].passing(allow_underclustered=False)
            assert k in runs_ks
            logger.detail(runs_ks[k].summarize())
            # Output each run's mutation rates and cluster proportions.
            for rank, run in enumerate(runs_k):
                write_mus(run, rank=rank, **path_kwargs)
                write_pis(run, rank=rank, **path_kwargs)
            logger.routine(f"Ended {num_runs_k} ({len(runs_k)} successful) "
                           f"run(s) of EM with {k} cluster(s)")
            # Check if this K is the best encounted so far. Use
            # allow_underclustered=True because the algorithm should
            # continue with the next K if this K is underclustered.
            if not (try_all_ks or k == find_best_k(runs_ks.values(),
                                                   allow_underclustered=True)):
                # The current k is not the best so far.
                break
        except Exception as error:
            logger.error(error)
            if not try_all_ks:
                # Break so that if clustering would fail for every K,
                # this part will not get stuck in a VERY long loop.
                break
    return runs_ks


@with_tmp_dir(pass_keep_tmp=False)
def cluster(dataset: MaskMutsDataset | JoinMaskMutsDataset, *,
            branch: str,
            tmp_dir: Path,
            min_clusters: int,
            max_clusters: int,
            try_all_ks: bool,
            write_all_ks: bool,
            max_em_runs: int,
            num_cpus: int,
            brotli_level: int,
            force: bool,
            cluster_pos_table: bool,
            cluster_abundance_table: bool,
            verify_times: bool,
            **kwargs):
    """ Cluster unique reads from one mask dataset. """
    # Check if the report file already exists.
    branches = path.add_branch(path.CLUSTER_STEP, branch, dataset.branches)
    report_file = ClusterReport.build_path({path.TOP: dataset.top,
                                            path.SAMPLE: dataset.sample,
                                            path.BRANCHES: branches,
                                            path.REF: dataset.ref,
                                            path.REG: dataset.region.name})
    if need_write(report_file, force):
        began = datetime.now()
        # Load the unique reads.
        tmp_clust_dir = path.buildpar(path.REG_DIR_SEGS,
                                      {path.TOP: tmp_dir,
                                       path.SAMPLE: dataset.sample,
                                       path.STEP: path.CLUSTER_STEP,
                                       path.BRANCHES: branches,
                                       path.REF: dataset.ref,
                                       path.REG: dataset.region.name})
        if dataset.min_mut_gap != 4:
            logger.warning("For clustering, it is highly recommended to use "
                           "the observer bias correction with min_mut_gap=4, "
                           f"but got min_mut_gap={dataset.min_mut_gap}")
        uniq_reads = UniqReads.from_dataset_contig(dataset, branch)
        # Run clustering for every number of clusters.
        if max_clusters >= 1:
            max_clusters_use = max_clusters
        elif max_clusters == 0:
            if try_all_ks:
                # Prevent accidentally forcing the use of a huge number
                # of clusters.
                raise ValueError("If using --try-all-ks, "
                                 "then must specify --max-clusters ≥ 1")
            # Set the maximum number of clusters to more than the number
            # of reads: effectively no limit.
            max_clusters_use = dataset.num_reads + 1
        else:
            raise ValueError(
                f"max_clusters must be ≥ 0, but got {max_clusters}"
            )
        runs_ks = run_ks(uniq_reads,
                         ks=range(min_clusters,
                                  max_clusters_use + 1),
                         try_all_ks=try_all_ks,
                         max_em_runs=max_em_runs,
                         num_cpus=num_cpus,
                         top=tmp_dir,
                         **kwargs)
        runs_ks_list = list(runs_ks.values())
        # Choose which numbers of clusters to write.
        if write_all_ks:
            write_ks = runs_ks_list
        elif (best_k := find_best_k(runs_ks_list)) >= 1:
            write_ks = [runs_ks[best_k]]
        else:
            logger.warning(f"No Ks passed filters for {dataset}: "
                           f"defaulting to ensemble average (K = 1)")
            write_ks = [runs_ks[1]]
        # Output the cluster memberships in batches of reads.
        batch_writer = ClusterBatchWriter(dataset,
                                          write_ks,
                                          brotli_level,
                                          tmp_dir,
                                          branch)
        batch_writer.write_batches()
        # Write the observed and expected counts for every best run.
        counts_dir = tmp_clust_dir.joinpath(path.CLUST_COUNTS_DIR)
        counts_dir.mkdir()
        write_obs_exp_counts(uniq_reads, runs_ks_list, counts_dir)
        # Summarize the runs in table and graph format.
        statistics_dir = tmp_clust_dir.joinpath(path.CLUST_STATS_DIR)
        statistics_dir.mkdir()
        write_summaries(runs_ks_list, statistics_dir)
        # Write the cluster report.
        ended = datetime.now()
        report = ClusterReport.from_clusters(runs_ks_list,
                                             uniq_reads,
                                             min_clusters=min_clusters,
                                             max_clusters=max_clusters,
                                             try_all_ks=try_all_ks,
                                             write_all_ks=write_all_ks,
                                             ks_written=batch_writer.ks_written,
                                             max_em_runs=max_em_runs,
                                             checksums=batch_writer.checksums,
                                             began=began,
                                             ended=ended,
                                             **kwargs)
        report_saved = report.save(tmp_dir)
        release_to_out(dataset.top, tmp_dir, report_saved.parent)
        # Write the tables.
        if cluster_pos_table or cluster_abundance_table:
            ClusterDatasetTabulator(
                dataset=ClusterMutsDataset(report_file,
                                           verify_times=verify_times),
                # count_pos must be True because counting positions is
                # required to adjust the counts, which are used by both
                # cluster_pos_table and cluster_abundance_table.
                count_pos=True,
                count_read=False,
                num_cpus=num_cpus,
            ).write_tables(pos=cluster_pos_table, clust=cluster_abundance_table)
    return report_file.parent
