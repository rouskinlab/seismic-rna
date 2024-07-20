from datetime import datetime
from logging import getLogger
from math import inf
from pathlib import Path

import numpy as np

from .compare import CompareRunsK, find_best_order, sort_replicate_runs
from .csv import write_log_counts, write_mus, write_props
from .em import EmClustering
from .report import ClusterReport
from .save import write_batches
from .uniq import UniqReads
from ..core.io import recast_file_path
from ..core.task import as_list_of_tuples, dispatch
from ..core.tmp import release_to_out
from ..core.types import get_max_uint
from ..core.write import need_write
from ..mask.data import load_mask_dataset
from ..mask.report import MaskReport

logger = getLogger(__name__)

SEED_DTYPE = np.uint32


def run_k(uniq_reads: UniqReads,
          k: int,
          n_runs: int, *,
          n_procs: int,
          **kwargs) -> list[EmClustering]:
    """ Run EM with a specific number of clusters. """
    logger.info(f"Began {n_runs} run(s) of EM with {k} cluster(s)")
    # Initialize one EmClustering object for each replicate run.
    replicates = [EmClustering(uniq_reads, k, **kwargs)
                  for _ in range(n_runs)]
    # Run independent replicates of the clustering algorithm.
    rng = np.random.default_rng()
    seeds = as_list_of_tuples(rng.integers(get_max_uint(SEED_DTYPE),
                                           size=n_runs,
                                           dtype=SEED_DTYPE))
    replicates = dispatch([rep.run for rep in replicates],
                          n_procs,
                          parallel=True,
                          pass_n_procs=False,
                          args=seeds)
    logger.info(f"Ended {n_runs} run(s) of EM with {k} cluster(s)")
    return sort_replicate_runs(replicates)


def run_ks(uniq_reads: UniqReads,
           min_clusters: int,
           max_clusters: int,
           n_runs: int, *,
           use_bic: bool = True,
           use_clust_corr: bool = False,
           prev_bic: float | None = None,
           min_iter: int,
           max_iter: int,
           em_thresh: float,
           n_procs: int,
           top: Path,
           **kwargs):
    """ Find the optimal order, from min_order to max_order. """
    if n_runs < 1:
        logger.warning(f"Number of runs must be ≥ 1: setting to 1")
        n_runs = 1
    path_kwargs = dict(top=top,
                       sample=uniq_reads.sample,
                       ref=uniq_reads.ref,
                       sect=uniq_reads.section.name)
    # Run clustering for each order from min_order to max_order.
    if min_clusters < 1:
        logger.warning(f"min_clusters must be ≥ 1, but got {min_clusters}; "
                       f"setting min_clusters to 1")
        min_clusters = 1
    if max_clusters < 0:
        logger.warning(f"max_clusters must be ≥ 0, but got {max_clusters}; "
                       f"setting max_clusters to 0")
        max_clusters = 0
    k = min_clusters
    clustering = True
    while clustering and (max_clusters == 0 or k <= max_clusters):
        # Cluster n_runs times with different starting points.
        runs = run_k(uniq_reads,
                     k,
                     n_runs=(n_runs if k > 1 else 1),
                     em_thresh=(em_thresh if k > 1 else inf),
                     min_iter=(min_iter * k if k > 1 else 2),
                     max_iter=(max_iter * k if k > 1 else 2),
                     n_procs=n_procs,
                     **kwargs)
        # Output tables of the mutation rates and cluster proportions
        # for every run.
        for rank, run in enumerate(runs):
            write_mus(run, rank=rank, **path_kwargs)
            write_props(run, rank=rank, **path_kwargs)
        # Compare all runs for this k.
        comparison = CompareRunsK(runs)
        if use_bic:
            # Compare the BIC for this k to lower ks, if any.
            best_bic = comparison.bic
            logger.debug(f"The best BIC with k={k} is {best_bic}")
            if prev_bic is not None:
                # Check if the BIC is better (smaller) than the BIC of
                # the previous k.
                if best_bic < prev_bic:
                    # If the BIC is less than the previous BIC, then
                    # this k is better than the previous k.
                    logger.debug(f"The BIC decreased from {prev_bic} "
                                 f"(k={k - 1}) to {best_bic} (k={k})")
                else:
                    # If the best BIC is ≥ the previous best BIC, then
                    # this k is worse than the previous k: stop.
                    logger.info(f"The BIC failed to decrease from {prev_bic} "
                                f"(k={k - 1}) to {best_bic} (k={k})")
                    clustering = False
            # Record the best BIC for comparison with the next k.
            prev_bic = best_bic
        yield comparison
        k += 1


def cluster(mask_report_file: Path,
            min_clusters: int,
            max_clusters: int,
            n_runs: int, *,
            n_procs: int,
            brotli_level: int,
            force: bool,
            tmp_dir: Path,
            **kwargs):
    """ Cluster unique reads from one mask dataset. """
    # Check if the cluster report file already exists.
    cluster_report_file = recast_file_path(mask_report_file,
                                           MaskReport,
                                           ClusterReport)
    if need_write(cluster_report_file, force):
        began = datetime.now()
        logger.info(f"Began clustering {mask_report_file}")
        # Load the unique reads.
        dataset = load_mask_dataset(mask_report_file)
        if dataset.min_mut_gap != 3:
            logger.warning("For clustering, it is highly recommended to use "
                           "the observer bias correction with min_mut_gap=3, "
                           f"but got min_mut_gap={dataset.min_mut_gap}")
        uniq_reads = UniqReads.from_dataset_contig(dataset)
        # Run clustering for every order.
        orders = list(run_ks(uniq_reads,
                             min_clusters=min_clusters,
                             max_clusters=max_clusters,
                             n_runs=n_runs,
                             n_procs=n_procs,
                             top=tmp_dir,
                             **kwargs))
        # Output the observed and expected counts for every best run.
        write_log_counts(orders,
                         top=tmp_dir,
                         sample=dataset.sample,
                         ref=dataset.ref,
                         sect=dataset.sect)
        # Output the cluster memberships in batches of reads.
        checksums = write_batches(dataset, orders, brotli_level, tmp_dir)
        ended = datetime.now()
        report = ClusterReport.from_clusters(orders,
                                             uniq_reads,
                                             max_clusters,
                                             n_runs,
                                             checksums=checksums,
                                             began=began,
                                             ended=ended,
                                             **kwargs)
        report_saved = report.save(tmp_dir)
        release_to_out(dataset.top, tmp_dir, report_saved.parent)
        order = find_best_order(orders)
        logger.info(f"Ended clustering {mask_report_file} to order {order}")
    return cluster_report_file

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
