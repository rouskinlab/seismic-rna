from datetime import datetime
from logging import getLogger
from math import inf
from pathlib import Path

import pandas as pd

from .compare import RunOrderResults, find_best_order, sort_replicate_runs
from .csv import write_log_counts, write_mus, write_props
from .em import EmClustering
from .io import ClustBatchIO
from .report import ClustReport
from .uniq import UniqReads
from ..core.parallel import dispatch
from ..core.write import need_write
from ..mask.data import MaskMerger

logger = getLogger(__name__)


def write_batches(dataset: MaskMerger,
                  ord_runs: dict[int, RunOrderResults],
                  brotli_level: int):
    """ Write the cluster memberships to batch files. """
    checksums = list()
    best_order = find_best_order(ord_runs)
    for batch_num in dataset.batch_nums:
        resps = pd.concat((ord_runs[order].best.get_resps(batch_num)
                           for order in range(1, best_order + 1)),
                          axis=1)
        batch = ClustBatchIO(sample=dataset.sample,
                             ref=dataset.ref,
                             sect=dataset.sect,
                             batch=batch_num,
                             resps=resps)
        _, checksum = batch.save(top=dataset.top,
                                 brotli_level=brotli_level,
                                 overwrite=True)
        checksums.append(checksum)
    return checksums


def run_order(uniq_reads: UniqReads,
              order: int,
              n_runs: int, *,
              min_iter: int,
              max_iter: int,
              conv_thresh: float,
              n_procs: int) -> list[EmClustering]:
    """ Run EM with a specific number of clusters. """
    logger.info(f"Began {n_runs} run(s) of EM with {order} cluster(s)")
    # Initialize one EmClustering object for each replicate run.
    replicates = [EmClustering(uniq_reads,
                               order,
                               conv_thresh=conv_thresh,
                               min_iter=min_iter,
                               max_iter=max_iter)
                  for _ in range(n_runs)]
    # Run independent replicates of the clustering algorithm.
    replicates = dispatch([rep.run for rep in replicates],
                          n_procs,
                          parallel=True,
                          pass_n_procs=False)
    logger.info(f"Ended {n_runs} run(s) of EM with {order} cluster(s)")
    return sort_replicate_runs(replicates)


def run_max_order(uniq_reads: UniqReads,
                  max_order: int,
                  n_runs: int, *,
                  min_iter: int,
                  max_iter: int,
                  conv_thresh: float,
                  n_procs: int,
                  top: Path):
    """
    Find the optimal order (i.e. number of clusters) for EM, up to max_order.
    """
    if n_runs < 1:
        logger.warning(f"Number of EM runs must be ≥ 1: setting to 1")
        n_runs = 1
    logger.info(f"Began clustering {uniq_reads} up to order {max_order} "
                f"with {n_runs} runs per order")
    path_kwargs = dict(top=top,
                       sample=uniq_reads.sample,
                       ref=uniq_reads.ref,
                       sect=uniq_reads.section.name)
    # For each order, keep the best run and a summary of all runs.
    results = dict()
    # Run clustering for each number of clusters, up to max_clusters.
    while len(results) < max_order:
        # Determine the current order for clustering.
        order = len(results) + 1
        # Run EM clustering n_runs times with different starting points.
        runs = run_order(uniq_reads,
                         order,
                         n_runs=(n_runs if order > 1 else 1),
                         conv_thresh=(conv_thresh if order > 1 else inf),
                         min_iter=(min_iter * order if order > 1 else 2),
                         max_iter=(max_iter * order if order > 1 else 2),
                         n_procs=n_procs)
        # Output tables of the mutation rates and cluster proportions
        # for every run.
        for rank, run in enumerate(runs):
            write_mus(run, rank=rank, **path_kwargs)
            write_props(run, rank=rank, **path_kwargs)
        # Compute a summary of all runs.
        results[order] = RunOrderResults(runs)
        # Compare the BIC for this order to lower orders, if order > 1.
        best_bic = results[order].best.bic
        logger.debug(f"The best BIC with {order} cluster(s) is {best_bic}")
        if order > 1:
            prev_bic = results[order - 1].best.bic
            # Check if the best BIC is better (smaller) than the best
            # BIC from clustering with one cluster fewer.
            if best_bic < prev_bic:
                # If the best BIC is < the previous best BIC, then the
                # current model is the best.
                logger.debug(f"The BIC decreased from {prev_bic} "
                             f"(order {order - 1}) to {best_bic} "
                             f"(order {order})")
            else:
                # If the best BIC is ≥ than the previous best BIC, then
                # this model is worse than the previous model. Stop.
                logger.info(f"The BIC failed to decrease from {prev_bic} "
                            f"(order {order - 1}) to {best_bic} "
                            f"(order {order}): stopping")
                break
    return results


def cluster(mask_report_file: Path,
            max_order: int,
            n_runs: int, *,
            min_muts: int,
            min_iter: int,
            max_iter: int,
            conv_thresh: float,
            brotli_level: int,
            n_procs: int,
            force: bool):
    """ Run all processes of clustering reads from one filter. """
    dataset = MaskMerger.load(mask_report_file)
    path_kwargs = dict(top=dataset.top,
                       sample=dataset.sample,
                       ref=dataset.ref,
                       sect=dataset.sect)
    # Check if the clustering report file already exists.
    cluster_report_file = ClustReport.build_path(**path_kwargs)
    if need_write(cluster_report_file, force):
        began = datetime.now()
        logger.info(f"Began clustering {dataset} up to order {max_order} "
                    f"cluster(s) and {n_runs} independent run(s) per order")
        # Get the unique reads.
        uniq_reads = UniqReads.from_dataset(dataset)
        # Run EM clustering for every number of clusters.
        results = run_max_order(uniq_reads,
                                max_order,
                                n_runs,
                                min_iter=min_iter,
                                max_iter=max_iter,
                                conv_thresh=conv_thresh,
                                n_procs=n_procs,
                                top=dataset.top)
        logger.info(
            f"Ended clustering {dataset}: {find_best_order(results)} clusters")
        # Output the observed and expected counts for every best run.
        write_log_counts(results, **path_kwargs)
        # Output the cluster memberships in batches of reads.
        checksums = write_batches(dataset, results, brotli_level)
        ended = datetime.now()
        report = ClustReport.from_clusters(results,
                                           uniq_reads,
                                           max_order,
                                           n_runs,
                                           min_iter=min_iter,
                                           max_iter=max_iter,
                                           conv_thresh=conv_thresh,
                                           checksums=checksums,
                                           began=began,
                                           ended=ended)
        report.save(top=dataset.top, overwrite=True)
    return cluster_report_file

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
