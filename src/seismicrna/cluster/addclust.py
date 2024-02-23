from datetime import datetime
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd
from click import command

from .compare import RunOrderResults, get_log_exp_obs_counts
from .compare import find_best_order
from .csv import get_count_path
from .data import ClusterMutsDataset
from .data import load_cluster_dataset
from .io import ClusterBatchIO
from .names import BIT_VECTOR_NAME, LOG_OBS_NAME
from .report import ClusterReport
from .uniq import UniqReads
from .write import run_orders
from ..core import path
from ..core.arg import (CMD_ADDCLUST,
                        docdef,
                        arg_input_path,
                        opt_max_clusters,
                        opt_brotli_level,
                        opt_parallel,
                        opt_max_procs)
from ..core.io import recast_file_path
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.report import (calc_dt_minutes,
                           Field,
                           ChecksumsF,
                           NumUniqReadKeptF,
                           TimeBeganF,
                           TimeTakenF,
                           MaxClustsF,
                           ClustNumRunsF,
                           MinIterClustF,
                           MaxIterClustF,
                           ClustConvThreshF,
                           ClustsConvF,
                           ClustsLogLikesF,
                           ClustsRMSDsF,
                           ClustsMeanRsF,
                           ClustsBicF,
                           NumClustsF)
from ..mask.data import load_mask_dataset
from ..mask.report import MaskReport

logger = getLogger(__name__)

BTYPE = ClusterBatchIO.btype()


def update_log_counts(new_orders: list[RunOrderResults],
                      top: Path,
                      sample: str,
                      ref: str,
                      sect: str):
    """ Update the expected log counts of unique bit vectors. """
    # Compute the log expected counts of the new orders.
    new_log_counts = get_log_exp_obs_counts(new_orders).sort_index()
    new_log_obs = new_log_counts.pop(LOG_OBS_NAME)
    # Build the path for the output file.
    file = get_count_path(top, sample, ref, sect)
    # Load the existing counts.
    original_log_counts = pd.read_csv(file, index_col=BIT_VECTOR_NAME)
    # Ensure the index and log observed counts match.
    if not new_log_counts.index.equals(original_log_counts.index):
        raise ValueError("Bit vectors differ between original "
                         f"({original_log_counts.index}) and new "
                         f"({new_log_counts.index}) orders")
    if not np.allclose(new_log_obs, original_log_counts[LOG_OBS_NAME]):
        raise ValueError("Log observed counts differ between original "
                         f"({original_log_counts[LOG_OBS_NAME]}) and new "
                         f"({new_log_obs}) orders")
    # Merge the original and new log counts.
    log_counts = pd.concat([original_log_counts, new_log_counts],
                           axis=1,
                           verify_integrity=True)
    # Write the updated log counts to the file.
    log_counts.to_csv(file)
    return file


def update_batches(dataset: ClusterMutsDataset,
                   new_orders: list[RunOrderResults],
                   brotli_level: int):
    """ Update the cluster memberships in batches. """
    checksums = list()
    for batch in dataset.iter_batches():
        # Merge the original responsibilities with the new ones.
        resps = pd.concat([batch.resps] + [runs.best.get_resps(batch.batch)
                                           for runs in new_orders],
                          axis=1,
                          verify_integrity=True)
        batch = ClusterBatchIO(sample=dataset.sample,
                               ref=dataset.ref,
                               sect=dataset.sect,
                               batch=batch.batch,
                               resps=resps)
        _, checksum = batch.save(top=dataset.top,
                                 brotli_level=brotli_level,
                                 force=True)
        checksums.append(checksum)
    return checksums


def update_field(report: ClusterReport,
                 field: Field,
                 new_orders: list[RunOrderResults],
                 attr: str):
    """ Merge the field from the original report with the field from the
    new orders. """
    original = report.get_field(field)
    if not isinstance(original, dict):
        raise TypeError(f"Expected dict, but got {type(original).__name__}")
    new = {runs.order: getattr(runs, attr) for runs in new_orders}
    if overlap := set(original) & set(new):
        raise ValueError(f"Original and new fields share keys: {overlap}")
    return original | new


def update_report(original_report: ClusterReport,
                  max_order: int,
                  best_order: int,
                  new_orders: list[RunOrderResults],
                  checksums: list[str],
                  began: datetime,
                  ended: datetime,
                  top: Path):
    new_report = ClusterReport(
        sample=original_report.sample,
        ref=original_report.ref,
        sect=original_report.sect,
        n_uniq_reads=original_report.get_field(NumUniqReadKeptF),
        max_order=max_order,
        num_runs=original_report.get_field(ClustNumRunsF),
        min_iter=original_report.get_field(MinIterClustF),
        max_iter=original_report.get_field(MaxIterClustF),
        conv_thresh=original_report.get_field(ClustConvThreshF),
        checksums={BTYPE: checksums},
        n_batches=len(checksums),
        converged=update_field(original_report,
                               ClustsConvF,
                               new_orders,
                               "converged"),
        log_likes=update_field(original_report,
                               ClustsLogLikesF,
                               new_orders,
                               "log_likes"),
        clusts_rmsds=update_field(original_report,
                                  ClustsRMSDsF,
                                  new_orders,
                                  "rmsds"),
        clusts_meanr=update_field(original_report,
                                  ClustsMeanRsF,
                                  new_orders,
                                  "meanr"),
        bic=update_field(original_report,
                         ClustsBicF,
                         new_orders,
                         "bic"),
        best_order=best_order,
        began=original_report.get_field(TimeBeganF),
        ended=ended,
        taken=(original_report.get_field(TimeTakenF)
               + calc_dt_minutes(began, ended)),
    )
    new_report.save(top, force=True)
    return new_report


def add_orders(report_file: Path,
               max_order: int, *,
               brotli_level: int,
               n_procs: int):
    """ Add orders to an existing report and dataset. """
    # Load the original cluster report.
    report = ClusterReport.load(report_file)
    original_max_order = report.get_field(MaxClustsF)
    if max_order > original_max_order:
        began = datetime.now()
        logger.info(
            f"Began adding clusters to {report_file} up to order {max_order}"
        )
        dataset = load_mask_dataset(recast_file_path(report_file,
                                                     ClusterReport,
                                                     MaskReport))
        original_best_order = report.get_field(NumClustsF)
        if original_best_order == original_max_order:
            # The original clustering stopped because the order reached
            # the maximum order, not because the BIC decreased.
            uniq_reads = UniqReads.from_dataset(dataset)
            # Add more clusters up to max_order.
            new_orders = list(run_orders(
                uniq_reads,
                original_max_order + 1,
                max_order,
                report.get_field(ClustNumRunsF),
                prev_bic=min(report.get_field(ClustsBicF).values()),
                min_iter=report.get_field(MinIterClustF),
                max_iter=report.get_field(MaxIterClustF),
                conv_thresh=report.get_field(ClustConvThreshF),
                n_procs=n_procs,
                top=dataset.top)
            )
            if new_orders:
                # One or more higher orders had better BICs: add them.
                best_order = find_best_order(new_orders)
                # Update the expected counts for each best run.
                update_log_counts(new_orders,
                                  top=dataset.top,
                                  sample=dataset.sample,
                                  ref=dataset.ref,
                                  sect=dataset.sect)
                # Output the cluster memberships in batches of reads.
                cluster_dataset = load_cluster_dataset(report_file)
                checksums = update_batches(cluster_dataset,
                                           new_orders,
                                           brotli_level)
            else:
                # No higher orders had better BICs.
                best_order = original_best_order
                checksums = report.get_field(ChecksumsF)[BTYPE]
        elif original_best_order < original_max_order:
            # The original clustering stopped because the BIC decreased
            # before reaching the maximum order, so no more clusters
            # should be added.
            new_orders = list()
            best_order = original_best_order
            checksums = report.get_field(ChecksumsF)[BTYPE]
        else:
            raise ValueError(
                f"The best order ({original_best_order}) exceeds the maximum "
                f"order ({original_max_order}) in {report_file}"
            )
        ended = datetime.now()
        update_report(report,
                      max_order,
                      best_order,
                      new_orders,
                      checksums,
                      began,
                      ended,
                      dataset.top)
        logger.info(
            f"Ended adding clusters to {report_file} up to order {max_order}"
        )
    else:
        logger.warning(f"New maximum order ({max_order}) is not greater than "
                       f"original ({original_max_order}): nothing to add")
    return report_file


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        max_clusters: int,
        brotli_level: int,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Add more clusters to a dataset that was already clustered. """
    # Find the cluster report files.
    report_files = path.find_files_chain(
        input_path, load_cluster_dataset.report_path_seg_types
    )
    # Cluster each mask dataset.
    return dispatch(add_orders,
                    max_procs,
                    parallel,
                    pass_n_procs=True,
                    args=as_list_of_tuples(report_files),
                    kwargs=dict(max_order=max_clusters,
                                brotli_level=brotli_level))


params = [
    # Input files
    arg_input_path,
    # Clustering options
    opt_max_clusters,
    # Compression
    opt_brotli_level,
    # Parallelization
    opt_max_procs,
    opt_parallel,
]


@command(CMD_ADDCLUST, params=params)
def cli(*args, **kwargs):
    """ Add more clusters to a dataset that was already clustered. """
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
