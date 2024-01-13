from datetime import datetime
from logging import getLogger
from pathlib import Path

from click import command

from .compare import RunOrderResults, find_best_order
from .csv import update_log_counts
from .data import load_cluster_dataset
from .io import ClustBatchIO
from .report import ClustReport
from .save import update_batches
from .uniq import UniqReads
from .write import run_orders
from ..core import path
from ..core.arg import (CMD_CLUSTUP,
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


def update_field(report: ClustReport,
                 field: Field,
                 orders: list[RunOrderResults],
                 attr: str):
    """ Merge the field from the original report with the field from the
    new orders. """
    original = report.get_field(field)
    if not isinstance(original, dict):
        raise TypeError(f"Expected dict, but got {type(original).__name__}")
    new = {runs.order: getattr(runs, attr) for runs in orders}
    if overlap := set(original) & set(new):
        raise ValueError(f"Original and new fields share keys: {overlap}")
    return original | new


def add_orders(cluster_report_file: Path,
               max_order: int, *,
               brotli_level: int,
               n_procs: int):
    """ Add orders to an existing report and dataset. """
    # Load the original cluster report.
    original_report = ClustReport.load(cluster_report_file)
    original_max_order = original_report.get_field(MaxClustsF)
    if max_order > original_max_order:
        new_began = datetime.now()
        # Determine clustering parameters from the report.
        n_runs = original_report.get_field(ClustNumRunsF)
        min_iter = original_report.get_field(MinIterClustF)
        max_iter = original_report.get_field(MaxIterClustF)
        conv_thresh = original_report.get_field(ClustConvThreshF)
        prev_bic = min(original_report.get_field(ClustsBicF).values())
        original_best_order = original_report.get_field(NumClustsF)
        original_began = original_report.get_field(TimeBeganF)
        taken = original_report.get_field(TimeTakenF)
        logger.info(f"Began adding clusters to {cluster_report_file} up to "
                    f"order {max_order}, with {n_runs} run(s) per order")
        # Load the unique reads.
        mask_dataset = load_mask_dataset(recast_file_path(cluster_report_file,
                                                          ClustReport,
                                                          MaskReport))
        uniq_reads = UniqReads.from_dataset(mask_dataset)
        # Run clustering for every new order.
        orders = list(run_orders(uniq_reads,
                                 original_max_order + 1,
                                 max_order,
                                 n_runs,
                                 prev_bic=prev_bic,
                                 min_iter=min_iter,
                                 max_iter=max_iter,
                                 conv_thresh=conv_thresh,
                                 n_procs=n_procs,
                                 top=mask_dataset.top))
        if orders:
            best_order = find_best_order(orders)
            # Update the observed and expected counts for each best run.
            try:
                update_log_counts(orders,
                                  top=mask_dataset.top,
                                  sample=mask_dataset.sample,
                                  ref=mask_dataset.ref,
                                  sect=mask_dataset.sect)
            except Exception as error:
                logger.error(f"Failed to update counts: {error}")
            # Output the cluster memberships in batches of reads.
            cluster_dataset = load_cluster_dataset(cluster_report_file)
            checksums = update_batches(cluster_dataset, orders, brotli_level)
            ended = datetime.now()
            new_report = ClustReport(
                sample=cluster_dataset.sample,
                ref=cluster_dataset.ref,
                sect=cluster_dataset.sect,
                n_uniq_reads=uniq_reads.num_uniq,
                max_order=max_order,
                num_runs=n_runs,
                min_iter=min_iter,
                max_iter=max_iter,
                conv_thresh=conv_thresh,
                checksums={ClustBatchIO.btype(): checksums},
                n_batches=len(checksums),
                converged=update_field(original_report,
                                       ClustsConvF,
                                       orders,
                                       "converged"),
                log_likes=update_field(original_report,
                                       ClustsLogLikesF,
                                       orders,
                                       "log_likes"),
                clusts_rmsds=update_field(original_report,
                                          ClustsRMSDsF,
                                          orders,
                                          "rmsds"),
                clusts_meanr=update_field(original_report,
                                          ClustsMeanRsF,
                                          orders,
                                          "meanr"),
                bic=update_field(original_report,
                                 ClustsBicF,
                                 orders,
                                 "bic"),
                best_order=best_order,
                began=original_began,
                ended=ended,
                taken=taken + calc_dt_minutes(new_began, ended),
            )
            new_report.save(cluster_dataset.top, force=True)
        else:
            best_order = original_best_order
        n_new = best_order - original_best_order
        logger.info(f"Ended adding {n_new} cluster(s) to {cluster_report_file}")
    else:
        logger.warning(f"New maximum order ({max_order}) is not greater than "
                       f"original ({original_max_order}): nothing to update")
    return cluster_report_file


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        max_clusters: int,
        brotli_level: int,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Add more clusters to a dataset that was already clustered. """
    if max_clusters == 0:
        # Exit immediately if the maximum number of clusters is 0.
        return list()
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


@command(CMD_CLUSTUP, params=params)
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
