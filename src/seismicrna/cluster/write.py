from functools import partial
from logging import getLogger
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pandas as pd

from .compare import (get_common_best_run_attr, get_log_exp_obs_counts,
                      find_best_order, RunOrderResults)
from .em import EmClustering
from .names import ORD_CLS_NAME, READ_NAME
from .report import ClustReport
from ..core import path
from ..core.files import digest_file
from ..mask.load import MaskLoader

logger = getLogger(__name__)

PRECISION = 6  # number of digits behind the decimal point


def get_table_path(out_dir: Path, sample: str, ref: str, sect: str,
                   table: str, k: int, run: int):
    """ Build a path for a table of clustering results. """
    return path.buildpar(path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg,
                         path.ClustTabSeg, top=out_dir, module=path.MOD_CLUST,
                         sample=sample, ref=ref, sect=sect, table=table, k=k,
                         run=run, ext=path.CSV_EXT)


def write_single_run_table(run: EmClustering,
                           rank: int,
                           output_func: Callable[[EmClustering], pd.DataFrame],
                           table: str):
    """ Write a DataFrame of one type of data from one independent run
    of EM clustering to a CSV file. """
    data = output_func(run)
    file = get_table_path(run.loader.out_dir, run.loader.sample, run.loader.ref,
                          run.loader.sect, table, run.order, rank)
    data.round(PRECISION).to_csv(file, header=True, index=True)
    logger.info(f"Wrote {table} of {run} to {file}")
    return file


write_props = partial(write_single_run_table,
                      output_func=EmClustering.output_props,
                      table=path.CLUST_PROP_RUN_TABLE)

write_mus = partial(write_single_run_table,
                    output_func=EmClustering.output_mus,
                    table=path.CLUST_MUS_RUN_TABLE)


def write_batch(resps: dict[int, pd.DataFrame], read_names: Sequence[str],
                out_dir: Path, sample: str, ref: str, sect: str, batch: int):
    """ Write the memberships of the reads in one batch to a file. """
    # Determine the path to the batch file.
    batch_file = ClustReport.build_batch_path(out_dir, batch,
                                              sample=sample, ref=ref,
                                              sect=sect, ext=path.CSVZIP_EXT)
    # Assemble the memberships of the selected reads into a DataFrame.
    batch_data = (((order, cluster),
                   np.round(cluster_resps.loc[read_names], PRECISION))
                  for order, ord_resps in resps.items()
                  for cluster, cluster_resps in ord_resps.items())
    batch_resps = pd.DataFrame.from_dict(dict(batch_data))
    # Rename the index.
    batch_resps.index.rename(READ_NAME, inplace=True)
    # Rename the column levels.
    batch_resps.columns.rename(ORD_CLS_NAME, inplace=True)
    # Write the memberships to the batch file.
    batch_resps.to_csv(batch_file, index=True)
    logger.debug(f"Wrote cluster memberships of batch {batch} to {batch_file}")
    # Return the checksum of the batch file.
    return digest_file(batch_file)


def write_batches(ord_runs: dict[int, RunOrderResults]):
    """ Write the cluster memberships to batch files. """
    # Get the data loader for the clustering runs.
    loader: MaskLoader = get_common_best_run_attr(ord_runs, "loader")
    # Compute cluster memberships up to and including the best order.
    resps = {order: ord_runs[order].best.output_resps()
             for order in range(1, find_best_order(ord_runs) + 1)}
    # Write the memberships of each batch of reads to a file and return
    # a list of the checksum of every file.
    return [write_batch(resps, read_names, loader.out_dir, loader.sample,
                        loader.ref, loader.sect, batch)
            for batch, read_names in enumerate(loader.iter_read_batches())]


def get_count_path(out_dir: Path, sample: str, ref: str, sect: str):
    """ Build a path for a table of bit vector counts. """
    return path.buildpar(path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg,
                         path.ClustCountSeg, top=out_dir, module=path.MOD_CLUST,
                         sample=sample, ref=ref, sect=sect, ext=path.CSVZIP_EXT)


def write_log_counts(ord_runs: dict[int, RunOrderResults]):
    """ Write the expected and observed log counts of unique bit vectors
    to a CSV file. """
    # Compute the log expected and observed counts.
    log_counts = get_log_exp_obs_counts(ord_runs)
    # Get the data loader for the clustering runs.
    loader = get_common_best_run_attr(ord_runs, "loader")
    # Build the path for the output file.
    file = get_count_path(loader.out_dir, loader.sample,
                          loader.ref, loader.sect)
    # Write the counts to the file.
    log_counts.to_csv(file)
    return file
