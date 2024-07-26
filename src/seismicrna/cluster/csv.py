from logging import getLogger
from pathlib import Path
from shutil import copy2
from typing import Callable

import pandas as pd

from .compare import EMRunsK, get_log_exp_obs_counts
from .em import EMRun
from .uniq import UniqReads
from ..core import path

logger = getLogger(__name__)

PRECISION = 6  # number of digits behind the decimal point


def get_table_path(top: Path,
                   sample: str,
                   ref: str,
                   sect: str,
                   table: str,
                   k: int,
                   run: int):
    """ Build a path for a table of clustering results. """
    return path.buildpar(*path.CLUST_TAB_SEGS,
                         top=top,
                         cmd=path.CMD_CLUST_DIR,
                         sample=sample,
                         ref=ref,
                         sect=sect,
                         table=table,
                         k=k,
                         run=run,
                         ext=path.CSV_EXT)


def write_single_run_table(run: EMRun,
                           top: Path,
                           sample: str,
                           ref: str,
                           sect: str,
                           rank: int, *,
                           output_func: Callable[[EMRun], pd.DataFrame],
                           table: str):
    """ Write a DataFrame of one type of data from one independent run
    of EM clustering to a CSV file. """
    data = output_func(run)
    file = get_table_path(top, sample, ref, sect, table, run.k, rank)
    data.round(PRECISION).to_csv(file, header=True, index=True)
    logger.info(f"Wrote {table} of {run} to {file}")
    return file


def write_props(run: EMRun,
                top: Path,
                sample: str,
                ref: str,
                sect: str,
                rank: int):
    return write_single_run_table(run,
                                  top,
                                  sample,
                                  ref,
                                  sect,
                                  rank,
                                  output_func=EMRun.get_props,
                                  table=path.CLUST_PROP_RUN_TABLE)


def write_mus(run: EMRun,
              top: Path,
              sample: str,
              ref: str,
              sect: str,
              rank: int):
    return write_single_run_table(run,
                                  top,
                                  sample,
                                  ref,
                                  sect,
                                  rank,
                                  output_func=EMRun.get_mus,
                                  table=path.CLUST_MUS_RUN_TABLE)


def get_count_path(top: Path, sample: str, ref: str, sect: str):
    """ Build a path for a table of bit vector counts. """
    return path.buildpar(*path.CLUST_COUNT_SEGS,
                         top=top,
                         cmd=path.CMD_CLUST_DIR,
                         sample=sample,
                         ref=ref,
                         sect=sect,
                         ext=path.CSVZIP_EXT)


def write_log_counts(uniq_reads: UniqReads,
                     ks: list[EMRunsK],
                     top: Path,
                     sample: str,
                     ref: str,
                     sect: str):
    """ Write the expected and observed log counts of unique bit vectors
    to a CSV file. """
    # Calculate the log expected and observed counts.
    log_counts = get_log_exp_obs_counts(uniq_reads, ks).sort_index()
    # Build the path for the output file.
    file = get_count_path(top, sample, ref, sect)
    # Write the log counts to the file.
    log_counts.to_csv(file)
    return file


def copy_single_run_table(to_dir: Path,
                          from_dir: Path,
                          sample: str,
                          ref: str,
                          sect: str,
                          table: str,
                          k: int,
                          run: int):
    """ Hard-link (if possible -- otherwise, copy) a table for a single
    run to a new directory. """
    src = get_table_path(from_dir, sample, ref, sect, table, k, run)
    dst = get_table_path(to_dir, sample, ref, sect, table, k, run)
    try:
        # First, attempt to hard-link the table.
        dst.hardlink_to(src)
    except OSError:
        # If that fails, then default to copying the table.
        copy2(src, dst)
    return dst


def copy_all_run_tables(to_dir: Path,
                        from_dir: Path,
                        sample: str,
                        ref: str,
                        sect: str,
                        max_k: int,
                        num_runs: int):
    return [copy_single_run_table(to_dir,
                                  from_dir,
                                  sample,
                                  ref,
                                  sect,
                                  table,
                                  k,
                                  run)
            for k in range(1, max_k + 1)
            for run in range(num_runs if k > 1 else 1)
            for table in (path.CLUST_MUS_RUN_TABLE, path.CLUST_PROP_RUN_TABLE)]

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
