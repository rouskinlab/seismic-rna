from pathlib import Path

from .em import EMRun
from ..core import path
from ..core.logs import logger


PRECISION = 6  # number of digits behind the decimal point


def get_table_path(top: Path,
                   sample: str,
                   ref: str,
                   reg: str,
                   table: str,
                   k: int,
                   run: int):
    """ Build a path for a table of clustering results. """
    return path.buildpar(*path.CLUST_TAB_SEGS,
                         top=top,
                         cmd=path.CMD_CLUST_DIR,
                         sample=sample,
                         ref=ref,
                         reg=reg,
                         table=table,
                         k=k,
                         run=run,
                         ext=path.CSV_EXT)


def write_single_run_table(run: EMRun,
                           top: Path,
                           sample: str,
                           ref: str,
                           reg: str,
                           rank: int, *,
                           attr: str,
                           table: str):
    """ Write a DataFrame of one type of data from one independent run
    of EM clustering to a CSV file. """
    data = getattr(run, attr)
    file = get_table_path(top, sample, ref, reg, table, run.k, rank)
    data.round(PRECISION).to_csv(file, header=True, index=True)
    logger.routine(f"Wrote {table} of {run} to {file}")
    return file


def write_pis(run: EMRun,
              top: Path,
              sample: str,
              ref: str,
              reg: str,
              rank: int):
    return write_single_run_table(run,
                                  top,
                                  sample,
                                  ref,
                                  reg,
                                  rank,
                                  attr="pis",
                                  table=path.CLUST_PARAM_PIS)


def write_mus(run: EMRun,
              top: Path,
              sample: str,
              ref: str,
              reg: str,
              rank: int):
    return write_single_run_table(run,
                                  top,
                                  sample,
                                  ref,
                                  reg,
                                  rank,
                                  attr="mus",
                                  table=path.CLUST_PARAM_MUS)

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
