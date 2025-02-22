from pathlib import Path

from .em import EMRun
from ..core import path
from ..core.logs import logger


PRECISION = 6  # number of digits behind the decimal point


def get_table_path(top: Path,
                   branches: dict[str, str],
                   sample: str,
                   ref: str,
                   reg: str,
                   table: str,
                   k: int,
                   run: int):
    """ Build a path for a table of clustering results. """
    return path.buildpar(path.CLUST_TAB_SEGS,
                         {path.TOP: top,
                          path.SAMPLE: sample,
                          path.STEP: path.CLUSTER_STEP,
                          path.BRANCHES: branches,
                          path.REF: ref,
                          path.REG: reg,
                          path.TABLE: table,
                          path.NCLUST: k,
                          path.RUN: run,
                          path.EXT: path.CSV_EXT})


def write_single_run_table(run: EMRun,
                           top: Path,
                           branches: dict[str, str],
                           sample: str,
                           ref: str,
                           reg: str,
                           rank: int, *,
                           attr: str,
                           table: str):
    """ Write a DataFrame of one type of data from one independent run
    of EM clustering to a CSV file. """
    data = getattr(run, attr)
    file = get_table_path(top, branches, sample, ref, reg, table, run.k, rank)
    data.round(PRECISION).to_csv(file, header=True, index=True)
    logger.action(f"Wrote {table} of {run} to {file}")
    return file


def write_pis(run: EMRun,
              top: Path,
              branches: dict[str, str],
              sample: str,
              ref: str,
              reg: str,
              rank: int):
    return write_single_run_table(run,
                                  top,
                                  branches,
                                  sample,
                                  ref,
                                  reg,
                                  rank,
                                  attr="pis",
                                  table=path.CLUST_PARAM_PIS)


def write_mus(run: EMRun,
              top: Path,
              branches: dict[str, str],
              sample: str,
              ref: str,
              reg: str,
              rank: int):
    return write_single_run_table(run,
                                  top,
                                  branches,
                                  sample,
                                  ref,
                                  reg,
                                  rank,
                                  attr="mus",
                                  table=path.CLUST_PARAM_MUS)
