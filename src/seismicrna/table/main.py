from itertools import chain
from logging import getLogger
from pathlib import Path

from click import command

from .write import write
from ..core import docdef, path
from ..core.cli import (opt_report, opt_table_cols,
                        opt_max_procs, opt_parallel, opt_rerun)
from ..core.parallel import dispatch

logger = getLogger(__name__)

params = [opt_report, opt_table_cols, opt_max_procs, opt_parallel, opt_rerun]


@command(path.MOD_TABLE, params=params)
def cli(*args, **kwargs):
    """ Tabulate per-base and per-read mutation counts from 'relate' and
    'mask', and mutation rates and read memberships from 'cluster'. """
    return run(*args, **kwargs)


@docdef.auto()
def run(report: tuple[str, ...], table_cols: str,
        max_procs: int, parallel: bool, **kwargs):
    """
    Run the table module.
    """
    files = list(map(Path, report))
    rels = path.find_files_chain(files, [path.RelateRepSeg])
    if rels:
        logger.debug(f"Found relate report files: {rels}")
    masks = path.find_files_chain(files, [path.MaskRepSeg])
    if masks:
        logger.debug(f"Found mask report files: {masks}")
    clusts = path.find_files_chain(files, [path.ClustRepSeg])
    if clusts:
        logger.debug(f"Found cluster report files: {clusts}")
    tasks = [(file, table_cols) for file in chain(rels, masks, clusts)]
    return list(chain(*dispatch(write, max_procs, parallel,
                                args=tasks, kwargs=kwargs,
                                pass_n_procs=False)))
