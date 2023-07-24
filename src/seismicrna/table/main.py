from collections import Counter
from itertools import chain
from logging import getLogger
from pathlib import Path
from typing import Iterable

from click import command

from .write import write
from ..core import docdef, path
from ..core.cli import (arg_input_file, opt_rels,
                        opt_max_procs, opt_parallel, opt_rerun)
from ..core.parallel import dispatch

logger = getLogger(__name__)

params = [arg_input_file, opt_rels, opt_max_procs, opt_parallel, opt_rerun]


@command(path.MOD_TABLE, params=params)
def cli(*args, **kwargs):
    """ Tabulate per-base and per-read mutation counts from 'relate' and
    'mask', and mutation rates and read memberships from 'cluster'. """
    return run(*args, **kwargs)


@docdef.auto()
def run(input_file: tuple[str, ...], rels: tuple[str, ...],
        max_procs: int, parallel: bool, **kwargs):
    """
    Run the table module.
    """
    report_files = list(map(Path, input_file))
    relate_reports = path.find_files_chain(report_files, [path.RelateRepSeg])
    if relate_reports:
        logger.debug(f"Found relate report files: {relate_reports}")
    mask_reports = path.find_files_chain(report_files, [path.MaskRepSeg])
    if mask_reports:
        logger.debug(f"Found mask report files: {mask_reports}")
    clust_reports = path.find_files_chain(report_files, [path.ClustRepSeg])
    if clust_reports:
        logger.debug(f"Found cluster report files: {clust_reports}")
    rels_str = join_rels(rels)
    tasks = [(file, rels_str) for file in chain(relate_reports,
                                                mask_reports,
                                                clust_reports)]
    return list(chain(*dispatch(write, max_procs, parallel,
                                args=tasks, kwargs=kwargs,
                                pass_n_procs=False)))


def join_rels(rels: Iterable[str]):
    """ Join relationships into one string. """
    counts = Counter(r for rel in rels for r in rel)
    if dups := [item for item, count in counts.items() if count > 1]:
        logger.warning(f"Got duplicate relationships: {dups}")
    return "".join(counts)
