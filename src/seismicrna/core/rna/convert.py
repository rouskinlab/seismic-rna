from logging import getLogger

from click import command

from .io import ct_to_db, db_to_ct
from .. import path
from ..arg import (CMD_CT2DB,
                   CMD_DB2CT,
                   arg_input_path,
                   opt_force,
                   opt_max_procs,
                   opt_parallel)
from ..run import run_func
from ..task import as_list_of_tuples, dispatch

logger = getLogger(__name__)


@run_func(logger.critical)
def run_ct_to_db(input_path: tuple[str, ...], *,
                 force: bool,
                 max_procs: int,
                 parallel: bool):
    """ Convert connectivity table (CT) to dot-bracket (DB) files. """
    ct_files = list(path.find_files_chain(input_path, [path.ConnectTableSeg]))
    return dispatch(ct_to_db,
                    max_procs,
                    parallel,
                    args=as_list_of_tuples(ct_files),
                    kwargs=dict(force=force),
                    pass_n_procs=False)


@run_func(logger.critical)
def run_db_to_ct(input_path: tuple[str, ...], *,
                 force: bool,
                 max_procs: int,
                 parallel: bool):
    """ Convert dot-bracket (DB) to connectivity table (CT) files. """
    db_files = list(path.find_files_chain(input_path, [path.DotBracketSeg]))
    return dispatch(db_to_ct,
                    max_procs,
                    parallel,
                    args=as_list_of_tuples(db_files),
                    kwargs=dict(force=force),
                    pass_n_procs=False)


params = [
    arg_input_path,
    opt_force,
    opt_max_procs,
    opt_parallel
]


@command(CMD_CT2DB, params=params)
def cli_ct2db(*args, **kwargs):
    """ Convert connectivity table (CT) to dot-bracket (DB) files. """
    run_ct_to_db(*args, **kwargs)


@command(CMD_DB2CT, params=params)
def cli_db2ct(*args, **kwargs):
    """ Convert dot-bracket (DB) to connectivity table (CT) files. """
    run_db_to_ct(*args, **kwargs)
