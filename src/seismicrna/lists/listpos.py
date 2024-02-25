from pathlib import Path

import pandas as pd
from click import command

from ..core import path
from ..core.arg import (CMD_LISTPOS,
                        docdef,
                        arg_input_path,
                        opt_complement,
                        opt_max_fmut_pos,
                        opt_force,
                        opt_max_procs,
                        opt_parallel)
from ..core.parallel import as_list_of_tuples, dispatch
from ..core.seq import POS_NAME
from ..core.write import need_write
from ..table.base import MUTAT_REL, PosTable
from ..table.load import find_pos_tables, load_pos_table

REF_NAME = "Reference"


def find_pos(table: PosTable,
             max_fmut_pos: float,
             complement: bool):
    # Initially select all unmasked positions.
    section = table.section.copy()
    positions = section.unmasked_int
    # Apply each filter.
    section.add_mask(
        "max_fmut_pos",
        positions[table.fetch_ratio(rel=MUTAT_REL,
                                    exclude_masked=True,
                                    squeeze=True)
                  > max_fmut_pos],
        complement=complement
    )
    return section.unmasked_int


def list_pos(table_file: Path, force: bool, **kwargs):
    """ List positions meeting specific criteria from the table. """
    table = load_pos_table(table_file)
    positions_file = path.buildpar(*path.SECT_DIR_SEGS,
                                   path.PosListSeg,
                                   top=table.top,
                                   sample=table.sample,
                                   cmd=path.CMD_LIST_DIR,
                                   ref=table.ref,
                                   sect=table.sect,
                                   ext=path.CSV_EXT)
    if need_write(positions_file, force):
        positions = pd.MultiIndex.from_product(
            [[table.ref], find_pos(table, **kwargs)],
            names=[REF_NAME, POS_NAME]
        )
        positions.to_frame(index=False).to_csv(positions_file, index=False)
    return positions_file


@docdef.auto()
def run(input_path: tuple[str, ...], *,
        max_fmut_pos,
        complement: bool,
        force: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ List positions meeting specific criteria from each table. """
    # Find the positional table files.
    pos_table_files = find_pos_tables(input_path)
    # List positions for each table.
    return dispatch(list_pos,
                    max_procs,
                    parallel,
                    pass_n_procs=False,
                    args=as_list_of_tuples(pos_table_files),
                    kwargs=dict(max_fmut_pos=max_fmut_pos,
                                complement=complement,
                                force=force))


params = [
    # Input files
    arg_input_path,
    # Selection
    opt_max_fmut_pos,
    opt_complement,
    # Effort
    opt_force,
    # Parallelization
    opt_max_procs,
    opt_parallel,
]


@command(CMD_LISTPOS, params=params)
def cli(*args, **kwargs):
    """ Add more clusters to a dataset that was already clustered. """
    return run(*args, **kwargs)
