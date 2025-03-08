from pathlib import Path
from typing import Iterable

from click import command

from .core.arg import (CMD_LIST,
                       arg_input_path,
                       opt_branch,
                       opt_min_ninfo_pos,
                       opt_max_fmut_pos,
                       opt_force,
                       opt_num_cpus)
from .core.lists import List, PositionList
from .core.run import run_func
from .core.table import MUTAT_REL, PositionTable, PositionTableLoader
from .core.task import dispatch
from .core.write import need_write
from .mask.lists import MaskPositionList
from .relate.lists import RelatePositionList


def find_pos(table: PositionTable,
             max_fmut_pos: float,
             complement: bool):
    # Initially select all unmasked positions.
    region = table.region.copy()
    positions = region.unmasked_int
    # Apply each filter.
    region.add_mask(
        "max_fmut_pos",
        positions[table.fetch_ratio(rel=MUTAT_REL,
                                    exclude_masked=True,
                                    squeeze=True)
                  > max_fmut_pos],
        complement=complement
    )
    return region.unmasked_int


def write_list(table: PositionTableLoader,
               list_type: type[List], *,
               branch: str,
               min_ninfo_pos: int,
               max_fmut_pos: float,
               force: bool):
    """ Write a List based on a Table. """
    list_file = list_type.get_path_from_table(table, branch)
    if need_write(list_file, force):
        kwargs = dict()
        if issubclass(list_type, PositionList):
            kwargs.update(min_ninfo_pos=min_ninfo_pos,
                          max_fmut_pos=max_fmut_pos)
        else:
            raise ValueError(list_type)
        new_list = list_type.from_table(table, branch=branch, **kwargs)
        new_file = new_list.save(table.top, force)
        assert new_file == list_file
    return list_file


def iter_tables(input_path: Iterable[str | Path], **kwargs):
    """ Iterate through all types of List and all Tables from which each
    type of List can be created. """
    if not isinstance(input_path, (tuple, list, set, dict)):
        # Make sure input_path is not an iterator that will be exhausted
        # on the first run-through.
        input_path = list(input_path)
    for list_type in [RelatePositionList,
                      MaskPositionList]:
        for table in list_type.get_table_type().load_tables(input_path,
                                                            **kwargs):
            yield table, list_type


@run_func(CMD_LIST)
def run(input_path: Iterable[str | Path], *,
        branch: str,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        force: bool,
        num_cpus: int) -> list[Path]:
    """ List positions to mask. """
    return dispatch(write_list,
                    num_cpus=num_cpus,
                    pass_num_cpus=False,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=list(iter_tables(input_path)),
                    kwargs=dict(branch=branch,
                                min_ninfo_pos=min_ninfo_pos,
                                max_fmut_pos=max_fmut_pos,
                                force=force))


params = [
    # Input and output files
    arg_input_path,
    opt_branch,
    # Selection
    opt_min_ninfo_pos,
    opt_max_fmut_pos,
    # Effort
    opt_force,
    # Parallelization
    opt_num_cpus,
]


@command(CMD_LIST, params=params)
def cli(*args, **kwargs):
    """ List positions to mask. """
    return run(*args, **kwargs)
