import gzip
from pathlib import Path

from click import command

from .core.arg import (CMD_MIGRATE,
                       arg_input_path,
                       opt_max_procs)
from .core.run import run_func
from .core.task import as_list_of_tuples, dispatch

JSON_INDENT = " " * 4
JSON_DELIM = f",\n{JSON_INDENT}"


def find_and_replace(file: Path, find: str, replace: str, count: int = 1):
    if file.suffix.endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open
    with open_func(file, "rt") as f:
        text = f.read()
    find_count = text.count(find)
    if 0 <= count != find_count:
        raise ValueError(f"{file} contains {repr(find)} {find_count} time(s)")
    text = text.replace(find, replace)
    with open_func(file, "wt") as f:
        f.write(text)
    return find_count


def migrate_relate_ref_dir(ref_dir: Path):
    relate_report_file = ref_dir.joinpath("relate-report.json")
    find_and_replace(relate_report_file,
                     '"Mark all ambiguous insertions and deletions":',
                     JSON_DELIM.join(
                         ['"Mark each insertion on the base to its 3\' (True) or 5\' (False) side": true',
                          f'"Mark all ambiguous insertions and deletions (indels)":']))
    find_and_replace(relate_report_file,
                     "qnames",
                     "readnames")
    return ref_dir


def migrate_mask_reg_dir(reg_dir: Path):
    mask_report_file = reg_dir.joinpath("mask-report.json")
    find_and_replace(mask_report_file,
                     "Section",
                     "Region",
                     count=3)
    find_and_replace(mask_report_file,
                     "section",
                     "region",
                     count=3)
    find_and_replace(mask_report_file,
                     "unambiguous",
                     "informative",
                     count=5)
    find_and_replace(mask_report_file,
                     '"Number of reads with too few bases covering the region":',
                     JSON_DELIM.join(['"Number of reads masked from a list": 0',
                                      '"Number of reads with too few bases covering the region":']))
    find_and_replace(mask_report_file,
                     '"Correct observer bias using a quick (typically linear time) heuristic":',
                     JSON_DELIM.join(['"Maximum number of iterations for masking (0 for no limit)": 1',
                                      '"Number of iterations until convergence (0 if not converged)": 1',
                                      '"Correct observer bias using a quick (typically linear time) heuristic":']))
    return reg_dir


def migrate_mask_ref_dir(ref_dir: Path, n_procs: int):
    dispatch(migrate_mask_reg_dir,
             max_procs=n_procs,
             pass_n_procs=False,
             args=as_list_of_tuples(ref_dir.iterdir()))
    return ref_dir


def migrate_cluster_reg_dir(reg_dir: Path):
    cluster_report_file = reg_dir.joinpath("cluster-report.json")
    find_and_replace(cluster_report_file,
                     "Section",
                     "Region")
    return reg_dir


def migrate_cluster_ref_dir(ref_dir: Path, n_procs: int):
    dispatch(migrate_cluster_reg_dir,
             max_procs=n_procs,
             pass_n_procs=False,
             args=as_list_of_tuples(ref_dir.iterdir()))
    return ref_dir


def migrate_table_reg_dir(reg_dir: Path):
    reg = reg_dir.name
    ref = reg_dir.parent.name
    sample_dir = reg_dir.parent.parent.parent
    for table_name in ["relate-per-pos.csv",
                       "relate-per-read.csv.gz",
                       "mask-per-pos.csv",
                       "mask-per-read.csv.gz",
                       "clust-per-pos.csv",
                       "clust-freq.csv"]:
        table_file = reg_dir.joinpath(table_name)
        if table_file.is_file():
            if "per" in table_name:
                find_and_replace(table_file,
                                 "Unambiguous",
                                 "Informative",
                                 count=-1)
            if table_name.startswith("relate"):
                cmd = "relate"
                to_dir = sample_dir.joinpath(cmd).joinpath(ref)
            elif table_name.startswith("mask"):
                cmd = "mask"
                to_dir = sample_dir.joinpath(cmd).joinpath(ref).joinpath(reg)
            elif table_name.startswith("clust"):
                cmd = "cluster"
                to_dir = sample_dir.joinpath(cmd).joinpath(ref).joinpath(reg)
            else:
                raise ValueError(table_name)
            if "per-pos" in table_name:
                table_type = "position"
            elif "per-read" in table_name:
                table_type = "read"
            elif "freq" in table_name:
                table_type = "abundance"
            else:
                raise ValueError(table_name)
            to_name = f"{cmd}-{table_type}-table{''.join(table_file.suffixes)}"
            to_file = to_dir.joinpath(to_name)
            table_file.rename(to_file)
    reg_dir.rmdir()
    return reg_dir


def migrate_table_ref_dir(ref_dir: Path, n_procs: int):
    dispatch(migrate_table_reg_dir,
             max_procs=n_procs,
             pass_n_procs=False,
             args=as_list_of_tuples(ref_dir.iterdir()))
    ref_dir.rmdir()
    return ref_dir


def migrate_fold_reg_dir(reg_dir: Path):
    for file in reg_dir.iterdir():
        if file.name.endswith("fold-report.json"):
            find_and_replace(file,
                             "Section",
                             "Region")
    return reg_dir


def migrate_fold_ref_dir(ref_dir: Path, n_procs: int):
    dispatch(migrate_fold_reg_dir,
             max_procs=n_procs,
             pass_n_procs=False,
             args=as_list_of_tuples(ref_dir.iterdir()))
    return ref_dir


def migrate_graph_file(graph_file: Path):
    if graph_file.suffix in [".html", ".svg"]:
        find_and_replace(graph_file, "Unambiguous", "Informative", count=-1)
        find_and_replace(graph_file, "masked", "filtered", count=-1)
        find_and_replace(graph_file, "' section '", "' region '")
    graph_name = graph_file.name
    if graph_name.split("_")[1] == "masked":
        graph_name = graph_name.replace("masked", "filtered", 1)
        graph_file = graph_file.rename(graph_file.parent.joinpath(graph_name))
    return graph_file


def migrate_graph_reg_dir(reg_dir: Path, n_procs: int):
    dispatch(migrate_graph_file,
             max_procs=n_procs,
             pass_n_procs=False,
             args=as_list_of_tuples(reg_dir.iterdir()))
    return reg_dir


def migrate_graph_ref_dir(ref_dir: Path, n_procs: int):
    dispatch(migrate_graph_reg_dir,
             max_procs=n_procs,
             pass_n_procs=True,
             args=as_list_of_tuples(ref_dir.iterdir()))
    return ref_dir


def migrate_sample_dir(sample_dir: Path, n_procs: int):
    relate_dir = sample_dir.joinpath("relate")
    if relate_dir.is_dir():
        dispatch(migrate_relate_ref_dir,
                 max_procs=n_procs,
                 pass_n_procs=False,
                 args=as_list_of_tuples(relate_dir.iterdir()))
    mask_dir = sample_dir.joinpath("mask")
    if mask_dir.is_dir():
        dispatch(migrate_mask_ref_dir,
                 max_procs=n_procs,
                 pass_n_procs=True,
                 args=as_list_of_tuples(mask_dir.iterdir()))
    cluster_dir = sample_dir.joinpath("cluster")
    if cluster_dir.is_dir():
        dispatch(migrate_cluster_ref_dir,
                 max_procs=n_procs,
                 pass_n_procs=True,
                 args=as_list_of_tuples(cluster_dir.iterdir()))
    table_dir = sample_dir.joinpath("table")
    if table_dir.is_dir():
        dispatch(migrate_table_ref_dir,
                 max_procs=n_procs,
                 pass_n_procs=True,
                 args=as_list_of_tuples(table_dir.iterdir()),
                 raise_on_error=True)
        try:
            table_dir.rmdir()
        except OSError:
            pass
    fold_dir = sample_dir.joinpath("fold")
    if fold_dir.is_dir():
        dispatch(migrate_fold_ref_dir,
                 max_procs=n_procs,
                 pass_n_procs=True,
                 args=as_list_of_tuples(fold_dir.iterdir()),
                 raise_on_error=True)
    graph_dir = sample_dir.joinpath("graph")
    if graph_dir.is_dir():
        dispatch(migrate_graph_ref_dir,
                 max_procs=n_procs,
                 pass_n_procs=True,
                 args=as_list_of_tuples(graph_dir.iterdir()),
                 raise_on_error=True)
    return sample_dir


def migrate_out_dir(out_dir: Path, n_procs: int):
    dispatch(migrate_sample_dir,
             max_procs=n_procs,
             pass_n_procs=True,
             args=as_list_of_tuples(out_dir.iterdir()))
    return out_dir


@run_func(CMD_MIGRATE)
def run(input_path: tuple[str, ...], *,
        max_procs: int) -> list[Path]:
    """ Migrate output directories from v0.21 to v0.22 """
    return dispatch(migrate_out_dir,
                    max_procs=max_procs,
                    pass_n_procs=True,
                    args=as_list_of_tuples(map(Path, input_path)))


params = [
    arg_input_path,
    opt_max_procs
]


@command(CMD_MIGRATE, params=params)
def cli(*args, **kwargs):
    """ Migrate output directories from v0.21 to v0.22 """
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
