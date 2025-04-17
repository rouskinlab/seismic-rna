from datetime import datetime
from itertools import chain
from pathlib import Path
from typing import Iterable

from click import command

from .report import FoldReport
from .rnastructure import fold, require_data_path
from .viennarna import rnafold
from ..cluster.data import ClusterPositionTableLoader
from ..core import path
from ..core.arg import (CMD_FOLD,
                        arg_input_path,
                        opt_branch,
                        opt_tmp_pfx,
                        opt_keep_tmp,
                        opt_fold_regions_file,
                        opt_fold_coords,
                        opt_fold_primers,
                        opt_fold_full,
                        opt_fold_vienna,
                        opt_fold_temp,
                        opt_fold_fpaired,
                        opt_fold_constraint,
                        opt_fold_commands,
                        opt_fold_md,
                        opt_fold_mfe,
                        opt_fold_max,
                        opt_fold_percent,
                        opt_verify_times,
                        opt_num_cpus,
                        opt_force,
                        optional_path,
                        extra_defaults)
from ..core.extern import (RNASTRUCTURE_FOLD_CMD,
                           require_dependency)
from ..core.io import calc_sha512_path
from ..core.rna import RNAProfile, ct_to_db
from ..core.run import run_func
from ..core.seq import DNA, RefRegions, RefSeqs, Region
from ..core.task import as_list_of_tuples, dispatch
from ..core.tmp import with_tmp_dir
from ..core.write import need_write
from ..mask.table import MaskPositionTableLoader


def load_foldable_tables(input_path: Iterable[str | Path], **kwargs):
    """ Load tables that can be folded. """
    # Since table_type.load_tables() will be called multiple times, make
    # sure it is not an exhaustible generator.
    paths = list(input_path)
    for table_type in [MaskPositionTableLoader, ClusterPositionTableLoader]:
        yield from table_type.load_tables(paths, **kwargs)


@with_tmp_dir(pass_keep_tmp=True)
def fold_region(rna: RNAProfile, *,
                out_dir: Path,
                tmp_dir: Path,
                branch: str,
                fold_vienna: bool,
                fold_constraint: Path | None,
                fold_commands: Path | None,
                fold_md: int,
                fold_mfe: bool,
                fold_max: int,
                fold_percent: float,
                force: bool,
                keep_tmp: bool,
                num_cpus: int):
    """ Fold a region of an RNA from one mutational profile. """
    branches = path.add_branch(path.FOLD_STEP, branch, rna.branches)
    report_file = FoldReport.build_path({path.TOP: out_dir,
                                         path.SAMPLE: rna.sample,
                                         path.BRANCHES: branches,
                                         path.REF: rna.ref,
                                         path.REG: rna.reg,
                                         path.PROFILE: rna.profile})
    if need_write(report_file, force):
        began = datetime.now()
        rna.to_varna_color_file(out_dir, branch)
        if fold_vienna:
            ct_file = rnafold(rna,
                              out_dir=out_dir,
                              tmp_dir=tmp_dir,
                              branch=branch,
                              fold_constraint=fold_constraint,
                              fold_commands=fold_commands,
                              fold_md=fold_md,
                              fold_mfe=fold_mfe,
                              fold_max=fold_max,
                              fold_percent=fold_percent,
                              keep_tmp=keep_tmp,
                              num_cpus=num_cpus)
        else:
            ct_file = fold(rna,
                           out_dir=out_dir,
                           tmp_dir=tmp_dir,
                           branch=branch,
                           fold_constraint=fold_constraint,
                           fold_md=fold_md,
                           fold_mfe=fold_mfe,
                           fold_max=fold_max,
                           fold_percent=fold_percent,
                           keep_tmp=keep_tmp,
                           num_cpus=num_cpus)
        ct_to_db(ct_file, force=True)
        ended = datetime.now()
        constraint_checksum = calc_sha512_path(fold_constraint) if fold_constraint else ''
        commands_checksum = calc_sha512_path(fold_commands) if fold_commands else ''
        report = FoldReport(branches=branches,
                            sample=rna.sample,
                            ref=rna.ref,
                            reg=rna.reg,
                            profile=rna.profile,
                            fold_vienna=fold_vienna,
                            fold_temp=rna.fold_temp,
                            fold_fpaired=rna.fold_fpaired,
                            fold_md=fold_md,
                            fold_mfe=fold_mfe,
                            fold_max=fold_max,
                            fold_percent=fold_percent,
                            constraint_checksum=constraint_checksum,
                            commands_checksum=commands_checksum,
                            began=began,
                            ended=ended)
        report.save(out_dir, force=True)
    return report_file


def fold_table(table: MaskPositionTableLoader | ClusterPositionTableLoader,
               regions: list[Region],
               fold_temp: float,
               fold_fpaired: float,
               num_cpus: int,
               **kwargs):
    """ Fold an RNA molecule from one table of reactivities. """
    return dispatch(fold_region,
                    num_cpus=num_cpus,
                    pass_num_cpus=True,
                    as_list=True,
                    ordered=False,
                    raise_on_error=False,
                    args=as_list_of_tuples(
                        table.iter_profiles(regions=regions,
                                            fold_temp=fold_temp,
                                            fold_fpaired=fold_fpaired)
                    ),
                    kwargs=dict(out_dir=table.top,
                                **kwargs))


@run_func(CMD_FOLD, extra_defaults=extra_defaults)
def run(input_path: Iterable[str | Path], *,
        branch: str,
        fold_coords: Iterable[tuple[str, int, int]],
        fold_primers: Iterable[tuple[str, DNA, DNA]],
        fold_regions_file: str | None,
        fold_full: bool,
        fold_vienna: bool,
        fold_temp: float,
        fold_fpaired: float,
        fold_constraint: str | None,
        fold_commands: str | None,
        fold_md: int,
        fold_mfe: bool,
        fold_max: int,
        fold_percent: float,
        tmp_pfx: str | Path,
        keep_tmp: bool,
        verify_times: bool,
        num_cpus: int,
        force: bool):
    """ Predict RNA secondary structures using mutation rates. """
    # Check for the dependencies and the DATAPATH environment variable.
    require_dependency(RNASTRUCTURE_FOLD_CMD, __name__)
    require_data_path()
    # List the tables.
    tables = list(load_foldable_tables(input_path, verify_times=verify_times))
    # Get the regions to fold for every reference sequence.
    ref_seqs = RefSeqs()
    for table in tables:
        ref_seqs.add(table.ref, table.refseq)
    fold_regions = RefRegions(ref_seqs,
                              regs_file=optional_path(fold_regions_file),
                              ends=fold_coords,
                              primers=fold_primers,
                              default_full=fold_full).dict
    # For each table whose reference had no regions defined, default to
    # the table's region.
    args = [(table, (fold_regions[table.ref]
                     if fold_regions[table.ref]
                     else [table.region]))
            for table in tables]
    # Fold the RNA profiles.
    return list(chain(*dispatch(
        fold_table,
        num_cpus=num_cpus,
        pass_num_cpus=True,
        as_list=False,
        ordered=False,
        raise_on_error=False,
        args=args,
        kwargs=dict(branch=branch,
                    tmp_pfx=tmp_pfx,
                    keep_tmp=keep_tmp,
                    fold_vienna=fold_vienna,
                    fold_temp=fold_temp,
                    fold_fpaired=fold_fpaired,
                    fold_constraint=optional_path(fold_constraint),
                    fold_commands=optional_path(fold_commands),
                    fold_md=fold_md,
                    fold_mfe=fold_mfe,
                    fold_max=fold_max,
                    fold_percent=fold_percent,
                    force=force)
    )))


params = [
    arg_input_path,
    opt_branch,
    opt_fold_regions_file,
    opt_fold_coords,
    opt_fold_primers,
    opt_fold_full,
    opt_fold_vienna,
    opt_fold_temp,
    opt_fold_fpaired,
    opt_fold_constraint,
    opt_fold_commands,
    opt_fold_md,
    opt_fold_mfe,
    opt_fold_max,
    opt_fold_percent,
    opt_tmp_pfx,
    opt_keep_tmp,
    opt_verify_times,
    opt_num_cpus,
    opt_force,
]


@command(CMD_FOLD, params=params)
def cli(*args, **kwargs):
    """ Predict RNA secondary structures using mutation rates. """
    return run(*args, **kwargs)
