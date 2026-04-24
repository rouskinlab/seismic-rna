from datetime import datetime
from itertools import chain
from pathlib import Path
from typing import Iterable

from click import command

from .profile import RNAFoldProfile
from .report import FoldReport
from .rnastructure import fold_shapeknots, require_data_path
from .viennarna import rnafold
from ..cluster.data import ClusterPositionTableLoader
from ..core import path
from ..core.arg import (CMD_FOLD,
                        FOLD_BACKEND_FOLD,
                        FOLD_BACKEND_SHAPEKNOTS,
                        FOLD_BACKEND_RNAFOLD,
                        arg_input_path,
                        opt_branch,
                        opt_tmp_pfx,
                        opt_keep_tmp,
                        opt_fold_regions_file,
                        opt_fold_coords,
                        opt_fold_primers,
                        opt_fold_full,
                        opt_fold_backend,
                        opt_fold_temp,
                        opt_fold_energy_method,
                        opt_shape_slope,
                        opt_shape_intercept,
                        opt_fold_quantile,
                        opt_fold_fpaired,
                        opt_fold_mu_eps,
                        opt_fold_constraint,
                        opt_fold_commands,
                        opt_fold_md,
                        opt_fold_mfe,
                        opt_fold_max,
                        opt_fold_percent,
                        opt_pseudoenergy_all,
                        opt_verify_times,
                        opt_num_cpus,
                        opt_force,
                        optional_path,
                        extra_defaults)
from ..core.error import IncompatibleOptionsError
from ..core.extern import (RNASTRUCTURE_FOLD_CMD,
                           RNASTRUCTURE_SHAPEKNOTS_CMD,
                           VIENNA_RNAFOLD_CMD,
                           VIENNA_RNASUBOPT_CMD,
                           require_dependency)
from ..core.io import calc_sha512_path
from ..core.rna import ct_to_db
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
def fold_region(rna: RNAFoldProfile, *,
                out_dir: Path,
                tmp_dir: Path,
                branch: str,
                fold_backend: str,
                shape_slope: float,
                shape_intercept: float,
                fold_constraint: Path | None,
                fold_commands: Path | None,
                fold_md: int,
                fold_mfe: bool,
                fold_max: int,
                fold_percent: float,
                pseudoenergy_all: bool,
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
        rna.write_varna_color_file(out_dir, branch)
        if fold_backend == FOLD_BACKEND_RNAFOLD:
            ct_file = rnafold(rna,
                              out_dir=out_dir,
                              tmp_dir=tmp_dir,
                              branch=branch,
                              shape_slope=shape_slope,
                              shape_intercept=shape_intercept,
                              fold_constraint=fold_constraint,
                              fold_commands=fold_commands,
                              fold_md=fold_md,
                              fold_mfe=fold_mfe,
                              fold_max=fold_max,
                              pseudoenergy_all=pseudoenergy_all,
                              keep_tmp=keep_tmp,
                              num_cpus=num_cpus)
        elif fold_backend in {FOLD_BACKEND_FOLD, FOLD_BACKEND_SHAPEKNOTS}:
            pseudoknots = fold_backend == FOLD_BACKEND_SHAPEKNOTS
            ct_file = fold_shapeknots(rna,
                                      out_dir=out_dir,
                                      tmp_dir=tmp_dir,
                                      branch=branch,
                                      pseudoknots=pseudoknots,
                                      shape_slope=shape_slope,
                                      shape_intercept=shape_intercept,
                                      fold_constraint=fold_constraint,
                                      fold_md=fold_md,
                                      fold_mfe=fold_mfe,
                                      fold_max=fold_max,
                                      fold_percent=fold_percent,
                                      keep_tmp=keep_tmp,
                                      num_cpus=num_cpus)
        else:
            raise ValueError(
                f"Invalid value for --fold-backend: {repr(fold_backend)}"
            )
        ct_to_db(ct_file, force=True)
        constraint_checksum = (calc_sha512_path(fold_constraint)
                               if fold_constraint else "")
        commands_checksum = (calc_sha512_path(fold_commands)
                             if fold_commands else "")
        ended = datetime.now()
        report = FoldReport(branches=branches,
                            sample=rna.sample,
                            ref=rna.ref,
                            reg=rna.reg,
                            profile=rna.profile,
                            fold_backend=fold_backend,
                            fold_energy_method=rna.fold_energy_method,
                            shape_slope=shape_slope,
                            shape_intercept=shape_intercept,
                            fold_quantile=rna.fold_quantile,
                            pseudoenergy_all=pseudoenergy_all,
                            fold_temp=rna.fold_temp,
                            fold_fpaired=rna.fold_fpaired,
                            fold_mu_eps=rna.mu_eps,
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
               fold_energy_method: str,
               fold_quantile: float,
               fold_fpaired: float,
               fold_mu_eps: float,
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
                        RNAFoldProfile.from_profile(
                            profile,
                            fold_temp=fold_temp,
                            fold_energy_method=fold_energy_method,
                            fold_quantile=fold_quantile,
                            fold_fpaired=fold_fpaired,
                            mu_eps=fold_mu_eps,
                        )
                        for profile in table.iter_profiles(regions=regions)
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
        fold_backend: str,
        fold_temp: float,
        fold_energy_method: str,
        shape_slope: float,
        shape_intercept: float,
        fold_quantile: float,
        fold_fpaired: float,
        fold_mu_eps: float,
        fold_constraint: str | None,
        fold_commands: str | None,
        fold_md: int,
        fold_mfe: bool,
        fold_max: int,
        fold_percent: float,
        pseudoenergy_all: bool,
        tmp_pfx: str | Path,
        keep_tmp: bool,
        verify_times: bool,
        num_cpus: int,
        force: bool):
    """ Predict RNA secondary structures using mutation rates. """
    # Check for dependencies.
    if fold_backend == FOLD_BACKEND_RNAFOLD:
        # Use ViennaRNA.
        require_dependency(VIENNA_RNAFOLD_CMD, __name__)
        require_dependency(VIENNA_RNASUBOPT_CMD, __name__)
    else:
        # Use RNAstructure.
        require_data_path()
        if not pseudoenergy_all:
            raise IncompatibleOptionsError(
                "--pseudoenergy-stacked requires "
                f"--fold-backend={FOLD_BACKEND_RNAFOLD}"
            )
        if fold_backend == FOLD_BACKEND_FOLD:
            require_dependency(RNASTRUCTURE_FOLD_CMD, __name__)
        elif fold_backend == FOLD_BACKEND_SHAPEKNOTS:
            require_dependency(RNASTRUCTURE_SHAPEKNOTS_CMD, __name__)
        else:
            raise ValueError(
                f"Invalid value for --fold-backend: {repr(fold_backend)}"
            )
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
                    fold_backend=fold_backend,
                    fold_temp=fold_temp,
                    fold_energy_method=fold_energy_method,
                    shape_slope=shape_slope,
                    shape_intercept=shape_intercept,
                    fold_quantile=fold_quantile,
                    fold_fpaired=fold_fpaired,
                    fold_mu_eps=fold_mu_eps,
                    fold_constraint=optional_path(fold_constraint),
                    fold_commands=optional_path(fold_commands),
                    fold_md=fold_md,
                    fold_mfe=fold_mfe,
                    fold_max=fold_max,
                    fold_percent=fold_percent,
                    pseudoenergy_all=pseudoenergy_all,
                    force=force)
    )))


params = [
    arg_input_path,
    opt_branch,
    opt_fold_regions_file,
    opt_fold_coords,
    opt_fold_primers,
    opt_fold_full,
    opt_fold_backend,
    opt_fold_temp,
    opt_fold_energy_method,
    opt_shape_slope,
    opt_shape_intercept,
    opt_fold_quantile,
    opt_fold_fpaired,
    opt_fold_mu_eps,
    opt_fold_constraint,
    opt_fold_commands,
    opt_fold_md,
    opt_fold_mfe,
    opt_fold_max,
    opt_fold_percent,
    opt_pseudoenergy_all,
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
