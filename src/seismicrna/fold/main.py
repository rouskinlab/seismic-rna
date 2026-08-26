from datetime import datetime
from itertools import chain
from pathlib import Path
from typing import Iterable

from click import command

from .profile import RNAFoldProfile
from .report import FoldReport
from .rnastructure import require_data_path, run_rnastructure
from .viennarna import run_rnafold
from ..cluster.data import ClusterPositionTableLoader
from ..core import path
from ..core.arg.cmd import CMD_FOLD
from ..core.arg.cli import (
    FOLD_BACKEND_AUTO,
    FOLD_BACKEND_RNASTRUCTURE,
    FOLD_BACKEND_VIENNARNA,
    FOLD_ENERGY_METHOD_AUTO,
    FOLD_ENERGY_METHOD_CORDERO,
    FOLD_ENERGY_METHOD_DEIGAN,
    FOLD_ENERGY_METHOD_EDDY,
    PROBE_DMS,
    arg_input_path,
    opt_branch,
    opt_tmp_pfx,
    opt_keep_tmp,
    opt_fold_dry_run,
    opt_fold_regions_file,
    opt_fold_coords,
    opt_fold_primers,
    opt_fold_table_region,
    opt_fold_backend,
    opt_pseudoknots,
    opt_fold_temp,
    opt_fold_energy_method,
    opt_deigan_slope,
    opt_deigan_intercept,
    opt_fold_quantile,
    opt_fold_fpaired,
    opt_fold_constraint,
    opt_fold_commands,
    opt_eddy_prior_paired_file,
    opt_eddy_prior_unpaired_file,
    opt_fold_md,
    opt_fold_mfe,
    opt_fold_max,
    opt_fold_percent,
    opt_fold_edelta,
    opt_fold_isolated,
    opt_verify_times,
    opt_num_cpus,
    opt_force,
    optional_path,
)
from ..core.error import IncompatibleOptionsError
from ..core.shell import (
    RNASTRUCTURE_FOLD_CMD,
    RNASTRUCTURE_SHAPEKNOTS_CMD,
    VIENNA_RNACOFOLD_CMD,
    VIENNA_RNAFOLD_CMD,
    VIENNA_RNASUBOPT_CMD,
)
from ..core.depend import require_dependency
from ..core.io.checksum import calc_sha512_path
from ..core.logs import logger
from ..core.rna.io import ct_to_db
from ..core.run import run_func
from ..core.seq.xna import DNA
from ..core.seq.region import RefRegions, Region
from ..core.seq.refs import RefSeqs
from ..core.task import as_list_of_tuples, dispatch
from ..core.tmp import with_tmp_dir
from ..core.write import need_write
from ..filter.table import FilterPositionTableLoader
from ..duplex.table import DuplexPositionTable, DuplexPositionTableLoader


def resolve_fold_backend(probe: str, fold_backend: str) -> str:
    if fold_backend == FOLD_BACKEND_AUTO:
        return (
            FOLD_BACKEND_RNASTRUCTURE if probe == PROBE_DMS else FOLD_BACKEND_VIENNARNA
        )
    return fold_backend


def resolve_fold_energy_method(probe: str, fold_energy_method: str) -> str:
    if fold_energy_method == FOLD_ENERGY_METHOD_AUTO:
        return (
            FOLD_ENERGY_METHOD_CORDERO
            if probe == PROBE_DMS
            else FOLD_ENERGY_METHOD_EDDY
        )
    return fold_energy_method


def check_fold_deps(fold_backend: str, pseudoknots: bool = False):
    if fold_backend == FOLD_BACKEND_VIENNARNA:
        if pseudoknots:
            raise ValueError(
                f"--pseudoknots is not supported with "
                f"--fold-backend={FOLD_BACKEND_VIENNARNA}; "
                f"use --fold-backend={FOLD_BACKEND_RNASTRUCTURE} instead"
            )
        require_dependency(VIENNA_RNAFOLD_CMD, __name__)
        require_dependency(VIENNA_RNASUBOPT_CMD, __name__)
    elif fold_backend == FOLD_BACKEND_RNASTRUCTURE:
        require_data_path()
        if pseudoknots:
            require_dependency(RNASTRUCTURE_SHAPEKNOTS_CMD, __name__)
        else:
            require_dependency(RNASTRUCTURE_FOLD_CMD, __name__)
    else:
        raise ValueError(f"Invalid value for --fold-backend: {repr(fold_backend)}")


def load_foldable_tables(input_path: Iterable[str | Path], **kwargs):
    """Load tables that can be folded."""
    # Since table_type.load_tables() will be called multiple times, make
    # sure it is not an exhaustible generator.
    paths = list(input_path)
    for table_type in [
        FilterPositionTableLoader,
        ClusterPositionTableLoader,
        DuplexPositionTableLoader,
    ]:
        yield from table_type.load_tables(paths, **kwargs)


@with_tmp_dir(pass_keep_tmp=True)
def fold_region(
    rna: RNAFoldProfile,
    *,
    out_dir: Path,
    tmp_dir: Path,
    branch: str,
    fold_dry_run: bool,
    fold_backend: str,
    pseudoknots: bool,
    fold_constraint: Path | None,
    fold_commands: Path | None,
    eddy_prior_paired_file: Path | None,
    eddy_prior_unpaired_file: Path | None,
    fold_md: int,
    fold_mfe: bool,
    fold_max: int,
    fold_percent: float,
    fold_edelta: float,
    fold_isolated: bool,
    force: bool,
    keep_tmp: bool,
    num_cpus: int,
):
    """Fold a region of an RNA from one mutational profile."""
    branches = path.add_branch(path.FOLD_STEP, branch, rna.branches)
    report_file = FoldReport.build_path(
        {
            path.TOP: out_dir,
            path.SAMPLE: rna.sample,
            path.BRANCHES: branches,
            path.REF: rna.ref,
            path.REG: rna.reg,
            path.PROFILE: rna.profile,
        }
    )
    if need_write(report_file, force):
        began = datetime.now()
        rna.write_varna_color_file(out_dir, branch)
        with logger.debug.single_context("folding {}", rna):
            fasta_tmp = rna.write_fasta(tmp_dir, branch)
            mus_file = rna.write_mus_file(tmp_dir, branch)
            ct_tmp = rna.get_ct_file(tmp_dir, branch)
            ct_out = rna.get_ct_file(out_dir, branch)
            vienna_tmp = None
            db_tmp = None
            try:
                if fold_backend == FOLD_BACKEND_VIENNARNA:
                    if pseudoknots:
                        raise IncompatibleOptionsError(
                            f"fold_backend={FOLD_BACKEND_VIENNARNA} "
                            f"is incompatible with pseudoknots={pseudoknots}"
                        )
                    vienna_tmp = rna.get_vienna_file(tmp_dir, branch)
                    db_tmp = rna.get_db_file(tmp_dir, branch)
                    run_rnafold(
                        fasta_tmp,
                        ct_tmp,
                        ct_out,
                        vienna_tmp,
                        db_tmp,
                        sp_data=mus_file,
                        sp_strategy=rna.rnafold_sp_strategy,
                        eddy_prior_paired_file=eddy_prior_paired_file,
                        eddy_prior_unpaired_file=eddy_prior_unpaired_file,
                        fold_constraint=fold_constraint,
                        fold_commands=fold_commands,
                        fold_temp_c=rna.fold_temp_c,
                        fold_isolated=fold_isolated,
                        fold_md=fold_md,
                        fold_max=fold_max,
                        fold_mfe=fold_mfe,
                        fold_edelta=fold_edelta,
                        end5=rna.region.end5,
                        num_cpus=num_cpus,
                        fold_dry_run=fold_dry_run,
                    )
                else:
                    rnastructure_shape_args = rna.get_rnastructure_shape_args(
                        tmp_dir, branch
                    )
                    run_rnastructure(
                        fasta_tmp,
                        ct_tmp,
                        ct_out,
                        pseudoknots=pseudoknots,
                        fold_temp_k=rna.fold_temp_k,
                        end5=rna.region.end5,
                        num_cpus=num_cpus,
                        fold_dry_run=fold_dry_run,
                        fold_constraint=fold_constraint,
                        fold_isolated=fold_isolated,
                        fold_md=fold_md,
                        fold_mfe=fold_mfe,
                        fold_max=fold_max,
                        fold_percent=fold_percent,
                        **rnastructure_shape_args,
                    )
            finally:
                if not keep_tmp:
                    fasta_tmp.unlink(missing_ok=True)
                    mus_file.unlink(missing_ok=True)
                    if ct_tmp != ct_out:
                        ct_tmp.unlink(missing_ok=True)
                    if db_tmp is not None:
                        db_tmp.unlink(missing_ok=True)
                    if vienna_tmp is not None:
                        vienna_tmp.unlink(missing_ok=True)
        ct_file = ct_out
        if not fold_dry_run:
            ct_to_db(ct_file, force=True)
        constraint_checksum = (
            calc_sha512_path(fold_constraint) if fold_constraint else ""
        )
        commands_checksum = calc_sha512_path(fold_commands) if fold_commands else ""
        eddy_prior_paired_file_checksum = (
            calc_sha512_path(eddy_prior_paired_file) if eddy_prior_paired_file else ""
        )
        eddy_prior_unpaired_file_checksum = (
            calc_sha512_path(eddy_prior_unpaired_file)
            if eddy_prior_unpaired_file
            else ""
        )
        ended = datetime.now()
        report = FoldReport(
            branches=branches,
            sample=rna.sample,
            ref=rna.ref,
            reg=rna.reg,
            profile=rna.profile,
            fold_dry_run=fold_dry_run,
            fold_backend=fold_backend,
            pseudoknots=pseudoknots,
            fold_energy_method=rna.fold_energy_method,
            deigan_slope=rna.deigan_slope,
            deigan_intercept=rna.deigan_intercept,
            fold_quantile=rna.fold_quantile,
            fold_isolated=fold_isolated,
            fold_temp=rna.fold_temp_c,
            fold_md=fold_md,
            fold_mfe=fold_mfe,
            fold_max=fold_max,
            fold_percent=fold_percent,
            fold_edelta=fold_edelta,
            constraint_checksum=constraint_checksum,
            commands_checksum=commands_checksum,
            eddy_prior_paired_file_checksum=eddy_prior_paired_file_checksum,
            eddy_prior_unpaired_file_checksum=eddy_prior_unpaired_file_checksum,
            began=began,
            ended=ended,
        )
        report.save(out_dir, force=True)
    return report_file


def fold_table(
    table: FilterPositionTableLoader | ClusterPositionTableLoader,
    regions: list[Region],
    fold_table_region: bool,
    fold_temp: float,
    fold_energy_method: str,
    fold_quantile: float,
    deigan_slope: float,
    deigan_intercept: float,
    eddy_prior_paired_file: Path | None,
    eddy_prior_unpaired_file: Path | None,
    num_cpus: int,
    keep_tmp: bool,
    fold_dry_run: bool,
    fold_backend: str,
    pseudoknots: bool,
    **kwargs,
):
    """Fold an RNA molecule from one table of reactivities."""
    probe = table._dataset.probe
    fold_backend = resolve_fold_backend(probe, fold_backend)
    fold_energy_method = resolve_fold_energy_method(probe, fold_energy_method)
    check_fold_deps(fold_backend, pseudoknots)
    return dispatch(
        fold_region,
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
                deigan_slope=deigan_slope,
                deigan_intercept=deigan_intercept,
            )
            for _, profile in table.iter_profiles(
                fold_table_region=fold_table_region, regions=regions
            )
        ),
        kwargs=dict(
            out_dir=table.top,
            keep_tmp=(keep_tmp or fold_dry_run),
            fold_dry_run=fold_dry_run,
            fold_backend=fold_backend,
            pseudoknots=pseudoknots,
            eddy_prior_paired_file=eddy_prior_paired_file,
            eddy_prior_unpaired_file=eddy_prior_unpaired_file,
            **kwargs,
        ),
    )


# Energy methods RNAcofold can apply to a duplex. Deigan reactivities go
# straight to --shape/--shapeMethod; Cordero (DMS) is applied as
# pseudo-mutation rates back-transformed to the Deigan form. RNAcofold
# has no modern --sp-strategy interface, so anything else (e.g. Eddy) is
# unsupported.
DUPLEX_ENERGY_METHODS = (FOLD_ENERGY_METHOD_DEIGAN, FOLD_ENERGY_METHOD_CORDERO)


def resolve_duplex_energy_method(probe: str, fold_energy_method: str):
    """Resolve the energy method for a duplex (folded by RNAcofold): auto
    -> Cordero for DMS, else Deigan. Any method RNAcofold cannot express
    as a SHAPE constraint raises, listing the compatible methods."""
    if fold_energy_method == FOLD_ENERGY_METHOD_AUTO:
        method = (
            FOLD_ENERGY_METHOD_CORDERO
            if probe == PROBE_DMS
            else FOLD_ENERGY_METHOD_DEIGAN
        )
    else:
        method = fold_energy_method
    if method not in DUPLEX_ENERGY_METHODS:
        raise IncompatibleOptionsError(
            f"RNAcofold cannot fold a duplex with the {method} energy method. "
            f"Cofolding supports only {' and '.join(DUPLEX_ENERGY_METHODS)}: "
            f"{FOLD_ENERGY_METHOD_DEIGAN} applies SHAPE reactivities directly "
            f"via --shapeMethod, and {FOLD_ENERGY_METHOD_CORDERO} (DMS) applies "
            "pseudo-energies; other strategies have no RNAcofold equivalent"
        )
    return method


@with_tmp_dir(pass_keep_tmp=True)
def fold_duplex_region(
    rna,
    out_dir: Path,
    *,
    tmp_dir: Path,
    branch: str,
    fold_dry_run: bool,
    fold_constraint: Path | None,
    fold_md: int,
    fold_isolated: bool,
    force: bool,
    keep_tmp: bool,
    num_cpus: int,
):
    """Fold a duplex profile with RNAcofold, writing a standard
    fold report so that graphing and drawing treat it like any fold."""
    from ..duplex.viennarna import run_rnacofold

    branches = path.add_branch(path.FOLD_STEP, branch, rna.branches)
    report_file = FoldReport.build_path(
        {
            path.TOP: out_dir,
            path.SAMPLE: rna.sample,
            path.BRANCHES: branches,
            path.REF: rna.ref,
            path.REG: rna.reg,
            path.PROFILE: rna.profile,
        }
    )
    if need_write(report_file, force):
        began = datetime.now()
        rna.write_varna_color_file(out_dir, branch)
        with logger.debug.single_context("cofolding {}", rna):
            fasta_tmp = rna.write_fasta(tmp_dir, branch)
            mus_file = rna.write_mus_file(tmp_dir, branch)
            ct_tmp = rna.get_ct_file(tmp_dir, branch)
            ct_out = rna.get_ct_file(out_dir, branch)
            vienna_tmp = rna.get_vienna_file(tmp_dir, branch)
            db_tmp = rna.get_db_file(tmp_dir, branch)
            try:
                run_rnacofold(
                    fasta_tmp,
                    ct_tmp,
                    ct_out,
                    vienna_tmp,
                    db_tmp,
                    shape_file=mus_file,
                    shape_method=rna.cofold_shape_method,
                    fold_constraint=fold_constraint,
                    fold_temp_c=rna.fold_temp_c,
                    fold_isolated=fold_isolated,
                    fold_md=fold_md,
                    end5=rna.region.end5,
                    fold_dry_run=fold_dry_run,
                )
            finally:
                if not keep_tmp:
                    fasta_tmp.unlink(missing_ok=True)
                    mus_file.unlink(missing_ok=True)
                    if ct_tmp != ct_out:
                        ct_tmp.unlink(missing_ok=True)
                    db_tmp.unlink(missing_ok=True)
                    vienna_tmp.unlink(missing_ok=True)
        if not fold_dry_run:
            ct_to_db(ct_out, force=True)
        constraint_checksum = (
            calc_sha512_path(fold_constraint) if fold_constraint else ""
        )
        report = FoldReport(
            branches=branches,
            sample=rna.sample,
            ref=rna.ref,
            reg=rna.reg,
            profile=rna.profile,
            fold_dry_run=fold_dry_run,
            fold_backend=FOLD_BACKEND_VIENNARNA,
            pseudoknots=False,
            fold_energy_method=rna.fold_energy_method,
            deigan_slope=rna.deigan_slope,
            deigan_intercept=rna.deigan_intercept,
            fold_quantile=rna.fold_quantile,
            fold_isolated=fold_isolated,
            fold_temp=rna.fold_temp_c,
            fold_md=fold_md,
            fold_mfe=True,
            fold_max=1,
            fold_percent=0.0,
            fold_edelta=0.0,
            constraint_checksum=constraint_checksum,
            commands_checksum="",
            eddy_prior_paired_file_checksum="",
            eddy_prior_unpaired_file_checksum="",
            began=began,
            ended=datetime.now(),
        )
        report.save(out_dir, force=True)
    return report_file


def fold_duplex_table(
    table,
    *,
    fold_table_region: bool,
    fold_temp: float,
    fold_energy_method: str,
    fold_quantile: float,
    fold_fpaired: float,
    deigan_slope: float,
    deigan_intercept: float,
    num_cpus: int,
    **kwargs,
):
    """Fold every duplex profile of a duplex table with RNAcofold."""
    from ..duplex.profile import RNACofoldProfile

    require_dependency(VIENNA_RNACOFOLD_CMD, __name__)
    # The probe is only needed to resolve the "auto" energy method; read
    # it lazily so an explicit method does not require the source dataset.
    probe = (
        table._dataset.probe if fold_energy_method == FOLD_ENERGY_METHOD_AUTO else None
    )
    fold_energy_method = resolve_duplex_energy_method(probe, fold_energy_method)
    duplexes = [
        (
            RNACofoldProfile.from_duplex(
                profile,
                table.cut,
                fold_temp=fold_temp,
                fold_energy_method=fold_energy_method,
                fold_quantile=fold_quantile,
                fold_fpaired=fold_fpaired,
                deigan_slope=deigan_slope,
                deigan_intercept=deigan_intercept,
            ),
            table.top,
        )
        for _, profile in table.iter_profiles(fold_table_region=fold_table_region)
    ]
    return dispatch(
        fold_duplex_region,
        num_cpus=num_cpus,
        pass_num_cpus=True,
        as_list=True,
        ordered=False,
        raise_on_error=False,
        args=duplexes,
        kwargs=kwargs,
    )


@run_func(CMD_FOLD)
def run(
    input_path: Iterable[str | Path],
    *,
    branch: str,
    fold_coords: Iterable[tuple[str, int, int]],
    fold_primers: Iterable[tuple[str, DNA, DNA]],
    fold_regions_file: str | None,
    fold_table_region: bool,
    fold_dry_run: bool,
    fold_backend: str,
    pseudoknots: bool,
    fold_energy_method: str,
    deigan_slope: float,
    deigan_intercept: float,
    fold_temp: float,
    fold_quantile: float,
    fold_fpaired: float,
    fold_constraint: str | None,
    fold_commands: str | None,
    eddy_prior_paired_file: str | None,
    eddy_prior_unpaired_file: str | None,
    fold_md: int,
    fold_mfe: bool,
    fold_max: int,
    fold_percent: float,
    fold_edelta: float,
    fold_isolated: bool,
    tmp_pfx: str | Path,
    keep_tmp: bool,
    verify_times: bool,
    num_cpus: int,
    force: bool,
):
    """Predict RNA secondary structures using mutation rates."""
    # Check for dependencies; deferred per table when backend is auto.
    if fold_backend != FOLD_BACKEND_AUTO:
        check_fold_deps(fold_backend, pseudoknots)
    # List the tables, separating duplex tables from single-molecule
    # tables since duplexes are folded with RNAcofold.
    tables = list(load_foldable_tables(input_path, verify_times=verify_times))
    duplex_tables = [
        table for table in tables if isinstance(table, DuplexPositionTable)
    ]
    tables = [table for table in tables if not isinstance(table, DuplexPositionTable)]
    # Fold the duplex tables with RNAcofold. If there are none, skip the
    # dispatch (which would otherwise warn that it was given no tasks).
    duplex_results = list()
    if duplex_tables:
        duplex_results.extend(
            chain(
                *dispatch(
                    fold_duplex_table,
                    num_cpus=num_cpus,
                    pass_num_cpus=True,
                    as_list=False,
                    ordered=False,
                    raise_on_error=False,
                    label="cofolding duplexes",
                    unit="table",
                    args=as_list_of_tuples(duplex_tables),
                    kwargs=dict(
                        fold_table_region=fold_table_region,
                        fold_temp=fold_temp,
                        fold_energy_method=fold_energy_method,
                        fold_quantile=fold_quantile,
                        fold_fpaired=fold_fpaired,
                        deigan_slope=deigan_slope,
                        deigan_intercept=deigan_intercept,
                        branch=branch,
                        tmp_pfx=tmp_pfx,
                        keep_tmp=(keep_tmp or fold_dry_run),
                        fold_dry_run=fold_dry_run,
                        fold_constraint=optional_path(fold_constraint),
                        fold_md=fold_md,
                        fold_isolated=fold_isolated,
                        force=force,
                    ),
                )
            )
        )
    # If the only inputs were duplex tables, there are no single-molecule
    # profiles to fold, so skip the dispatch (which would otherwise warn
    # that it was given no tasks).
    if not tables:
        return duplex_results
    # Get the regions to fold for every single-molecule reference.
    ref_seqs = RefSeqs()
    for table in tables:
        ref_seqs.add(table.ref, table.refseq)
    fold_regions = RefRegions(
        ref_seqs,
        regs_file=optional_path(fold_regions_file),
        ends=fold_coords,
        primers=fold_primers,
        default_full=(not fold_table_region),
    ).dict
    # For each table whose reference had no regions defined, default to
    # the table's region.
    args = [
        (table, (region if (region := fold_regions[table.ref]) else [table.region]))
        for table in tables
    ]
    # Fold the RNA profiles.
    return duplex_results + list(
        chain(
            *dispatch(
                fold_table,
                num_cpus=num_cpus,
                pass_num_cpus=True,
                as_list=False,
                ordered=False,
                raise_on_error=False,
                label="folding tables",
                unit="table",
                args=args,
                kwargs=dict(
                    fold_table_region=fold_table_region,
                    branch=branch,
                    tmp_pfx=tmp_pfx,
                    keep_tmp=keep_tmp,
                    fold_dry_run=fold_dry_run,
                    fold_backend=fold_backend,
                    pseudoknots=pseudoknots,
                    fold_temp=fold_temp,
                    fold_energy_method=fold_energy_method,
                    deigan_slope=deigan_slope,
                    deigan_intercept=deigan_intercept,
                    fold_quantile=fold_quantile,
                    fold_constraint=optional_path(fold_constraint),
                    fold_commands=optional_path(fold_commands),
                    eddy_prior_paired_file=optional_path(eddy_prior_paired_file),
                    eddy_prior_unpaired_file=optional_path(eddy_prior_unpaired_file),
                    fold_md=fold_md,
                    fold_mfe=fold_mfe,
                    fold_max=fold_max,
                    fold_percent=fold_percent,
                    fold_edelta=fold_edelta,
                    fold_isolated=fold_isolated,
                    force=force,
                ),
            )
        )
    )


params = [
    arg_input_path,
    opt_branch,
    opt_fold_regions_file,
    opt_fold_coords,
    opt_fold_primers,
    opt_fold_table_region,
    opt_fold_dry_run,
    opt_fold_backend,
    opt_pseudoknots,
    opt_fold_energy_method,
    opt_fold_temp,
    opt_deigan_slope,
    opt_deigan_intercept,
    opt_fold_quantile,
    opt_fold_fpaired,
    opt_fold_constraint,
    opt_fold_commands,
    opt_eddy_prior_paired_file,
    opt_eddy_prior_unpaired_file,
    opt_fold_md,
    opt_fold_mfe,
    opt_fold_max,
    opt_fold_percent,
    opt_fold_edelta,
    opt_fold_isolated,
    opt_tmp_pfx,
    opt_keep_tmp,
    opt_verify_times,
    opt_num_cpus,
    opt_force,
]


@command(CMD_FOLD, params=params)
def cli(*args, **kwargs):
    """Predict RNA secondary structures using mutation rates."""
    return run(*args, **kwargs)
