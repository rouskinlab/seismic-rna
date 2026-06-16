import os
from pathlib import Path
from typing import Iterable

from click import command

from .ref import get_fasta_path
from ..core import path
from ..core.arg import (
    FOLD_BACKEND_VIENNARNA,
    arg_fasta,
    opt_fold_backend,
    opt_pseudoknots,
    opt_sim_dir,
    opt_tmp_pfx,
    opt_profile_name,
    opt_probe,
    opt_fold_regions_file,
    opt_fold_coords,
    opt_fold_primers,
    opt_fold_constraint,
    opt_fold_temp,
    opt_fold_md,
    opt_fold_mfe,
    opt_fold_max,
    opt_fold_min,
    opt_fold_percent,
    opt_fold_edelta,
    opt_keep_tmp,
    opt_force,
    opt_num_cpus,
    optional_path,
)
from ..core.rna import from_ct
from ..core.logs import logger
from ..core.run import run_func
from ..core.seq import DNA, RefRegions, Region, parse_fasta, write_fasta
from ..core.task import as_list_of_tuples, dispatch
from ..core.write import need_write
from ..fold.main import check_fold_deps, resolve_fold_backend
from ..fold.profile import celsius_to_kelvin, guess_temperature_to_celsius
from ..fold.rnastructure import run_rnastructure
from ..fold.viennarna import run_rnafold

COMMAND = __name__.split(os.path.extsep)[-1]


def get_ct_path(top: Path, region: Region, profile: str):
    """Get the path of a connectivity table (CT) file."""
    return path.buildpar(
        path.CT_FILE_LAST_SEGS,
        {
            path.TOP: top.joinpath(path.SIM_PARAM_DIR),
            path.REF: region.ref,
            path.REG: region.name,
            path.PROFILE: profile,
            path.EXT: path.CT_EXT,
        },
    )


def fold_region(
    region: Region,
    *,
    sim_dir: Path,
    tmp_dir: Path,
    profile_name: str,
    fold_backend: str,
    pseudoknots: bool,
    fold_constraint: Path | None,
    fold_temp: float,
    fold_md: int,
    fold_mfe: bool,
    fold_max: int,
    fold_min: int,
    fold_percent: float,
    fold_edelta: float,
    keep_tmp: bool,
    force: bool,
    num_cpus: int,
):
    """
    Predict RNA secondary structures for one region using RNAstructure.

    Parameters
    ----------
    region: Region
        Sequence region to fold.
    sim_dir: Path
        Simulation output directory; CT file is written under its
        parameter subdirectory.
    tmp_dir: Path
        Directory for temporary FASTA and intermediate CT files.
    profile_name: str
        Profile label embedded in the output CT file path.
    fold_constraint: Path | None
        Path to a folding constraint file; None for no constraints.
    fold_temp: float
        Folding temperature; interpreted as Celsius if in the typical
        physiological range, otherwise as Kelvin.
    fold_md: int
        Maximum distance between paired bases (0 for no limit).
    fold_mfe: bool
        Whether to predict only the minimum free energy structure.
    fold_max: int
        Maximum number of structures to predict.
    fold_percent: float
        Maximum percent energy difference from the MFE structure.
    fold_edelta: float
        Maximum absolute energy difference (kcal/mol) from the MFE for
        suboptimal structures (RNAsubopt only).
    keep_tmp: bool
        Whether to retain temporary files after folding.
    force: bool
        Whether to overwrite an existing output CT file.
    num_cpus: int
        Number of CPU cores to use.

    Returns
    -------
    Path
        Path of the written CT file.
    """
    ct_sim = get_ct_path(sim_dir, region, profile_name)
    if need_write(ct_sim, force):
        fasta_tmp = get_fasta_path(tmp_dir, region.ref)
        ct_tmp = get_ct_path(tmp_dir, region, profile_name)
        vienna_tmp = (
            ct_tmp.with_suffix(path.VIENNA_EXT)
            if fold_backend == FOLD_BACKEND_VIENNARNA
            else None
        )
        db_tmp = (
            ct_tmp.with_suffix(path.DB_EXT)
            if fold_backend == FOLD_BACKEND_VIENNARNA
            else None
        )
        try:
            write_fasta(fasta_tmp, [(region.ref, region.seq.tr())], force=force)
            fold_temp_c = guess_temperature_to_celsius(fold_temp)
            if fold_backend == FOLD_BACKEND_VIENNARNA:
                run_rnafold(
                    fasta_tmp,
                    ct_tmp,
                    ct_sim,
                    vienna_tmp,
                    db_tmp,
                    sp_data=None,
                    sp_strategy=None,
                    eddy_prior_paired_file=None,
                    eddy_prior_unpaired_file=None,
                    fold_constraint=fold_constraint,
                    fold_commands=None,
                    fold_temp_c=fold_temp_c,
                    fold_isolated=False,
                    fold_md=fold_md,
                    fold_max=fold_max,
                    fold_mfe=fold_mfe,
                    fold_edelta=fold_edelta,
                    end5=region.end5,
                    num_cpus=num_cpus,
                )
            else:
                run_rnastructure(
                    fasta_tmp,
                    ct_tmp,
                    ct_sim,
                    pseudoknots=pseudoknots,
                    fold_temp_k=celsius_to_kelvin(fold_temp_c),
                    dms_file=None,
                    shape_file=None,
                    deigan_slope=None,
                    deigan_intercept=None,
                    fold_constraint=fold_constraint,
                    fold_isolated=False,
                    fold_md=fold_md,
                    fold_mfe=fold_mfe,
                    fold_max=fold_max,
                    fold_percent=fold_percent,
                    end5=region.end5,
                    num_cpus=num_cpus,
                )
        finally:
            if not keep_tmp:
                fasta_tmp.unlink(missing_ok=True)
                if vienna_tmp is not None:
                    vienna_tmp.unlink(missing_ok=True)
                    db_tmp.unlink(missing_ok=True)
                if ct_tmp != ct_sim:
                    ct_tmp.unlink(missing_ok=True)
    if not fold_mfe:
        num_structures = sum(1 for _ in from_ct(ct_sim))
        if num_structures < fold_min:
            ct_sim.unlink(missing_ok=True)
            raise ValueError(
                f"Folded {num_structures} structure(s) for "
                f"{region.ref}/{region.name} but need >= {fold_min}"
            )
    return ct_sim


@run_func(COMMAND, with_tmp=True, pass_keep_tmp=True)
def run(
    fasta: str | Path,
    *,
    sim_dir: str | Path,
    profile_name: str,
    probe: str,
    fold_backend: str,
    pseudoknots: bool,
    fold_coords: Iterable[tuple[str, int, int]],
    fold_primers: Iterable[tuple[str, DNA, DNA]],
    fold_regions_file: str | None,
    fold_constraint: str | None,
    fold_temp: float,
    fold_md: int,
    fold_mfe: bool,
    fold_max: int,
    fold_min: int,
    fold_percent: float,
    fold_edelta: float,
    keep_tmp: bool,
    tmp_dir: Path,
    force: bool,
    num_cpus: int,
):
    """
    Fold regions of a reference FASTA file and write CT files.

    Parameters
    ----------
    fasta: str | Path
        Path to the reference FASTA file.
    sim_dir: str | Path
        Simulation output directory for writing CT parameter files.
    profile_name: str
        Profile label embedded in each output CT file path.
    fold_coords: Iterable[tuple[str, int, int]]
        Explicit (ref, end5, end3) coordinate tuples defining regions.
    fold_primers: Iterable[tuple[str, DNA, DNA]]
        Primer sequences used to define region boundaries.
    fold_regions_file: str | None
        Path to a file listing regions to fold; None to use coords/primers.
    fold_constraint: str | None
        Path to a folding constraint file; None for no constraints.
    fold_temp: float
        Folding temperature (Celsius or Kelvin).
    fold_md: int
        Maximum pairing distance (0 for no limit).
    fold_mfe: bool
        Whether to predict only the minimum free energy structure.
    fold_max: int
        Maximum number of structures per region.
    fold_min: int
        Minimum number of structures required per region.
    fold_percent: float
        Maximum percent energy difference from the MFE structure.
    fold_edelta: float
        Maximum absolute energy difference (kcal/mol) from the MFE for
        suboptimal structures (RNAsubopt only).
    keep_tmp: bool
        Whether to retain temporary files.
    tmp_dir: Path
        Directory for temporary files (injected by `run_func`).
    force: bool
        Whether to overwrite existing CT files.
    num_cpus: int
        Number of CPU cores to use.

    Returns
    -------
    list[Path]
        Paths of all written CT files.
    """
    fold_backend = resolve_fold_backend(probe, fold_backend)
    check_fold_deps(fold_backend, pseudoknots)
    if not fold_mfe and fold_min > fold_max:
        logger.warning(
            "fold_min ({}) > fold_max ({}): setting fold_min = fold_max",
            fold_min,
            fold_max,
        )
        fold_min = fold_max
    # List the regions.
    regions = RefRegions(
        parse_fasta(Path(fasta), DNA),
        regs_file=(Path(fold_regions_file) if fold_regions_file else None),
        ends=fold_coords,
        primers=fold_primers,
    )
    sim_dir = Path(sim_dir)
    sim_dir.mkdir(parents=True, exist_ok=True)
    return dispatch(
        fold_region,
        num_cpus=num_cpus,
        pass_num_cpus=True,
        as_list=True,
        ordered=False,
        raise_on_error=False,
        args=as_list_of_tuples(regions.regions),
        kwargs=dict(
            sim_dir=sim_dir,
            tmp_dir=tmp_dir,
            profile_name=profile_name,
            fold_backend=fold_backend,
            pseudoknots=pseudoknots,
            fold_constraint=optional_path(fold_constraint),
            fold_temp=fold_temp,
            fold_md=fold_md,
            fold_mfe=fold_mfe,
            fold_max=fold_max,
            fold_min=fold_min,
            fold_percent=fold_percent,
            fold_edelta=fold_edelta,
            keep_tmp=keep_tmp,
            force=force,
        ),
    )


params = [
    arg_fasta,
    opt_sim_dir,
    opt_tmp_pfx,
    opt_profile_name,
    opt_probe,
    opt_fold_backend,
    opt_pseudoknots,
    opt_fold_regions_file,
    opt_fold_coords,
    opt_fold_primers,
    opt_fold_constraint,
    opt_fold_temp,
    opt_fold_md,
    opt_fold_mfe,
    opt_fold_max,
    opt_fold_min,
    opt_fold_percent,
    opt_fold_edelta,
    opt_keep_tmp,
    opt_force,
    opt_num_cpus,
]


@command(COMMAND, params=params)
def cli(*args, **kwargs):
    """Simulate secondary structure(s) a reference sequence."""
    run(*args, **kwargs)
