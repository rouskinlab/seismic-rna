""" Wrapper around RNAstructure from the Mathews Lab at the University
of Rochester: https://rna.urmc.rochester.edu/RNAstructure.html
"""

import os
import re
from pathlib import Path
from shutil import which

from .dryrun import dry_run
from .profile import ZERO_CELSIUS
from ..core.arg import opt_fold_temp
from ..core.error import IncompatibleValuesError
from ..core.extern import (RNASTRUCTURE_FOLD_CMD,
                           RNASTRUCTURE_FOLD_SMP_CMD,
                           RNASTRUCTURE_SHAPEKNOTS_CMD,
                           args_to_cmd,
                           run_cmd)
from ..core.logs import logger
from ..core.rna import renumber_ct
from ..core.write import need_write, write_mode

ENERGY_UNIT = "kcal/mol"
FOLD_SMP_NUM_THREADS = "OMP_NUM_THREADS"
DATAPATH = "DATAPATH"
DATAPATH_FILES = """
autodetect.txt
average_ensemble_defect.model
coaxial.dat
coaxstack.dat
dangle.dat
data_assemble_training_Multifind_predict_ensemble_z_final_svmformat.model
description.txt
design.DNA.Helices.dat
design.DNA.Loops.dat
design.RNA.Helices.dat
design.RNA.Loops.dat
dists
dna.coaxial.dg
dna.coaxial.dh
dna.coaxstack.dg
dna.coaxstack.dh
dna.dangle.dg
dna.dangle.dh
dna.dynalignmiscloop.dg
dna.dynalignmiscloop.dh
dna.hexaloop.dg
dna.hexaloop.dh
dna.int11.dg
dna.int11.dh
dna.int21.dg
dna.int21.dh
dna.int22.dg
dna.int22.dh
dna.loop.dg
dna.loop.dh
dna.miscloop.dg
dna.miscloop.dh
dna.specification.dat
dna.stack.dg
dna.stack.dh
dna.tloop.dg
dna.tloop.dh
dna.triloop.dg
dna.triloop.dh
dna.tstack.dg
dna.tstack.dh
dna.tstackcoax.dg
dna.tstackcoax.dh
dna.tstackh.dg
dna.tstackh.dh
dna.tstacki.dg
dna.tstacki.dh
dna.tstacki1n.dg
dna.tstacki1n.dh
dna.tstacki23.dg
dna.tstacki23.dh
dna.tstackm.dg
dna.tstackm.dh
dnacoaxial.dat
dnacoaxstack.dat
dnadangle.dat
dnahexaloop.dat
dnaint11.dat
dnaint21.dat
dnaint22.dat
dnaloop.dat
dnamiscloop.dat
dnastack.dat
dnatloop.dat
dnatriloop.dat
dnatstack.dat
dnatstackcoax.dat
dnatstackh.dat
dnatstacki.dat
dnatstacki1n.dat
dnatstacki23.dat
dnatstackm.dat
dynalignmiscloop.dat
fam_hmm_pars.dat
helix.dat
helixdr.dat
hexaloop.dat
hexaloop.dh
int11.dat
int21.dat
int22.dat
loop.dat
miscloop.dat
new_training_z_ave.scale.model
new_training_z_std.scale.model
pseudconst.dat
rna.coaxial.dg
rna.coaxial.dh
rna.coaxstack.dg
rna.coaxstack.dh
rna.cov.dg
rna.cov.dh
rna.dangle.dg
rna.dangle.dh
rna.dynalignmiscloop.dg
rna.hexaloop.dg
rna.hexaloop.dh
rna.int11.dg
rna.int11.dh
rna.int21.dg
rna.int21.dh
rna.int22.dg
rna.int22.dh
rna.loop.dg
rna.loop.dh
rna.miscloop.dg
rna.miscloop.dh
rna.param_map.dg
rna.specification.dat
rna.stack.dg
rna.stack.dh
rna.tloop.dg
rna.tloop.dh
rna.triloop.dg
rna.triloop.dh
rna.tstack.dg
rna.tstack.dh
rna.tstackcoax.dg
rna.tstackcoax.dh
rna.tstackh.dg
rna.tstackh.dh
rna.tstacki.dg
rna.tstacki.dh
rna.tstacki1n.dg
rna.tstacki1n.dh
rna.tstacki23.dg
rna.tstacki23.dh
rna.tstackm.dg
rna.tstackm.dh
rsample
stack.dat
stack.dh
stack.ds
stackdr.dat
stackdr.dh
stackdr.ds
std_ensemble_defect.model
tloop.dat
triloop.dat
tstack.dat
tstackcoax.dat
tstackh.dat
tstacki.dat
tstacki1n.dat
tstacki23.dat
tstackm.dat
"""


def check_data_path(data_path: str | Path | None = None) -> Path:
    """ Confirm the DATAPATH environment variable indicates the correct
    directory. """
    if data_path is None:
        data_path = os.environ.get(DATAPATH)
    if data_path is None:
        raise OSError(f"The {DATAPATH} environment variable is not set")
    if not isinstance(data_path, Path):
        data_path = Path(data_path)
    # Check if the path indicated by DATAPATH exists on the file system.
    if not data_path.is_dir():
        raise FileNotFoundError(f"{data_path} is not a directory")
    # Check if all expected files in the DATAPATH directory exist.
    extant_files = set(os.listdir(data_path))
    for expected_file in DATAPATH_FILES.strip().split():
        if expected_file not in extant_files:
            raise FileNotFoundError(f"{data_path} is missing the required "
                                    f"file {repr(expected_file)}")
    return data_path


def _guess_data_path_conda():
    """ Guess the DATAPATH if RNAstructure was installed with Conda. """
    fold_path = which(RNASTRUCTURE_FOLD_CMD)
    if fold_path is None:
        raise OSError(
            f"RNAstructure not seem to be installed: {RNASTRUCTURE_FOLD_CMD}"
        )
    fold_path = Path(fold_path)
    env_dir = fold_path.parent.parent
    data_path = env_dir.joinpath("share", "rnastructure", "data_tables")
    if not data_path.is_dir():
        raise OSError("It seems RNAstructure is not installed with Conda: "
                      f"{data_path} does not exist")
    check_data_path(data_path)
    logger.detail(f"Successfully guessed {DATAPATH}: {data_path}")
    return data_path


def _guess_data_path_manual():
    """ Guess the DATAPATH if RNAstructure was installed manually
    (e.g. by downloading from the Mathews Lab website). """
    fold_path = which(RNASTRUCTURE_FOLD_CMD)
    if fold_path is None:
        raise OSError(
            f"RNAstructure not seem to be installed: {RNASTRUCTURE_FOLD_CMD}"
        )
    fold_path = Path(fold_path)
    data_path = fold_path.parent.parent.joinpath("data_tables")
    check_data_path(data_path)
    logger.detail(f"Successfully guessed {DATAPATH}: {data_path}")
    return data_path


def guess_data_path():
    """ Guess the DATAPATH. """
    errors = list()
    try:
        return check_data_path()
    except OSError as error:
        errors.append(error)
        logger.warning(f"The {DATAPATH} environment variable is not valid; "
                       f"attempting to guess it")
    for guess_func in [_guess_data_path_conda,
                       _guess_data_path_manual]:
        try:
            return guess_func()
        except OSError as error:
            errors.append(error)
    raise OSError("\n".join(f" -> {error}" for error in errors))


def require_data_path():
    """ Return an error message if the DATAPATH is not valid. """
    try:
        data_path = guess_data_path()
    except OSError as error:
        raise OSError(
            f"RNAstructure requires an environment variable called {DATAPATH} "
            f"to point to the directory in which its thermodynamic parameters "
            f"are located, but\n{error}\nFor more information, please refer to "
            f"https://rna.urmc.rochester.edu/Text/Thermodynamics.html"
        ) from None
    # Set the DATAPATH environment variable if it is not already set.
    data_path_str = str(data_path)
    if os.environ.get(DATAPATH) != data_path_str:
        os.environ[DATAPATH] = data_path_str
    return data_path


def make_rnastructure_cmd(fasta_file: Path,
                          ct_file: Path, *,
                          pseudoknots: bool,
                          fold_constraint: Path | None,
                          dms_file: Path | None,
                          shape_file: Path | None,
                          deigan_intercept: float | None,
                          deigan_slope: float | None,
                          fold_temp_k: float | None,
                          fold_isolated: bool,
                          fold_md: int,
                          fold_mfe: bool,
                          fold_max: int,
                          fold_percent: float,
                          num_cpus: int = 1):
    """ Make a command for 'Fold', 'Fold-smp', or 'ShapeKnots'. """
    if pseudoknots:
        if num_cpus > 1:
            logger.warning(
                f"ShapeKnots cannot use {num_cpus} threads; defaulting to 1"
            )
        args = [RNASTRUCTURE_SHAPEKNOTS_CMD]
    else:
        if num_cpus > 1:
            # Fold with multiple threads using the Fold-smp program.
            args = [RNASTRUCTURE_FOLD_SMP_CMD]
            os.environ[FOLD_SMP_NUM_THREADS] = str(num_cpus)
        else:
            # Fold with one thread using the Fold program.
            args = [RNASTRUCTURE_FOLD_CMD]
    if fold_constraint is not None:
        # File of constraints.
        args.extend(["--constraint", fold_constraint])
    if dms_file is not None:
        # File of DMS reactivities.
        if shape_file is not None:
            raise IncompatibleValuesError("Cannot give both DMS and SHAPE files")
        args.extend(["--DMS", dms_file])
    if shape_file is not None:
        # File of SHAPE reactivities.
        args.extend(["--SHAPE", shape_file])
    if deigan_intercept is not None:
        # SHAPE intercept parameter (kcal/mol).
        args.extend(["--SHAPEintercept", deigan_intercept])
    if deigan_slope is not None:
        # SHAPE slope parameter (kcal/mol).
        args.extend(["--SHAPEslope", deigan_slope])
    if fold_temp_k is not None:
        # Temperature of folding (Kelvin).
        if pseudoknots:
            default_temp_k = opt_fold_temp.default + ZERO_CELSIUS
            if abs(fold_temp_k - default_temp_k) > 0.01:
                logger.warning(
                    f"ShapeKnots cannot fold at {fold_temp_k} K; "
                    f"defaulting to {default_temp_k} K"
                )
        else:
            args.extend(["--temperature", fold_temp_k])
    if fold_isolated:
        # Allow isolated pairs.
        if pseudoknots:
            logger.warning(
                "ShapeKnots does not support --isolated; "
                "isolated pairs cannot be allowed with this backend"
            )
        else:
            args.append("--isolated")
    if fold_md > 0:
        # Maximum distance between paired bases.
        args.extend(["--maxdistance", fold_md])
    if fold_mfe:
        # Predict only the minimum free energy structure.
        if pseudoknots:
            # ShapeKnots has no --MFE flag; request a single structure instead.
            args.extend(["--maximum", 1])
        else:
            args.append("--MFE")
    else:
        if fold_max > 0:
            # Maximum number of structures.
            args.extend(["--maximum", fold_max])
        if fold_percent > 0.:
            # Maximum % difference between free energies of structures.
            args.extend(["--percent", fold_percent])
    # Input and output files.
    args.extend([fasta_file, ct_file])
    return args_to_cmd(args)


def run_rnastructure(fasta_tmp: Path,
                     ct_tmp: Path,
                     ct_out: Path, *,
                     pseudoknots: bool,
                     fold_temp_k: float | None,
                     dms_file: Path | None,
                     shape_file: Path | None,
                     deigan_slope: float | None,
                     deigan_intercept: float | None,
                     fold_constraint: Path | None,
                     fold_isolated: bool,
                     fold_md: int,
                     fold_mfe: bool,
                     fold_max: int,
                     fold_percent: float,
                     end5: int,
                     num_cpus: int,
                     fold_dry_run: bool = False):
    """ Run Fold/ShapeKnots on pre-built paths, retitle, and renumber. """
    if not pseudoknots and num_cpus > 1:
        parallel_options = [True, False]
    else:
        parallel_options = [False]
    fold_cmds = [
        make_rnastructure_cmd(
            fasta_tmp,
            ct_tmp,
            pseudoknots=pseudoknots,
            fold_temp_k=fold_temp_k,
            dms_file=dms_file,
            shape_file=shape_file,
            deigan_slope=deigan_slope,
            deigan_intercept=deigan_intercept,
            fold_constraint=fold_constraint,
            fold_isolated=fold_isolated,
            fold_md=fold_md,
            fold_mfe=fold_mfe,
            fold_max=fold_max,
            fold_percent=fold_percent,
            num_cpus=(num_cpus if parallel else 1),
        )
        for parallel in parallel_options
    ]
    if fold_dry_run:
        dry_run(fold_cmds, ct_tmp.parent)
    else:
        fold_cmd = fold_cmds.pop(0)
        while True:
            try:
                run_cmd(fold_cmd)
            except RuntimeError as error:
                if fold_cmds:
                    logger.warning(error)
                    fold_cmd = fold_cmds.pop(0)
                else:
                    raise
            else:
                break
        retitle_ct(ct_tmp, ct_tmp, force=True)
        renumber_ct(ct_tmp, ct_out, end5, force=True)


class RNAStructureConnectivityTableTitleLineFormatError(ValueError):
    """ Error in the format of a CT title line from RNAStructure. """


class ConnectivityTableAlreadyRetitledError(RuntimeError):
    """ A CT file was already retitled. """


def parse_rnastructure_ct_title(line: str):
    """ Parse a title in a CT file from RNAstructure, in this format:

    {length}  ENERGY = {energy}  {ref}

    where {length} is the number of positions in the structure, {ref} is
    the name of the reference, and {energy} is the predicted free energy
    of folding.
    Also handle the edge case when RNAstructure predicts no base pairs
    (and thus does not write the free energy) by returning 0.

    Parameters
    ----------
    line: str
        Line containing the title of the structure.

    Returns
    -------
    tuple[int, float, str]
        Tuple of number of positions in the structure, predicted free
        energy of folding, and name of the reference sequence.
    """
    # Parse the line assuming it contains an energy term.
    if m := re.match(r"\s*([0-9]+)\s+ENERGY = (-?[0-9.]+)\s+(\S+)", line):
        length, energy, ref = m.groups()
    else:
        # If that failed, then check if the line was already retitled.
        try:
            parse_energy(line)
        except (ValueError, TypeError):
            # The line was not retitled: parse it assuming it has no
            # energy term (which happens if no base pairs exist).
            if m := re.match(r"\s*([0-9]+)\s+(\S+)", line):
                length, ref = m.groups()
            else:
                # The line violated the basic length-and-title format.
                raise RNAStructureConnectivityTableTitleLineFormatError(line)
            logger.warning("CT line contains no energy term (probably because "
                           f"no base pairs were predicted): {repr(line)}")
            energy = 0.
        else:
            raise ConnectivityTableAlreadyRetitledError(line)
    return int(length), float(energy), ref


def format_retitled_ct_line(length: int, ref: str, uniqid: int, energy: float):
    """ Format a new CT title line including unique identifiers:

    {length}    {ref} #{uniqid}: {energy}

    where {length} is the number of positions in the structure (required
    for all CT files), {ref} is the name of the reference, {uniqid} is
    the unique identifier, and {energy} is the free energy of folding.

    Parameters
    ----------
    length: int
        Number of positions in the structure.
    uniqid: int
        Unique identifier (non-negative integer).
    ref: str
        Name of the reference.
    energy: float
        Free energy of folding (kcal/mol).

    Returns
    -------
    str
        Formatted CT title line.
    """
    return f"{length}\t{ref} #{uniqid}: {energy} {ENERGY_UNIT}\n"


def retitle_ct(ct_input: Path, ct_output: Path, force: bool = False):
    """ Retitle the structures in a CT file produced by RNAstructure.

    The default titles follow this format:

    ENERGY = {energy}  {reference}

    where {reference} is the name of the reference sequence and {energy}
    is the predicted free energy of folding.

    The major problem with this format is that structures can have equal
    predicted free energies, so the titles of the structures can repeat,
    which would cause some functions (e.g. graphing ROC curves) to fail.

    This function assigns a unique integer to each structure (starting
    with 0 for the minimum free energy and continuing upwards), which
    ensures that no two structures have identical titles.

    Parameters
    ----------
    ct_input: Path
        Path of the CT file to retitle.
    ct_output: Path
        Path of the CT file to which to write the retitled information.
    force: bool = False
        Overwrite the output CT file if it already exists.
    """
    if need_write(ct_output, force):
        # Read all lines from the input file.
        lines = list()
        with open(ct_input) as f:
            uniqid = 0
            while title_line := f.readline():
                # Parse and reformat the title line.
                n, energy, ref = parse_rnastructure_ct_title(title_line)
                lines.append(format_retitled_ct_line(n, ref, uniqid, energy))
                # Add the lines that encode the structure.
                for _ in range(n):
                    lines.append(f.readline())
                uniqid += 1
        # Write the reformatted lines to the output file.
        text = "".join(lines)
        with open(ct_output, write_mode(force=True)) as f:
            f.write(text)
        logger.routine(f"Retitled CT file {ct_input}"
                       + (f" to {ct_output}"
                          if ct_input != ct_output
                          else ""))


def parse_energy(line: str):
    """ Parse the predicted free energy of folding from a line in format

    {length}    {ref} #{uniqid}: {energy}

    where {length} is the number of positions in the structure (required
    for all CT files), {ref} is the name of the reference, {uniqid} is
    the unique identifier, and {energy} is the free energy of folding.

    Parameters
    ----------
    line: str
        Line from which to parse the energy.

    Returns
    -------
    float
        Free energy of folding.
    """
    _, energy = line.split(": ")
    value, unit = energy.split()
    if unit != ENERGY_UNIT:
        raise ValueError(f"Expected energy to have units of {ENERGY_UNIT}, "
                         f"but got {repr(unit)}")
    return float(value)
