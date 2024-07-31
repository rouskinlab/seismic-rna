"""
Struct -- RNAstructure Module

Wrapper around RNAstructure from the Mathews Lab at U of Rochester:
https://rna.urmc.rochester.edu/RNAstructure.html
"""

import os
import re
from logging import getLogger
from pathlib import Path
from shutil import which

from ..core.arg import docdef
from ..core.extern import (RNASTRUCTURE_FOLD_CMD,
                           RNASTRUCTURE_FOLD_SMP_CMD,
                           args_to_cmd,
                           run_cmd)
from ..core.rna import RNAProfile, renumber_ct
from ..core.write import need_write, write_mode

logger = getLogger(__name__)

ENERGY_UNIT = "kcal/mol"
FOLD_SMP_NUM_THREADS = "OMP_NUM_THREADS"
DATAPATH = "DATAPATH"
DATAPATH_FILES = """
autodetect.txt
average_ensemble_defect.model
b-test.coaxial.dg
b-test.coaxial.dh
b-test.coaxstack.dg
b-test.coaxstack.dh
b-test.dangle.dg
b-test.dangle.dh
b-test.dynalignmiscloop.dg
b-test.hexaloop.dg
b-test.hexaloop.dh
b-test.int11.dg
b-test.int11.dh
b-test.int21.dg
b-test.int21.dh
b-test.int22.dg
b-test.int22.dh
b-test.loop.dg
b-test.miscloop.dg
b-test.specification.dat
b-test.stack.dg
b-test.stack.dh
b-test.tloop.dg
b-test.tloop.dh
b-test.triloop.dg
b-test.triloop.dh
b-test.tstack.dg
b-test.tstack.dh
b-test.tstackcoax.dg
b-test.tstackcoax.dh
b-test.tstackh.dg
b-test.tstackh.dh
b-test.tstacki.dg
b-test.tstacki.dh
b-test.tstacki1n.dg
b-test.tstacki1n.dh
b-test.tstacki23.dg
b-test.tstacki23.dh
b-test.tstackm.dg
b-test.tstackm.dh
coaxial.dat
coaxial.dh
coaxstack.dat
coaxstack.dh
dangle.dat
dangle.dh
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
dnacoaxial.dh
dnacoaxstack.dat
dnacoaxstack.dh
dnadangle.dat
dnadangle.dh
dnadynalignmiscloop.dat
dnadynalignmiscloop.dh
dnahexaloop.dat
dnahexaloop.dh
dnaint11.dat
dnaint11.dh
dnaint21.dat
dnaint21.dh
dnaint22.dat
dnaint22.dh
dnaloop.dat
dnaloop.dh
dnamiscloop.dat
dnamiscloop.dh
dnastack.dat
dnastack.dh
dnatloop.dat
dnatloop.dh
dnatriloop.dat
dnatriloop.dh
dnatstack.dat
dnatstack.dh
dnatstackcoax.dat
dnatstackcoax.dh
dnatstackh.dat
dnatstackh.dh
dnatstacki.dat
dnatstacki.dh
dnatstacki1n.dat
dnatstacki1n.dh
dnatstacki23.dat
dnatstacki23.dh
dnatstackm.dat
dnatstackm.dh
dynalignmiscloop.dat
fam_hmm_pars.dat
helix.dat
helixdr.dat
hexaloop.dat
hexaloop.dh
int11.dat
int11.dh
int21.dat
int21.dh
int22-exp.dh
int22.dat
int22.dh
loop.dat
loop.dh
miscloop.dat
miscloop.dh
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
tloop.dh
triloop.dat
triloop.dh
tstack.dat
tstack.dh
tstackcoax.dat
tstackcoax.dh
tstackh.dat
tstackh.dh
tstacki.dat
tstacki.dh
tstacki1n.dat
tstacki1n.dh
tstacki23.dat
tstacki23.dh
tstackm.dat
tstackm.dh
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
    logger.debug(f"Successfully guessed {DATAPATH}: {data_path}")
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
    logger.debug(f"Successfully guessed {DATAPATH}: {data_path}")
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
    for attempt in [_guess_data_path_conda, _guess_data_path_manual]:
        try:
            return attempt()
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


def make_fold_cmd(fasta_file: Path,
                  ct_file: Path, *,
                  dms_file: Path | None,
                  fold_constraint: Path | None,
                  fold_temp: float,
                  fold_md: int,
                  fold_mfe: bool,
                  fold_max: int,
                  fold_percent: float,
                  n_procs: int = 1):
    if n_procs > 1:
        # Fold with multiple threads using the Fold-smp program.
        cmd = [RNASTRUCTURE_FOLD_SMP_CMD]
        os.environ[FOLD_SMP_NUM_THREADS] = str(n_procs)
    else:
        # Fold with one thread using the Fold program.
        cmd = [RNASTRUCTURE_FOLD_CMD]
    if dms_file is not None:
        # File of DMS reactivities.
        cmd.extend(["--DMS", dms_file])
    if fold_constraint is not None:
        # File of constraints.
        cmd.extend(["--constraint", fold_constraint])
    # Temperature of folding (Kelvin).
    cmd.extend(["--temperature", fold_temp])
    if fold_md > 0:
        # Maximum distance between paired bases.
        cmd.extend(["--maxdistance", fold_md])
    if fold_mfe:
        # Predict only the minimum free energy structure.
        cmd.append("--MFE")
    else:
        # Maximum number of structures.
        cmd.extend(["--maximum", fold_max])
        # Maximum % difference between free energies of structures.
        cmd.extend(["--percent", fold_percent])
    # Input and output files.
    cmd.append(fasta_file)
    cmd.append(ct_file)
    return cmd


@docdef.auto()
def fold(rna: RNAProfile, *,
         fold_temp: float,
         fold_constraint: Path | None = None,
         fold_md: int,
         fold_mfe: bool,
         fold_max: int,
         fold_percent: float,
         out_dir: Path,
         tmp_dir: Path,
         keep_tmp: bool,
         n_procs: int):
    """ Run the 'Fold' or 'Fold-smp' program of RNAstructure. """
    ct_out = rna.get_ct_file(out_dir)
    # Temporary FASTA file for the RNA.
    fasta_tmp = rna.to_fasta(tmp_dir)
    # Path of the temporary CT file.
    ct_tmp = rna.get_ct_file(tmp_dir)
    # DMS reactivities file for the RNA.
    dms_file = rna.to_dms(tmp_dir)
    try:
        # Run the command.
        fold_cmds = {
            smp: args_to_cmd(make_fold_cmd(fasta_tmp,
                                           ct_tmp,
                                           dms_file=dms_file,
                                           fold_constraint=fold_constraint,
                                           fold_temp=fold_temp,
                                           fold_md=fold_md,
                                           fold_mfe=fold_mfe,
                                           fold_max=fold_max,
                                           fold_percent=fold_percent,
                                           n_procs=(n_procs if smp else 1)))
            for smp in [True, False]
        }
        try:
            run_cmd(fold_cmds[True])
        except RuntimeError as error:
            logger.warning(
                f"Unable to fold using {RNASTRUCTURE_FOLD_SMP_CMD}:\n{error}"
            )
            run_cmd(fold_cmds[False])
        # Reformat the CT file title lines so that each is unique.
        retitle_ct(ct_tmp, ct_tmp, force=True)
        # Renumber the CT file so that it has the same numbering scheme
        # as the section, rather than always starting at 1, the latter
        # of which is always output by the Fold program.
        renumber_ct(ct_tmp, ct_out, rna.section.end5, force=True)
    finally:
        if not keep_tmp:
            # Delete the temporary files.
            fasta_tmp.unlink(missing_ok=True)
            dms_file.unlink(missing_ok=True)
            if ct_tmp != ct_out:
                ct_tmp.unlink(missing_ok=True)
    logger.info(f"Predicted structure of {rna} to {ct_out}")
    return ct_out


def parse_rnastructure_ct_title(line: str):
    """ Parse a title in a CT file from RNAstructure, in this format:

    {length}  ENERGY = {energy}  {ref}

    where {length} is the number of positions in the structure, {ref} is
    the name of the reference, and {energy} is the predicted free energy
    of folding.
    Also handle the edge case when RNAstructure predicts no base pairs
    (and thus does not write the free energy) by returning NaN.

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
        # If that failed, then parse the line assuming it does not.
        if m := re.match(r"\s*([0-9]+)\s+(\S+)", line):
            length, ref = m.groups()
        else:
            # The line violated the basic length-and-title format.
            raise ValueError(f"Failed to parse CT title line: {repr(line)}")
        logger.warning("CT line contains no energy term (probably because no "
                       f"base pairs were predicted): {repr(line)}")
        energy = "nan"
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
        logger.info(f"Retitled CT file {ct_input}"
                    + (f" to {ct_output}" if ct_input != ct_output else ""))


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
