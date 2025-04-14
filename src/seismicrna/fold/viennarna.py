""" Wrapper around ViennaRNA from Lorenz and Hofacker at the
University of Vienna: https://www.tbi.univie.ac.at/RNA/
"""

import os
import re
import math
import numpy as np
import pandas as pd
from pathlib import Path

from ..core.arg import docdef
from ..core.extern import (VIENNA_RNAFOLD_CMD,
                           VIENNA_RNASUBOPT_CMD,
                           cmds_to_redirect_in,
                           cmds_to_redirect_out,
                           args_to_cmd,
                           run_cmd)
from ..core.logs import logger
from ..core.rna import RNAProfile, renumber_ct, db_to_ct
from ..core.write import need_write, write_mode

from .parameterize_archiveii import calc_rnastructure_defaults

ENERGY_UNIT = "kcal/mol"
RNAFOLD_NUM_THREADS = "OMP_NUM_THREADS"


def make_fold_cmd(fasta_file: Path,
                  vienna_file: Path, *,
                  fold_delta_e: float,
                  dms_file: Path | None,
                  shape_method: str | None,
                  fold_constraint: Path | None,
                  fold_commands: Path | None,
                  fold_temp: float,
                  fold_md: int,
                  fold_mfe: bool,
                  num_cpus: int = 1,
                  **kwargs):
    if fold_mfe:
        cmd = [VIENNA_RNAFOLD_CMD]
    else:
        cmd = [VIENNA_RNASUBOPT_CMD, "-e", fold_delta_e]
        os.environ[RNAFOLD_NUM_THREADS] = str(num_cpus)
    cmd.append("--noLP")
    if dms_file is not None:
        # File of DMS reactivities.
        assert shape_method is not None
        cmd.extend(["--shape", dms_file,
                    "--shapeMethod", shape_method])
    if fold_constraint is not None:
        # File of constraints.
        cmd.extend(["--constraint", fold_constraint])
    if fold_commands is not None:
        # File of commands.
        cmd.extend(["--commands", fold_commands])
    # Temperature of folding (Kelvin).
    cmd.extend(["--temp", fold_temp])
    if fold_md > 0:
        # Maximum distance between paired bases.
        cmd.extend(["--maxBPspan", fold_md])
    # Input and output files.
    cmd = args_to_cmd(cmd)
    cmd = cmds_to_redirect_in([cmd, str(fasta_file)])
    cmd = cmds_to_redirect_out([cmd, str(vienna_file)])
    return cmd


@docdef.auto()
def rnafold(rna: RNAProfile, *,
         branch: str,
         fold_temp: float,
         fold_constraint: Path | None = None,
         fold_md: int,
         fold_mfe: bool,
         fold_max: int,
         fold_percent: float,
         out_dir: Path,
         tmp_dir: Path,
         keep_tmp: bool,
         num_cpus: int,
         **kwargs):
    """ Run the 'RNAFold' or 'RNAsubopt' program of ViennaRNA. """
    logger.routine(f"Began folding {rna}")
    ct_out = rna.get_ct_file(out_dir, branch)
    # Temporary FASTA file for the RNA.
    fasta_tmp = rna.to_fasta(tmp_dir, branch)
    # Path of the temporary CT file.
    ct_tmp = rna.get_ct_file(tmp_dir, branch)
    # Path of the temporary vienna file.
    vienna_tmp = rna.get_vienna_file(tmp_dir, branch)
    # DMS reactivities file for the RNA.
    dms_file = rna.to_dms(tmp_dir, branch)

    #TODO Reimplement builtin method
    dms_data = pd.read_table(dms_file, header=None, index_col=0, names=["Position", "Reactivity"])
    dms = dms_data["Reactivity"].to_numpy()
    _, _, pseudoenergies = calc_rnastructure_defaults(dms)
    dms_data["Cordero"] = [pseudoenergies[i] for i in range(len(dms_data.index))]
    b = min(dms_data["Cordero"])
    m = (max(dms_data["Cordero"]) - b)/math.log(2)
    shape_method = f"Dm{m}b{b}"
    dms_data["Scaled Reactivity"] = (np.expm1((dms_data["Cordero"] - b) / m))
    dms_data["Scaled Reactivity"].to_csv(dms_file, sep='\t', header=False)
    dms_data["Deigan_Transformed"] = (m * np.log1p(dms_data["Scaled Reactivity"])) + b
    assert np.allclose(dms_data["Cordero"], dms_data["Deigan_Transformed"])
    fold_commands = None # TODO: enable command files
    try:
        # Run the command.
        fold_cmds = {
            smp: make_fold_cmd(fasta_tmp,
                               vienna_tmp,
                               fold_delta_e=5, #TODO set default
                               dms_file=dms_file,
                               shape_method=shape_method,
                               fold_constraint=fold_constraint,
                               fold_commands=fold_commands,
                               fold_temp=fold_temp-273.15,
                               fold_md=fold_md,
                               fold_mfe=fold_mfe,
                               fold_max=fold_max,
                               fold_percent=fold_percent,
                               num_cpus=(num_cpus if smp else 1))
            for smp in [True, False]
        }
        try:
            run_cmd(fold_cmds[True])
        except RuntimeError as error:
            logger.warning(error)
            run_cmd(fold_cmds[False])
        db_to_ct(vienna_tmp, force=True)
        # Reformat the CT file title lines so that each is unique.
        retitle_ct(ct_tmp, ct_tmp, force=True)
        # Renumber the CT file so that it has the same numbering scheme
        # as the region, rather than always starting at 1, the latter
        # of which is always output by the Fold program.
        renumber_ct(ct_tmp, ct_out, rna.region.end5, force=True)
    finally:
        if not keep_tmp:
            # Delete the temporary files.
            fasta_tmp.unlink(missing_ok=True)
            dms_file.unlink(missing_ok=True)
            if ct_tmp != ct_out:
                ct_tmp.unlink(missing_ok=True)
    logger.routine(f"Ended folding {rna}")
    return ct_out


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
