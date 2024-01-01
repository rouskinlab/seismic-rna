"""
Struct -- RNAstructure Module

Wrapper around RNAstructure from the Mathews Lab at U of Rochester:
https://rna.urmc.rochester.edu/RNAstructure.html
"""

import re
from logging import getLogger
from pathlib import Path

from ..core import path
from ..core.extern import (RNASTRUCTURE_CT2DOT_CMD,
                           RNASTRUCTURE_DOT2CT_CMD,
                           RNASTRUCTURE_FOLD_CMD,
                           args_to_cmd,
                           run_cmd)
from ..core.rna import RNAProfile, renumber_ct
from ..core.write import need_write, write_mode

logger = getLogger(__name__)


def fold(rna: RNAProfile, *,
         fold_temp: float,
         fold_constraint: Path | None,
         fold_md: int,
         fold_mfe: bool,
         fold_max: int,
         fold_percent: float,
         out_dir: Path,
         temp_dir: Path,
         keep_temp: bool,
         force: bool):
    """ Run the 'Fold' program of RNAstructure. """
    ct_file = rna.get_ct_file(out_dir)
    if need_write(ct_file, force):
        cmd = [RNASTRUCTURE_FOLD_CMD, "--temperature", fold_temp]
        # Constraints.
        if fold_constraint is not None:
            cmd.extend(["--constraint", fold_constraint])
        # Maximum distance between paired bases.
        if fold_md > 0:
            cmd.extend(["--maxdistance", fold_md])
        # Predict only the minimum free energy structure.
        if fold_mfe:
            cmd.append("--MFE")
        # Maximum number of structures.
        cmd.extend(["--maximum", fold_max])
        # Maximum % difference between free energies of structures.
        cmd.extend(["--percent", fold_percent])
        # DMS reactivities file for the RNA.
        cmd.extend(["--DMS", dms_file := rna.to_dms(temp_dir)])
        # Temporary FASTA file for the RNA.
        cmd.append(fasta := rna.to_fasta(temp_dir))
        # Path of the temporary CT file.
        cmd.append(ct_temp := rna.get_ct_file(temp_dir))
        try:
            # Run the command.
            run_cmd(args_to_cmd(cmd))
            # Reformat the CT file title lines so that each is unique.
            retitle_ct_structures(ct_temp, ct_temp, force=True)
            # Renumber the CT file so that it has the same numbering
            # scheme as the section, rather than always starting at 1,
            # the latter of which is always output by the Fold program.
            renumber_ct(ct_temp, ct_file, rna.section.end5, force)
        finally:
            if not keep_temp:
                # Delete the temporary files.
                fasta.unlink(missing_ok=True)
                dms_file.unlink(missing_ok=True)
                ct_temp.unlink(missing_ok=True)
    return ct_file


def ct2dot(ct_file: Path, number: int | str = "all"):
    """ Make a dot-bracket (DB) file of a connectivity-table (CT) file.

    Parameters
    ----------
    ct_file: pathlib.Path
        Path to the CT file.
    number: int | str = "all"
        Number of the structure to convert, or "all" to convert all.

    Returns
    -------
    pathlib.Path
        Path to the DB file.
    """
    db_file = ct_file.with_suffix(path.DB_EXT)
    cmd = [RNASTRUCTURE_CT2DOT_CMD, ct_file, number, db_file]
    run_cmd(args_to_cmd(cmd))
    return db_file


def dot2ct(db_file: Path):
    """ Make a connectivity-table (CT) file of a dot-bracket (DB) file.

    Parameters
    ----------
    db_file: pathlib.Path
        Path to the DB file.

    Returns
    -------
    pathlib.Path
        Path to the CT file.
    """
    ct_file = db_file.with_suffix(path.CT_EXT)
    cmd = [RNASTRUCTURE_DOT2CT_CMD, db_file, ct_file]
    run_cmd(args_to_cmd(cmd))
    return ct_file


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
        Free energy of folding.

    Returns
    -------
    str
        Formatted CT title line.
    """
    return f"{length}\t{ref} #{uniqid}: {energy}\n"


def retitle_ct_structures(ct_input: Path, ct_output: Path, force: bool = False):
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
        with open(ct_output, write_mode(force)) as f:
            f.write(text)


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
    _, energy = line.split(":")
    return float(energy)

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
