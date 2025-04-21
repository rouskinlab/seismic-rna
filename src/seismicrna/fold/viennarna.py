""" Wrapper around ViennaRNA from Lorenz and Hofacker at the
University of Vienna: https://www.tbi.univie.ac.at/RNA/
"""

import os
import numpy as np
import pandas as pd
import shutil

from pathlib import Path

from .rnastructure import retitle_ct

from .profile import RNAFoldProfile
from ..core.arg import docdef
from ..core.extern import (VIENNA_RNAFOLD_CMD,
                           VIENNA_RNASUBOPT_CMD,
                           cmds_to_redirect_in,
                           cmds_to_redirect_out,
                           args_to_cmd,
                           run_cmd)
from ..core.logs import logger
from ..core.rna import renumber_ct, db_to_ct
from ..core.write import need_write, write_mode

ENERGY_UNIT = "kcal/mol"
RNAFOLD_NUM_THREADS = "OMP_NUM_THREADS"


ZERO_CELSIUS = 273.15  # Kelvin


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
    # Temperature of folding (convert Kelvin to Celsius).
    cmd.extend(["--temp", fold_temp - ZERO_CELSIUS])
    if fold_md > 0:
        # Maximum distance between paired bases.
        cmd.extend(["--maxBPspan", fold_md])
    # Input and output files.
    cmd = args_to_cmd(cmd)
    cmd = cmds_to_redirect_in([cmd, str(fasta_file)])
    cmd = cmds_to_redirect_out([cmd, str(vienna_file)])
    return cmd


def calc_bp_pseudoenergy(seq_len: int,
                         pseudoenergies: pd.Series,
                         out_file: Path):
    """
    Identify (i, j) pairs, with j >= i + 3, where the sum of pseudoenergies is non-zero.

    The output file will have lines formatted as ViennaRNA commands:
      "E i j 1 <bp_pseudoenergy>"

    Parameters:
      seq_len: The total length of the sequence. Positions are assumed to be numbered 1..seq_len.
      pseudoenergies: A pandas Series where each key is a position and its value is the energy.
                      If a position is missing, it is assumed to have energy 0.
      out_file: Path to the file where the results will be written.
    """
    # Create a full array of energies with default 0, for positions 1..seq_len.
    energies_full = np.zeros(seq_len, dtype=np.float64)
    # The Series index uses 1-indexing; fill in available energies.
    indices = np.array(pseudoenergies.index)
    energies_full[indices - 1] = pseudoenergies.values

    outputs = []  # Collect valid rows from all valid offsets.

    # Loop over allowed offsets (j - i) starting from 3 (i.e., j >= i + 3).
    for offset in range(3, seq_len):
        # Generate starting positions: 1, 2, …, seq_len - offset.
        i_vals = np.arange(1, seq_len - offset + 1)
        j_vals = i_vals + offset

        # Compute energy sums for each (i, j) pair using vectorized slicing.
        energy_sum = energies_full[i_vals - 1] + energies_full[j_vals - 1]

        # Filter pairs that have a positive energy sum.
        valid = energy_sum > 0
        if valid.any():
            valid_i = i_vals[valid]
            valid_j = j_vals[valid]
            valid_sum = energy_sum[valid]
            # Prepare data as columns: "E i j 1 valid_sum".
            data = np.column_stack((
                valid_i,
                valid_j,
                np.ones(valid_i.shape[0], dtype=int),
                valid_sum
            ))
            outputs.append(data)

    # If there are valid pairs, write them to file.
    if outputs:
        result = np.vstack(outputs)

        # Determine mode: append if file exists, else write.
        mode = "a" if out_file.exists() else "w"

        # If appending, check if we need to add a starting newline.
        need_newline = False
        if mode == "a" and out_file.stat().st_size > 0:
            with open(out_file, "rb") as f:
                # Move to the last character of the file.
                f.seek(-1, 2)
                last_char = f.read(1)
            if last_char != b'\n':
                need_newline = True

        with open(out_file, mode) as f:
            if need_newline:
                f.write("\n")
            # Write the result with the format: "E i j 1 <energy>".
            np.savetxt(f, result, fmt="E %d %d %d %.6f")

    return out_file


@docdef.auto()
def rnafold(rna: RNAFoldProfile, *,
            branch: str,
            fold_constraint: Path | None = None,
            fold_commands: Path | None = None,
            fold_md: int,
            fold_mfe: bool,
            fold_max: int,
            fold_percent: float,
            out_dir: Path,
            tmp_dir: Path,
            keep_tmp: bool,
            num_cpus: int):
    """ Run the 'RNAFold' or 'RNAsubopt' program of ViennaRNA. """
    logger.routine(f"Began folding {rna}")
    ct_out = rna.get_ct_file(out_dir, branch)
    # Temporary FASTA file for the RNA.
    fasta_tmp = rna.write_fasta(tmp_dir, branch)
    # Path of the temporary CT file.
    ct_tmp = rna.get_ct_file(tmp_dir, branch)
    # Path of the temporary DB file.
    db_tmp = rna.get_db_file(tmp_dir, branch)
    # Path of the temporary vienna file.
    vienna_tmp = rna.get_vienna_file(tmp_dir, branch)
    # Path of the temporary command file.
    command_tmp = rna.get_command_file(tmp_dir, branch)
    # DMS reactivities file for the RNA.
    dms_file = rna.to_pseudomus(tmp_dir, branch)

    apply_all_paired = False

    if fold_commands:
        shutil.copy2(fold_commands, command_tmp)
    if apply_all_paired:
        dms_file = None
        command_file = calc_bp_pseudoenergy(len(rna.seq), rna.pseudoenergies, command_tmp)
    else:
        command_file = None
    try:
        # Run the command.
        fold_cmds = {
            smp: make_fold_cmd(fasta_tmp,
                               vienna_tmp,
                               fold_delta_e=5,  # TODO set default
                               dms_file=dms_file,
                               shape_method=rna.shape_method,
                               fold_constraint=fold_constraint,
                               fold_commands=command_file,
                               fold_temp=rna.fold_temp,
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
        # Extract energy values from vienna file to DB file.
        extract_energies(vienna_tmp, db_tmp, force=True)
        # Convert DB file to CT
        db_to_ct(db_tmp, force=True)
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


def extract_energies(vienna_input: Path, db_output: Path, force: bool = False):
    """ Extract the free energies from a vienna file and prepend them to the file name.

    The title will follow this format:

    ENERGY = {energy}  {reference}

    where {reference} is the name of the reference sequence and {energy}
    is the predicted free energy of folding.

    Parameters
    ----------
    vienna_input: Path
        Path of the vienna file to extract energies from.
    db_output: Path
        Path of the DB file to which to write the extracted information.
    force: bool = False
        Overwrite the output DB file if it already exists.
    """
    if need_write(db_output, force):
        # Read all lines from the input file.
        lines = list()
        with open(vienna_input) as f:
            while title_line := f.readline():
                seq_line = f.readline()
                struct_line = f.readline()
                if " " not in struct_line:
                    logger.error(("No energy value could be "
                                  f"parsed from the vienna line {struct_line}"))
                struct_line_parts = struct_line.split(" ")
                struct_line = struct_line_parts[0] + "\n"
                energy = struct_line_parts[-1].strip("\n").strip("()")
                title_line = f">ENERGY = {energy} {title_line.strip('>')}"
                lines.extend([title_line, seq_line, struct_line])
        # Write the reformatted lines to the output file.
        text = "".join(lines)
        with open(db_output, write_mode(force=True)) as f:
            f.write(text)
        logger.routine(f"Energies extracted from file {vienna_input}"
                       + f" to {db_output}")
