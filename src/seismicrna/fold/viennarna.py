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
                           cmds_to_pipe,
                           cmds_to_series,
                           cmds_to_redirect_in,
                           cmds_to_redirect_out,
                           args_to_cmd,
                           run_cmd)
from ..core.logs import logger
from ..core.path import VIENNA_EXT, VIENNA_SUBOPT_EXT
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
    os.environ[RNAFOLD_NUM_THREADS] = str(num_cpus)
    final_cmds = list()
    cmds = [[VIENNA_RNAFOLD_CMD]]
    if not fold_mfe:
        cmds.append([VIENNA_RNASUBOPT_CMD])
    for cmd in cmds:
        cmd.append("--noLP")
        if VIENNA_RNAFOLD_CMD in cmd:
            cmd.append("--noPS")
        else:
            cmd.extend(["--sorted", "--en-only"])
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
        vienna_suffix = VIENNA_SUBOPT_EXT if VIENNA_RNASUBOPT_CMD in cmd else VIENNA_EXT
        num_structs = 5 #TODO: Fixme
        cmd = cmds_to_pipe([cmd, f"head -n {num_structs+1}"])
        cmd = cmds_to_redirect_out([cmd, str(vienna_file.with_suffix(vienna_suffix))])
        final_cmds.append(cmd)
    final_cmd = cmds_to_series(final_cmds)
    return final_cmd


def calc_bp_pseudoenergy(seq_len: int,
                         pseudoenergies: pd.Series,
                         out_file: Path):
    """
    Identify (i, j) pairs, with j >= i + 3, where the sum of pseudoenergies is non-zero,
    and write them in ViennaRNA command format:

      "E i j 1 <bp_pseudoenergy>"

    Parameters:
      seq_len: The total length of the sequence. Positions are assumed to be numbered 1..seq_len.
      pseudoenergies: A pandas Series of pseudoenergies, indexed with multiindex (Position, Base).
      out_file: Path to the file where the results will be written.
    """
    # Extract Position level and drop NaNs
    pos_idx = pseudoenergies.index.get_level_values("Position").astype(int)
    energy_vals = pseudoenergies.values
    mask = ~pd.isna(energy_vals)

    # Full-length energy array, default zero
    energies_full = np.zeros(seq_len, dtype=np.float64)
    energies_full[pos_idx[mask] - 1] = energy_vals[mask]

    # Collect valid (i, j) pairs
    parts = []
    for offset in range(3, seq_len):
        i_vals = np.arange(1, seq_len - offset + 1)
        j_vals = i_vals + offset
        sums = energies_full[i_vals - 1] + energies_full[j_vals - 1]
        valid = np.abs(sums) > 0
        if valid.any():
            block = np.column_stack([
                i_vals[valid],
                j_vals[valid],
                np.ones(valid.sum(), dtype=int),
                sums[valid]
            ])
            parts.append(block)

    if parts:
        result = np.vstack(parts)
        with open(out_file, 'w') as f:
            np.savetxt(f, result, fmt="E %d %d %d %.6f")

    return out_file


def append_or_copy(source: Path, dest: Path):
    if dest.exists():
        with source.open('rb') as fsrc, dest.open('ab') as fdst:
            shutil.copyfileobj(fsrc, fdst)
    else:
        shutil.copy2(source, dest)


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

    apply_all_paired = True

    if apply_all_paired:
        command_file = calc_bp_pseudoenergy(len(rna.seq), rna.pseudoenergies, command_tmp)
        probe_file = None
    else:
        command_file = None
        probe_file = dms_file

    if fold_commands:
        append_or_copy(fold_commands, command_tmp)

    try:
        # Run the command.
        fold_cmd = make_fold_cmd(fasta_tmp,
                                 vienna_tmp,
                                 fold_delta_e=5,  # TODO set default
                                 dms_file=probe_file,
                                 shape_method=rna.shape_method,
                                 fold_constraint=fold_constraint,
                                 fold_commands=command_file,
                                 fold_temp=rna.fold_temp,
                                 fold_md=fold_md,
                                 fold_mfe=fold_mfe,
                                 fold_max=fold_max,
                                 fold_percent=fold_percent,
                                 num_cpus=num_cpus)
        run_cmd(fold_cmd)
        # Extract energy values from vienna file to DB file.
        extract_energies(vienna_tmp, db_tmp, force=True)
        if not fold_mfe:
            get_subopt(vienna_tmp.with_suffix(VIENNA_SUBOPT_EXT), db_tmp)
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


def get_subopt(subopt_out: Path,
               db_target: Path):
    """ Extract suboptimal structures from the output of RNAsubopt and add them to a vienna file. """
    lines = list()
    with open(subopt_out) as f:
        while line := f.readline():
            if line.startswith(">"):
                seq_title = line.split(" ")[0]
                seq_line = f.readline()
                continue
            assert line.count(" ") == 1
            struct, energy = line.split(" ")
            struct = struct + "\n"
            energy = energy.strip("\n").strip("()")
            title_line = f">ENERGY = {energy} {seq_title.strip('>')}\n"
            lines.extend([title_line, struct])
    text = "".join(lines)
    with open(db_target, mode="a") as f:
            f.write(text)
    logger.routine(f"Suboptimal structures from {subopt_out} written"
                       + f" to {db_target}")


def extract_energies(vienna_input: Path, db_output: Path, force: bool = False):
    """ Extract the free energies from a vienna file and prepend them to the reference name.

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
