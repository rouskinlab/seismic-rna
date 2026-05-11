""" Wrapper around ViennaRNA from Lorenz and Hofacker at the
University of Vienna: https://www.tbi.univie.ac.at/RNA/
"""

import os
from pathlib import Path

from .dryrun import dry_run
from .profile import RNAFoldProfile
from .rnastructure import retitle_ct
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


def make_rnafold_cmd(fasta_file: Path,
                     vienna_file: Path, *,
                     sp_data: Path | None,
                     sp_strategy: str | None,
                     fold_constraint: Path | None,
                     fold_commands: Path | None,
                     fold_temp_c: float,
                     fold_isolated: bool,
                     fold_md: int,
                     fold_max: int,
                     fold_mfe: bool,
                     num_cpus: int = 1):
    """ Build the shell command to run RNAfold (and optionally RNAsubopt).

    Parameters
    ----------
    fasta_file: Path
        Input FASTA file containing the RNA sequence.
    vienna_file: Path
        Path prefix for the output vienna file (suffix is determined
        automatically based on ``fold_mfe``).
    sp_data: Path or None
        File of per-position reactivity data for soft-constraints; None
        disables soft constraints.
    sp_strategy: str or None
        Soft-constraint strategy passed to ``--sp-strategy``; None
        omits the flag.
    fold_constraint: Path or None
        Hard-constraint file passed to ``--constraint``; None omits the
        flag.
    fold_commands: Path or None
        Commands file passed to ``--commands``; None omits the flag.
    fold_temp_c: float
        Folding temperature in degrees Celsius.
    fold_isolated: bool
        If True, allow isolated base pairs; if False, pass ``--noLP``.
    fold_md: int
        Maximum base-pair span in nucleotides; 0 disables the limit.
    fold_max: int
        Maximum number of structures to keep (passed to ``head``).
    fold_mfe: bool
        If True, run only RNAfold (MFE structure); if False, also run
        RNAsubopt for suboptimal structures.
    num_cpus: int, optional
        Number of threads for RNAfold (default 1).

    Returns
    -------
    str
        A shell command string ready to be executed.
    """
    os.environ[RNAFOLD_NUM_THREADS] = str(num_cpus)
    final_cmds = list()
    cmds = [[VIENNA_RNAFOLD_CMD]]
    if not fold_mfe:
        cmds.append([VIENNA_RNASUBOPT_CMD])
    for cmd in cmds:
        if VIENNA_RNAFOLD_CMD in cmd:
            cmd.append("--noPS")
        else:
            cmd.extend(["--sorted", "--en-only"])
        if sp_data is not None:
            # File of reactivities.
            cmd.extend(["--sp-data", sp_data])
            if sp_strategy is not None:
                cmd.extend(["--sp-strategy", sp_strategy])
        if fold_constraint is not None:
            # File of constraints.
            cmd.extend(["--constraint", fold_constraint])
        if fold_commands is not None:
            # File of commands.
            cmd.extend(["--commands", fold_commands])
        # Temperature of folding (Celsius).
        cmd.extend(["--temp", fold_temp_c])
        if not fold_isolated:
            # Forbid isolated pairs.
            cmd.append("--noLP")
        if fold_md > 0:
            # Maximum distance between paired bases.
            cmd.extend(["--maxBPspan", fold_md])
        # Input and output files.
        cmd = args_to_cmd(cmd)
        cmd = cmds_to_redirect_in([cmd, str(fasta_file)])
        vienna_suffix = VIENNA_SUBOPT_EXT if VIENNA_RNASUBOPT_CMD in cmd else VIENNA_EXT
        cmd = cmds_to_pipe([cmd, f"head -n {fold_max + 1}"])
        cmd = cmds_to_redirect_out([cmd, str(vienna_file.with_suffix(vienna_suffix))])
        final_cmds.append(cmd)
    final_cmd = cmds_to_series(final_cmds)
    return final_cmd


def rnafold(rna: RNAFoldProfile, *,
            branch: str,
            out_dir: Path,
            tmp_dir: Path,
            keep_tmp: bool,
            fold_dry_run: bool,
            fold_mfe: bool,
            **kwargs):
    """ Run the 'RNAFold' or 'RNAsubopt' program of ViennaRNA. """
    logger.routine(f"Began folding {rna}")
    # Output CT file.
    ct_out = rna.get_ct_file(out_dir, branch)
    # Temporary FASTA file for the RNA.
    fasta_tmp = rna.write_fasta(tmp_dir, branch)
    # Path of the temporary CT file.
    ct_tmp = rna.get_ct_file(tmp_dir, branch)
    # Path of the temporary DB file.
    db_tmp = rna.get_db_file(tmp_dir, branch)
    # Path of the temporary vienna file.
    vienna_tmp = rna.get_vienna_file(tmp_dir, branch)
    # DMS reactivities file for the RNA.
    mus_file = rna.write_mus_file(tmp_dir, branch)
    try:
        # Generate the command.
        fold_cmd = make_rnafold_cmd(fasta_tmp,
                                    vienna_tmp,
                                    sp_data=mus_file,
                                    sp_strategy=rna.rnafold_sp_strategy,
                                    fold_temp_c=rna.fold_temp_c,
                                    fold_mfe=fold_mfe,
                                    **kwargs)
        if fold_dry_run:
            dry_run([fold_cmd], ct_tmp.parent)
        else:
            # Run the command.
            run_cmd(fold_cmd)
            # Extract energy values from vienna file to DB file.
            extract_energies(vienna_tmp, db_tmp, force=True)
            if not fold_mfe:
                get_subopt(vienna_tmp.with_suffix(VIENNA_SUBOPT_EXT), db_tmp)
            # Convert DB file to CT
            db_to_ct(db_tmp, force=True)
            # Reformat the CT file title lines so that each is unique.
            retitle_ct(ct_tmp, ct_tmp, force=True)
            # Renumber the CT file so that it has the same numbering
            # scheme as the region, rather than always starting at 1,
            # the latter of which is always output by the Fold program.
            renumber_ct(ct_tmp, ct_out, rna.region.end5, force=True)
    finally:
        if not keep_tmp:
            # Delete the temporary files.
            fasta_tmp.unlink(missing_ok=True)
            if ct_tmp != ct_out:
                ct_tmp.unlink(missing_ok=True)
            db_tmp.unlink(missing_ok=True)
            vienna_tmp.unlink(missing_ok=True)
            mus_file.unlink(missing_ok=True)
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
                f.readline()
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
