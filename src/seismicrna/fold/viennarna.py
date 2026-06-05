"""Wrapper around ViennaRNA from Lorenz and Hofacker at the
University of Vienna: https://www.tbi.univie.ac.at/RNA/
"""

import os
from pathlib import Path

from .dryrun import dry_run
from .rnastructure import retitle_ct
from ..core.extern import (
    VIENNA_RNAFOLD_CMD,
    VIENNA_RNASUBOPT_CMD,
    cmds_to_pipe,
    cmds_to_redirect_in,
    cmds_to_redirect_out,
    args_to_cmd,
    run_cmd,
)
from ..core.logs import logger
from ..core.rna import renumber_ct, db_to_ct
from ..core.write import need_write, write_mode

ENERGY_UNIT = "kcal/mol"
RNAFOLD_NUM_THREADS = "OMP_NUM_THREADS"


def make_rnafold_cmd(
    fasta_file: Path,
    vienna_file: Path,
    *,
    sp_data: Path | None,
    sp_strategy: str | None,
    eddy_prior_paired_file: Path | None,
    eddy_prior_unpaired_file: Path | None,
    fold_constraint: Path | None,
    fold_commands: Path | None,
    fold_temp_c: float,
    fold_isolated: bool,
    fold_md: int,
    fold_max: int,
    fold_mfe: bool,
    fold_edelta: float,
    num_cpus: int = 1,
):
    """Build the shell command to run either RNAfold or RNAsubopt.

    Parameters
    ----------
    fasta_file: Path
        Input FASTA file containing the RNA sequence.
    vienna_file: Path
        Output path for the intermediate vienna file.
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
        If True, run only RNAfold (MFE structure); if False, run
        RNAsubopt for suboptimal structures.
    fold_edelta: float
        Maximum absolute energy difference (kcal/mol) from the MFE for
        suboptimal structures; passed to RNAsubopt as ``-e``.
    num_cpus: int, optional
        Number of threads for RNAfold (default 1).

    Returns
    -------
    str
        A shell command string ready to be executed.
    """
    os.environ[RNAFOLD_NUM_THREADS] = str(num_cpus)
    if fold_mfe:
        # Use RNAfold to get the MFE structure.
        args = [VIENNA_RNAFOLD_CMD, "--noPS"]
        # RNAfold: 3 lines per structure (title + sequence + structure).
        head_n = fold_max * 3
    else:
        # Use RNAsubopt to get the MFE and suboptimal structures.
        args = [
            VIENNA_RNASUBOPT_CMD,
            "--deltaEnergy",
            fold_edelta,
            "--sorted",
            "--en-only",
        ]
        # RNAsubopt: 2 header lines (title + sequence) + fold_max structure lines.
        head_n = fold_max + 2
    if sp_data is not None:
        # File of reactivities.
        args.extend(["--sp-data", sp_data])
        if sp_strategy is not None:
            args.extend(["--sp-strategy", sp_strategy])
    if eddy_prior_paired_file is not None:
        args.extend(["--sp-data", eddy_prior_paired_file, "--sp-strategy", "Pp"])
    if eddy_prior_unpaired_file is not None:
        args.extend(["--sp-data", eddy_prior_unpaired_file, "--sp-strategy", "Pu"])
    if fold_constraint is not None:
        # File of constraints.
        args.extend(["--constraint", fold_constraint])
    if fold_commands is not None:
        # File of commands.
        args.extend(["--commands", fold_commands])
    # Temperature of folding (Celsius).
    args.extend(["--temp", fold_temp_c])
    if not fold_isolated:
        # Forbid isolated pairs.
        args.append("--noLP")
    if fold_md > 0:
        # Maximum distance between paired bases.
        args.extend(["--maxBPspan", fold_md])
    # Input and output files.
    cmd = args_to_cmd(args)
    cmd = cmds_to_redirect_in([cmd, str(fasta_file)])
    cmd = cmds_to_pipe([cmd, f"head -n {head_n}"])
    cmd = cmds_to_redirect_out([cmd, str(vienna_file)])
    return cmd


def run_rnafold(
    fasta_tmp: Path,
    ct_tmp: Path,
    ct_out: Path,
    vienna_tmp: Path,
    db_tmp: Path,
    *,
    sp_data: Path | None,
    sp_strategy: str | None,
    eddy_prior_paired_file: Path | None,
    eddy_prior_unpaired_file: Path | None,
    fold_constraint: Path | None,
    fold_commands: Path | None,
    fold_temp_c: float,
    fold_isolated: bool,
    fold_md: int,
    fold_max: int,
    fold_mfe: bool,
    fold_edelta: float,
    end5: int,
    num_cpus: int,
    fold_dry_run: bool = False,
):
    """Run RNAfold/RNAsubopt on pre-built paths, convert to CT, retitle, and
    renumber."""
    fold_cmd = make_rnafold_cmd(
        fasta_tmp,
        vienna_tmp,
        sp_data=sp_data,
        sp_strategy=sp_strategy,
        eddy_prior_paired_file=eddy_prior_paired_file,
        eddy_prior_unpaired_file=eddy_prior_unpaired_file,
        fold_constraint=fold_constraint,
        fold_commands=fold_commands,
        fold_temp_c=fold_temp_c,
        fold_isolated=fold_isolated,
        fold_md=fold_md,
        fold_max=fold_max,
        fold_mfe=fold_mfe,
        fold_edelta=fold_edelta,
        num_cpus=num_cpus,
    )
    if fold_dry_run:
        dry_run([fold_cmd], ct_tmp.parent)
    else:
        run_cmd(fold_cmd)
        if fold_mfe:
            extract_energies(vienna_tmp, db_tmp, force=True)
        else:
            get_subopt(vienna_tmp, db_tmp, force=True)
        db_to_ct(db_tmp, force=True)
        retitle_ct(ct_tmp, ct_tmp, force=True)
        renumber_ct(ct_tmp, ct_out, end5, force=True)


def get_subopt(subopt_out: Path, db_target: Path, force: bool = False):
    """Extract suboptimal structures from the output of RNAsubopt and write them
    to a DB file.

    RNAsubopt output format (per sequence):

    .. code-block:: text

        >NAME [N]          <- title; N = energy window in 0.01 kcal/mol units
        SEQ  MIN_E  MAX_E  <- sequence followed by energy-range values
        STRUCT  ENERGY     <- one line per suboptimal structure
        ...

    All fields on lines 2+ are whitespace-separated (3 spaces in practice).
    """
    if need_write(db_target, force):
        lines = list()
        seq_title = None
        seq = None
        first_struct = True
        with open(subopt_out) as f:
            while line := f.readline():
                if line.startswith(">"):
                    # ">seqname [N]" → keep only ">seqname"
                    seq_title = line.split()[0]
                    # "SEQUENCE   MIN_E   MAX_E" → keep only the sequence
                    seq = f.readline().split()[0]
                    first_struct = True
                    continue
                parts = line.split()
                if len(parts) != 2:
                    continue
                struct, energy = parts
                title_line = f">ENERGY = {energy} {seq_title.strip('>')}\n"
                if first_struct:
                    lines.extend([title_line, seq + "\n", struct + "\n"])
                    first_struct = False
                else:
                    lines.extend([title_line, struct + "\n"])
        text = "".join(lines)
        with open(db_target, mode=write_mode(force=True)) as f:
            f.write(text)
        logger.debug(
            "Suboptimal structures from {} written to {}", subopt_out, db_target
        )


def extract_energies(vienna_input: Path, db_output: Path, force: bool = False):
    """Extract the free energies from a vienna file and prepend them to the reference
    name.

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
                    logger.error(
                        "No energy value could be parsed from the vienna line {}",
                        struct_line,
                    )
                struct_line_parts = struct_line.split(" ")
                struct_line = struct_line_parts[0] + "\n"
                energy = struct_line_parts[-1].strip("\n").strip("()")
                title_line = f">ENERGY = {energy} {title_line.strip('>')}"
                lines.extend([title_line, seq_line, struct_line])
        # Write the reformatted lines to the output file.
        text = "".join(lines)
        with open(db_output, write_mode(force=True)) as f:
            f.write(text)
        logger.debug(
            "Energies extracted from file {} to {}", vienna_input, db_output
        )
