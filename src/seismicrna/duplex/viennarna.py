"""Wrapper around RNAcofold from the ViennaRNA package by Lorenz and
Hofacker at the University of Vienna: https://www.tbi.univie.ac.at/RNA/

RNAcofold predicts the secondary structure of a duplex of two RNA
strands. Its input and output separate the two strands with an
ampersand (``&``); this module removes the ampersand so that the duplex
is represented as a single, continuously numbered sequence, exactly like
any structure produced by ``seismic fold``.
"""

from pathlib import Path

from .dataset import LINKER, LINKER_LENGTH
from .profile import STRAND_SEP
from ..fold.dryrun import dry_run
from ..fold.rnastructure import retitle_ct
from ..core.shell import (
    VIENNA_RNACOFOLD_CMD,
    cmds_to_redirect_in,
    cmds_to_redirect_out,
    args_to_cmd,
    run_cmd,
)
from ..core.logs import logger
from ..core.rna.io import renumber_ct, db_to_ct
from ..core.write import need_write, write_mode


def make_rnacofold_cmd(
    fasta_file: Path,
    vienna_file: Path,
    *,
    shape_file: Path | None,
    shape_method: str | None,
    fold_constraint: Path | None,
    fold_temp_c: float,
    fold_isolated: bool,
    fold_md: int,
):
    """Build the shell command to run RNAcofold.

    RNAcofold offers only the older SHAPE interface (``--shape`` /
    ``--shapeMethod``), not the newer ``--sp-data`` / ``--sp-strategy``
    interface of RNAfold, so reactivities are passed accordingly.

    Parameters
    ----------
    fasta_file: Path
        Input FASTA file containing the two strands separated by ``&``.
    vienna_file: Path
        Output path for the intermediate vienna file.
    shape_file: Path or None
        File of per-position reactivity data passed to ``--shape``; None
        disables soft constraints.
    shape_method: str or None
        SHAPE incorporation method passed to ``--shapeMethod`` (e.g.
        ``Dm1.8b-0.6`` for Deigan); None omits the flag.
    fold_constraint: Path or None
        Hard-constraint file passed to ``--constraint``; None omits the
        flag.
    fold_temp_c: float
        Folding temperature in degrees Celsius.
    fold_isolated: bool
        If True, allow isolated base pairs; if False, pass ``--noLP``.
    fold_md: int
        Maximum base-pair span in nucleotides; 0 disables the limit.

    Returns
    -------
    str
        A shell command string ready to be executed.
    """
    args = [VIENNA_RNACOFOLD_CMD, "--noPS"]
    if shape_file is not None:
        # File of reactivities.
        args.extend(["--shape", shape_file])
        if shape_method is not None:
            args.extend(["--shapeMethod", shape_method])
    if fold_constraint is not None:
        # File of constraints.
        args.extend(["--constraint", fold_constraint])
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
    cmd = cmds_to_redirect_out([cmd, str(vienna_file)])
    return cmd


def run_rnacofold(
    fasta_tmp: Path,
    ct_tmp: Path,
    ct_out: Path,
    vienna_tmp: Path,
    db_tmp: Path,
    *,
    shape_file: Path | None,
    shape_method: str | None,
    fold_constraint: Path | None,
    fold_temp_c: float,
    fold_isolated: bool,
    fold_md: int,
    end5: int,
    fold_dry_run: bool = False,
):
    """Run RNAcofold on pre-built paths, convert to CT, retitle, and
    renumber."""
    fold_cmd = make_rnacofold_cmd(
        fasta_tmp,
        vienna_tmp,
        shape_file=shape_file,
        shape_method=shape_method,
        fold_constraint=fold_constraint,
        fold_temp_c=fold_temp_c,
        fold_isolated=fold_isolated,
        fold_md=fold_md,
    )
    if fold_dry_run:
        dry_run([fold_cmd], ct_tmp.parent)
    else:
        run_cmd(fold_cmd)
        extract_duplex(vienna_tmp, db_tmp, force=True)
        db_to_ct(db_tmp, force=True)
        retitle_ct(ct_tmp, ct_tmp, force=True)
        renumber_ct(ct_tmp, ct_out, end5, force=True)


def extract_duplex(vienna_input: Path, db_output: Path, force: bool = False):
    """Convert the output of RNAcofold to a dot-bracket (DB) file.

    RNAcofold prints, for each duplex, a title line, the two strands
    joined by ``&``, and the structure of the two strands joined by
    ``&`` followed by the free energy in parentheses:

    .. code-block:: text

        >NAME
        SEQ5&SEQ3
        STRUCT5&STRUCT3 (ENERGY)

    This function replaces the ``&`` with a short unpaired linker in
    both the sequence and the structure, so that the duplex becomes a
    single molecule in which each inter-strand base pair encloses the
    linker as a loop, and prepends the free energy to the title
    (matching ``seismic fold``):

    .. code-block:: text

        >ENERGY = {energy} {NAME}
        SEQ5{linker}SEQ3
        STRUCT5{....}STRUCT3

    Parameters
    ----------
    vienna_input: Path
        Path of the RNAcofold output file to convert.
    db_output: Path
        Path of the DB file to which to write the converted structure.
    force: bool = False
        Overwrite the output DB file if it already exists.
    """
    if need_write(db_output, force):
        lines = list()
        with open(vienna_input) as f:
            while title_line := f.readline():
                seq_line = f.readline()
                struct_line = f.readline()
                if " " not in struct_line:
                    logger.error(
                        "No energy value could be parsed from the RNAcofold line {}",
                        struct_line,
                    )
                # Split the structure from the trailing "(energy)" term.
                # The dot-bracket structure has no spaces, but RNAcofold
                # right-aligns the energy in its field (e.g. "( -0.40)"),
                # so split on the first space, not the last.
                struct, _, energy = struct_line.strip().partition(" ")
                energy = energy.strip().strip("()")
                # Replace the strand break with the linker, turning each
                # inter-strand base pair into a loop of unpaired bases.
                seq = seq_line.strip().replace(STRAND_SEP, LINKER)
                struct = struct.strip().replace(STRAND_SEP, "." * LINKER_LENGTH)
                title_line = f">ENERGY = {energy} {title_line.strip('>')}"
                lines.extend([title_line, f"{seq}\n", f"{struct}\n"])
        text = "".join(lines)
        with open(db_output, write_mode(force=True)) as f:
            f.write(text)
        logger.debug("Duplex structure from {} written to {}", vienna_input, db_output)
