"""
Struct -- RNAstructure Module

Wrapper around RNAstructure from the Mathews Lab at U of Rochester:
https://rna.urmc.rochester.edu/RNAstructure.html
"""

from logging import getLogger
from pathlib import Path

from ..core import path
from ..core.rna import RnaProfile
from ..core.shell import run_cmd, RNASTRUCTURE_FOLD_CMD

logger = getLogger(__name__)


def fold(rna: RnaProfile, *,
         out_dir: Path, temp_dir: Path, save_temp: bool, rerun: bool):
    """ Run the 'Fold' program of RNAstructure. """
    ct_file = rna.ct_file(out_dir)
    if rerun or not ct_file.is_file():
        cmd = [RNASTRUCTURE_FOLD_CMD]
        # Write the DMS reactivities file for the RNA.
        cmd.extend(["--DMS", rna.to_dms(temp_dir)])
        # Write a temporary FASTA file for the RNA.
        cmd.append(fasta := rna.to_fasta(temp_dir))
        try:
            # Get the path of the output CT file.
            cmd.append(ct_file)
            # Run the command.
            run_cmd(cmd, check_is_after=[ct_file])
        finally:
            if not save_temp:
                # Delete the temporary files.
                fasta.unlink(missing_ok=True)
    else:
        logger.warning(f"File exists: {ct_file}")
    return ct_file


def ct2dot(ct_file: Path, number: int | str = "all"):
    """ Convert a CT file to a DOT file. """
    dot_file = ct_file.with_suffix(path.DOT_EXT)
    cmd = ["ct2dot", ct_file, number, dot_file]
    run_cmd(cmd, check_is_before=[ct_file], check_is_after=[dot_file])
    return dot_file


def dot2ct(dot_file: Path):
    """ Convert a DOT file to a CT file. """
    ct_file = dot_file.with_suffix(path.CT_EXT)
    cmd = ["dot2ct", dot_file, ct_file]
    run_cmd(cmd, check_is_before=[dot_file], check_is_after=[ct_file])
    return ct_file
