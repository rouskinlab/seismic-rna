from itertools import chain, filterfalse
import logging
from pathlib import Path
import shlex
import subprocess
from typing import Any, Sequence

logger = logging.getLogger(__name__)

# Commands for external applications
BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
FASTQC_CMD = "fastqc"
RNASTRUCTURE_FOLD_CMD = "Fold"
SAMTOOLS_CMD = "samtools"
WHICH_CMD = "which"


# Command utility functions

def run_cmd(args: list[Any], *,
            check_is_before: Sequence[Path] = (),
            check_no_before: Sequence[Path] = (),
            check_is_after: Sequence[Path] = (),
            check_no_after: Sequence[Path] = (),
            check_created: Sequence[Path] = (),
            check_deleted: Sequence[Path] = ()):
    """ Run a command via subprocess.run(), with logging. """
    # Check created and deleted are shortcuts: expand them.
    if check_created:
        # Created files must exist after and not before.
        check_no_before = list(set(chain(check_no_before, check_created)))
        check_is_after = list(set(chain(check_is_after, check_created)))
    if check_deleted:
        # Deleted files must exist before and not after.
        check_is_before = list(set(chain(check_is_before, check_deleted)))
        check_no_after = list(set(chain(check_no_after, check_deleted)))
    # Use shlex to place quotes around arguments containing whitespace.
    cmd = shlex.join(map(str, args))
    # Check if any required input files are missing.
    if missing := list(filterfalse(Path.exists, check_is_before)):
        raise FileNotFoundError(f"Missing input files: {missing}")
    # Check if any expected output files already exist.
    if exists := list(filter(Path.exists, check_no_before)):
        raise FileExistsError(f"Existing output files: {exists}")
    # Log the command with which the process was run.
    logger.debug(f"Shell $ {cmd}")
    # Run the process and capture the output.
    process = subprocess.run(cmd, check=True, shell=True, capture_output=True)
    # Log the output of the process.
    log_process(process)
    # Check if any expected output files are missing.
    if missing := list(filterfalse(Path.exists, check_is_after)):
        raise FileNotFoundError(f"Missing output files: {missing}")
    # Check if any expected deleted files still exist.
    if exists := list(filter(Path.exists, check_no_after)):
        raise FileExistsError(f"Existing input files: {exists}")
    return process


def log_process(process: subprocess.CompletedProcess):
    """ Log the output and error messages of a process. """
    if process.stdout:
        logger.info(f"STDOUT of {process.args}:\n{process.stdout.decode()}")
    if process.stderr:
        logger.info(f"STDERR of {process.args}:\n{process.stderr.decode()}")
