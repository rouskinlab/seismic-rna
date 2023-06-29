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
SAMTOOLS_CMD = "samtools"


# Command utility functions

def run_cmd(args: list[Any],
            check_is_before: Sequence[Path] = (),
            check_no_before: Sequence[Path] = (),
            check_is_after: Sequence[Path] = (),
            check_no_after: Sequence[Path] = (),
            check_created: Sequence[Path] = (),
            check_deleted: Sequence[Path] = ()):
    """ Run a command via subprocess.run(), with logging. """
    verify_str = ", ".join(map(str, check_created))
    # Use shlex to place quotes around arguments containing whitespace.
    cmd = shlex.join(map(str, args))
    # Check if any required input files are missing.
    if missing := list(filterfalse(Path.exists,
                                   chain(check_is_before, check_deleted))):
        raise FileNotFoundError(f"Missing input files: {missing}")
    # Check if any expected output files already exist.
    if exists := list(filter(Path.exists,
                             chain(check_no_before, check_created))):
        raise FileExistsError(f"Existing output files: {exists}")
    # Log the command with which the process was run.
    logger.debug(f"Shell $ {cmd}")
    # Run the process and capture the output.
    process = subprocess.run(cmd, check=True, shell=True, capture_output=True)
    # Log the output of the process.
    log_process(process)
    # Check if any expected output files are missing.
    if missing := list(filterfalse(Path.exists,
                                   chain(check_is_after, check_created))):
        raise FileNotFoundError(f"Missing output files: {missing}")
    # Check if any expected deleted files still exist.
    if exists := list(filter(Path.exists,
                             chain(check_no_after, check_deleted))):
        raise FileExistsError(f"Existing input files: {exists}")
    return process


def log_process(process: subprocess.CompletedProcess):
    """ Log the output and error messages of a process. """
    if process.stdout:
        logger.info(f"STDOUT of {process.args}:\n{process.stdout.decode()}")
    if process.stderr:
        logger.info(f"STDERR of {process.args}:\n{process.stderr.decode()}")
