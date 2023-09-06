import logging
import shlex
import subprocess
from functools import wraps
from pathlib import Path
from subprocess import CompletedProcess
from typing import Any, Callable

logger = logging.getLogger(__name__)

# Commands for external applications
BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
ECHO_CMD = "echo"
FASTQC_CMD = "fastqc"
RNASTRUCTURE_FOLD_CMD = "Fold"
GREP_CMD = "grep"
GUNZIP_CMD = "gunzip"
SAMTOOLS_CMD = "samtools"
WORD_COUNT_CMD = "wc"
WHICH_CMD = "which"


def args_to_cmd(args: list[Any]):
    """ Join a list of arguments into a command with shlex. """
    return shlex.join(map(str, args))


def cmds_to_pipe(cmds: list[str]):
    """ Join commands into a pipeline. """
    return " | ".join(cmds)


def cmds_to_series(cmds: list[str]):
    """ Run commands in series. """
    return " ; ".join(cmds)


def cmds_to_subshell(cmds: list[str]):
    """ Run commands in a subshell. """
    return f"( {cmds_to_series(cmds)} )"


def run_cmd(cmd: str):
    """ Run a command via subprocess.run(), with logging. """
    # Log the command with which the process was run.
    logger.debug(f"Running command via the shell:\n{cmd}")
    # Run the process and capture the output.
    process = subprocess.run(cmd, check=True, shell=True, capture_output=True)
    # Log the output of the process.
    if process.stdout:
        logger.debug(f"STDOUT of {process.args}:\n{process.stdout.decode()}")
    if process.stderr:
        logger.debug(f"STDERR of {process.args}:\n{process.stderr.decode()}")
    return process


def iopaths(has_ipath: bool = True, has_opath: bool = True):
    """ Given a function that takes an input path, output path, both, or
    neither, and potentially other positional/keyword arguments, return
    another function that takes input and output paths as its first two
    arguments and calls the original function with the right paths. """

    def decorator(func: Callable[[Any, Any, Any], str]):

        if has_ipath:
            if has_opath:
                def with_iopaths(ipath, opath, *args, **kwargs):
                    return func(ipath, opath, *args, **kwargs)

            else:
                def with_iopaths(ipath, __, *args, **kwargs):
                    return func(ipath, *args, **kwargs)
        else:
            if has_opath:
                def with_iopaths(_, opath, *args, **kwargs):
                    return func(opath, *args, **kwargs)

            else:
                def with_iopaths(_, __, *args, **kwargs):
                    return func(*args, **kwargs)

        return wraps(func)(with_iopaths)

    return decorator


class ShellCommand(object):
    """ Command that can be run in the shell. """

    def __init__(self, action: str,
                 mkcmd: Callable[[Any, Any, Any], str],
                 parse: Callable[[CompletedProcess], Any] | None = None,
                 ipath: bool = True,
                 opath: bool = True):
        self._make_command = iopaths(ipath, opath)(mkcmd)
        self._parse_output = parse
        self._action = action

    def _format_action(self, ipath, opath):
        action = self._action
        if ipath:
            action = f"{action} {ipath}"
        if opath:
            action = f"{action} to {opath}"
        return action

    def __call__(self,
                 ipath: Path | Any | None = None,
                 opath: Path | Any | None = None, /,
                 **kwargs):
        if opath:
            # Make the parent directory of the output file, if it does
            # not already exist.
            opath.parent.mkdir(parents=True, exist_ok=True)
        action = self._format_action(ipath, opath)
        logger.info(f"Began {action}")
        # Generate and run the command.
        process = run_cmd(self._make_command(ipath, opath, **kwargs))
        logger.info(f"Ended {action}")
        return self._parse_output(process) if self._parse_output else process
