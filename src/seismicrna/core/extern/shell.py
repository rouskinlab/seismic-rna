import logging
import shlex
from functools import wraps
from pathlib import Path
from subprocess import CompletedProcess, run
from typing import Any, Callable

logger = logging.getLogger(__name__)

# Commands for external applications
BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
ECHO_CMD = "echo"
FASTQC_CMD = "fastqc"
RNASTRUCTURE_FOLD_CMD = "Fold"
RNASTRUCTURE_FOLD_SMP_CMD = "Fold-smp"
GUNZIP_CMD = "gunzip"
SAMTOOLS_CMD = "samtools"
WORD_COUNT_CMD = "wc"


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


def run_cmd(cmd: str, text: bool | None = True):
    """ Run a command via subprocess.run(), with logging. """
    # Log the command with which the process was run.
    logger.debug(f"Running command via the shell:\n{cmd}")
    # Run the process and capture the output.
    process = run(cmd,
                  shell=True,
                  capture_output=text is not None,
                  text=text)
    # Format a message depending on whether the process passed.
    passed = process.returncode == 0
    status = "PASSED" if passed else f"FAILED with code {process.returncode}"
    message = "\n".join([f"Shell command {status}:\n{cmd}\n",
                         f"STDOUT:\n{process.stdout}\n",
                         f"STDERR:\n{process.stderr}\n"])
    if not passed:
        raise RuntimeError(message)
    logger.debug(message)
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
