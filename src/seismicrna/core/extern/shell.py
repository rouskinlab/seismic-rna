import shlex
from functools import wraps
from pathlib import Path
from subprocess import CompletedProcess, run
from typing import Any, Callable

from ..logs import logger

# Commands for external applications
ECHO_CMD = "echo"
WORD_COUNT_CMD = "wc"
GUNZIP_CMD = "gunzip"
FASTP_CMD = "fastp"
BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
SAMTOOLS_CMD = "samtools"
RNASTRUCTURE_FOLD_CMD = "Fold"
RNASTRUCTURE_FOLD_SMP_CMD = "Fold-smp"
JAVA_CMD = "java"
JAR_CMD = "-jar"
JGO_CMD = "jgo"


class ShellCommandFailedError(RuntimeError):
    """ A command failed that was run through the shell. """


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
    logger.action(f"Began running shell command:\n{cmd}")
    # Run the process and capture the output.
    process = run(cmd,
                  shell=True,
                  capture_output=text is not None,
                  text=text)
    failed = process.returncode != 0
    result = f"FAILED (code {process.returncode})" if failed else "SUCCEEDED"
    message = "\n".join([f"Shell command {result}:\n{cmd}\n",
                         f"STDOUT:\n{process.stdout}\n",
                         f"STDERR:\n{process.stderr}\n"])
    if failed:
        raise ShellCommandFailedError(message)
    logger.detail(message)
    logger.action("Ended running shell command")
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

    def __init__(self,
                 action: str,
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
        logger.routine(f"Began {action}")
        # Generate and run the command.
        process = run_cmd(self._make_command(ipath, opath, **kwargs))
        logger.routine(f"Ended {action}")
        if self._parse_output:
            logger.routine(f"Began parsing output of {action}")
            output = self._parse_output(process)
            logger.routine(f"Ended parsing output of {action}")
            return output
        return process
