import logging
import shlex
import subprocess
from pathlib import Path
from typing import Any, Callable

logger = logging.getLogger(__name__)

# Commands for external applications
BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
FASTQC_CMD = "fastqc"
GUNZIP_CMD = "gunzip"
RNASTRUCTURE_FOLD_CMD = "Fold"
SAMTOOLS_CMD = "samtools"
WORD_COUNT_CMD = "wc"
WHICH_CMD = "which"


def log_process(process: subprocess.CompletedProcess):
    """ Log the output and error messages of a process. """
    if process.stdout:
        logger.info(f"STDOUT of {process.args}:\n{process.stdout.decode()}")
    if process.stderr:
        logger.info(f"STDERR of {process.args}:\n{process.stderr.decode()}")


def make_cmd(args: list):
    """ Join a list of arguments into a command with shlex. """
    return shlex.join(map(str, args))


def join_cmds(args: list[list], sep: str = " ; "):
    """ Given a list of commands, each of which is a list of arguments,
    join each list of arguments into a command, then join the commands
    into a single string with the `sep` string between each command. """
    return sep.join(map(make_cmd, args))


def make_pipeline_cmd(args: list[list]):
    return join_cmds(args, " | ")


def run_cmd(args: list, cmd_func: Callable[[list], str] = make_cmd):
    """ Run a command via subprocess.run(), with logging. """
    # Use shlex to place quotes around arguments containing whitespace.
    cmd = cmd_func(args)
    # Log the command with which the process was run.
    logger.debug(f"Shell $ {cmd}")
    # Run the process and capture the output.
    process = subprocess.run(cmd, check=True, shell=True, capture_output=True)
    # Log the output of the process.
    log_process(process)
    return process


class PipelineStep(object):
    """ Run a step that accepts input and produces output. """

    def __init__(self,
                 get_args_func: Callable[[Any, Any, Any], list],
                 make_cmd_func: Callable[[Any], str],
                 action: str = "processing"):
        self._get_args_func = get_args_func
        self._make_cmd_func = make_cmd_func
        self._action = action

    def __call__(self,
                 in_file: Path | Any | None,
                 out_file: Path | Any | None,
                 *args, **kwargs):
        action = self._action
        if in_file:
            action = f"{action} {in_file}"
        if out_file:
            action = f"{action} to {out_file}"
            # Make the parent directory of the output file, if it does
            # not already exist.
            out_file.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Began {action}")
        # Generate and run the command.
        process = run_cmd(self._get_args_func(in_file,
                                              out_file,
                                              *args, **kwargs),
                          self._make_cmd_func)
        logger.info(f"Ended {action}")
        return process


class ParsedPipelineStep(PipelineStep):
    """ Like PipelineStep, but instead of returning a CompletedProcess
    when called, parse its stdout/stderr with the given parser function
    and return the result. """
    
    def __init__(self,
                 get_args_func: Callable[[Any, Any, Any], list],
                 make_cmd_func: Callable[[Any], str],
                 parser_func: Callable[[subprocess.CompletedProcess], Any],
                 action: str = "processing"):
        super().__init__(get_args_func, make_cmd_func, action)
        self._parser = parser_func

    def __call__(self, *args, **kwargs):
        return self._parser(super().__call__(*args, **kwargs))
