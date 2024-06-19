from concurrent.futures import Future, ProcessPoolExecutor
from itertools import chain, filterfalse, repeat
from logging import getLogger
from typing import Any, Callable, Iterable

from .logs import exc_info, get_config, set_config

logger = getLogger(__name__)


def get_num_parallel(n_tasks: int,
                     max_procs: int,
                     parallel: bool,
                     hybrid: bool = False) -> tuple[int, int]:
    """
    Determine how to parallelize the tasks.

    Parameters
    ----------
    n_tasks: int
        Number of tasks to parallelize. Must be ≥ 1.
    max_procs: int
        Maximum number of processes to run at one time. Must be ≥ 1.
    parallel: bool
        Whether multiple tasks may be run in parallel. If False, then
        the number of tasks to run in parallel is set to 1, but the
        number of processes to run for each task may be > 1.
    hybrid: bool = False
        Whether to allow both multiple tasks to run in parallel and,
        at the same, each task to run multiple processes in parallel.

    Returns
    -------
    tuple[int, int]
        - Number of tasks to run in parallel. Always ≥ 1.
        - Number of processes to run for each task. Always ≥ 1.
    """
    if n_tasks >= 1 and max_procs >= 1:
        # This function only works if there is at least one task to
        # parallelize, at least one process is allowed, and parallel
        # is a valid option.
        if parallel:
            # Multiple tasks may be run in parallel. The number of tasks
            # run in parallel cannot exceed 1) the total number of tasks
            # and 2) the user-specified maximum number of processes.
            n_tasks_parallel = min(n_tasks, max_procs)
        else:
            # Otherwise, only one task at a time can be run.
            n_tasks_parallel = 1
        if n_tasks_parallel == 1 or hybrid:
            # Each individual task can be run by multiple processes in
            # parallel, as long as either 1) multiple tasks are not run
            # simultaneously in parallel (i.e. n_tasks_parallel == 1)
            # or 2) the calling function sets hybrid=True, which lets
            # multiple tasks run in parallel and each run with multiple
            # processes. Only the alignment module can simultaneously
            # run multiple tasks and multiple processes for each task
            # because its two most computation-heavy processes (cutadapt
            # and bowtie2) come with their own parallelization abilities
            # that can work independently of Python's multiprocessing
            # module. However, the other modules (e.g. vectoring) are
            # parallelized using the multiprocessing module, which does
            # not support "nesting" parallelization in multiple layers.
            # Because n_tasks_parallel is either 1 or the smallest of
            # n_tasks and n_procs (both of which are ≥ 1), it must be
            # that 1 ≤ n_tasks_parallel ≤ n_procs, and therefore that
            # 1 ≤ n_procs / n_tasks_parallel ≤ n_procs, so the
            # integer quotient must be a valid number of processes.
            n_procs_per_task = max_procs // n_tasks_parallel
        else:
            # Otherwise, only one process can work on each task.
            n_procs_per_task = 1
    else:
        logger.warning("Defaulting to 1 process due to invalid number of "
                       f"tasks ({n_tasks}) and/or processes ({max_procs}).")
        n_tasks_parallel = 1
        n_procs_per_task = 1
    return n_tasks_parallel, n_procs_per_task


def fmt_func_args(func: Callable, *args, **kwargs):
    """ Format the name and arguments of a function as a string. """
    fargs = ", ".join(chain(map(repr, args),
                            (f"{kw}={repr(arg)}"
                             for kw, arg in kwargs.items())))
    return f"{func.__name__}({fargs})"


class Task(object):
    """ Wrap a parallelizable task in a try-except block so that if it
    fails, it just returns `None` rather than crashing the other tasks
    being run in parallel. """

    def __init__(self, func: Callable):
        self._func = func
        self._config = get_config()

    def __call__(self, *args, **kwargs):
        """ Call the task's function in a try-except block, return the
        result if it succeeds, and return None otherwise. """
        if get_config() != self._config:
            # Tasks running in parallel may not have the same top logger
            # as the parent process (this seems to be system dependent).
            # If not, then this task's top logger must be configured to
            # match the configuration of the parent process.
            set_config(*self._config)
        task = fmt_func_args(self._func, *args, **kwargs)
        try:
            logger.debug(f"Began task {task}")
            result = self._func(*args, **kwargs)
        except Exception as error:
            logger.error(f"Failed task {task}:\n{error}\n", exc_info=exc_info())
        else:
            logger.debug(f"Ended task {task}:\n{result}\n")
            return result


def dispatch(funcs: list[Callable] | Callable,
             max_procs: int, parallel: bool, *,
             hybrid: bool = False,
             pass_n_procs: bool = True,
             args: list[tuple] | tuple = (),
             kwargs: dict[str, Any] | None = None):
    """
    Run one or more tasks in series or in parallel, depending on the
    number of tasks, the maximum number of processes, and whether tasks
    are allowed to be run in parallel.

    Parameters
    ----------
    funcs: list[Callable] | Callable
        The function(s) to run. Can be a list of functions or a single
        function that is not in a list. If a single function, then if
        `args` is a tuple, it is called once with that tuple as its
        positional arguments; and if `args` is a list of tuples, it is
        called for each tuple of positional arguments in `args`.
    max_procs: int
        See docstring for `get_num_parallel`.
    parallel: bool
        See docstring for `get_num_parallel`.
    hybrid: bool = False
        See docstring for `get_num_parallel`.
    pass_n_procs: bool = True
        Whether to pass the number of processes to the function as the
        keyword argument `n_procs`.
    args: list[tuple] | tuple = ()
        Positional arguments to pass to each function in `funcs`. Can be
        a list of tuples of positional arguments or a single tuple that
        is not in a list. If a single tuple, then each function receives
        `args` as positional arguments. If a list, then `args` must be
        the same length as `funcs`; each function `funcs[i]` receives
        `args[i]` as positional arguments.
    kwargs: dict[str, Any] | None = None
        Keyword arguments to pass to every function call.

    Returns
    -------
    list
        List of the return value of each run.
    """
    # Default to an empty dict if kwargs is not given.
    if kwargs is None:
        kwargs = dict()
    if callable(funcs):
        if isinstance(args, tuple):
            # If args is a tuple, make it the sole element of a list.
            args = [args]
        else:
            # Ensure that every item in args is a tuple.
            if nontuple := list(filterfalse(lambda x: isinstance(x, tuple),
                                            args)):
                raise TypeError(f"Got non-tuple args: {nontuple}")
        # If a function is given rather than an iterable of functions,
        # then put the function in a list whose length equal that of the
        # list of arguments.
        funcs = list(repeat(funcs, len(args)))
    else:
        # Ensure that every item in funcs is actually callable.
        if uncallable := list(filterfalse(callable, funcs)):
            raise TypeError(f"Got uncallable funcs: {uncallable}")
        if isinstance(args, tuple):
            # If args is a tuple, repeat it once for each function.
            args = list(repeat(args, len(funcs)))
    # Ensure that numbers of functions and argument tuples match.
    if (n_tasks := len(funcs)) != len(args):
        raise ValueError(f"Got {len(funcs)} funcs but {len(args)} args")
    if n_tasks == 0:
        # No tasks to run: return.
        logger.warning("No tasks were given to dispatch")
        return list()
    # Determine how to parallelize each task.
    n_tasks_parallel, n_procs_per_task = get_num_parallel(n_tasks,
                                                          max_procs,
                                                          parallel,
                                                          hybrid=hybrid)
    if pass_n_procs:
        # Add the number of processes as a keyword argument.
        kwargs = {**kwargs, "n_procs": n_procs_per_task}
    if n_tasks_parallel > 1:
        # Run the tasks in parallel.
        with ProcessPoolExecutor(max_workers=n_tasks_parallel) as pool:
            logger.info(f"Opened pool of {n_tasks_parallel} processes")
            # Initialize an empty list of tasks to run.
            tasks: list[Future] = list()
            for func, task_args in zip(funcs, args, strict=True):
                # Create a new task and submit it to the process pool.
                task = Task(func)
                tasks.append(pool.submit(task, *task_args, **kwargs))
            # Run all the tasks in parallel and collect the results as
            # they become available.
            logger.info(f"Waiting for {n_tasks} tasks to finish")
            results = [task.result() for task in tasks]
        logger.info(f"Closed pool of {n_tasks_parallel} processes")
    else:
        # Run the tasks in series.
        logger.info(f"Began running {n_tasks} task(s) in series")
        # Initialize an empty list of results from the tasks.
        results = list()
        for func, task_args in zip(funcs, args, strict=True):
            # Create a new task, run it in the current process, and add
            # its result to the list of results.
            task = Task(func)
            results.append(task(*task_args, **kwargs))
        logger.info(f"Ended running {n_tasks} task(s) in series")
    # Remove any failed runs (None values) from results.
    results = [result for result in results if result is not None]
    n_pass = len(results)
    n_fail = n_tasks - n_pass
    if n_fail:
        p_fail = n_fail / n_tasks * 100.
        logger.warning(
            f"Failed {n_fail} of {n_tasks} task(s) ({round(p_fail, 1)} %)"
        )
    else:
        logger.info(f"All {n_tasks} task(s) completed successfully")
    return results


def as_list_of_tuples(args: Iterable[Any]):
    """ Given an iterable of arguments, return a list of 1-item tuples,
    each containing one of the given arguments. This function is useful
    for creating a list of tuples to pass to the `args` parameter of
    `dispatch`. """
    return [(arg,) for arg in args]

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
