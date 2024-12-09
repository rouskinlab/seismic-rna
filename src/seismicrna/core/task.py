from concurrent.futures import Future, ProcessPoolExecutor
from itertools import chain, filterfalse, repeat
from typing import Any, Callable, Iterable

from .logs import logger, get_config, set_config


def calc_pool_size(num_tasks: int, max_procs: int):
    """ Calculate the size of a process pool.

    Parameters
    ----------
    num_tasks: int
        Number of tasks to parallelize. Must be ≥ 1.
    max_procs: int
        Maximum number of processes to run at one time. Must be ≥ 1.

    Returns
    -------
    tuple[int, int]
        - Number of tasks to run in parallel. Always ≥ 1.
        - Number of processes to run for each task. Always ≥ 1.
    """
    logger.detail("Began calculating pool size")
    if max_procs < 1:
        logger.warning(f"max_procs must be ≥ 1, but got {max_procs}; "
                       f"defaulting to 1")
        max_procs = 1
    if num_tasks < 1:
        logger.warning(f"num_tasks must be ≥ 1, but got {num_tasks}; "
                       f"defaulting to 1")
        num_tasks = 1
    # The number of tasks that can run simultaneously is the smallest
    # of (a) the number of tasks and (b) the number of processors minus
    # one, since one processor must be reserved for the parent process
    # that is managing the process pool.
    max_child_procs = max(max_procs - 1, 1)
    max_simultaneous = min(num_tasks, max_child_procs)
    if max_simultaneous > 1:
        # Parallelize the tasks, controlled by the parent process, and
        # distribute the child processors evenly among the pooled tasks.
        pool_size = max_simultaneous
        num_procs_per_task = max_child_procs // pool_size
    else:
        # Run tasks serially; each task runs in the same process as the
        # parent and can thus have all processors.
        pool_size = 1
        num_procs_per_task = max_procs
    logger.detail(f"Ended calculating pool size: {pool_size} "
                  f"({num_procs_per_task} processors per task)")
    return pool_size, num_procs_per_task


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
            # Tasks running in parallel may not have the same logger as
            # the parent process (this seems to be system-dependent).
            # If not, then this task's top logger must be configured to
            # match the configuration of the parent process.
            set_config(*self._config)
            close_file_stream = True
        else:
            close_file_stream = False
        task = fmt_func_args(self._func, *args, **kwargs)
        try:
            logger.task(f"Began task {task}")
            result = self._func(*args, **kwargs)
        except Exception as error:
            logger.error(error)
        else:
            logger.task(f"Ended task {task}:\n{result}\n")
            return result
        finally:
            if close_file_stream and logger.file_stream is not None:
                # If the logger's configuration needed to be set, then
                # it is not the same logger as for the parent process.
                # That means that it is using a separate file stream,
                # which should be closed explicitly to free up file
                # resources when this task finishes.
                logger.file_stream.close()


def dispatch(funcs: list[Callable] | Callable,
             max_procs: int,
             pass_n_procs: bool = True,
             raise_on_error: bool = False,
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
        Maximum number of processes to run at one time. Must be ≥ 1.
    pass_n_procs: bool
        Whether to pass the number of processes to the function as the
        keyword argument `n_procs`.
    raise_on_error: bool
        Whether to raise an error if any tasks fail (if False, only log
        a warning message).
    args: list[tuple] | tuple
        Positional arguments to pass to each function in `funcs`. Can be
        a list of tuples of positional arguments or a single tuple that
        is not in a list. If a single tuple, then each function receives
        `args` as positional arguments. If a list, then `args` must be
        the same length as `funcs`; each function `funcs[i]` receives
        `args[i]` as positional arguments.
    kwargs: dict[str, Any] | None
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
    pool_size, n_procs_per_task = calc_pool_size(n_tasks, max_procs)
    if pass_n_procs:
        # Add the number of processes as a keyword argument.
        kwargs = {**kwargs, "n_procs": n_procs_per_task}
    if pool_size > 1:
        # Run the tasks in parallel.
        with ProcessPoolExecutor(max_workers=pool_size) as pool:
            logger.task(f"Opened pool of {pool_size} processes")
            # Initialize an empty list of tasks to run.
            tasks: list[Future] = list()
            for func, task_args in zip(funcs, args, strict=True):
                # Create a new task and submit it to the process pool.
                task = Task(func)
                tasks.append(pool.submit(task, *task_args, **kwargs))
            # Run all the tasks in parallel and collect the results as
            # they become available.
            logger.task(f"Waiting for {n_tasks} tasks to finish")
            results = [task.result() for task in tasks]
        logger.task(f"Closed pool of {pool_size} processes")
    else:
        # Run the tasks in series.
        logger.task(f"Began running {n_tasks} task(s) in series")
        # Initialize an empty list of results from the tasks.
        results = list()
        for func, task_args in zip(funcs, args, strict=True):
            # Create a new task, run it in the current process, and add
            # its result to the list of results.
            task = Task(func)
            results.append(task(*task_args, **kwargs))
        logger.task(f"Ended running {n_tasks} task(s) in series")
    # Remove any failed runs (None values) from results.
    results = [result for result in results if result is not None]
    n_pass = len(results)
    n_fail = n_tasks - n_pass
    if n_fail:
        p_fail = n_fail / n_tasks * 100.
        message = f"Failed {n_fail} of {n_tasks} task(s) ({round(p_fail, 1)} %)"
        if raise_on_error:
            raise RuntimeError(message)
        logger.warning(message)
    else:
        logger.task(f"All {n_tasks} task(s) completed successfully")
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
