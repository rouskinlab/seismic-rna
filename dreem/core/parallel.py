from concurrent.futures import Future, ProcessPoolExecutor
from functools import wraps
from itertools import filterfalse, repeat
from logging import getLogger
import os
from pathlib import Path
from shutil import rmtree
from typing import Any, Callable, Iterable

logger = getLogger(__name__)

# DREEM lock directory (for function lock_output)
LOCK_DIR = ".dreem-lock"


def lock_temp_dir(run: Callable):
    @wraps(run)
    def wrapper(*args, temp_dir: str | Path, save_temp: bool, **kwargs):
        lock_error = (f"The directory {temp_dir} is currently being used by "
                      f"another run of DREEM. If possible, use a different "
                      f"temporary directory that is not in use. If another "
                      f"run of DREEM crashed and left this directory locked, "
                      f"then you may maually delete {temp_dir}.")
        # Determine whether the temporary directory and the lock exist.
        lock = os.path.join(temp_dir, LOCK_DIR)
        try:
            os.mkdir(lock)
        except FileExistsError:
            # The lock already exists, which means another instance of
            # DREEM is using this temporary directory.
            raise SystemExit(lock_error)
        except FileNotFoundError:
            # The temporary directory does not exist yet, so create it
            # along with a lock.
            try:
                os.makedirs(lock, exist_ok=False)
            except FileExistsError:
                # If this error happens, it is due to a very unlikely
                # race condition wherein another run of DREEM raises a
                # FileNotFoundError from the step os.mkdir(lock), then
                # this run of DREEM does the same, then the first run
                # creates the directory with the step os.makedirs(lock),
                # and then this run tries to do the same thing but fails
                # because the directory was created moments before.
                raise SystemExit(lock_error)
            temp_dir_existed_before = False
            logger.debug(f"Created and locked temporary directory: {temp_dir}")
        else:
            # The temporary directory had existed, but the lock had not.
            temp_dir_existed_before = True
            logger.debug(f"Locked temporary directory: {temp_dir}")
        # The lock now exists, and any other run of DREEM that tires to
        # use the same lock will exit before it can use the temporary
        # directory or delete the lock. Thus, this run must delete the
        # lock when it exits, regardless of the circumstances.
        try:
            if temp_dir_existed_before and not save_temp:
                raise SystemExit(f"The directory {temp_dir} existed before "
                                 f"this run of DREEM was launched, and may "
                                 f"thus contain arbitrary files. Because DREEM "
                                 f"would delete this temporary directory after "
                                 f"finishing, causing unintentional loss of "
                                 f"any and all data in this directory, DREEM "
                                 f"is exiting as a precaution. Either specify "
                                 f"a temporary directory that does not exist "
                                 f"yet or use the option --save-temp (CLI) / "
                                 f"save_temp=True (API) to tell DREEM to not "
                                 f"delete the temporary directory on exiting.")
            # Run the wrapped function and capture its result.
            res = run(*args, temp_dir=temp_dir, save_temp=save_temp, **kwargs)
            if not save_temp:
                # If the run completes successfully and the temporary
                # directory should not be saved, then delete it.
                rmtree(temp_dir, ignore_errors=True)
                logger.debug(f"Deleted temporary directory: {temp_dir}")
            return res
        finally:
            # Always ensure that the temporary directory is unlocked
            # upon exiting.
            try:
                os.rmdir(lock)
                logger.debug(f"Unlocked temporary directory: {temp_dir}")
            except FileNotFoundError:
                pass

    # Return the decorator.
    return wrapper


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
    fargs = ", ".join(map(repr, args))
    fkwargs = ", ".join(f"{kw}={repr(arg)}" for kw, arg in kwargs.items())
    return f"{func.__name__}({fargs}, {fkwargs})"


class Task(object):
    """ Wrap a parallelizable task in a try-except block so that if it
    fails, it just returns `None` rather than crashing the other tasks
    being run in parallel. """
    def __init__(self, func: Callable):
        self._func = func

    def __call__(self, *args, **kwargs):
        """ Call the task's function in a try-except block, return the
        result if it succeeds, and return `None` otherwise. """
        task = fmt_func_args(self._func, *args, **kwargs)
        logger.debug(f"Began task {task}")
        try:
            result = self._func(*args, **kwargs)
        except Exception as error:
            logger.error(f"Failed task {task}: {error}", exc_info=True)
        else:
            logger.debug(f"Ended task {task}: {result}")
            return result


def dispatch(funcs: list[Callable] | Callable,
             max_procs: int, parallel: bool, *,
             hybrid: bool = False,
             pass_n_procs: bool = True,
             drop_failed: bool = True,
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
        See docstring for ```get_num_parallel```.
    parallel: bool
        See docstring for ```get_num_parallel```.
    hybrid: bool = False
        See docstring for ```get_num_parallel```.
    pass_n_procs: bool = True
        Whether to pass the number of processes to the function as the
        keyword argument `n_procs`.
    drop_failed: bool = True
        If True, remove any failed runs from the results, where failure
        is indicated by the run returning a value of `None`.
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
                tasks.append(pool.submit(Task(func), *task_args, **kwargs))
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
            results.append(Task(func)(*task_args, **kwargs))
        logger.info(f"Ended running {n_tasks} task(s) in series")
    if drop_failed:
        # Remove any failed runs (None values) from results.
        results = [result for result in results if result is not None]
        n_pass = len(results)
        n_fail = n_tasks - n_pass
        logger.info(
            f"Tasks passed: {n_pass: >6d} ({100 * n_pass / n_tasks: >6.2f} %)")
        logger.info(
            f"Tasks failed: {n_fail: >6d} ({100 * n_fail / n_tasks: >6.2f} %)")
    return results


def as_list_of_tuples(args: Iterable[Any]):
    """ Given an iterable of arguments, return a list of 1-item tuples,
    each containing one of the given arguments. This function is useful
    for creating a list of tuples to pass to the `args` parameter of
    `dispatch`. """
    return [(arg,) for arg in args]
