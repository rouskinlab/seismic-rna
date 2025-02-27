from concurrent.futures import ProcessPoolExecutor
from inspect import getmodule
from itertools import filterfalse, repeat
from typing import Any, Callable, Iterable

from .logs import logger, get_config, set_config
from .validate import require_equal


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
    return pool_size, num_procs_per_task


class Task(object):
    """ Wrap a parallelizable task in a try-except block so that if it
    fails, it just returns `None` rather than crashing the other tasks
    being run in parallel. """

    def __init__(self, func: Callable):
        self._func = func
        self._config = get_config()

    @property
    def name(self):
        return f"{getmodule(self._func).__name__}.{self._func.__name__}"

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
        try:
            logger.task(f"Began {self.name}")
            result = self._func(*args, **kwargs)
        except Exception as error:
            logger.error(error)
        else:
            logger.task(f"Ended {self.name}")
            return result
        finally:
            if close_file_stream and logger.file_stream is not None:
                # If the logger's configuration needed to be set, then
                # it is not the same logger as for the parent process.
                # That means that it is using a separate file stream,
                # which should be closed explicitly to free up file
                # resources when this task finishes.
                logger.file_stream.close()


def dispatch(funcs: Callable | list[Callable],
             max_procs: int,
             pass_n_procs: bool = True,
             raise_on_error: bool = False,
             args: tuple | Iterable[tuple] = (),
             kwargs: dict[str, Any] | None = None):
    """
    Run one or more tasks in series or in parallel, depending on the
    number of tasks, the maximum number of processes, and whether tasks
    are allowed to be run in parallel.

    Parameters
    ----------
    funcs: Callable | list[Callable]
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
    args: tuple | Iterable[tuple]
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
            args = list(args)
            # Ensure that every item in args is a tuple.
            nontuple = list(filterfalse(lambda x: isinstance(x, tuple), args))
            if nontuple:
                raise TypeError(f"Got non-tuple args: {nontuple}")
        # If a function is given rather than an iterable of functions,
        # then put the function in a list whose length equal that of the
        # list of arguments.
        funcs = list(repeat(funcs, len(args)))
    else:
        # Ensure that every item in funcs is actually callable.
        uncallable = list(filterfalse(callable, funcs))
        if uncallable:
            raise TypeError(f"Got uncallable funcs: {uncallable}")
        if isinstance(args, tuple):
            # If args is a tuple, repeat it once for each function.
            args = list(repeat(args, len(funcs)))
    # Ensure that numbers of functions and argument tuples match.
    num_tasks = len(funcs)
    require_equal("len(funcs)", num_tasks, len(args), "len(args)")
    if num_tasks == 0:
        # No tasks to run: return.
        logger.task("No tasks were given to dispatch")
        return list()
    # Determine how to parallelize each task.
    pool_size, n_procs_per_task = calc_pool_size(num_tasks, max_procs)
    if pass_n_procs:
        # Add the number of processes as a keyword argument.
        kwargs = {**kwargs, "n_procs": n_procs_per_task}
        logger.detail(f"Calculated size of process pool: {pool_size}, "
                      f"each with {n_procs_per_task} processor(s)")
    else:
        logger.detail(f"Calculated size of process pool: {pool_size}")
    if pool_size > 1:
        # Run the tasks in parallel.
        with ProcessPoolExecutor(max_workers=pool_size) as pool:
            logger.task(f"Opened pool of {pool_size} processes")
            # Create a Future for each task and submit it to the pool.
            futures = [pool.submit(Task(func), *task_args, **kwargs)
                       for func, task_args in zip(funcs, args, strict=True)]
            # Run all the tasks in parallel and collect the results as
            # they become available.
            logger.task(f"Waiting for {num_tasks} tasks to finish")
            results = [future.result() for future in futures]
        logger.task(f"Closed pool of {pool_size} processes")
    else:
        # Run the tasks in series.
        logger.task(f"Began running {num_tasks} task(s) in series")
        # Initialize an empty list of results from the tasks.
        results = list()
        for func, task_args in zip(funcs, args, strict=True):
            # Create a new task, run it in the current process, and add
            # its result to the list of results.
            task = Task(func)
            results.append(task(*task_args, **kwargs))
        logger.task(f"Ended running {num_tasks} task(s) in series")
    # Remove any failed runs (None values) from results.
    results = [result for result in results if result is not None]
    num_pass = len(results)
    num_fail = num_tasks - num_pass
    if num_fail:
        message = f"Failed {num_fail} of {num_tasks} task(s)"
        if raise_on_error:
            raise RuntimeError(message)
        logger.warning(message)
    else:
        logger.task(f"All {num_tasks} task(s) completed successfully")
    return results


def as_list_of_tuples(args: Iterable[Any]):
    """ Given an iterable of arguments, return a list of 1-item tuples,
    each containing one of the given arguments. This function is useful
    for creating a list of tuples to pass to the `args` parameter of
    `dispatch`. """
    return [(arg,) for arg in args]
