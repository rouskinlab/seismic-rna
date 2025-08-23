from concurrent.futures import ProcessPoolExecutor, as_completed
from inspect import getmodule
from itertools import filterfalse, repeat
from typing import Any, Callable, Iterable

from .logs import Level, logger, get_config, set_config
from .validate import require_equal


def calc_pool_size(num_tasks: int, num_cpus: int):
    """ Calculate the size of a process pool.

    Parameters
    ----------
    num_tasks: int
        Number of tasks to parallelize. Must be ≥ 1.
    num_cpus: int
        Number of CPUs available. Must be ≥ 1.

    Returns
    -------
    tuple[int, int]
        - Size of the pool (number of concurrent tasks). Always ≥ 1.
        - Number of CPUs for each task in the pool. Always ≥ 1.
    """
    if num_cpus < 1:
        logger.warning(f"num_cpus must be ≥ 1, but got {num_cpus}; "
                       f"defaulting to 1")
        num_cpus = 1
    if num_tasks < 1:
        logger.warning(f"num_tasks must be ≥ 1, but got {num_tasks}; "
                       f"defaulting to 1")
        num_tasks = 1
    # The number of tasks that can run concurrently is the smallest of
    # (a) the number of tasks and (b) the number of processors.
    num_cpus_for_tasks = max(num_cpus, 1)
    num_simultaneous_tasks = min(num_tasks, num_cpus_for_tasks)
    if num_simultaneous_tasks > 1:
        # Parallelize the tasks, controlled by the parent process, and
        # distribute the child processors evenly among the pooled tasks.
        pool_size = num_simultaneous_tasks
        num_cpus_per_task = num_cpus_for_tasks // pool_size
    else:
        # Run tasks serially; each task runs in the same process as the
        # parent and can thus have all processors.
        pool_size = 1
        num_cpus_per_task = num_cpus
    return pool_size, num_cpus_per_task


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
        verbosity = self._config.verbosity
        description_items = list()
        if verbosity >= Level.ACTION:
            description_items.extend(map(repr, args))
        if verbosity >= Level.ROUTINE:
            description_items.extend(f"{k}={repr(v)}"
                                     for k, v in kwargs.items())
        description = f"{self.name}({', '.join(description_items)})"
        try:
            logger.task(f"Began {description}")
            result = self._func(*args, **kwargs)
            logger.task(f"Ended {description}")
            return result
        finally:
            if close_file_stream and logger.file_stream is not None:
                # If the logger's configuration needed to be set, then
                # it is not the same logger as for the parent process.
                # That means that it is using a separate file stream,
                # which should be closed explicitly to free up file
                # resources when this task finishes.
                logger.file_stream.close()


def _dispatch(funcs: Callable | list[Callable], *,
              num_cpus: int,
              pass_num_cpus: bool,
              ordered: bool,
              raise_on_error: bool,
              args: tuple | Iterable[tuple] = (),
              kwargs: dict[str, Any] | None = None):
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
    pool_size, num_cpus_per_task = calc_pool_size(num_tasks, num_cpus)
    if pass_num_cpus:
        # Add the number of processes as a keyword argument.
        kwargs = {**kwargs, "num_cpus": num_cpus_per_task}
        logger.detail(f"Calculated size of process pool: {pool_size}, "
                      f"each with {num_cpus_per_task} processor(s)")
    else:
        logger.detail(f"Calculated size of process pool: {pool_size}")
    # Run the tasks.
    num_failed = 0
    if pool_size > 1:
        # Run the tasks in parallel.
        with ProcessPoolExecutor(max_workers=pool_size) as pool:
            logger.task(f"Opened pool of {pool_size} processes")
            # Create and submit a Future for each task.
            futures = [pool.submit(Task(func), *task_args, **kwargs)
                       for func, task_args in zip(funcs, args, strict=True)]
            logger.task(f"Waiting for {num_tasks} tasks to finish")
            for future in (futures if ordered else as_completed(futures)):
                try:
                    yield future.result()
                except Exception as error:
                    if raise_on_error:
                        raise error
                    else:
                        logger.error(error)
                        num_failed += 1
        logger.task(f"Closed pool of {pool_size} processes")
    else:
        # Run the tasks in series.
        logger.task(f"Began running {num_tasks} task(s) in series")
        for func, task_args in zip(funcs, args, strict=True):
            try:
                task = Task(func)
                yield task(*task_args, **kwargs)
            except Exception as error:
                if raise_on_error:
                    raise error
                else:
                    logger.error(error)
                    num_failed += 1
        logger.task(f"Ended running {num_tasks} task(s) in series")
    if num_failed:
        message = f"Failed {num_failed} of {num_tasks} task(s)"
        if raise_on_error:
            raise RuntimeError(message)
        logger.warning(message)
    else:
        logger.task(f"All {num_tasks} task(s) completed successfully")


def dispatch(funcs: Callable | list[Callable], *,
             num_cpus: int,
             pass_num_cpus: bool,
             as_list: bool,
             ordered: bool,
             raise_on_error: bool,
             args: tuple | Iterable[tuple] = (),
             kwargs: dict[str, Any] | None = None):
    """ Run one or more tasks in series or in parallel, depending on the
    number of tasks and the maximum number of CPUs.

    Parameters
    ----------
    funcs: Callable | list[Callable]
        The function(s) to run. Can be a list of functions or a single
        function that is not in a list. If a single function, then if
        `args` is a tuple, it is called once with that tuple as its
        positional arguments; and if `args` is a list of tuples, it is
        called for each tuple of positional arguments in `args`.
    num_cpus: int
        Number of CPUs available. Must be ≥ 1.
    pass_num_cpus: bool
        Pass the number of processes to the function(s) in `funcs` as
        the keyword argument `num_cpus`.
    as_list: bool
        Return results as a list (if True) or an iterator (if False).
    ordered: bool
        Return results in the same order as they were given in `funcs`
        and/or `args` (if True) or in order of completion (if False).
    raise_on_error: bool
        If any task fails, then raise the exception that it raises (if
        True) or log that exception as an error (if False).
    args: tuple | Iterable[tuple]
        Positional arguments to pass to each function in `funcs`. Can be
        a list of tuples of positional arguments or a single tuple that
        is not in a list. If a single tuple, then each function receives
        `args` as positional arguments. If a list, then `args` must be
        the same length as `funcs`; each function `funcs[i]` receives
        `args[i]` as positional arguments.
    kwargs: dict[str, Any] | None
        Keyword arguments to pass to every function call.
    """
    results = _dispatch(funcs,
                        num_cpus=num_cpus,
                        pass_num_cpus=pass_num_cpus,
                        ordered=ordered,
                        raise_on_error=raise_on_error,
                        args=args,
                        kwargs=kwargs)
    return list(results) if as_list else iter(results)


def as_list_of_tuples(args: Iterable[Any]):
    """ Given an iterable of arguments, return a list of 1-item tuples,
    each containing one of the given arguments. This function is useful
    for creating a list of tuples to pass to the `args` parameter of
    `dispatch`. """
    return [(arg,) for arg in args]
