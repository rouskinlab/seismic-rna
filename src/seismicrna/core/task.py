from concurrent.futures import (
    FIRST_COMPLETED,
    Future,
    ProcessPoolExecutor,
    as_completed,
    wait,
)
from inspect import getmodule
from itertools import filterfalse, repeat
from os import getpid
from typing import Any, Callable, Iterable

from .logs import Level, logger, get_config, set_config
from .progress import DEFAULT_UNIT, ProgressBar, console_has_bars, drawing, refresh_all
from .progress import set_shared_config as set_progress_shared_config
from .validate import require_equal

# Seconds to wait for a task to finish before redrawing the progress bar.
# A task in another process erases the bar to write a log message, so this is
# how long the bar can be missing; redrawing more often than this would keep
# it on screen slightly more of the time but write to the console constantly
# throughout a long task.
REFRESH_INTERVAL = 0.05


def calc_pool_size(num_tasks: int, num_cpus: int, force_serial: bool = False):
    """Calculate the size of a process pool.

    Parameters
    ----------
    num_tasks: int
        Number of tasks to parallelize. Must be ≥ 1.
    num_cpus: int
        Number of CPUs available. Must be ≥ 1.
    force_serial: bool
        If False, tasks are run in parallel (reserving one processor for the
        parent process and distributing the remaining processors among tasks).
        If True, tasks run serially, but each task can still use the remaining
        processors (i.e. max_procs - 1, at minimum 1).

    Returns
    -------
    tuple[int, int]
        - Size of the pool (number of concurrent tasks). Always ≥ 1.
        - Number of CPUs for each task in the pool. Always ≥ 1.
    """
    if num_cpus < 1:
        logger.warning("num_cpus must be ≥ 1, but got {}; defaulting to 1", num_cpus)
        num_cpus = 1
    if num_tasks < 1:
        logger.warning("num_tasks must be ≥ 1, but got {}; defaulting to 1", num_tasks)
        num_tasks = 1
    # The number of tasks that can run concurrently is the smallest of
    # (a) the number of tasks and (b) the number of processors.
    num_cpus_for_tasks = max(num_cpus, 1)
    num_simultaneous_tasks = min(num_tasks, num_cpus_for_tasks)
    if num_simultaneous_tasks > 1 and not force_serial:
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
    """Wrap a parallelizable task in a try-except block so that if it
    fails, it just returns `None` rather than crashing the other tasks
    being run in parallel."""

    def __init__(self, func: Callable):
        self._func = func
        self._config = get_config()
        # Record the process and the open logging contexts of the process
        # that creates this task, so that if the task runs in a child
        # process, its log can begin one level deeper than its parent's.
        self._pid = getpid()
        self._context_levels = list(logger.context_levels)
        # Record whether any process is drawing progress bars, so that if
        # the task runs in another process, that process erases them before
        # writing a log message over them. This must be true of the console,
        # not of this process: a task that dispatches tasks of its own draws
        # nothing itself, but the process that began the command still does.
        self._console_has_bars = console_has_bars()

    @property
    def name(self):
        return f"{getmodule(self._func).__name__}.{self._func.__name__}"

    def __call__(self, *args, **kwargs):
        """Call the task's function in a try-except block, return the
        result if it succeeds, and return None otherwise."""
        if get_config() != self._config:
            # Tasks running in parallel may not have the same logger as
            # the parent process (this seems to be system-dependent).
            # If not, then this task's top logger must be configured to
            # match the configuration of the parent process.
            set_config(*self._config)
            close_file_stream = True
        else:
            close_file_stream = False
        pid = getpid()
        if pid != self._pid:
            # This task is running in a child process, whose logger starts
            # with no contexts. Reproduce the parent's contexts and add one
            # always-visible level (the lowest level, shown at any verbosity)
            # so the child's log nests one level deeper than its parent's.
            logger.context_levels = list(self._context_levels) + [min(Level)]
            # This task shares the console with its parent process, which
            # is the only process that may draw a progress bar on it.
            set_progress_shared_config(self._console_has_bars)
            logger.debug("Process {} is a child of process {}", pid, self._pid)
        try:
            return self._func(*args, **kwargs)
        finally:
            if close_file_stream and logger.file_stream is not None:
                # If the logger's configuration needed to be set, then
                # it is not the same logger as for the parent process.
                # That means that it is using a separate file stream,
                # which should be closed explicitly to free up file
                # resources when this task finishes.
                logger.file_stream.close()


def _iter_finished(futures: list[Future], ordered: bool, bar: ProgressBar):
    """Yield each future as it finishes, redrawing the progress bar while
    waiting for the next one.

    A task that runs in another process erases the progress bars before it
    writes a log message, since only the process that owns them can draw
    them; redrawing them here, in that process, brings them back below the
    message, no matter how many messages the tasks write."""
    if not drawing():
        # No bar is on the console, so there is nothing to draw again while
        # waiting: just take each task as it finishes, without waking up
        # periodically. This asks about the console rather than about this
        # task's own bar, because these tasks may write messages over the
        # bars of the tasks that contain them even when they have none.
        yield from (futures if ordered else as_completed(futures))
        return
    pending = list(futures)
    while pending:
        # If the results must be returned in order, then wait only for the
        # first future, even if later ones finish before it.
        done, _ = wait(
            pending[:1] if ordered else pending,
            timeout=REFRESH_INTERVAL,
            return_when=FIRST_COMPLETED,
        )
        # Draw the whole block again, not just this task's bar: a task in
        # another process erases every line of it to write a message. This
        # happens on every pass, not only while waiting, so that the bars
        # come back promptly even when tasks are finishing in a burst.
        refresh_all()
        if not done:
            continue
        for future in list(pending):
            if future not in done:
                if ordered:
                    break
                continue
            pending.remove(future)
            yield future


def _dispatch(
    funcs: Callable | list[Callable],
    *,
    num_cpus: int,
    pass_num_cpus: bool,
    ordered: bool,
    raise_on_error: bool,
    force_serial: bool = False,
    label: str = "",
    unit: str = DEFAULT_UNIT,
    args: tuple | Iterable[tuple] = (),
    kwargs: dict[str, Any] | None = None,
):
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
        logger.warning("No tasks were given to dispatch")
        return list()
    # Determine how to parallelize each task.
    pool_size, num_cpus_per_task = calc_pool_size(
        num_tasks, num_cpus, force_serial=force_serial
    )
    logger.trace(
        "Run {} task(s) simultaneously, each using {} processor(s)",
        pool_size,
        num_cpus_per_task,
    )
    if pass_num_cpus:
        # Add the number of processes as a keyword argument.
        kwargs = {**kwargs, "num_cpus": num_cpus_per_task}
    # Run the tasks, counting each one that finishes (whether or not it
    # succeeds) on a progress bar; if these tasks have no label, then no bar
    # is drawn. The bar is advanced explicitly on both paths rather than in a
    # finally block, which would also count the GeneratorExit raised if this
    # generator is abandoned.
    with ProgressBar(label, num_tasks, unit) as bar:
        if pool_size > 1:
            # Run the tasks in parallel.
            with ProcessPoolExecutor(max_workers=pool_size) as pool:
                logger.trace("Opened process pool with {} processors", pool_size)
                # Create and submit a Future for each task.
                futures = [
                    pool.submit(Task(func), *task_args, **kwargs)
                    for func, task_args in zip(funcs, args, strict=True)
                ]
                for future in _iter_finished(futures, ordered, bar):
                    try:
                        result = future.result()
                    except Exception as error:
                        bar.tick()
                        if raise_on_error:
                            raise error
                        else:
                            logger.error(error)
                    else:
                        bar.tick()
                        yield result
            logger.trace("Closed process pool with {} processors", pool_size)
        else:
            # Run the tasks in series.
            for func, task_args in zip(funcs, args, strict=True):
                try:
                    task = Task(func)
                    result = task(*task_args, **kwargs)
                except Exception as error:
                    bar.tick()
                    if raise_on_error:
                        raise error
                    else:
                        logger.error(error)
                else:
                    bar.tick()
                    yield result


def dispatch(
    funcs: Callable | list[Callable],
    *,
    num_cpus: int,
    pass_num_cpus: bool,
    as_list: bool,
    ordered: bool,
    raise_on_error: bool,
    force_serial: bool = False,
    label: str = "",
    unit: str = DEFAULT_UNIT,
    args: tuple | Iterable[tuple] = (),
    kwargs: dict[str, Any] | None = None,
):
    """Run one or more tasks in series or in parallel, depending on the
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
    force_serial: bool
        Run the tasks in series even if multiple CPUs are available
        (e.g. because each task parallelizes its own work and would
        otherwise spawn a nested pool of processes).
    label: str
        Name of these tasks on the progress bar. If empty (the default),
        then these tasks get no progress bar, and whichever task contains
        them keeps the bar. Name only tasks that a user would want to
        follow, which means roughly one task per sample or per file: the
        parts of one sample's analysis (batches of reads, regions, runs
        of an algorithm) begin and end too quickly to follow, and naming
        them would replace the bar many times per second. Note also that
        with `ordered` set to True, the progress bar advances in bursts,
        since each task's result is awaited in the order given, not the
        order finished.
    unit: str
        What one task is, named on the progress bar beside the numbers so
        that they read as e.g. "sample 2/3" rather than a bare "2/3".
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
    results = _dispatch(
        funcs,
        num_cpus=num_cpus,
        pass_num_cpus=pass_num_cpus,
        ordered=ordered,
        raise_on_error=raise_on_error,
        force_serial=force_serial,
        label=label,
        unit=unit,
        args=args,
        kwargs=kwargs,
    )
    return list(results) if as_list else iter(results)


def as_list_of_tuples(args: Iterable[Any]):
    """Given an iterable of arguments, return a list of 1-item tuples,
    each containing one of the given arguments. This function is useful
    for creating a list of tuples to pass to the `args` parameter of
    `dispatch`."""
    return [(arg,) for arg in args]
