import atexit
import os
import warnings
from collections import namedtuple
from functools import wraps
from os import getpid
from typing import Callable

from tqdm import tqdm

from . import logs
from .logs import DEFAULT_VERBOSITY, Level, logger

# Stop tqdm from starting its monitor thread, whose only purpose is to redraw
# bars that have not been updated recently. SEISMIC-RNA redraws them itself
# while waiting for tasks (see core.task), and a background thread would be a
# hazard here: tasks run in processes forked from this one, and forking from a
# process whose threads hold locks (such as the lock on standard error) can
# deadlock the child.
tqdm.monitor_interval = 0

# Progress bars are opt-in: the command line interface turns them on, while
# the Python API leaves them off, as it has no terminal of its own.
DEFAULT_ENABLED = False

# Verbosity below which progress bars are hidden. The point of --quiet is to
# produce little or no output on the terminal, which a progress bar would
# defeat.
MIN_VERBOSITY = Level.INFO

# Color of the bar itself, used only if log messages are also colored.
BAR_COLOR = "green"

# Number of characters in the gauge, and in the description that precedes it.
# Both are fixed so that the gauge stays in the same columns of the terminal
# as tasks begin and end: a gauge that slid left and right, or that stretched
# to whatever room the description left it, would be hard to read at a glance.
BAR_WIDTH = 14
DESC_WIDTH = 21

# Mark that the description was too long to fit and was cut short.
DESC_CUT = "…"

# Layout of one progress bar, which tqdm renders on exactly one line. The
# items are named after the gauge ("sample 2/3" rather than a bare "2/3") so
# that the numbers cannot be read as counting whatever the task is doing.
BAR_FORMAT = (
    f"{{desc}} |{{bar:{BAR_WIDTH}}}| "
    "{unit} {n_fmt}/{total_fmt} [{elapsed}<{remaining}]"
)


def format_done(elapsed: str) -> str:
    """Layout of a task that has finished.

    It shows the time the task took rather than the time still to come, of
    which there is none. That time is written into the layout rather than
    left for tqdm to work out, because the bar of a finished task is drawn
    again every time a message is logged, and a time worked out then would be
    the time since the task began, which keeps rising long after it ended."""
    return (
        f"{{desc}} |{{bar:{BAR_WIDTH}}}| {{unit}} {{n_fmt}}/{{total_fmt}} [{elapsed}]"
    )


# What a task counts, if it does not say.
DEFAULT_UNIT = "item"

# Return to the start of the line and erase from there to the end of the
# screen, which is where the progress bar is drawn. A process that shares the
# console with the process drawing the bar writes this before every message,
# since it cannot draw the bar again afterwards: erasing the bar is the only
# way to keep a message from being written over it.
ERASE_BARS = "\r\033[J"


def format_desc(label: str) -> str:
    """Format the description of a progress bar.

    The label is padded, or cut short, to a fixed width, so that the gauge
    after it never moves and the gauges of nested tasks line up in a column.
    Nothing precedes the label: a bar needs no level, as every bar is a bar,
    and no process ID, as every bar belongs to the process running the
    command."""
    if len(label) > DESC_WIDTH:
        label = f"{label[: DESC_WIDTH - len(DESC_CUT)]}{DESC_CUT}"
    return f"{label:<{DESC_WIDTH}}"


# Whether progress bars are currently enabled, the ID of the process that
# enabled them (see _owns_console), and whether another process is drawing
# progress bars on this process's console (see _shares_console).
_enabled = DEFAULT_ENABLED
_pid = getpid()
_shared = False

# Bars being drawn in this process, in the order their tasks began. Each one
# occupies its own line of the console, the outermost task on top; tqdm
# assigns and frees those lines as bars open and close.
_open_bars: list["ProgressBar"] = list()

# Rows to leave for the log above the block of bars, and how many rows to
# assume the console has if it cannot be measured.
LOG_ROWS = 4
DEFAULT_ROWS = 24


def _max_bars() -> int:
    """How many bars may stay on the console at once.

    Every task that finishes keeps its bar, so that the console ends up
    showing what each task did and how long it took. The block may not fill
    the console, though, or there would be no room for the log; once it is
    full, the bar of each further task goes into the log instead.

    The console is measured directly rather than with shutil, which measures
    standard output: the bars are on standard error, which is a terminal even
    when standard output has been redirected somewhere else."""
    try:
        rows = os.get_terminal_size(logs.stderr.fileno()).lines
    except Exception:
        rows = DEFAULT_ROWS
    return max(rows - LOG_ROWS, 1)


def _commit_bars():
    """Leave the block of bars on the console for good.

    The cursor sits at the start of the first bar, where a task in another
    process expects it, so it is moved below the last bar; otherwise the next
    thing written to the console, such as the shell prompt, would overwrite
    the record of the run."""
    if _owns_console() and _open_bars:
        try:
            logs.stderr.write("\n" * len(_open_bars))
            logs.stderr.flush()
        except Exception:
            pass
    _open_bars.clear()


def drawing() -> bool:
    """Whether any bar is on the console, and so needs drawing again.

    This asks about the console, not about one task: a task that draws no bar
    of its own still runs while the tasks containing it have bars on screen,
    and those bars must be drawn again as its own tasks write messages over
    them."""
    return bool(_open_bars) and _owns_console()


def _live_bars():
    """Bars whose tasks are still running."""
    return [bar for bar in _open_bars if not bar.done]


def refresh_all():
    """Draw every bar again, the outermost task first.

    A task in another process erases the whole block of bars before writing a
    log message, since it cannot know how many lines the block has; drawing
    only the innermost bar again would leave the tasks that contain it blank.

    The whole block is drawn while holding the lock that tqdm takes for each
    bar, which a task in another process takes to write a message: drawing
    one bar at a time would let a message land between two of them, in the
    middle of the block."""
    with tqdm.get_lock():
        for bar in _open_bars:
            bar.refresh()


def _owns_console() -> bool:
    """Whether this process may draw progress bars on the console.

    A process that is forked while running a task inherits both the console
    hook and the progress bars of its parent; drawing from the child would
    overwrite the output of the parent, so only the process that enabled the
    progress bars ever draws them."""
    return _enabled and getpid() == _pid


def _shares_console() -> bool:
    """Whether another process is drawing progress bars on this console.

    This is checked only after _owns_console has returned False, so progress
    bars being enabled here means that this process was forked from the one
    that enabled them (and thus shares its console) but has not yet been told
    so by its task, which happens when the task starts running."""
    return _shared or _enabled


# How Python showed warnings before core.progress took over, so that it can
# be put back when progress bars are disabled.
_show_warning = None


def _write_warning(message, category, filename, lineno, file=None, line=None):
    """Write a warning above the progress bars.

    Python writes warnings straight to standard error, which would put them
    through the bars, since they do not go through the logger."""
    write_console(warnings.formatwarning(message, category, filename, lineno, line))


def console_has_bars() -> bool:
    """Whether any process is drawing progress bars on this console.

    A task passes this to the tasks that it dispatches in turn, so that a
    task nested any number of processes deep still erases the bars before it
    writes a message, even though only the process that began the command
    draws them."""
    return _enabled or _shared


def _console_is_terminal() -> bool:
    """Whether the console (stderr) is a terminal."""
    try:
        return logs.stderr.isatty()
    except (AttributeError, ValueError):
        # The stream has no isatty method, or has been detached or closed.
        return False


def write_console(body: str):
    """Write the body of a log message above any open progress bars.

    This function is installed as core.logs.console_hook while progress bars
    are enabled, so that messages scroll up the terminal while the bars stay
    pinned to the bottom. If drawing fails for any reason, progress bars are
    disabled for the rest of the run and the message is written directly, so
    that drawing can never break logging. Nothing here may log a message,
    which would recurse back into this function."""
    try:
        if _owns_console():
            # tqdm.write clears the bars, writes the body, and redraws the
            # bars; the body already ends with a line break.
            tqdm.write(body, file=logs.stderr, end="")
            return
        if _shares_console():
            # Another process is drawing the bar, so it cannot be redrawn
            # here; erase it and let that process draw it again. Hold the
            # lock that tqdm takes whenever it draws: this process inherited
            # that lock when it was forked from the process that draws the
            # bar, so holding it keeps the two from writing at the same time
            # and interleaving a message with a bar.
            with tqdm.get_lock():
                logs.stderr.write(f"{ERASE_BARS}{body}")
                logs.stderr.flush()
            return
    except Exception:
        erase_config()
    logs.stderr.write(body)


class ProgressBar(object):
    """Progress of one task, shown on its own line at the bottom of the
    console.

    A progress bar is a context manager that finishes itself when its block
    ends, whether it ends normally, because of an exception, or because a
    generator containing it was closed. Every method does nothing if progress
    bars are disabled, so a caller never needs to check first.

    Tasks nest: a task of the workflow dispatches tasks of its own, each of
    which may dispatch more. Each task that draws a bar gets its own line, so
    that the step of the workflow stays visible above the sample it is
    working on. Their descriptions are padded to the same width, so the
    gauges line up in a column.

    A task with no name (an empty label) gets no bar: every method still
    works, but nothing is drawn and no line is taken. This is how the parts
    of one sample's analysis stay out of the way of the sample-level progress
    that a user actually follows.

    When a task finishes, its bar is held on its line until another task
    takes that line or the command ends, rather than blanking the line in the
    moment between one task finishing and the next beginning.

    Examples
    --------
    >>> with ProgressBar("aligning", 3) as bar:
    ...     for item in range(3):
    ...         bar.tick()
    """

    __slots__ = ["label", "total", "count", "unit", "done", "_shown", "_bar"]

    def __init__(self, label: str, total: int, unit: str = DEFAULT_UNIT):
        self.label = label
        self.total = total
        self.unit = unit
        self.count = 0
        self.done = False
        # Draw a bar only for a task that is named, and only in the process
        # that draws bars.
        self._shown = bool(label) and _owns_console()
        if self._shown:
            # Take the next line below the tasks already on the console, so
            # that the block grows downward in the order the tasks began.
            _open_bars.append(self)
        self._bar = tqdm(
            desc=self._format_desc(),
            total=total,
            initial=self.number,
            unit=unit,
            file=logs.stderr,
            leave=False,
            dynamic_ncols=True,
            disable=not self._shown,
            colour=BAR_COLOR if logger.log_color else None,
            bar_format=BAR_FORMAT,
        )

    @property
    def shown(self) -> bool:
        """Whether this bar is being drawn on the console."""
        return self._shown

    @property
    def number(self) -> int:
        """Number of the item being worked on, counting the first item as 1.

        The gauge tracks this number rather than the number of items that
        have finished, so that it agrees with the number shown and is full
        while the last item is running (immediately, for a task of one item).
        Once every item has finished, this stays at the number of the last
        item, since no further item begins."""
        return min(self.count + 1, self.total)

    def _format_desc(self) -> str:
        """Format this bar's description, which names the task; what it is
        counting is named beside the numbers, after the gauge."""
        return format_desc(self.label)

    def _erase(self):
        """Erase this bar from the console for good."""
        try:
            self._bar.disable = False
            self._bar.close()
        except Exception:
            pass

    def _freeze(self):
        """Stop the time on this bar at the time the task took, so that it
        does not keep rising every time the bar is drawn again."""
        try:
            elapsed = self._bar.format_interval(self._bar.format_dict["elapsed"])
            self._bar.bar_format = format_done(elapsed)
        except Exception:
            pass

    def _graduate(self):
        """Erase this bar from the block and write it into the log, where it
        stays as a record of what the task did and how long it took."""
        self._freeze()
        try:
            done = self._bar.format_meter(**self._bar.format_dict)
        except Exception:
            done = ""
        self._erase()
        if done:
            write_console(f"{done}\n")

    def set_label(self, label: str):
        """Rename the bar, e.g. when a multi-step task moves to the next
        step, and draw it again."""
        self.label = label
        self._bar.set_description_str(self._format_desc(), refresh=True)

    def tick(self, n: int = 1, label: str = ""):
        """Advance by `n` items, first renaming the bar if `label` is
        given."""
        number = self.number
        self.count = min(self.count + n, self.total)
        if label:
            self.label = label
            self._bar.set_description_str(self._format_desc(), refresh=False)
        self._bar.update(self.number - number)

    def update(self, count: int):
        """Set the number of items completed to an absolute value."""
        self.tick(count - self.count)

    def refresh(self):
        """Draw the bar again without changing it.

        A task running in another process erases the bar before writing a log
        message, since it cannot draw the bar again itself, so the process
        that owns the bar redraws it while waiting for its tasks."""
        self._bar.refresh()

    def close(self):
        """Finish this task, holding its bar on its line until another task
        takes the line or the command ends.

        Doing so must never raise an exception, since a bar can be closed
        while another exception (including GeneratorExit) is propagating, and
        closing a bar twice must be harmless, since a bar can be closed both
        by its own block and by erase_config."""
        if self.done or not self._shown:
            self.done = True
            return
        self.done = True
        if len(_open_bars) <= _max_bars():
            # Keep this finished task on its line, so that the console ends
            # up showing every task and how long it took.
            self._freeze()
            self.refresh()
            return
        # The block is full, so this task goes into the log instead, and
        # gives its line back for the next task to use.
        try:
            _open_bars.remove(self)
        except ValueError:
            return
        self._graduate()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


class ConsoleStream(object):
    """File-like object that writes above any open progress bars.

    It exists for code that writes to the console without going through the
    logger, in particular unittest, which writes directly to standard error.
    Text is buffered until a line is complete, because the progress bar is
    drawn immediately below the last complete line, so writing a partial line
    would leave the cursor in the middle of it. If no progress bar is being
    drawn, text is passed straight through and nothing is buffered."""

    __slots__ = ["_buffer"]

    def __init__(self):
        self._buffer = ""

    def write(self, text: str):
        """Write text, holding back any characters after the last line
        break until the line they belong to is complete."""
        if not _owns_console():
            logs.stderr.write(text)
            return len(text)
        lines = f"{self._buffer}{text}".split("\n")
        # The text after the last line break belongs to a line that is not
        # yet complete, so keep buffering it.
        self._buffer = lines.pop()
        for line in lines:
            write_console(f"{line}\n")
        return len(text)

    def flush(self):
        """Flush the console.

        Any partial line stays buffered: it is written as soon as the line is
        complete, or when this stream is closed."""
        if not _owns_console():
            logs.stderr.flush()

    def close(self):
        """Write any buffered partial line, ending it with a line break."""
        if self._buffer:
            write_console(f"{self._buffer}\n")
            self._buffer = ""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def close_bars(func: Callable):
    """After the function exits, close any progress bars it left open.

    A progress bar normally closes itself at the end of its own block, but a
    bar opened inside a generator (as dispatch does when `as_list` is False)
    closes only once that generator is closed or garbage collected, which
    never happens if the caller abandons it partway through. Closing the
    leftover bars when a command ends keeps an abandoned bar from lingering
    on the console for the rest of the run."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        opened = {id(bar) for bar in _open_bars}
        try:
            return func(*args, **kwargs)
        finally:
            for bar in _live_bars():
                if id(bar) not in opened:
                    bar.close()

    return wrapper


ProgressConfig = namedtuple("ProgressConfig", ["enabled"])


def erase_config():
    """Erase the existing progress bar configuration.

    Every task still running is finished, and the whole block of bars is left
    on the console as a record of the run, with the cursor moved below it so
    that whatever is written next does not overwrite it. The hook is removed
    from core.logs, so that log messages are written directly to standard
    error again."""
    global _enabled, _shared, _show_warning
    for bar in _live_bars():
        bar.close()
    _commit_bars()
    _enabled = False
    _shared = False
    logs.console_hook = None
    if _show_warning is not None:
        warnings.showwarning = _show_warning
        _show_warning = None


def set_config(enabled: bool = DEFAULT_ENABLED, verbosity: int = DEFAULT_VERBOSITY):
    """Configure progress bars.

    Parameters
    ----------
    enabled: bool = False
        Whether to draw progress bars. They are drawn only if this is True,
        the verbosity is at least MIN_VERBOSITY, and the console (standard
        error) is a terminal; otherwise, this parameter has no effect.
        Enabling progress bars installs a hook in core.logs so that log
        messages are written above the bars rather than through them.
    verbosity: int = 0
        Verbosity of the log, as in core.logs.set_config. This verbosity is
        given explicitly rather than read from the logger because a command
        may change the verbosity of the logger while it runs (as the test
        command does) without meaning to hide its own progress bars.
    """
    global _enabled, _pid, _show_warning
    erase_config()
    if enabled and verbosity >= MIN_VERBOSITY and _console_is_terminal():
        _enabled = True
        _pid = getpid()
        logs.console_hook = write_console
        _show_warning = warnings.showwarning
        warnings.showwarning = _write_warning


def set_shared_config(shared: bool = False):
    """Configure this process to share the console with the process that is
    drawing the progress bars.

    A task that runs in another process must not draw progress bars of its
    own, but its log messages would otherwise be written over the bars that
    its parent process is drawing. Sharing the console makes every message
    from this process erase the bars first, so that its parent can draw them
    again where they belong: below the message. This must be configured
    explicitly rather than inferred, because a task process may be forked
    from its parent (inheriting its configuration) or started fresh.

    Parameters
    ----------
    shared: bool = False
        Whether another process is drawing progress bars on this console.
    """
    global _enabled, _shared
    _enabled = False
    _shared = shared
    # Forget, but do not close, any bars inherited from the parent process:
    # closing them would draw over the console that the parent owns.
    _open_bars.clear()
    logs.console_hook = write_console if shared else None


def get_config():
    """Get the configuration parameters of the progress bars."""
    return ProgressConfig(enabled=_enabled)


def restore_config(func: Callable):
    """After the function exits, restore the progress bar configuration that
    was in place before the function ran."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        config = get_config()
        try:
            return func(*args, **kwargs)
        finally:
            set_config(**config._asdict())

    return wrapper


# Leave the console clean if SEISMIC-RNA exits while progress bars are open.
atexit.register(erase_config)
