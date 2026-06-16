import os
import shutil
import textwrap
from collections import namedtuple
from datetime import datetime
from enum import IntEnum
from functools import cache, wraps
from pathlib import Path
from sys import stderr
from traceback import format_exception
from typing import Callable, Optional, TextIO


class Level(IntEnum):
    """Level of a logging message."""

    ERROR = -2
    WARNING = -1
    INFO = 0
    DEBUG = 1
    TRACE = 2


DEFAULT_COLOR = True
DEFAULT_EXIT_ON_ERROR = True
DEFAULT_VERBOSITY = Level.INFO

# Number of spaces by which to indent each level of nesting depth.
INDENT = 2 * " "


def format_sample_reference_region(sample: str, ref: str, region: str = ""):
    formatted = f"sample {repr(sample)} over reference {repr(ref)}"
    if region:
        return f"{formatted} region {repr(region)}"
    return formatted


class FileStream(object):
    """Log to a file, opened lazily on first write."""

    __slots__ = ["file_path", "_file"]

    def __init__(self, file_path: str | Path):
        self.file_path = Path(file_path)
        self._file = None

    @property
    def stream(self) -> TextIO:
        if self._file is None:
            self.file_path.parent.mkdir(exist_ok=True, parents=True)
            self._file = open(self.file_path, "a")
        return self._file

    def close(self):
        """Close the file stream."""
        if self._file is not None:
            self._file.close()
            self._file = None


def wrap_console(content: str, first_indent: str, cont_indent: str, width: int):
    """Wrap content to at most `width` columns.

    The first visual line is prefixed with `first_indent` and every other
    visual line (whether produced by wrapping or by a line break already in
    `content`) with `cont_indent`. Words longer than the available width are
    not broken, so the wrapped lines may occasionally exceed `width`."""
    lines = []
    for i, line in enumerate(content.split("\n")):
        lead = first_indent if i == 0 else cont_indent
        if len(lead) + len(line) <= width:
            # The line already fits; keep it (and its whitespace) intact.
            lines.append(lead + line)
        else:
            lines.extend(
                textwrap.wrap(
                    line,
                    width=width,
                    initial_indent=lead,
                    subsequent_indent=cont_indent,
                    break_long_words=False,
                    break_on_hyphens=False,
                )
                or [lead]
            )
    return "\n".join(lines)


def format_body(level: Level, content: str, depth: int) -> str:
    """Format a log message body, shared by console and file output.

    The level and process ID form a fixed-width prefix so that columns stay
    aligned; the message text is indented by nesting depth and wrapped to one
    less than the terminal width, with continuation lines aligned under the
    message column."""
    prefix = f"{level.name.capitalize():<7} {str(os.getpid())[:7]:>7} "
    indent = INDENT * depth
    first_indent = prefix + indent
    cont_indent = " " * len(prefix) + indent
    width = max(shutil.get_terminal_size().columns - 1, len(cont_indent) + 1)
    return wrap_console(content, first_indent, cont_indent, width) + "\n"


class AnsiCode(object):
    """Format text with ANSI codes."""

    # Control codes.
    START = "\033["  # Indicate the start of an ANSI format code.
    END = "m"  # Indicate the end of an ANSI format code.
    RESET = 0  # Reset any existing formatting.

    @classmethod
    @cache
    def format_color(cls, color: int):
        """Make a format string for one 256-color code."""
        if not 0 <= color < 256:
            raise ValueError(f"Invalid ANSI 256-color code: {color}")
        # 38 means set the foreground color (i.e. of the text itself).
        # 5 means use 256-color mode.
        # The next number is the color in 256-color mode.
        return f"{cls.START}38;5;{color}{cls.END}"

    @classmethod
    @cache
    def format(cls, code: int):
        """Make a format string for one ANSI code."""
        return f"{cls.START}{code}{cls.END}"

    @classmethod
    @cache
    def reset(cls):
        """Convenience function to end formatting."""
        return cls.format(cls.RESET)


# Level-specific color formatting.
# The color of each code can be displayed with this Bash command:
# for i in {0..255}; do echo -ne "\033[38;5;${i}m  ${i} "; done
LEVEL_COLORS = {
    Level.ERROR: AnsiCode.format_color(160),
    Level.WARNING: AnsiCode.format_color(214),
    Level.INFO: AnsiCode.format_color(28),
    Level.DEBUG: AnsiCode.format_color(69),
    Level.TRACE: AnsiCode.format_color(247),
}


def color_wrap(level: Level, body: str) -> str:
    """Bracket a formatted body with the ANSI color for the given level."""
    try:
        return f"{LEVEL_COLORS[level]}{body}{AnsiCode.reset()}"
    except KeyError:
        return body


class LoggingContext(object):
    """Context manager that delimits a task in the log.

    On entry, it logs ``Began {name}`` and pushes its level onto the logger's
    context stack so that messages logged inside the block are indented. On
    normal exit, it pops the context back off and logs ``Ended {name}``. If the
    block exits because of an exception, the context is still popped, but the
    ``Ended`` message is skipped (the task did not finish successfully)."""

    __slots__ = ["_logger", "_level", "_template", "_began_ended", "_args", "_kwargs"]

    def __init__(
        self,
        logger: "Logger",
        level: Level,
        template: object,
        began_ended: bool,
        args: tuple,
        kwargs: dict,
    ):
        self._logger = logger
        self._level = level
        self._template = template
        self._began_ended = began_ended
        self._args = args
        self._kwargs = kwargs

    def __enter__(self):
        if self._began_ended:
            template = f"Began {self._template}"
        else:
            template = self._template
        self._logger._log(self._level, template, self._args, self._kwargs)
        self._logger.context_levels.append(self._level)
        return self._logger

    def __exit__(self, exc_type, exc_value, traceback):
        # Guard against an empty stack so that logging can never crash the
        # program, even if the context stack was erased while this block ran.
        try:
            self._logger.context_levels.pop()
        except IndexError:
            self._logger.warning(
                "The logging context was erased in the middle of a context block: "
                "this indicates a bug in or misuse of the logging system."
            )
        if self._began_ended:
            if exc_type is not None:
                self._logger.warning(
                    f"FAILED {self._template}", self._args, self._kwargs
                )
            else:
                self._logger._log(
                    self._level, f"Ended {self._template}", self._args, self._kwargs
                )


class LevelLogger(object):
    """Log messages at one logging level.

    Calling the object logs its template at this level, formatting it with the
    given positional and keyword arguments via ``str.format`` only if the level
    is visible at the current verbosity. The ``begin`` method returns a context
    manager that brackets a task with ``Began``/``Ended`` messages and indents
    the messages logged inside it."""

    __slots__ = ["_logger", "_level"]

    def __init__(self, logger: "Logger", level: Level):
        self._logger = logger
        self._level = level

    def __call__(self, template: object, *args, **kwargs):
        self._logger._log(self._level, template, args, kwargs)

    def single_context(self, template: object, *args, **kwargs):
        """Context manager that logs a message at the beginning and indents
        messages logged inside the block."""
        return LoggingContext(self._logger, self._level, template, False, args, kwargs)

    def double_context(self, template: object, *args, **kwargs):
        """Context manager that logs a message at the beginning and the end
        and indents messages logged inside the block."""
        return LoggingContext(self._logger, self._level, template, True, args, kwargs)


class Logger(object):
    """Log messages to the console and to files."""

    def __init__(self):
        self.verbosity = DEFAULT_VERBOSITY
        self.log_color = DEFAULT_COLOR
        self.console_enabled = False
        self.file_stream = None
        self.exit_on_error = DEFAULT_EXIT_ON_ERROR
        # Levels of the currently open task contexts, used to indent messages.
        self.context_levels = []
        # Create one callable LevelLogger per logging level, e.g. so that
        # self.trace(template, *args) logs at the TRACE level.
        self.error = LevelLogger(self, Level.ERROR)
        self.warning = LevelLogger(self, Level.WARNING)
        self.info = LevelLogger(self, Level.INFO)
        self.debug = LevelLogger(self, Level.DEBUG)
        self.trace = LevelLogger(self, Level.TRACE)

    def _build_content(self, template: object, args: tuple, kwargs: dict) -> str:
        """Build the final message string from a template and optional args."""
        if isinstance(template, BaseException):
            return "".join(format_exception(template))
        if isinstance(template, str):
            s = template
        else:
            s = str(template)
        if args or kwargs:
            return s.format(*args, **kwargs)
        return s

    def _log(self, level: Level, template: object, args: tuple, kwargs: dict):
        """Log a message, formatting the template only if the level is visible."""
        if level <= Level.ERROR and self.exit_on_error:
            if isinstance(template, BaseException):
                raise template
            raise RuntimeError(self._build_content(template, args, kwargs))
        # Short-circuit if the verbosity is less than the logging level
        # or no logging outputs are enabled.
        if level > self.verbosity or not (self.console_enabled or self.file_stream):
            return
        content = self._build_content(template, args, kwargs)
        # Indent by the contexts that are visible at the current verbosity.
        depth = sum(1 for lvl in self.context_levels if lvl <= self.verbosity)
        body = format_body(level, content, depth)
        if self.console_enabled:
            if self.log_color:
                body = color_wrap(level, body)
            stderr.write(body)
        if self.file_stream is not None:
            timestamp = datetime.now().strftime(r"%Y-%m-%d %H:%M:%S.%f")
            self.file_stream.stream.write(f"{timestamp}  {body}")


logger = Logger()


LoggerConfig = namedtuple(
    "LoggerConfig", ["verbosity", "log_file_path", "log_color", "exit_on_error"]
)


def erase_config():
    """Erase the existing logger configuration.

    The context stack is deliberately left untouched: it is runtime state
    managed by the task contexts (and reset per worker by Task), not part of
    the configuration. Clearing it here would corrupt the stack if logging is
    reconfigured while a task context is open (e.g. the test command, which
    reconfigures logging inside its own ``begin`` block)."""
    logger.verbosity = DEFAULT_VERBOSITY
    logger.log_color = DEFAULT_COLOR
    logger.console_enabled = False
    if logger.file_stream is not None:
        logger.file_stream.close()
    logger.file_stream = None
    logger.exit_on_error = DEFAULT_EXIT_ON_ERROR


def set_config(
    verbosity: int = DEFAULT_VERBOSITY,
    log_file_path: str | Path | None = None,
    log_color: bool = DEFAULT_COLOR,
    exit_on_error: bool = DEFAULT_EXIT_ON_ERROR,
):
    """Configure the main logger with handlers and verbosity.

    Parameters
    ----------
    verbosity: int = 0
        Verbosity level; messages at or below this level are shown on the
        console and written to the log file.
    log_file_path: str | Path | None
        Path of the log file; logging to a file is disabled if None.
        The file records messages at the same verbosity as the console,
        with each line prefixed by a timestamp.
    log_color: bool = True
        Whether to use ANSI color codes in console output.
    exit_on_error: bool
        Whether to raise an exception (or exit) when an error is logged.
    """
    erase_config()
    logger.verbosity = verbosity
    logger.log_color = log_color
    logger.console_enabled = True
    if log_file_path is not None:
        logger.file_stream = FileStream(log_file_path)
    logger.exit_on_error = exit_on_error


def get_config():
    """Get the configuration parameters of the logger."""
    return LoggerConfig(
        verbosity=logger.verbosity,
        log_file_path=(
            logger.file_stream.file_path if logger.file_stream is not None else None
        ),
        log_color=logger.log_color,
        exit_on_error=logger.exit_on_error,
    )


def log_exceptions(default: Optional[Callable]):
    """If any exception occurs, catch it and return the default."""

    def decorator(func: Callable):

        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as error:
                logger.error(error)
                return default() if default is not None else None

        return wrapper

    return decorator


def restore_config(func: Callable):
    """After the function exits, restore the logging configuration (and the
    context stack) that was in place before the function ran."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        config = get_config()
        context_levels = list(logger.context_levels)
        try:
            return func(*args, **kwargs)
        finally:
            set_config(**config._asdict())
            logger.context_levels = context_levels

    return wrapper


# Initialize logging so that it will run by default if using the API.
set_config()
