import os
import shutil
import textwrap
from abc import ABC, abstractmethod
from collections import namedtuple
from datetime import datetime
from enum import IntEnum
from functools import cache, wraps
from pathlib import Path
from sys import stderr
from traceback import format_exception
from typing import Callable, Iterable, Optional, TextIO


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
FILE_VERBOSITY = Level.TRACE

# Number of spaces by which to indent each level of nesting depth.
INDENT = "  "


class Message(object):
    """Message with a logging level and enclosing nesting contexts."""

    __slots__ = ["level", "content", "context_levels"]

    def __init__(
        self, level: Level, content: object, context_levels: Iterable[Level] = ()
    ):
        self.level = level
        self.content = content
        # Levels of the task contexts enclosing this message, from outermost
        # to innermost; snapshotted so that later changes do not affect it.
        self.context_levels = tuple(context_levels)

    def depth(self, verbosity: int):
        """Indentation depth at the given verbosity: the number of enclosing
        contexts that are themselves visible at that verbosity. Contexts whose
        level is above the verbosity (i.e. hidden) add no indentation, so
        visible messages never jump several levels at once."""
        return sum(1 for level in self.context_levels if level <= verbosity)

    def __str__(self):
        content = self.content
        if isinstance(content, BaseException):
            content = "".join(format_exception(content))
        elif not isinstance(content, str):
            content = str(content)
        return content


class Filterer(object):
    """Filter messages before logging."""

    __slots__ = ["verbosity"]

    def __init__(self, verbosity: int):
        self.verbosity = verbosity

    def __call__(self, message: Message):
        return message.level <= self.verbosity


class Formatter(object):
    """Format messages before logging."""

    __slots__ = ["formatter"]

    def __init__(self, formatter: Callable[[Message, int], str]):
        self.formatter = formatter

    def __call__(self, message: Message, depth: int):
        return self.formatter(message, depth)


class Stream(ABC):
    """Log to a stream, such as to the console or to a file."""

    __slots__ = ["filterer", "formatter"]

    def __init__(self, filterer: Filterer, formatter: Formatter):
        self.filterer = filterer
        self.formatter = formatter

    @property
    @abstractmethod
    def stream(self) -> TextIO:
        """Text stream to which messages will be logged after filtering
        and formatting."""

    def log(self, message: Message):
        """Log a message to the stream."""
        if self.filterer(message):
            # Indent by the contexts that are visible at this stream's own
            # verbosity, so each stream's indentation matches what it shows.
            depth = message.depth(self.filterer.verbosity)
            self.stream.write(self.formatter(message, depth))


class ConsoleStream(Stream):
    """Log to the console's stderr stream."""

    @property
    def stream(self):
        return stderr


class FileStream(Stream):
    """Log to a file."""

    __slots__ = ["file_path", "_file"]

    def __init__(self, file_path: str | Path, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.file_path = Path(file_path)
        self._file = None

    @property
    def stream(self):
        if self._file is None:
            # Create and open the file if it does not already exist.
            self.file_path.parent.mkdir(exist_ok=True, parents=True)
            self._file = open(self.file_path, "a")
        return self._file

    def close(self):
        """Close the file stream."""
        if self._file is not None:
            self._file.close()
            self._file = None


def indent_content(content: str, indent: str):
    """Indent the first line and every continuation line of content."""
    return indent + content.replace("\n", "\n" + indent)


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


def format_console_plain(message: Message, depth: int = 0):
    """Format a message to log on the console without color."""
    # The level and process ID form a fixed-width prefix so that columns
    # stay aligned; the message text is indented by nesting depth and wrapped
    # to one less than the terminal width, with continuation lines aligned
    # under the message column.
    prefix = f"{message.level.name: <8}{os.getpid()}  "
    indent = INDENT * depth
    first_indent = prefix + indent
    cont_indent = " " * len(prefix) + indent
    width = max(shutil.get_terminal_size().columns - 1, len(cont_indent) + 1)
    return wrap_console(str(message), first_indent, cont_indent, width) + "\n"


class AnsiCode(object):
    """Format text with ANSI codes."""

    # Control codes.
    START = "\033["  # Indicate the start of an ANSI format code.
    END = "m"  # Indicate the end of an ANSI format code.
    RESET = 0  # Reset any existing formatting.
    # Format codes.
    BOLD = 1

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
# The color of each code can be displayed in a bash terminal as follows:
#   for i in {0..255}; do
#       echo -ne "\033[38;5;${i}m  ${i} "
#   done
LEVEL_COLORS = {
    Level.ERROR: AnsiCode.format_color(160),
    Level.WARNING: AnsiCode.format_color(214),
    Level.INFO: AnsiCode.format_color(28),
    Level.DEBUG: AnsiCode.format_color(69),
    Level.TRACE: AnsiCode.format_color(247),
}


def format_console_color(message: Message, depth: int = 0):
    """Format a message to log on the console with color."""
    # Get the ANSI format codes based on the message's logging level.
    fmt = LEVEL_COLORS.get(message.level, AnsiCode.reset())
    # Wrap the formatted text with ANSI format codes.
    return "".join([fmt, format_console_plain(message, depth), AnsiCode.reset()])


def format_logfile(message: Message, depth: int = 0):
    """Format a message to write into the log file."""
    timestamp = datetime.now().strftime("on %Y-%m-%d at %H:%M:%S.%f")
    indent = INDENT * depth
    body = indent_content(str(message), indent)
    return f"LOGMSG> {message.level.name} {os.getpid()} {timestamp}\n{body}\n\n"


class TaskContext(object):
    """Context manager that delimits a task in the log.

    On entry, it logs ``Began {name}`` and pushes its level onto the logger's
    context stack so that messages logged inside the block are indented. On
    normal exit, it pops the context back off and logs ``Ended {name}``. If the
    block exits because of an exception, the context is still popped, but the
    ``Ended`` message is skipped (the task did not finish successfully)."""

    __slots__ = ["_logger", "_level", "_name"]

    def __init__(self, logger: "Logger", level: Level, name: object):
        self._logger = logger
        self._level = level
        self._name = name

    def __enter__(self):
        self._logger._log(self._level, f"Began {self._name}")
        self._logger.context_levels.append(self._level)
        return self._logger

    def __exit__(self, exc_type, exc_value, traceback):
        self._logger.context_levels.pop()
        if exc_type is None:
            self._logger._log(self._level, f"Ended {self._name}")
        return False


class LevelLogger(object):
    """Log messages at one logging level.

    Calling the object logs its argument verbatim at this level. The
    ``begin`` method returns a context manager that brackets a task with
    ``Began``/``Ended`` messages and indents the messages logged inside it."""

    __slots__ = ["_logger", "_level"]

    def __init__(self, logger: "Logger", level: Level):
        self._logger = logger
        self._level = level

    def __call__(self, content: object):
        self._logger._log(self._level, content)

    def begin(self, name: object):
        """Context manager that logs ``Began {name}``/``Ended {name}`` and
        indents the messages logged inside the block."""
        return TaskContext(self._logger, self._level, name)


class Logger(object):
    """Log messages to the console and to files."""

    __slots__ = ["console_stream", "file_stream", "exit_on_error", "context_levels"] + [
        level.name.lower() for level in Level
    ]

    def __init__(
        self,
        console_stream: ConsoleStream | None = None,
        file_stream: FileStream | None = None,
        exit_on_error: bool = DEFAULT_EXIT_ON_ERROR,
    ):
        self.console_stream = console_stream
        self.file_stream = file_stream
        self.exit_on_error = exit_on_error
        # Levels of the currently open task contexts, used to indent messages.
        self.context_levels = []
        # Create one callable LevelLogger per logging level, e.g. so that
        # self.trace(content) logs content at the TRACE level.
        for level in Level:
            setattr(self, level.name.lower(), LevelLogger(self, level))

    def _log(self, level: Level, content: object):
        """Create and log a message to the stream(s)."""
        message = Message(level, content, self.context_levels)
        if level == Level.ERROR and self.exit_on_error:
            if isinstance(content, BaseException):
                raise content
            raise RuntimeError(str(message))
        if self.console_stream is not None:
            self.console_stream.log(message)
        if self.file_stream is not None:
            self.file_stream.log(message)


logger = Logger()


LoggerConfig = namedtuple(
    "LoggerConfig", ["verbosity", "log_file_path", "log_color", "exit_on_error"]
)


def erase_config():
    """Erase the existing logger configuration."""
    logger.console_stream = None
    if logger.file_stream is not None:
        logger.file_stream.close()
    logger.file_stream = None
    logger.exit_on_error = DEFAULT_EXIT_ON_ERROR
    # Clear the context stack so that no indentation leaks across runs and
    # so that each parallel worker starts with no enclosing contexts.
    logger.context_levels = []


def set_config(
    verbosity: int = 0,
    log_file_path: str | Path | None = None,
    log_color: bool = True,
    exit_on_error: bool = DEFAULT_EXIT_ON_ERROR,
):
    """Configure the main logger with handlers and verbosity.

    Parameters
    ----------
    verbosity: int = 0
        Verbosity level; messages at or below this level are shown.
    log_file_path: str | Path | None
        Path of the log file; logging to a file is disabled if None.
    log_color: bool = True
        Whether to use ANSI color codes in console output.
    exit_on_error: bool
        Whether to raise an exception (or exit) when an error is logged.
    """
    # Erase any existing configuration.
    erase_config()
    # Set up logger.
    logger.console_stream = ConsoleStream(
        Filterer(verbosity),
        Formatter(format_console_color if log_color else format_console_plain),
    )
    if log_file_path is not None:
        logger.file_stream = FileStream(
            log_file_path, Filterer(FILE_VERBOSITY), Formatter(format_logfile)
        )
    logger.exit_on_error = exit_on_error


def get_config():
    """Get the configuration parameters of a logger."""
    if logger.console_stream is not None:
        verbosity = logger.console_stream.filterer.verbosity
        log_color = logger.console_stream.formatter.formatter is format_console_color
    else:
        verbosity = DEFAULT_VERBOSITY
        log_color = DEFAULT_COLOR
    if logger.file_stream is not None:
        log_file_path = logger.file_stream.file_path
    else:
        log_file_path = None
    return LoggerConfig(
        verbosity=verbosity,
        log_file_path=log_file_path,
        log_color=log_color,
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
    """After the function exits, restore the logging configuration that
    was in place before the function ran."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        config = get_config()
        try:
            return func(*args, **kwargs)
        finally:
            set_config(**config._asdict())

    return wrapper


# Initialize logging so that it will run by default if using the API.
set_config()
