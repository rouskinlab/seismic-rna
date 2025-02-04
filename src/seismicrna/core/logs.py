from abc import ABC, abstractmethod
from collections import namedtuple
from datetime import datetime
from enum import IntEnum
from functools import cache, wraps
from pathlib import Path
from sys import stderr
from traceback import format_exception, format_exception_only
from typing import Callable, Optional, TextIO


class Level(IntEnum):
    """ Level of a logging message. """
    FATAL = -3
    ERROR = -2
    WARNING = -1
    STATUS = 0
    TASK = 1
    ACTION = 2
    ROUTINE = 3
    DETAIL = 4


DEFAULT_COLOR = True
DEFAULT_EXIT_ON_ERROR = False
DEFAULT_VERBOSITY = Level.STATUS
FILE_VERBOSITY = Level.DETAIL
EXC_INFO_VERBOSITY = Level.TASK


class Message(object):
    """ Message with a logging level. """
    __slots__ = ["level", "content"]

    def __init__(self, level: Level, content: object):
        self.level = level
        self.content = content

    def __str__(self):
        content = self.content
        if isinstance(content, BaseException):
            if exc_info():
                formatter = format_exception
            else:
                formatter = format_exception_only
            content = "".join(formatter(content))
        elif not isinstance(content, str):
            content = str(content)
        return content


class Filterer(object):
    """ Filter messages before logging. """
    __slots__ = ["verbosity"]

    def __init__(self, verbosity: int):
        self.verbosity = verbosity

    def __call__(self, message: Message):
        return message.level <= self.verbosity


class Formatter(object):
    """ Filter messages before logging. """
    __slots__ = ["formatter"]

    def __init__(self, formatter: Callable[[Message], str]):
        self.formatter = formatter

    def __call__(self, message: Message):
        return self.formatter(message)


class Stream(ABC):
    """ Log to a stream, such as to the console or to a file. """
    __slots__ = ["filterer", "formatter"]

    def __init__(self, filterer: Filterer, formatter: Formatter):
        self.filterer = filterer
        self.formatter = formatter

    @property
    @abstractmethod
    def stream(self) -> TextIO:
        """ Text stream to which messages will be logged after filtering
        and formating. """

    def log(self, message: Message):
        """ Log a message to the stream. """
        if self.filterer(message):
            self.stream.write(self.formatter(message))


class ConsoleStream(Stream):
    """ Log to the console's stderr stream. """

    @property
    def stream(self):
        return stderr


class FileStream(Stream):
    """ Log to a file. """
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
        """ Close the file stream. """
        if self._file is not None:
            self._file.close()
            self._file = None


def format_console_plain(message: Message):
    """ Format a message to log on the console without color. """
    return f"{message.level.name: <8}{message}\n"


class AnsiCode(object):
    """ Format text with ANSI codes. """
    # Control codes.
    START = "\033["  # Indicate the start of an ANSI format code.
    END = "m"  # Indicate the end of an ANSI format code.
    RESET = 0  # Reset any existing formatting.
    # Format codes.
    BOLD = 1

    @classmethod
    @cache
    def format_color(cls, color: int):
        """ Make a format string for one 256-color code. """
        if not 0 <= color < 256:
            raise ValueError(f"Invalid ANSI 256-color code: {color}")
        # 38 means set the foreground color (i.e. of the text itself).
        # 5 means use 256-color mode.
        # The next number is the color in 256-color mode.
        return f"{cls.START}38;5;{color}{cls.END}"

    @classmethod
    @cache
    def format(cls, code: int):
        """ Make a format string for one ANSI code. """
        return f"{cls.START}{code}{cls.END}"

    @classmethod
    @cache
    def reset(cls):
        """ Convenience function to end formatting. """
        return cls.format(cls.RESET)


# Level-specific color formatting.
# The color of each code can be displayed in a bash terminal as follows:
#   for i in {0..255}; do
#       echo -ne "\033[38;5;${i}m  ${i} "
#   done
LEVEL_COLORS = {
    Level.FATAL: "".join([AnsiCode.format_color(198),
                          AnsiCode.format(AnsiCode.BOLD)]),
    Level.ERROR: AnsiCode.format_color(160),
    Level.WARNING: AnsiCode.format_color(214),
    Level.STATUS: AnsiCode.format_color(28),
    Level.TASK: AnsiCode.format_color(38),
    Level.ACTION: AnsiCode.format_color(69),
    Level.ROUTINE: AnsiCode.format_color(147),
    Level.DETAIL: AnsiCode.format_color(247),
}


def format_console_color(message: Message):
    """ Format a message to log on the console with color. """
    # Get the ANSI format codes based on the message's logging level.
    fmt = LEVEL_COLORS.get(message.level, AnsiCode.reset())
    # Wrap the formatted text with ANSI format codes.
    return "".join([fmt, format_console_plain(message), AnsiCode.reset()])


def format_logfile(message: Message):
    """ Format a message to write into the log file. """
    timestamp = datetime.now().strftime("on %Y-%m-%d at %H:%M:%S.%f")
    return f"LOGMSG> {message.level.name} {timestamp}\n{message}\n\n"


class Logger(object):
    """ Log messages to the console and to files. """
    __slots__ = ["console_stream", "file_stream", "exit_on_error"]

    def __init__(self,
                 console_stream: ConsoleStream | None = None,
                 file_stream: FileStream | None = None,
                 exit_on_error: bool = DEFAULT_EXIT_ON_ERROR):
        self.console_stream = console_stream
        self.file_stream = file_stream
        self.exit_on_error = exit_on_error

    def _log(self, level: Level, content: object):
        """ Create and log a message to the stream(s). """
        message = Message(level, content)
        if level <= Level.ERROR and self.exit_on_error:
            if isinstance(content, BaseException):
                raise content
            raise RuntimeError(str(message))
        if self.console_stream is not None:
            self.console_stream.log(message)
        if self.file_stream is not None:
            self.file_stream.log(message)

    def fatal(self, content: object):
        self._log(Level.FATAL, content)

    def error(self, content: object):
        self._log(Level.ERROR, content)

    def warning(self, content: object):
        self._log(Level.WARNING, content)

    def status(self, content: object):
        self._log(Level.STATUS, content)

    def task(self, content: object):
        self._log(Level.TASK, content)

    def action(self, content: object):
        self._log(Level.ACTION, content)

    def routine(self, content: object):
        self._log(Level.ROUTINE, content)

    def detail(self, content: object):
        self._log(Level.DETAIL, content)


logger = Logger()


LoggerConfig = namedtuple("LoggerConfig",
                          ["verbosity",
                           "log_file_path",
                           "log_color",
                           "exit_on_error"])


def erase_config():
    """ Erase the existing logger configuration. """
    logger.console_stream = None
    logger.file_stream = None
    logger.exit_on_error = DEFAULT_EXIT_ON_ERROR


def set_config(verbosity: int = 0,
               log_file_path: str | Path | None = None,
               log_color: bool = True,
               exit_on_error: bool = DEFAULT_EXIT_ON_ERROR):
    """ Configure the main logger with handlers and verbosity. """
    # Erase any existing configuration.
    erase_config()
    # Set up logger.
    logger.console_stream = ConsoleStream(Filterer(verbosity),
                                          Formatter(format_console_color
                                                    if log_color
                                                    else format_console_plain))
    if log_file_path is not None:
        logger.file_stream = FileStream(log_file_path,
                                        Filterer(FILE_VERBOSITY),
                                        Formatter(format_logfile))
    logger.exit_on_error = exit_on_error


def get_config():
    """ Get the configuration parameters of a logger. """
    if logger.console_stream is not None:
        verbosity = logger.console_stream.filterer.verbosity
        log_color = (logger.console_stream.formatter.formatter
                     is format_console_color)
    else:
        verbosity = DEFAULT_VERBOSITY
        log_color = DEFAULT_COLOR
    if logger.file_stream is not None:
        log_file_path = logger.file_stream.file_path
    else:
        log_file_path = None
    return LoggerConfig(verbosity=verbosity,
                        log_file_path=log_file_path,
                        log_color=log_color,
                        exit_on_error=logger.exit_on_error)


def exc_info():
    """ Whether to log exception information. """
    return get_config().verbosity >= EXC_INFO_VERBOSITY


def log_exceptions(default: Optional[Callable]):
    """ If any exception occurs, catch it and return the default. """

    def decorator(func: Callable):

        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as error:
                logger.fatal(error)
                return default() if default is not None else None

        return wrapper

    return decorator


def restore_config(func: Callable):
    """ After the function exits, restore the logging configuration that
    was in place before the function ran. """

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
