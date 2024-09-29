from collections import namedtuple
from datetime import datetime
from enum import IntEnum
from functools import cache, wraps
from sys import stderr
from traceback import format_exception, format_exception_only
from typing import Callable, Optional, TextIO


class Level(IntEnum):
    """ Level of a logging message. """
    FATAL = -3
    ERROR = -2
    WARNING = -1
    STATUS = 0
    PROCESS = 1
    ROUTINE = 2
    DETAIL = 3


class Message(object):
    """ Message with a logging level. """
    __slots__ = ["level", "content", "args", "kwargs"]

    def __init__(self, level: Level, content: object, *args, **kwargs):
        self.level = level
        self.content = content
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        content = self.content
        if isinstance(content, BaseException):
            content = "".join(format_exception(content)
                              if exc_info()
                              else format_exception_only(content))
        else:
            if not isinstance(content, str):
                content = str(content)
            if self.args or self.kwargs:
                content = content.format(*self.args, **self.kwargs)
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


class Stream(object):
    """ Direct logging output to one stream, such as to the console or
    to a file. """
    __slots__ = ["stream", "filterer", "formatter"]

    def __init__(self,
                 stream: TextIO,
                 filterer: Filterer,
                 formatter: Formatter):
        self.stream = stream
        self.filterer = filterer
        self.formatter = formatter

    def log(self, message: Message):
        """ Log a message to the stream. """
        if self.filterer(message):
            self.stream.write(self.formatter(message))


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
# The color of each code can be visualized in a terminal as follows:
#   for i in {0..255}; do
#       echo -ne "\033[38;5;${i}m  ${i} "
#   done
LEVEL_COLORS = {
    Level.FATAL: "".join([AnsiCode.format_color(198),
                          AnsiCode.format(AnsiCode.BOLD)]),
    Level.ERROR: AnsiCode.format_color(160),
    Level.WARNING: AnsiCode.format_color(214),
    Level.STATUS: AnsiCode.format_color(28),
    Level.PROCESS: AnsiCode.format_color(75),
    Level.ROUTINE: AnsiCode.format_color(99),
    Level.DETAIL: AnsiCode.format_color(242),
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
    __slots__ = ["console_stream", "file_stream", "raise_on_error"]

    def __init__(self,
                 console_stream: Stream | None = None,
                 file_stream: Stream | None = None,
                 raise_on_error: bool = False):
        self.console_stream = console_stream
        self.file_stream = file_stream
        self.raise_on_error = raise_on_error

    def _log_message(self, message: Message):
        """ Log a message to the stream(s). """
        if self.file_stream is not None:
            self.file_stream.log(message)
        if self.console_stream is not None:
            self.console_stream.log(message)

    def _log_object(self, level: Level, content: object, *args, **kwargs):
        """ Create and log a message to the stream(s). """
        message = Message(level, content, *args, **kwargs)
        if self.raise_on_error and level <= Level.ERROR:
            if isinstance(content, BaseException):
                raise content
            raise RuntimeError(str(message))
        self._log_message(message)

    def fatal(self, text: object, *args, **kwargs):
        self._log_object(Level.FATAL, text, *args, **kwargs)

    def error(self, text: object, *args, **kwargs):
        self._log_object(Level.ERROR, text, *args, **kwargs)

    def warning(self, text: object, *args, **kwargs):
        self._log_object(Level.WARNING, text, *args, **kwargs)

    def status(self, text: object, *args, **kwargs):
        self._log_object(Level.STATUS, text, *args, **kwargs)

    def process(self, text: object, *args, **kwargs):
        self._log_object(Level.PROCESS, text, *args, **kwargs)

    def routine(self, text: object, *args, **kwargs):
        self._log_object(Level.ROUTINE, text, *args, **kwargs)

    def detail(self, text: object, *args, **kwargs):
        self._log_object(Level.DETAIL, text, *args, **kwargs)


logger = Logger()

DEFAULT_COLOR = True
DEFAULT_RAISE = False
DEFAULT_VERBOSITY = Level.STATUS
FILE_VERBOSITY = Level.ROUTINE
EXC_INFO_VERBOSITY = Level.ROUTINE


def erase_config():
    """ Erase the existing logger configuration. """
    logger.console_stream = None
    logger.file_stream = None
    logger.raise_on_error = DEFAULT_RAISE


def set_config(verbosity: int = 0,
               log_file: TextIO | None = None,
               log_color: bool = True,
               raise_on_error: bool = DEFAULT_RAISE):
    """ Configure the main logger with handlers and verbosity. """
    # Erase any existing configuration.
    erase_config()
    # Set up logger.
    logger.console_stream = Stream(stderr,
                                   Filterer(verbosity),
                                   Formatter(format_console_color
                                             if log_color
                                             else format_console_plain))
    if log_file is not None:
        logger.file_stream = Stream(log_file,
                                    Filterer(FILE_VERBOSITY),
                                    Formatter(format_logfile))
    logger.raise_on_error = raise_on_error


LoggerConfig = namedtuple("LoggerConfig",
                          ["verbosity",
                           "log_file",
                           "log_color",
                           "raise_on_error"])


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
        log_file = logger.file_stream.stream
    else:
        log_file = None
    return LoggerConfig(verbosity=verbosity,
                        log_file=log_file,
                        log_color=log_color,
                        raise_on_error=logger.raise_on_error)


def exc_info():
    """ Whether to log exception information. """
    return get_config().verbosity >= EXC_INFO_VERBOSITY


def log_exceptions(logging_method: Callable, default: Optional[Callable]):
    """ If any exception occurs, catch it and return an empty list. """
    try:
        method_logger = getattr(logging_method, "__self__")
    except AttributeError:
        raise TypeError("logging_method is not an instance method")
    if method_logger is not logger:
        raise TypeError(f"logging_method must be a method of {logger}, "
                        f"but it is a method of {method_logger}")
    if logging_method not in [method_logger.fatal,
                              method_logger.error,
                              method_logger.warning]:
        raise ValueError("logging_method must be fatal, error, or warning, "
                         f"but got {logging_method}")

    def decorator(func: Callable):

        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as error:
                logging_method(error)
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
