"""
Core -- Logging Module

Purpose
-------
Central manager of logging.
"""

import logging
from collections import namedtuple
from functools import cache, wraps
from typing import Callable, Optional

MAX_VERBOSE = 2
MAX_QUIET = 2
FILE_MSG_FORMAT = "LOGMSG>\t%(asctime)s\t%(name)s\t%(levelname)s\n%(message)s\n"
STREAM_MSG_FORMAT = "%(levelname)s\t%(message)s"
LEVELS = {(2, 0): logging.DEBUG,
          (1, 0): logging.INFO,
          (0, 0): logging.WARNING,
          (0, 1): logging.ERROR,
          (0, 2): logging.CRITICAL}
VERBOSITIES = {level: verbosity for verbosity, level in LEVELS.items()}
DEFAULT_LEVEL = logging.WARNING
DEFAULT_RAISE = False


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
    def fmt_color(cls, color: int):
        """ Make a format string for one 256-color code. """
        if not 0 <= color < 256:
            raise ValueError(f"Invalid ANSI 256-color code: {color}")
        # 38 means set the foreground color (i.e. of the text itself).
        # 5 means use 256-color mode.
        # The next number is the color in 256-color mode.
        return f"{cls.START}38;5;{color}{cls.END}"

    @classmethod
    @cache
    def fmt(cls, code: int):
        """ Make a format string for one ANSI code. """
        return f"{cls.START}{code}{cls.END}"

    @classmethod
    @cache
    def reset(cls):
        """ Convenience function to end formatting. """
        return cls.fmt(cls.RESET)


class ColorFormatter(logging.Formatter):
    # The color of each code can be visualized in a terminal as follows:
    # for i in {0..255}; do
    #     echo -ne "\033[38;5;${i}m  ${i} "
    # done
    ansi_codes = {
        logging.DEBUG: AnsiCode.fmt_color(244),
        logging.INFO: AnsiCode.fmt_color(75),
        logging.WARNING: AnsiCode.fmt_color(214),
        logging.ERROR: AnsiCode.fmt_color(160),
        logging.CRITICAL: "".join([AnsiCode.fmt_color(201),
                                   AnsiCode.fmt(AnsiCode.BOLD)])
    }

    def format(self, record: logging.LogRecord) -> str:
        """ Log the message in color by adding an ANSI color escape code
        to the beginning and a color stopping code to the end. """
        # Get the ANSI format codes based on the record's logging level.
        fmt = self.ansi_codes.get(record.levelno, AnsiCode.reset())
        # Wrap the formatted text with ANSI format codes.
        return "".join([fmt, super().format(record), AnsiCode.reset()])


class RaisableLogger(logging.Logger):

    @staticmethod
    def _handle_error(log_func: Callable, msg: object, *args, **kwargs):
        """ Handle logging an error message. """
        if get_config().raise_on_error:
            if isinstance(msg, BaseException):
                raise msg
            raise RuntimeError(msg)
        log_func(msg, *args, **kwargs)

    def __init__(self, *args, raise_on_error: bool = DEFAULT_RAISE, **kwargs):
        super().__init__(*args, **kwargs)
        self.raise_on_error = raise_on_error

    def error(self, msg: object, *args, **kwargs):
        self._handle_error(super().error, msg, *args, **kwargs)

    def critical(self, msg: object, *args, **kwargs):
        self._handle_error(super().critical, msg, *args, **kwargs)


logging.setLoggerClass(RaisableLogger)


def get_top_logger() -> RaisableLogger:
    """ Return the top-level logger. """
    if __name__ != (expect_name := "seismicrna.core.logs"):
        raise ValueError(f"Expected {__file__} named {repr(expect_name)}, "
                         f"but got {repr(__name__)}")
    top_logger_name = __name__.split(".")[0]
    return logging.getLogger(top_logger_name)


def get_verbosity(verbose: int = 0, quiet: int = 0):
    """ Get the logging level based on the verbose and quiet arguments.

    Parameters
    ----------
    verbose: int [0, 2]
        0 (): Log only warnings and errors
        1 (-v): Also log status updates
        2 (-vv): Also log detailed information (useful for debugging)
    quiet: int [0, 2]
        0 (): Suppress only status updates and detailed information
        1 (-q): Also suppress warnings
        2 (-qq): Also suppress non-critical error messages (discouraged)

    Giving both `verbose` and `quiet` flags causes the verbosity
    to default to `verbose=0`, `quiet=0`.
    """
    logger = get_top_logger()
    # Limit verbose and quiet to 2.
    if verbose > MAX_VERBOSE:
        logger.warning(f"Setting 'verbose' to {MAX_VERBOSE} (got {verbose})")
        verbose = MAX_VERBOSE
    if quiet > MAX_QUIET:
        logger.warning(f"Setting 'quiet' to {MAX_QUIET} (got {quiet})")
        quiet = MAX_QUIET
    # Set logging level based on verbose and quiet.
    try:
        return LEVELS[verbose, quiet]
    except KeyError:
        logger.warning(f"Invalid options: verbose={verbose}, "
                       f"quiet={quiet}. Setting both to 0")
        return get_verbosity()


def erase_config():
    """ Reset the logging configuration to the defaults. """
    logger = get_top_logger()
    logger.setLevel(DEFAULT_LEVEL)
    logger.raise_on_error = DEFAULT_RAISE
    # Need to use logger.handlers.copy() because logger.handlers will be
    # modified by logger.removeHandler(); iterating over logger.handlers
    # itself can therefore fail to remove all handlers.
    for handler in logger.handlers.copy():
        if isinstance(handler, logging.FileHandler):
            handler.close()
        logger.removeHandler(handler)
    if logger.handlers:
        raise RuntimeError(f"Failed to remove all handlers from {logger}; "
                           f"remaining handlers are {logger.handlers}")


def set_config(verbose: int = 0,
               quiet: int = 0,
               log_file: str | None = None,
               log_color: bool = True,
               raise_on_error: bool = DEFAULT_RAISE):
    """ Configure the main logger with handlers and verbosity. """
    # Erase any existing configuration.
    erase_config()
    # Set up logger.
    logger = get_top_logger()
    logger.setLevel(get_verbosity(verbose=MAX_VERBOSE))
    logger.raise_on_error = raise_on_error
    # Add stream handler.
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(get_verbosity(verbose, quiet))
    stream_handler.setFormatter(ColorFormatter(STREAM_MSG_FORMAT) if log_color
                                else logging.Formatter(STREAM_MSG_FORMAT))
    logger.addHandler(stream_handler)
    # Add file handler.
    if log_file is not None:
        file_handler = logging.FileHandler(log_file, "a")
        file_handler.setLevel(get_verbosity(verbose=MAX_VERBOSE))
        file_handler.setFormatter(logging.Formatter(FILE_MSG_FORMAT))
        logger.addHandler(file_handler)


LoggerConfig = namedtuple("LoggerConfig",
                          ["verbose",
                           "quiet",
                           "log_file",
                           "log_color",
                           "raise_on_error"])


def get_config():
    """ Get the configuration parameters of a logger. """
    logger = get_top_logger()
    verbose = 0
    quiet = 0
    log_file = None
    log_color = False
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            log_file = handler.baseFilename
        elif isinstance(handler, logging.StreamHandler):
            verbose, quiet = VERBOSITIES.get(handler.level, (verbose, quiet))
            if isinstance(handler.formatter, ColorFormatter):
                log_color = True
    return LoggerConfig(verbose=verbose,
                        quiet=quiet,
                        log_file=log_file,
                        log_color=log_color,
                        raise_on_error=logger.raise_on_error)


def exc_info():
    """ Whether to log exception information. """
    return get_config().verbose == MAX_VERBOSE


def log_exceptions(logging_method: Callable, default: Optional[Callable]):
    """ If any exception occurs, catch it and return an empty list. """
    try:
        logger = getattr(logging_method, "__self__")
    except AttributeError:
        raise TypeError("logging_method is not an instance method")
    if not isinstance(logger, logging.Logger):
        raise TypeError("logging_method must be instance method of a Logger, "
                        f"but got {logging_method}")
    if logging_method not in [logger.critical, logger.error, logger.warning]:
        raise ValueError("logging_method must be critical, error, or warning, "
                         f"but got {logging_method}")

    def decorator(func: Callable):

        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as error:
                logging_method(error, exc_info=exc_info())
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
