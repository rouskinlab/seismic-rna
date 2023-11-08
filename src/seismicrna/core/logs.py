"""
Core -- Logging Module

Purpose
-------
Central manager of logging.
"""

import logging
from functools import cache

MAX_VERBOSE = 2
MAX_QUIET = 2
FILE_MSG_FORMAT = "LOGMSG>\t%(asctime)s\t%(name)s\t%(levelname)s\n%(message)s\n"
STREAM_MSG_FORMAT = "%(levelname)s\t%(message)s"
LEVELS = {(2, 0): logging.DEBUG,
          (1, 0): logging.INFO,
          (0, 0): logging.WARNING,
          (0, 1): logging.ERROR,
          (0, 2): logging.CRITICAL}


class AnsiCode(object):
    """ Format text with ANSI codes. """
    END = 0
    BOLD = 1
    ULINE = 4
    RED = 91
    GREEN = 92
    YELLOW = 93
    BLUE = 94
    PURPLE = 95
    CYAN = 96
    CODES = END, BOLD, ULINE, RED, GREEN, YELLOW, BLUE, PURPLE, CYAN

    @classmethod
    @cache
    def fmt(cls, code: int):
        """ Format one color code into text. """
        if code not in cls.CODES:
            raise ValueError(f"Invalid ANSI color code: {code}")
        return f"\033[{code}m"

    @classmethod
    def end(cls):
        """ Convenience function to end formatting. """
        return cls.fmt(cls.END)

    @classmethod
    def wrap(cls, text: str, *codes: int, end: bool = True):
        """ Wrap text with ANSI color code(s). """
        return f"{''.join(map(cls.fmt, codes))}{text}{cls.end() if end else ''}"


class ColorFormatter(logging.Formatter):
    ansi_codes = {
        logging.DEBUG: (AnsiCode.BLUE,),
        logging.INFO: (AnsiCode.CYAN,),
        logging.WARNING: (AnsiCode.YELLOW,),
        logging.ERROR: (AnsiCode.RED,),
        logging.CRITICAL: (AnsiCode.PURPLE, AnsiCode.BOLD),
    }

    def format(self, record: logging.LogRecord) -> str:
        """ Log the message in color by adding an ANSI color escape code
        to the beginning and a color stopping code to the end. """
        # Get the ANSI format codes based on the record's logging level.
        # Wrap the formatted text with ANSI format codes.
        return AnsiCode.wrap(super().format(record),
                             *self.ansi_codes.get(record.levelno,
                                                  (AnsiCode.END,)))


def get_top_logger():
    """ Return the top-level logger. """
    if __name__ != (expect_name := "seismicrna.core.logs"):
        raise ValueError(
            f"{__file__} is named '{__name__}' (expected '{expect_name}')")
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

    # Limit verbose and quiet to 2.
    if verbose > MAX_VERBOSE:
        logging.warning(f"Setting 'verbose' to {MAX_VERBOSE} (got {verbose})")
        verbose = MAX_VERBOSE
    if quiet > MAX_QUIET:
        logging.warning(f"Setting 'quiet' to {MAX_QUIET} (got {quiet})")
        quiet = MAX_QUIET

    # Set logging level based on verbose and quiet.
    try:
        return LEVELS[verbose, quiet]
    except KeyError:
        get_top_logger().warning(f"Invalid options: verbose={verbose}, "
                                 f"quiet={quiet}. Setting both to 0")
        return get_verbosity()


def config(verbose: int,
           quiet: int,
           log_file: str | None = None,
           log_color: bool = True):
    """ Configure the main logger with handlers and verbosity. """
    # Set up logger.
    logger = get_top_logger()
    logger.setLevel(get_verbosity(verbose=MAX_VERBOSE))
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

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
