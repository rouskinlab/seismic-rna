from logging import getLogger
from pathlib import Path

logger = getLogger(__name__)


def need_write(file: Path, force: bool = False, warn: bool = True):
    """ Determine whether a file must be written.

    Parameters
    ----------
    file: Path
        File for which to check the need for writing.
    force: bool = False
        Force the file to be written, even if it already exists.
    warn: bool = True
        If the file does not need to be written, then log a warning.

    Returns
    -------
    bool
        Whether the file must be written.
    """
    if force or not file.is_file():
        return True
    if warn:
        logger.warning(f"File exists: {file}")
    return False


def write_mode(force: bool = False, binary: bool = False):
    """ Get the mode in which to open a file for writing.

    Parameters
    ----------
    force: bool = False
        Force the file to be written, truncating the file if it exists.
        If False and the file exists, a FileExistsError will be raised.
    binary: bool = False
        Write the file in binary mode instead of text mode.

    Returns
    -------
    str
        The mode argument for the builtin function `open()`.
    """
    return "".join(["w" if force else "x", "b" if binary else ""])
