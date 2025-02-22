from pathlib import Path

from .logs import logger


def need_write(query: str | Path, force: bool = False, warn: bool = True):
    """ Determine whether a file/directory must be written.

    Parameters
    ----------
    query: str | Path
        File or directory for which to check the need for writing.
    force: bool = False
        Force the query to be written, even if it already exists.
    warn: bool = True
        If the query does not need to be written, then log a warning.

    Returns
    -------
    bool
        Whether the file must be written.
    """
    if not isinstance(query, Path):
        query = Path(query)
    if force or not query.exists():
        return True
    if warn:
        logger.warning(f"{query} exists: use --force to overwrite")
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
