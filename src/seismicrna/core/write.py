from pathlib import Path


def need_write(file: Path, force: bool = False):
    """ Determine whether a file must be written.

    Parameters
    ----------
    file: Path
        File for which to check the need for writing.
    force: bool = False
        Force the file to be written, even if it already exists.

    Returns
    -------
    bool
        Whether the file must be written.
    """
    return force or not file.is_file()


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
