from hashlib import sha512
from pathlib import Path

from ..validate import require_isinstance


class BadChecksumError(ValueError):
    """ A file or piece of data has the wrong checksum. """


def calc_sha512_bytes(data: bytes):
    """ Calculate the SHA-512 checksum of a bytestring.

    Parameters
    ----------
    data: bytes
        Data of which to calculate the checksum

    Returns
    -------
    SHA-512 checksum as a hexadecimal string.
    """
    require_isinstance("data", data, bytes)
    return sha512(data).hexdigest()


def calc_sha512_path(path: str | Path, chunk_size: int = 8192):
    """ Calculate the SHA-512 checksum of a file.

    Parameters
    ----------
    path: str | Path
        File path
    chunk_size: int
        Size of chunks to read when processing a file

    Returns
    -------
    SHA-512 checksum as a hexadecimal string.
    """
    require_isinstance("path", path, (str, Path))
    hasher = sha512()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()
