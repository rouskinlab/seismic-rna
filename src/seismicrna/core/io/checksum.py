from hashlib import sha512
from pathlib import Path


class BadChecksumError(ValueError):
    """ A file or piece of data has the wrong checksum. """


def calc_sha512_bytes(data: bytes):
    """
    Calculate the SHA512 checksum of a bytestring.

    Parameters
    ----------
    data:
        Bytes object.

    Returns
    -------
    SHA512 checksum as a hexadecimal string.
    """
    if not isinstance(data, bytes):
        raise ValueError(f"Expected a bytes object but got {type(data)}")

    return sha512(data).hexdigest()


def calc_sha512_path(path: str | Path, chunk_size=8192):
    """
    Calculate the SHA512 checksum of a file.

    Parameters
    ----------
    path:
        File path as a string/Path object.
    chunk_size: 
        Size of chunks to read when processing a file.

    Returns
    -------
    SHA512 checksum as a hexadecimal string.
    """
    if not isinstance(path, (str, Path)):
        raise ValueError(f"Expected a str or pathlib.Path but got {type(path)}")

    hasher = sha512()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)

    return hasher.hexdigest()
