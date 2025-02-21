from hashlib import sha512
from pathlib import Path

class BadChecksumError(ValueError):
    """ A file or piece of data has the wrong checksum. """


def calc_sha512_digest(target: str | Path | bytes, chunk_size=8192):
    """
    Calculate the SHA512 checksum of a file (given a path) or a bytestring.
    
    Parameters
    ----------
    target: 
        File path as a string/Path object, or a bytes object.
    chunk_size: 
        Size of chunks to read when processing a file.
    
    Returns
    -------
    SHA512 checksum as a hexadecimal string.
    """
    hasher = sha512()
    if isinstance(target, bytes):
        # Directly hash the bytes if input is a bytes object
        hasher.update(target)
    elif isinstance(target, (str, Path)):
        # Assume it's a file path and read in chunks
        with open(target, "rb") as f:
            while chunk := f.read(chunk_size):
                hasher.update(chunk)
    else:
        raise ValueError("Target must be a file path or a bytes object.")

    return hasher.hexdigest()