from hashlib import md5
from logging import getLogger
from pathlib import Path

logger = getLogger(__name__)


def digest_file(file_path: Path) -> str:
    """
    Compute the checksum of a file.

    Parameters
    ----------
    file_path: Path
        Path of the file on which to compute the checksum. Can be
        any type that the open() function recognizes as a path.

    Returns
    -------
    str
        Checksum of the file (in hexadecimal)
    """
    with open(file_path, "rb") as f:
        digest = md5(f.read()).hexdigest()
    logger.debug(f"Computed MD5 digest of {file_path}: {digest}")
    return digest
