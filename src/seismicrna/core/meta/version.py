"""
Version information for SEISMIC-RNA
"""

from logging import getLogger
from typing import Any

logger = getLogger(__name__)

__version__ = "0.9.0"

VERSION_DELIM = '.'


def parse_version(version: str = __version__):
    try:
        major, minor, patch = map(int, version.split(VERSION_DELIM))
    except ValueError:
        raise ValueError(f"Malformatted version number: {version}")
    return major, minor, patch


MAJOR, MINOR, PATCH = parse_version()


def format_version(major: int = MAJOR, minor: int = MINOR, patch: int = PATCH):
    return VERSION_DELIM.join(map(str, (major, minor, patch)))


def compatible_versions(version1: str, version2: str = __version__):
    match = 0
    for x1, x2 in zip(parse_version(version1),
                      parse_version(version2),
                      strict=True):
        if match == 0:
            if x1 == x2:
                match = x1
            else:
                return False
    return True


def check_compatibility(version: str,
                        what: Any = "An item",
                        strict: bool = True):
    if not compatible_versions(version):
        state = "is not" if strict else "may not be"
        msg = (f"{what} is labeled with version {version} of SEISMIC-RNA, "
               f"which {state} compatible with your version ({__version__})")
        if strict:
            raise RuntimeError(msg)
        logger.warning(msg)
