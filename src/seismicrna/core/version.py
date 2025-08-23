import re


__version__ = "0.24.3"


def parse_version(version: str = __version__):
    """ Major and minor versions, patch, and pre-release tag. """
    match = re.match(r"^([0-9]+)[.]([0-9]+)[.]([0-9]+)([a-z]*[0-9]*)$", version)
    if not match:
        raise ValueError(f"Malformatted version: {repr(version)}")
    major, minor, patch = map(int, match.groups()[:3])
    prtag = match.groups()[3]
    return major, minor, patch, prtag


MAJOR, MINOR, PATCH, PRTAG = parse_version()


def format_version(major: int = MAJOR,
                   minor: int = MINOR,
                   patch: int = PATCH,
                   prtag: str = PRTAG):
    return f"{major}.{minor}.{patch}{prtag}"
