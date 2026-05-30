import glob as _glob

from click import Path

_WILDCARD_CHARS = frozenset("*?[")


def _has_wildcards(value: str) -> bool:
    return any(c in value for c in _WILDCARD_CHARS)


class GlobPath(Path):
    """A click.Path that expands glob patterns via Python's glob module."""

    name = "glob_path"

    def convert(self, value, param, ctx):
        s = value if isinstance(value, str) else str(value)
        if _has_wildcards(s):
            matches = sorted(_glob.glob(s, recursive=True))
            if not matches:
                # Raise an error if no files match a glob pattern to be
                # consistent with other options that fail if the path
                # does not exist and to prevent a silent error if a typo
                # in a glob pattern causes it to match no files.
                raise FileNotFoundError(f"No files matched {repr(s)}")
            return tuple(Path.convert(self, m, param, ctx) for m in matches)
        return Path.convert(self, value, param, ctx)


def flatten_glob_results(ctx, param, value):
    """Flatten a multiple=True value where entries may be tuples of paths."""
    if value is None:
        return value
    flat = []
    for entry in value:
        if isinstance(entry, tuple):
            flat.extend(entry)
        else:
            flat.append(entry)
    return tuple(flat)


def expand_ct_pos_5(ctx, param, value):
    """Pair each glob-expanded path with its associated integer position."""
    if value is None:
        return value
    flat = []
    for entry in value:
        path, pos = entry
        if isinstance(path, tuple):
            for p in path:
                flat.append((p, pos))
        else:
            flat.append((path, pos))
    return tuple(flat)
