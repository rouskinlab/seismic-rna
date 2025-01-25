from pathlib import Path
from typing import Iterable

from ..core import path
from ..core.rna import from_ct


def find_ct_files(files: Iterable[str | Path]):
    """ Yield a file for each given file/directory of a table. """
    yield from path.find_files_chain(map(Path, files), [path.ConnectTableSeg])


def load_ct_structs(files: Iterable[str | Path]):
    """ Yield an RNA structure generator for each CT file. """
    for file in find_ct_files(files):
        yield from_ct(file)
