from __future__ import annotations
from abc import ABC
from functools import cached_property

from ..core import path
from ..core.io.brickle import RegBrickleIO
from ..core.io.file import RegFileIO
from ..core.seq.xna import DNA


class DuplexFile(path.HasRegFilePath, ABC):
    @classmethod
    def get_step(cls):
        return path.DUPLEX_STEP


class DuplexIO(DuplexFile, RegFileIO, ABC):
    pass


class DuplexRefseqIO(RegBrickleIO, DuplexFile, RegFileIO):
    """The fused reference sequence of a duplex, stored like any other
    reference sequence (a brickle file referenced by checksum), in the
    duplex's region directory alongside its table and report."""

    @classmethod
    def get_file_seg_type(cls):
        return path.RefseqFileSeg

    def __init__(self, *args, refseq: DNA, **kwargs):
        super().__init__(*args, **kwargs)
        self._s = refseq.compress()

    @cached_property
    def refseq(self):
        return self._s.decompress()
