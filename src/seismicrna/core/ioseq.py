from functools import cached_property

from . import path
from .clicmd import CMD_REL
from .iofile import SavedRef, SavedBrickle
from .seq import DNA


class SavedRefseq(SavedRef, SavedBrickle):

    @classmethod
    def file_seg_type(cls):
        return path.RefseqFileSeg

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: CMD_REL}

    def __init__(self, *args, refseq: DNA, **kwargs):
        super().__init__(*args, **kwargs)
        self._s = refseq.compress()

    @cached_property
    def refseq(self):
        return self._s.decompress()
