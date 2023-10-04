from functools import cached_property

from . import path
from .output import PickleOutput, RefOutput
from .seq import DNA


class RefseqFile(PickleOutput, RefOutput):

    @classmethod
    def file_seg_type(cls):
        return path.RefseqFileSeg

    def __init__(self, *args, refseq: DNA, **kwargs):
        super().__init__(*args, **kwargs)
        self._s = refseq.compress()

    @cached_property
    def refseq(self):
        return self._s.decompress()
