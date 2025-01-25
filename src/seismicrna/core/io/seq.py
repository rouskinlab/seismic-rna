from functools import cached_property

from .file import RefIO, BrickleIO
from .. import path
from ..seq import DNA


class RefseqIO(RefIO, BrickleIO):

    @classmethod
    def file_seg_type(cls):
        return path.RefseqFileSeg

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: path.RELATE_STEP}

    def __init__(self, *args, refseq: DNA, **kwargs):
        super().__init__(*args, **kwargs)
        self._s = refseq.compress()

    @cached_property
    def refseq(self):
        return self._s.decompress()
