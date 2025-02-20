import warnings
from functools import cached_property

from .brickle import BrickleIO
from .file import RefFileIO
from .. import path
from ..seq import DNA


class RefseqIO(RefFileIO, BrickleIO):

    @classmethod
    def get_file_seg_type(cls):
        return path.RefseqFileSeg

    @classmethod
    def get_step(cls):
        return path.RELATE_STEP

    def __init__(self, *args, refseq: DNA, **kwargs):
        warnings.warn(
            "RefseqIO is deprecated and will be removed in version 0.25. "
            "Run seismic migrate to convert output files to the new format "
            "and suppress this warning.",
            FutureWarning
        )
        super().__init__(*args, **kwargs)
        self._s = refseq.compress()

    @cached_property
    def refseq(self):
        return self._s.decompress()
