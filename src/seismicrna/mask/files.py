from abc import ABC

from ..core import path
from ..core.clicmd import CMD_MASK
from ..core.batch import MaskBatch
from ..core.iobatch import SavedBatch
from ..core.iofile import SavedSect


class SavedMask(SavedSect, ABC):

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: CMD_MASK}


class SavedMaskBatch(SavedBatch, SavedMask, MaskBatch):

    @classmethod
    def file_seg_type(cls):
        return path.MaskBatSeg
