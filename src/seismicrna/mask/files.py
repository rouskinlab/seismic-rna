from abc import ABC

from .batch import MaskBatch
from ..core import path
from ..core.cmd import CMD_MASK
from ..core.files import PickleSectFile


class MaskFile(PickleSectFile, ABC):

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: CMD_MASK}


class MaskBatchFile(MaskBatch, MaskFile):

    @classmethod
    def file_seg_type(cls):
        return path.MaskBatSeg
