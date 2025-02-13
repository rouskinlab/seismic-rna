from abc import ABC

from .batch import MaskReadBatch
from ..core import path
from ..core.io import ReadBatchIO, RegFileIO, RegBrickleIO


class MaskFile(path.HasRegFilePath, ABC):

    @classmethod
    def get_step(cls):
        return path.MASK_STEP


class MaskIO(MaskFile, RegFileIO, ABC):
    pass


class MaskBatchIO(MaskReadBatch, ReadBatchIO, RegBrickleIO, MaskIO):

    @classmethod
    def get_file_seg_type(cls):
        return path.MaskBatSeg
