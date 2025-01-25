from abc import ABC

from .batch import MaskReadBatch
from ..core import path
from ..core.io import ReadBatchIO, RegIO


class MaskIO(RegIO, ABC):

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: path.MASK_STEP}


class MaskBatchIO(ReadBatchIO, MaskIO, MaskReadBatch):

    @classmethod
    def file_seg_type(cls):
        return path.MaskBatSeg
