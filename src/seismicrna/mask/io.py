from abc import ABC

from .batch import MaskReadBatch
from ..core import path
from ..core.io import ReadBatchIO, SectIO


class MaskIO(SectIO, ABC):

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: path.CMD_MSK_DIR}


class MaskBatchIO(ReadBatchIO, MaskIO, MaskReadBatch):

    @classmethod
    def file_seg_type(cls):
        return path.MaskBatSeg
