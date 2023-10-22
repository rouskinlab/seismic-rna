from abc import ABC

from .batch import ClustReadBatch
from ..core import path
from ..core.io import ReadBatchIO, SectIO


class ClustIO(SectIO, ABC):

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: path.CMD_CLS_DIR}


class ClustBatchIO(ReadBatchIO, ClustIO, ClustReadBatch):

    @classmethod
    def file_seg_type(cls):
        return path.ClustBatSeg
