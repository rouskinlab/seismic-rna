from abc import ABC

from .batch import ClusterReadBatch
from ..core import path
from ..core.io import ReadBatchIO, SectIO


class ClusterIO(SectIO, ABC):

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: path.CMD_CLS_DIR}


class ClusterBatchIO(ReadBatchIO, ClusterIO, ClusterReadBatch):

    @classmethod
    def file_seg_type(cls):
        return path.ClustBatSeg
