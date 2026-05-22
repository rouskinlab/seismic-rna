from abc import ABC

from .batch import FilterReadBatch
from ..core import path
from ..core.io import ReadBatchIO, RegFileIO, RegBrickleIO


class FilterFile(path.HasRegFilePath, ABC):

    @classmethod
    def get_step(cls):
        return path.FILTER_STEP


class FilterIO(FilterFile, RegFileIO, ABC):
    pass


class FilterBatchIO(FilterReadBatch, ReadBatchIO, RegBrickleIO, FilterIO):

    @classmethod
    def get_file_seg_type(cls):
        return path.FilterBatSeg
