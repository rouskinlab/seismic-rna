from abc import ABC

import numpy as np

from .batch import FilterMutsBatch, FilterReadBatch
from ..core import path
from ..core.io import ReadBatchIO, RegFileIO, RegBrickleIO
from ..core.seq import Region


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

    def __init__(self,
                 *args,
                 region: Region | None = None,
                 seg_end5s: np.ndarray | None = None,
                 seg_end3s: np.ndarray | None = None,
                 muts: dict | None = None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.region = region
        self.seg_end5s = seg_end5s
        self.seg_end3s = seg_end3s
        self.muts = muts

    def to_muts_batch(self) -> FilterMutsBatch:
        """ Build a FilterMutsBatch from stored self-contained data. """
        return FilterMutsBatch(batch=self.batch,
                               read_nums=self.read_nums,
                               region=self.region,
                               seg_end5s=self.seg_end5s,
                               seg_end3s=self.seg_end3s,
                               muts=self.muts,
                               sanitize=False)
