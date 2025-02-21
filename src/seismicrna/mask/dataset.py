from abc import ABC
from functools import cached_property

import numpy as np

from .batch import MaskMutsBatch, apply_mask
from .io import MaskBatchIO
from .report import MaskReport, JoinMaskReport
from ..core.dataset import (LoadedDataset,
                            LoadFunction,
                            MergedUnbiasDataset,
                            MultistepDataset,
                            UnbiasDataset)
from ..core.join import (BATCH_NUM,
                         READ_NUMS,
                         SEG_END5S,
                         SEG_END3S,
                         MUTS,
                         JoinMutsDataset)
from ..core.rel import RelPattern
from ..core.report import (CountMutsF,
                           CountRefsF,
                           MinMutGapF,
                           PosKeptF,
                           RefF,
                           RegF,
                           End5F,
                           End3F,
                           QuickUnbiasF,
                           QuickUnbiasThreshF,
                           JoinedClustersF)
from ..core.seq import Region
from ..relate.batch import RelateMutsBatch
from ..relate.dataset import AverageDataset, load_relate_dataset


class MaskDataset(AverageDataset, ABC):
    """ Dataset of masked data. """


class MaskReadDataset(MaskDataset, LoadedDataset, UnbiasDataset):
    """ Load batches of masked data. """

    @classmethod
    def get_report_type(cls):
        return MaskReport

    @classmethod
    def get_batch_type(cls):
        return MaskBatchIO

    @property
    def min_mut_gap(self):
        return self.report.get_field(MinMutGapF)

    @property
    def quick_unbias(self):
        return self.report.get_field(QuickUnbiasF)

    @property
    def quick_unbias_thresh(self):
        return self.report.get_field(QuickUnbiasThreshF)

    @property
    def pos_kept(self):
        """ Positions kept after masking. """
        return self.report.get_field(PosKeptF)

    @cached_property
    def pattern(self):
        return RelPattern(self.report.get_field(CountMutsF),
                          self.report.get_field(CountRefsF))


class MaskMutsDataset(MaskDataset, MultistepDataset, UnbiasDataset):
    """ Chain mutation data with masked reads. """

    MASK_NAME = "mask"

    @classmethod
    def get_dataset1_load_func(cls):
        return load_relate_dataset

    @classmethod
    def get_dataset2_type(cls):
        return MaskReadDataset

    @property
    def pattern(self):
        return self.dataset2.pattern

    @property
    def min_mut_gap(self):
        return getattr(self.dataset2, "min_mut_gap")

    @property
    def quick_unbias(self):
        return getattr(self.dataset2, "quick_unbias")

    @property
    def quick_unbias_thresh(self):
        return getattr(self.dataset2, "quick_unbias_thresh")

    @cached_property
    def region(self):
        region = Region(ref=self.report.get_field(RefF),
                        seq=self.refseq,
                        name=self.report.get_field(RegF),
                        end5=self.report.get_field(End5F),
                        end3=self.report.get_field(End3F))
        region.add_mask(self.MASK_NAME,
                        getattr(self.dataset2, "pos_kept"),
                        complement=True)
        return region

    def _integrate(self, batch1: RelateMutsBatch, batch2: MaskBatchIO):
        if self.masked_read_nums is not None:
            read_nums = np.setdiff1d(batch2.read_nums,
                                     self.masked_read_nums.get(batch2.batch),
                                     assume_unique=True)
        else:
            read_nums = batch2.read_nums
        return apply_mask(batch1,
                          read_nums,
                          self.region,
                          sanitize=False)


class JoinMaskMutsDataset(MaskDataset, JoinMutsDataset, MergedUnbiasDataset):

    @classmethod
    def get_report_type(cls):
        return JoinMaskReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_mask_dataset

    @classmethod
    def get_batch_type(cls):
        return MaskMutsBatch

    @classmethod
    def name_batch_attrs(cls):
        return [BATCH_NUM, READ_NUMS, SEG_END5S, SEG_END3S, MUTS]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        clusts = self.report.get_field(JoinedClustersF, missing_ok=True)
        if clusts is not None:
            raise TypeError(f"{self} has no clusters, but got {clusts}")


load_mask_dataset = LoadFunction(MaskMutsDataset, JoinMaskMutsDataset)
