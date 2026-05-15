from abc import ABC
from functools import cached_property

import numpy as np

from .batch import MaskMutsBatch, apply_mask
from ..core.arg import MUT_COLLISIONS_MERGE
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
                           MutCollisionsF,
                           PosKeptF,
                           ProbeF,
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
    def mut_collisions(self):
        return self.report.get_field(MutCollisionsF)

    @property
    def quick_unbias(self):
        return self.report.get_field(QuickUnbiasF)

    @property
    def quick_unbias_thresh(self):
        return self.report.get_field(QuickUnbiasThreshF)

    @property
    def probe(self):
        return self.report.get_field(ProbeF)

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
    def mut_collisions(self):
        return getattr(self.dataset2, "mut_collisions")

    @property
    def quick_unbias(self):
        return getattr(self.dataset2, "quick_unbias")

    @property
    def quick_unbias_thresh(self):
        return getattr(self.dataset2, "quick_unbias_thresh")

    @property
    def probe(self):
        return getattr(self.dataset2, "probe")

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
        """ Combine a relate batch with a mask batch into a MaskMutsBatch.

        Parameters
        ----------
        batch1: RelateMutsBatch
            Batch of mutation data from the relate step.
        batch2: MaskBatchIO
            Batch of read numbers that passed the mask filters.

        Returns
        -------
        MaskMutsBatch
            A batch containing only the reads and positions that pass
            the mask, clipped to the dataset's region.
        """
        if self.masked_read_nums is not None:
            read_nums = np.setdiff1d(batch2.read_nums,
                                     self.masked_read_nums.get(batch2.batch),
                                     assume_unique=True)
        else:
            read_nums = batch2.read_nums
        masked_batch = apply_mask(batch1,
                                  read_nums,
                                  self.region,
                                  sanitize=False)
        if self.min_mut_gap > 0 and self.mut_collisions == MUT_COLLISIONS_MERGE:
            return MaskMutsBatch(
                batch=masked_batch.batch,
                read_nums=masked_batch.read_nums,
                region=masked_batch.region,
                seg_end5s=masked_batch.seg_end5s,
                seg_end3s=masked_batch.seg_end3s,
                muts=masked_batch.merge_close_muts(self.pattern,
                                                   self.min_mut_gap),
            )
        return masked_batch


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
