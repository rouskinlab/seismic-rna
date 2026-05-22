from abc import ABC
from functools import cached_property

import numpy as np

from .batch import FilterMutsBatch, apply_filters
from ..core.arg import MUT_COLLISIONS_MERGE
from .io import FilterBatchIO
from .report import FilterReport, JoinFilterReport
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
from ..idmut.batch import IDmutMutsBatch
from ..idmut.dataset import AverageDataset, load_idmut_dataset


class FilterDataset(AverageDataset, ABC):
    """ Dataset of filtered data. """


class FilterReadDataset(FilterDataset, LoadedDataset, UnbiasDataset):
    """ Load batches of filtered data. """

    @classmethod
    def get_report_type(cls):
        return FilterReport

    @classmethod
    def get_batch_type(cls):
        return FilterBatchIO

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
        """ Positions kept after filtering. """
        return self.report.get_field(PosKeptF)

    @cached_property
    def pattern(self):
        return RelPattern(self.report.get_field(CountMutsF),
                          self.report.get_field(CountRefsF))


class FilterMutsDataset(FilterDataset, MultistepDataset, UnbiasDataset):
    """ Chain mutation data with filtered reads. """

    FILTER_NAME = "filter"

    @classmethod
    def get_dataset1_load_func(cls):
        return load_idmut_dataset

    @classmethod
    def get_dataset2_type(cls):
        return FilterReadDataset

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
        region.add_mask(self.FILTER_NAME,
                        getattr(self.dataset2, "pos_kept"),
                        complement=True)
        return region

    def _integrate(self, batch1: IDmutMutsBatch, batch2: FilterBatchIO):
        """ Combine an idmut batch with a filter batch into a FilterMutsBatch.

        Parameters
        ----------
        batch1: IDmutMutsBatch
            Batch of mutation data from the IDmut step.
        batch2: FilterBatchIO
            Batch of read numbers that passed the filter filters.

        Returns
        -------
        FilterMutsBatch
            A batch containing only the reads and positions that pass
            the filter, clipped to the dataset's region.
        """
        if self.masked_read_nums is not None:
            read_nums = np.setdiff1d(batch2.read_nums,
                                     self.masked_read_nums.get(batch2.batch),
                                     assume_unique=True)
        else:
            read_nums = batch2.read_nums
        filtered_batch = apply_filters(batch1,
                                  read_nums,
                                  self.region,
                                  sanitize=False)
        if self.min_mut_gap > 0 and self.mut_collisions == MUT_COLLISIONS_MERGE:
            return FilterMutsBatch(
                batch=filtered_batch.batch,
                read_nums=filtered_batch.read_nums,
                region=filtered_batch.region,
                seg_end5s=filtered_batch.seg_end5s,
                seg_end3s=filtered_batch.seg_end3s,
                muts=filtered_batch.merge_close_muts(self.pattern,
                                                   self.min_mut_gap),
            )
        return filtered_batch


class JoinFilterMutsDataset(FilterDataset, JoinMutsDataset, MergedUnbiasDataset):

    @classmethod
    def get_report_type(cls):
        return JoinFilterReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_filter_dataset

    @classmethod
    def get_batch_type(cls):
        return FilterMutsBatch

    @classmethod
    def name_batch_attrs(cls):
        return [BATCH_NUM, READ_NUMS, SEG_END5S, SEG_END3S, MUTS]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        clusts = self.report.get_field(JoinedClustersF, missing_ok=True)
        if clusts is not None:
            raise TypeError(f"{self} has no clusters, but got {clusts}")


load_filter_dataset = LoadFunction(FilterMutsDataset, JoinFilterMutsDataset)
