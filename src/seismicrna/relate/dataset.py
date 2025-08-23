from abc import ABC
from functools import cached_property

from .io import ReadNamesBatchIO, RelateBatchIO, RefseqIO
from .report import RelateReport, PoolReport
from ..core import path
from ..core.dataset import (Dataset,
                            LoadedDataset,
                            MergedRegionDataset,
                            MutsDataset,
                            LoadFunction,
                            TallDataset)
from ..core.header import NO_K, NO_KS
from ..core.report import RefseqChecksumF
from ..core.seq import FULL_NAME, Region


class AverageDataset(Dataset, ABC):
    """ Dataset of population average data. """

    @property
    def ks(self):
        return NO_KS

    @property
    def best_k(self):
        return NO_K


class RelateDataset(AverageDataset, ABC):
    """ Dataset of relationships. """


class RelateMutsDataset(RelateDataset, LoadedDataset, MutsDataset):
    """ Dataset of mutations from the Relate step. """

    @classmethod
    def get_batch_type(cls):
        return RelateBatchIO

    @classmethod
    def get_report_type(cls):
        return RelateReport

    @cached_property
    def refseq(self):
        return RefseqIO.load(
            RefseqIO.build_path({path.TOP: self.top,
                                 path.SAMPLE: self.sample,
                                 path.BRANCHES: self.branches,
                                 path.REF: self.ref}),
            checksum=self.report.get_field(RefseqChecksumF)
        ).refseq

    @cached_property
    def region(self):
        return Region(ref=self.ref,
                      seq=self.refseq,
                      end5=1,
                      end3=len(self.refseq),
                      name=FULL_NAME)

    @property
    def pattern(self):
        return None

    @cached_property
    def paired(self):
        """ Whether the reads are paired-end. """
        if self.num_batches == 0:
            return False
        return self.get_batch(0).num_segments == 2

    def get_batch(self, batch: int):
        # Load the saved batch, which is a RelateBatchIO instance.
        relate_batch = super().get_batch(batch)
        # Generate a RelateRegionMutsBatch from that batch.
        return relate_batch.to_region_batch(self.region)


class PoolDataset(TallDataset, ABC):
    """ Pooled dataset of relationships. """


class PoolMutsDataset(RelateDataset,
                      PoolDataset,
                      MutsDataset,
                      MergedRegionDataset):
    """ Load pooled batches of relationships. """

    @classmethod
    def get_report_type(cls):
        return PoolReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_relate_dataset

    @cached_property
    def region(self):
        return self._get_common_attr("region")


class NamesDataset(AverageDataset, ABC):

    @classmethod
    def kind(cls):
        return "names"


class ReadNamesDataset(NamesDataset, LoadedDataset):
    """ Dataset of read names from the Relate step. """

    @classmethod
    def get_batch_type(cls):
        return ReadNamesBatchIO

    @classmethod
    def get_report_type(cls):
        return RelateReport

    @property
    def pattern(self):
        return None


class PoolReadNamesDataset(NamesDataset, PoolDataset):
    """ Pooled Dataset of read names. """

    @classmethod
    def get_report_type(cls):
        return PoolReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_read_names_dataset


load_relate_dataset = LoadFunction(RelateMutsDataset, PoolMutsDataset)
load_read_names_dataset = LoadFunction(ReadNamesDataset, PoolReadNamesDataset)
