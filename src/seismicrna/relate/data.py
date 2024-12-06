from functools import cached_property

from .batch import RelateBatch
from .io import ReadNamesBatchIO, RelateBatchIO
from .report import RelateReport, PoolReport
from ..core.data import (LoadedDataset,
                         LoadedMutsDataset,
                         LoadFunction,
                         TallDataset,
                         TallMutsDataset)


class RelateDataset(LoadedMutsDataset):
    """ Dataset of mutations from the Relate step. """

    @classmethod
    def get_batch_type(cls):
        return RelateBatchIO

    @classmethod
    def get_report_type(cls):
        return RelateReport

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
        relate_batch = super().get_batch(batch)
        return RelateBatch(batch=relate_batch.batch,
                           seg_end5s=relate_batch.seg_end5s,
                           seg_end3s=relate_batch.seg_end3s,
                           muts=relate_batch.muts,
                           region=self.region,
                           sanitize=False)


class PoolDataset(TallMutsDataset):
    """ Load pooled batches of relationships. """

    @classmethod
    def get_report_type(cls):
        return PoolReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_relate_dataset


class ReadNamesDataset(LoadedDataset):
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


class PoolReadNamesDataset(TallDataset):
    """ Pooled Dataset of read names. """

    @classmethod
    def get_report_type(cls):
        return PoolReport

    @classmethod
    def get_dataset_load_func(cls):
        return load_read_names_dataset


load_relate_dataset = LoadFunction(RelateDataset, PoolDataset)
load_read_names_dataset = LoadFunction(ReadNamesDataset, PoolReadNamesDataset)

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
