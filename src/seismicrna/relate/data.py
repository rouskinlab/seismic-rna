from .batch import RelateRefseqBatch
from .io import RelateBatchIO
from .report import RelateReport
from ..core.data import LoadedMutsDataset, PooledDataset


class RelateLoader(LoadedMutsDataset):
    """ Load batches of relation vectors. """

    @classmethod
    def get_batch_type(cls):
        return RelateBatchIO

    @classmethod
    def get_report_type(cls):
        return RelateReport

    @property
    def pattern(self):
        return None

    def get_batch(self, batch: int):
        relate_batch = super().get_batch(batch)
        # Add the reference sequence to the batch.
        if relate_batch.max_pos != len(self.refseq):
            raise ValueError(f"Reference sequence is {len(self.refseq)} nt, "
                             f"but {relate_batch} has {relate_batch.max_pos}")
        return RelateRefseqBatch(refseq=self.refseq,
                                 batch=relate_batch.batch,
                                 muts=relate_batch.muts,
                                 end5s=relate_batch.end5s,
                                 mid5s=relate_batch.mid5s,
                                 mid3s=relate_batch.mid3s,
                                 end3s=relate_batch.end3s,
                                 sanitize=False)


class PooledRelateLoader(PooledDataset):

    @classmethod
    def get_dataset_type(cls):
        return RelateLoader

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
