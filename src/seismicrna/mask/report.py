from .io import MaskIO, MaskBatchIO
from ..core import path
from ..core.report import (BatchedReport,
                           SampleF,
                           RefF,
                           SectF,
                           End5F,
                           End3F,
                           CountMutsF,
                           CountRefsF,
                           ExclGUF,
                           ExclPolyAF,
                           ExclUserPosF,
                           DiscontigF,
                           NumDiscontigF,
                           MinNInfoPosF,
                           MaxFMutPosF,
                           NumPosInitF,
                           NumPosCutGUF,
                           NumPosCutPolyAF,
                           NumPosCutListF,
                           NumPosCutLoInfoF,
                           NumPosCutHiMutF,
                           NumPosKeptF,
                           PosCutGUF,
                           PosCutPolyAF,
                           PosCutListF,
                           PosCutLoInfoF,
                           PosCutHiMutF,
                           PosKeptF,
                           MinNCovReadF,
                           MinFInfoReadF,
                           MaxFMutReadF,
                           MinMutGapF,
                           NumReadsInitF,
                           NumReadsLoNCovF,
                           NumReadsLoInfoF,
                           NumReadsHiMutF,
                           NumReadsCloseMutF,
                           NumReadsKeptF,
                           QuickUnbiasF,
                           QuickUnbiasThreshF)


class MaskReport(BatchedReport, MaskIO):

    @classmethod
    def file_seg_type(cls):
        return path.MaskRepSeg

    @classmethod
    def _batch_types(cls):
        return MaskBatchIO,

    @classmethod
    def fields(cls):
        return [
            # Sample, reference, and section information.
            SampleF,
            RefF,
            SectF,
            End5F,
            End3F,
            # Types of mutations and matches to count.
            CountMutsF,
            CountRefsF,
            # Position filtering parameters.
            ExclGUF,
            ExclPolyAF,
            ExclUserPosF,
            MinNInfoPosF,
            MaxFMutPosF,
            # Position filtering results.
            NumPosInitF,
            NumPosCutGUF,
            NumPosCutPolyAF,
            NumPosCutListF,
            NumPosCutLoInfoF,
            NumPosCutHiMutF,
            NumPosKeptF,
            PosCutGUF,
            PosCutPolyAF,
            PosCutListF,
            PosCutLoInfoF,
            PosCutHiMutF,
            PosKeptF,
            # Read filtering parameters.
            MinNCovReadF,
            DiscontigF,
            MinFInfoReadF,
            MaxFMutReadF,
            MinMutGapF,
            # Read filtering results.
            NumReadsInitF,
            NumReadsLoNCovF,
            NumDiscontigF,
            NumReadsLoInfoF,
            NumReadsHiMutF,
            NumReadsCloseMutF,
            NumReadsKeptF,
            # Observer bias correction.
            QuickUnbiasF,
            QuickUnbiasThreshF,
        ] + super().fields()

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.CMD_MASK_DIR}

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
