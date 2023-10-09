from .batch import MaskBatch
from .files import MaskFile
from ..core import path
from ..core.cmd import CMD_MASK
from ..core.report import (BatchReport,
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
                           MinInfoPosF,
                           MaxMutPosF,
                           NumPosInitF,
                           NumPosCutGUF,
                           NumPosCutPolyAF,
                           NumPosCutUserF,
                           NumPosCutLoInfoF,
                           NumPosCutHiMutF,
                           NumPosKeptF,
                           PosCutGUF,
                           PosCutPolyAF,
                           PosCutUserF,
                           PosCutLoInfoF,
                           PosCutHiMutF,
                           PosKeptF,
                           MinInfoReadF,
                           MaxMutReadF,
                           MinMutGapF,
                           NumReadsInitF,
                           NumReadsLoInfoF,
                           NumReadsHiMutF,
                           NumReadsCloseMutF,
                           NumReadsKeptF)


class MaskReport(BatchReport, MaskFile):

    @classmethod
    def file_seg_type(cls):
        return path.MaskRepSeg

    @classmethod
    def _batch_types(cls):
        return MaskBatch,

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
            MinInfoPosF,
            MaxMutPosF,
            # Position filtering results.
            NumPosInitF,
            NumPosCutGUF,
            NumPosCutPolyAF,
            NumPosCutUserF,
            NumPosCutLoInfoF,
            NumPosCutHiMutF,
            NumPosKeptF,
            PosCutGUF,
            PosCutPolyAF,
            PosCutUserF,
            PosCutLoInfoF,
            PosCutHiMutF,
            PosKeptF,
            # Read filtering parameters.
            MinInfoReadF,
            MaxMutReadF,
            MinMutGapF,
            # Read filtering results.
            NumReadsInitF,
            NumReadsLoInfoF,
            NumReadsHiMutF,
            NumReadsCloseMutF,
            NumReadsKeptF,
        ] + super().fields()

    @classmethod
    def path_segs(cls):
        return (path.SampSeg,
                path.CmdSeg,
                path.RefSeg,
                path.SectSeg,
                path.MaskRepSeg)

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: CMD_MASK}

    @classmethod
    def get_batch_seg(cls):
        return path.MaskBatSeg

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
