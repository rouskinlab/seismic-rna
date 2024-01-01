from __future__ import annotations

from ..core import path
from ..core.report import (Report,
                           SampleF,
                           RefF,
                           SectF,
                           ProfileF,
                           Quantile,
                           FoldTempF,
                           FoldMaxDistF,
                           FoldMinFreeEnergyF,
                           FoldMaxStructsF,
                           FoldPercent)


class FoldReport(Report):

    @classmethod
    def fields(cls):
        return [SampleF,
                RefF,
                SectF,
                ProfileF,
                Quantile,
                FoldTempF,
                FoldMaxDistF,
                FoldMinFreeEnergyF,
                FoldMaxStructsF,
                FoldPercent] + super().fields()

    @classmethod
    def file_seg_type(cls):
        return path.FoldRepSeg

    @classmethod
    def dir_seg_types(cls):
        return path.SampSeg, path.CmdSeg, path.RefSeg, path.SectSeg

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.CMD_FOLD_DIR}

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
