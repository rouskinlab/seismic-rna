from __future__ import annotations

from ..core import path
from ..core.io import (Report,
                       SampleF,
                       IsDemultF,
                       IsPairedEndF,
                       PhredEncF,
                       UseFastqcF,
                       UseCutadaptF,
                       CutadaptQ1,
                       CutadaptQ2,
                       CutadaptG1,
                       CutadaptA1,
                       CutadaptG2,
                       CutadaptA2,
                       CutadaptOverlap,
                       CutadaptErrors,
                       CutadaptIndels,
                       CutadaptNextSeq,
                       CutadaptNoTrimmed,
                       CutadaptNoUntrimmed,
                       CutadaptMinLength,
                       Bowtie2Local,
                       Bowtie2Discord,
                       Bowtie2Mixed,
                       Bowtie2Dovetail,
                       Bowtie2Contain,
                       Bowtie2Unal,
                       Bowtie2ScoreMin,
                       Bowtie2MinLength,
                       Bowtie2MaxLength,
                       Bowtie2GBar,
                       Bowtie2SeedLength,
                       Bowtie2SeedInterval,
                       Bowtie2ExtTries,
                       Bowtie2Reseed,
                       Bowtie2Dpad,
                       Bowtie2Orient,
                       MinMapQual,
                       ReadsInit,
                       ReadsTrim,
                       ReadsAlign,
                       ReadsDedup,
                       ReadsRefs)


class AlignReport(Report):

    @classmethod
    def fields(cls):
        return [SampleF,
                IsDemultF,
                IsPairedEndF,
                PhredEncF,
                UseFastqcF,
                UseCutadaptF,
                CutadaptQ1,
                CutadaptQ2,
                CutadaptG1,
                CutadaptA1,
                CutadaptG2,
                CutadaptA2,
                CutadaptOverlap,
                CutadaptErrors,
                CutadaptIndels,
                CutadaptNextSeq,
                CutadaptNoTrimmed,
                CutadaptNoUntrimmed,
                CutadaptMinLength,
                Bowtie2Local,
                Bowtie2Discord,
                Bowtie2Mixed,
                Bowtie2Dovetail,
                Bowtie2Contain,
                Bowtie2Unal,
                Bowtie2ScoreMin,
                Bowtie2MinLength,
                Bowtie2MaxLength,
                Bowtie2GBar,
                Bowtie2SeedLength,
                Bowtie2SeedInterval,
                Bowtie2ExtTries,
                Bowtie2Reseed,
                Bowtie2Dpad,
                Bowtie2Orient,
                MinMapQual,
                ReadsInit,
                ReadsTrim,
                ReadsAlign,
                ReadsDedup,
                ReadsRefs] + super().fields()

    @classmethod
    def dir_seg_types(cls):
        return path.SampSeg, path.CmdSeg

    @classmethod
    def file_seg_type(cls):
        return path.AlignRepSeg

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.CMD_ALN_DIR}

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
