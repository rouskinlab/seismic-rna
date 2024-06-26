from __future__ import annotations

from abc import ABC, abstractmethod

from ..core import path
from ..core.report import (Report,
                           SampleF,
                           RefF,
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
                           Bowtie2Un,
                           Bowtie2ScoreMin,
                           Bowtie2MinLengthF,
                           Bowtie2MaxLengthF,
                           Bowtie2GBarF,
                           Bowtie2SeedLength,
                           Bowtie2SeedInterval,
                           Bowtie2ExtTries,
                           Bowtie2Reseed,
                           Bowtie2Dpad,
                           Bowtie2Orient,
                           SepStrandsF,
                           F1R2PlusF,
                           MinusLabelF,
                           MinMapQualF,
                           MinReadsF,
                           AlignReadsInitF,
                           ReadsTrimF,
                           ReadsAlignF,
                           ReadsDedupF,
                           ReadsRefsF)


class AlignReport(Report, ABC):

    @classmethod
    @abstractmethod
    def fields(cls):
        return [IsDemultF,
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
                Bowtie2ScoreMin,
                Bowtie2MinLengthF,
                Bowtie2MaxLengthF,
                Bowtie2GBarF,
                Bowtie2SeedLength,
                Bowtie2SeedInterval,
                Bowtie2ExtTries,
                Bowtie2Reseed,
                Bowtie2Dpad,
                Bowtie2Orient,
                Bowtie2Un,
                MinMapQualF,
                SepStrandsF,
                F1R2PlusF,
                MinusLabelF,
                MinReadsF,
                AlignReadsInitF,
                ReadsTrimF,
                ReadsAlignF,
                ReadsDedupF,
                ReadsRefsF] + super().fields()

    @classmethod
    def dir_seg_types(cls):
        return path.SampSeg, path.CmdSeg

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.CMD_ALIGN_DIR}


class AlignSampleReport(AlignReport):

    @classmethod
    def fields(cls):
        return [SampleF] + super().fields()

    @classmethod
    def file_seg_type(cls):
        return path.AlignSampleRepSeg

    def __init__(self, ref: str | None = None, **kwargs):
        if ref is not None:
            raise TypeError(f"Got an unexpected reference name: {repr(ref)}")
        super().__init__(demultiplexed=False, **kwargs)


class AlignRefReport(AlignReport):

    @classmethod
    def fields(cls):
        return [SampleF, RefF] + super().fields()

    @classmethod
    def file_seg_type(cls):
        return path.AlignRefRepSeg

    def __init__(self, ref: str, **kwargs):
        if ref is None:
            raise TypeError(f"Expected a reference name, but got {repr(ref)}")
        super().__init__(ref=ref, demultiplexed=True, **kwargs)

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
