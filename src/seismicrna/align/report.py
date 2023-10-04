from __future__ import annotations
from logging import getLogger

from ..core import path, report
from ..core.cmd import CMD_ALIGN
from ..core.report import Report

logger = getLogger(__name__)


class AlignReport(Report):

    @classmethod
    def field_names(cls):
        return (report.SampleF.key,
                report.IsDemultF.key,
                report.IsPairedEndF.key,
                report.PhredEncF.key,
                report.UseFastqcF.key,
                report.UseCutadaptF.key,
                report.CutadaptQ1.key,
                report.CutadaptQ2.key,
                report.CutadaptG1.key,
                report.CutadaptA1.key,
                report.CutadaptG2.key,
                report.CutadaptA2.key,
                report.CutadaptOverlap.key,
                report.CutadaptErrors.key,
                report.CutadaptIndels.key,
                report.CutadaptNextSeq.key,
                report.CutadaptNoTrimmed.key,
                report.CutadaptNoUntrimmed.key,
                report.CutadaptMinLength.key,
                report.Bowtie2Local.key,
                report.Bowtie2Discord.key,
                report.Bowtie2Mixed.key,
                report.Bowtie2Dovetail.key,
                report.Bowtie2Contain.key,
                report.Bowtie2Unal.key,
                report.Bowtie2ScoreMin.key,
                report.Bowtie2MinLength.key,
                report.Bowtie2MaxLength.key,
                report.Bowtie2GBar.key,
                report.Bowtie2SeedLength.key,
                report.Bowtie2SeedInterval.key,
                report.Bowtie2ExtTries.key,
                report.Bowtie2Reseed.key,
                report.Bowtie2Dpad.key,
                report.Bowtie2Orient.key,
                report.MinMapQual.key,
                report.ReadsInit.key,
                report.ReadsTrim.key,
                report.ReadsAlign.key,
                report.ReadsDedup.key,
                report.ReadsRefs.key)

    @classmethod
    def dir_seg_types(cls):
        return path.SampSeg, path.CmdSeg

    @classmethod
    def file_seg_type(cls):
        return path.AlignRepSeg

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: CMD_ALIGN}

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
