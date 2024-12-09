from __future__ import annotations

from functools import cache
from pathlib import Path

from .io import RelateIO, ReadNamesBatchIO, RelateBatchIO
from ..core import path
from ..core.io import RefIO
from ..core.report import (Report,
                           BatchedRefseqReport,
                           RefF,
                           SampleF,
                           PooledSamplesF,
                           MinMapQualF,
                           MinReadsF,
                           NumReadsXamF,
                           NumReadsRelF,
                           Insert3F,
                           AmbindelF,
                           OverhangsF,
                           MinPhredF,
                           PhredEncF,
                           ClipEnd5F,
                           ClipEnd3F)

BATCH_INDEX_COL = "Read Name"


class RelateReport(BatchedRefseqReport, RelateIO):

    @classmethod
    def fields(cls):
        return [SampleF,
                RefF,
                MinMapQualF,
                PhredEncF,
                MinPhredF,
                Insert3F,
                AmbindelF,
                OverhangsF,
                ClipEnd5F,
                ClipEnd3F,
                MinReadsF,
                NumReadsXamF,
                NumReadsRelF] + super().fields()

    @classmethod
    def file_seg_type(cls):
        return path.RelateRepSeg

    @classmethod
    def _batch_types(cls):
        return ReadNamesBatchIO, RelateBatchIO

    def refseq_file(self, top: Path):
        return refseq_file_path(top,
                                self.get_field(SampleF),
                                self.get_field(RefF))


class PoolReport(Report, RefIO):

    @classmethod
    def file_seg_type(cls):
        return path.RelateRepSeg

    @classmethod
    def fields(cls):
        return [
            # Sample and reference.
            SampleF,
            RefF,
            # Pooled samples.
            PooledSamplesF,
        ] + super().fields()

    @classmethod
    def path_segs(cls):
        return (path.SampSeg,
                path.CmdSeg,
                path.RefSeg,
                path.RelateRepSeg)

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.CMD_REL_DIR}


@cache
def refseq_file_seg_types():
    return RelateReport.seg_types()[:-1] + (path.RefseqFileSeg,)


@cache
def refseq_file_auto_fields():
    return {**RelateReport.auto_fields(), path.EXT: path.BROTLI_PICKLE_EXT}


def refseq_file_path(top: Path, sample: str, ref: str):
    return path.build(*refseq_file_seg_types(),
                      **refseq_file_auto_fields(),
                      top=top,
                      sample=sample,
                      ref=ref)

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
