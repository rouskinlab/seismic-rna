from __future__ import annotations

from functools import cache
from logging import getLogger
from pathlib import Path

from .files import QnamesBatchFile, RelateBatchFile, RelateFile
from ..core import path
from ..core.cmd import CMD_REL
from ..core.report import (RefBatchReport,
                           calc_speed,
                           calc_taken,
                           RefF,
                           SampleF,
                           NumReadsRel,
                           TimeBeganF,
                           TimeEndedF,
                           TimeTakenF,
                           SpeedF)

logger = getLogger(__name__)

BATCH_INDEX_COL = "Read Name"


class RelateReport(RefBatchReport, RelateFile):

    @classmethod
    def fields(cls):
        return [
            SampleF,
            RefF,
            NumReadsRel,
            TimeBeganF,
            TimeEndedF,
            TimeTakenF,
            SpeedF,
        ] + super().fields()

    @classmethod
    def file_seg_type(cls):
        return path.RelateRepSeg

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: CMD_REL}

    @classmethod
    def _batch_types(cls):
        return QnamesBatchFile, RelateBatchFile

    def __init__(self, *, taken=calc_taken, speed=calc_speed, **kwargs):
        # Note that the named keyword arguments must come after **kwargs
        # because they are calculated using the values of the arguments
        # in **kwargs. If **kwargs was given last, those values would be
        # undefined when the named keyword arguments would be computed.
        super().__init__(**kwargs, taken=taken, speed=speed)

    def refseq_file(self, top: Path):
        return refseq_file_path(top,
                                self.get_field(SampleF),
                                self.get_field(RefF))


@cache
def refseq_file_seg_types():
    return RelateReport.seg_types()[:-1] + (path.RefseqFileSeg,)


@cache
def refseq_file_auto_fields():
    return {**RelateReport.auto_fields(), path.EXT: path.PICKLE_BROTLI_EXT}


def refseq_file_path(top: Path, sample: str, ref: str):
    return path.build(*refseq_file_seg_types(),
                      **refseq_file_auto_fields(),
                      top=top,
                      sample=sample,
                      ref=ref)

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
