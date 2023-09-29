from __future__ import annotations
from logging import getLogger

from .batch import QnamesBatch, RelateBatch
from .output import RelateOutput
from ..core import path
from ..core.cmd import CMD_REL
from ..core.report import BatchReport, calc_seqlen, calc_time_taken, calc_speed

logger = getLogger(__name__)

BATCH_INDEX_COL = "Read Name"


class RelateReport(BatchReport, RelateOutput):

    @classmethod
    def field_names(cls) -> tuple[str, ...]:
        return ("sample",
                "ref",
                "seq",
                "length",
                "n_reads_rel") + super().field_names() + ("began",
                                                          "ended",
                                                          "taken",
                                                          "speed")

    @classmethod
    def file_seg_type(cls):
        return path.RelateRepSeg

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: CMD_REL}

    @classmethod
    def _batch_types(cls):
        return QnamesBatch, RelateBatch

    def __init__(self, *,
                 length=calc_seqlen,
                 taken=calc_time_taken,
                 speed=calc_speed,
                 **kwargs):
        # Note that the named keyword arguments must come after **kwargs
        # because they are calculated using the values of the arguments
        # in **kwargs. If **kwargs was given last, those values would be
        # undefined when the named keyword arguments would be computed.
        super().__init__(**kwargs, length=length, taken=taken, speed=speed)

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
