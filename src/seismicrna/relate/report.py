from __future__ import annotations
from logging import getLogger

from ..core import path
from ..core.cmd import CMD_REL
from ..core.report import (BatchReport, calc_seqlen, calc_time_taken,
                           calc_n_batches, calc_speed)

logger = getLogger(__name__)


BATCH_INDEX_COL = "Read Name"


class RelateReport(BatchReport):
    __slots__ = ("sample", "ref", "seq", "length",
                 "n_reads_rel_pass", "n_reads_rel_fail",
                 "checksums", "n_batches",
                 "began", "ended", "taken", "speed")

    def __init__(self, *,
                 length=calc_seqlen,
                 n_batches=calc_n_batches,
                 taken=calc_time_taken,
                 speed=calc_speed,
                 **kwargs):
        # Note that the named keyword arguments must come after **kwargs
        # because they are calculated using the values of the arguments
        # in **kwargs. If **kwargs was given last, those values would be
        # undefined when the named keyword arguments would be computed.
        super().__init__(**kwargs, length=length, n_batches=n_batches,
                         taken=taken, speed=speed)

    @classmethod
    def path_segs(cls):
        return path.SampSeg, path.CmdSeg, path.RefSeg, path.RelateRepSeg

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: CMD_REL}

    @classmethod
    def get_batch_seg(cls):
        return path.RelateBatSeg

    @classmethod
    def default_batch_ext(cls):
        return path.PARQ_EXTS[0]

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
