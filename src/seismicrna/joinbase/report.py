from ..core import path
from ..core.io import SectIO
from ..core.report import (Report,
                           SampleF,
                           RefF,
                           SectF,
                           JoinedSectionsF,
                           JoinedClustersF)


class JoinMaskReport(Report, SectIO):

    @classmethod
    def file_seg_type(cls):
        return path.MaskRepSeg

    @classmethod
    def fields(cls):
        return [
            # Sample and reference.
            SampleF,
            RefF,
            SectF,
            # Joined data.
            JoinedSectionsF,
        ] + super().fields()

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.CMD_MASK_DIR}


class JoinClusterReport(Report, SectIO):

    @classmethod
    def file_seg_type(cls):
        return path.ClustRepSeg

    @classmethod
    def fields(cls):
        return [
            # Sample and reference.
            SampleF,
            RefF,
            SectF,
            # Joined data.
            JoinedSectionsF,
            JoinedClustersF,
        ] + super().fields()

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.CMD_CLUST_DIR}

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
