from abc import ABC
from pathlib import Path

import pandas as pd

from .batch import ClusterReadBatch
from .compare import EMRunsK
from ..core import path
from ..core.header import ClustHeader
from ..core.io import ReadBatchIO, SectIO
from ..mask.data import MaskMutsDataset


class ClusterIO(SectIO, ABC):

    @classmethod
    def auto_fields(cls):
        return super().auto_fields() | {path.CMD: path.CMD_CLUST_DIR}


class ClusterBatchIO(ReadBatchIO, ClusterIO, ClusterReadBatch):

    @classmethod
    def file_seg_type(cls):
        return path.ClustBatSeg


def write_batches(dataset: MaskMutsDataset,
                  ks: list[EMRunsK],
                  brotli_level: int,
                  top: Path):
    """ Write the cluster memberships to batch files. """
    checksums = list()
    read_nums = dict()

    def get_read_nums(num: int):
        if (nums := read_nums.get(num)) is not None:
            return nums
        for batch_ in dataset.iter_batches():
            read_nums[batch_.batch] = batch_.read_nums
        return read_nums[num]

    # Filter the numbers of clusters to keep only those with at least
    # one successful run.
    ks = [runs for runs in ks if runs.best is not None]
    # Write each cluster batch.
    for batch_num in dataset.batch_nums:
        resps = [runs.best.get_resps(batch_num) for runs in ks]
        if resps:
            resps = pd.concat(resps, axis=1)
        else:
            resps = pd.DataFrame(index=get_read_nums(batch_num),
                                 columns=ClustHeader(ks=[]).index)
        batch = ClusterBatchIO(sample=dataset.sample,
                               ref=dataset.ref,
                               sect=dataset.sect,
                               batch=batch_num,
                               resps=resps)
        _, checksum = batch.save(top, brotli_level=brotli_level)
        checksums.append(checksum)
    ks_written = [runs.k for runs in ks]
    return checksums, ks_written

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
