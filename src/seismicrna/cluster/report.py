from datetime import datetime

from .compare import RunOrderResults, find_best_order
from .io import ClusterIO, ClusterBatchIO
from .uniq import UniqReads
from ..core import path
from ..core.report import (BatchedReport,
                           SampleF,
                           RefF,
                           SectF,
                           NumUniqReadKeptF,
                           MaxClustsF,
                           ClustNumRunsF,
                           MinIterClustF,
                           MaxIterClustF,
                           ClustConvThreshF,
                           ClustsConvF,
                           ClustsLogLikesF,
                           ClustsRMSDsF,
                           ClustsMeanRsF,
                           ClustsBicF,
                           NumClustsF)


class ClusterReport(BatchedReport, ClusterIO):

    @classmethod
    def file_seg_type(cls):
        return path.ClustRepSeg

    @classmethod
    def _batch_types(cls):
        return ClusterBatchIO,

    @classmethod
    def fields(cls):
        return [
            # Sample, reference, and section information.
            SampleF,
            RefF,
            SectF,
            NumUniqReadKeptF,
            # Clustering parameters.
            MaxClustsF,
            ClustNumRunsF,
            MinIterClustF,
            MaxIterClustF,
            ClustConvThreshF,
            # Clustering results.
            ClustsConvF,
            ClustsLogLikesF,
            ClustsRMSDsF,
            ClustsMeanRsF,
            ClustsBicF,
            NumClustsF,
        ] + super().fields()

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.CMD_CLUST_DIR}

    @classmethod
    def from_clusters(cls,
                      orders: list[RunOrderResults],
                      uniq_reads: UniqReads,
                      max_order: int,
                      num_runs: int, *,
                      min_iter: int,
                      max_iter: int,
                      em_thresh: float,
                      checksums: list[str],
                      began: datetime,
                      ended: datetime):
        """ Create a ClusterReport from EmClustering objects. """
        return cls(sample=uniq_reads.sample,
                   ref=uniq_reads.ref,
                   sect=uniq_reads.section.name,
                   n_uniq_reads=uniq_reads.num_uniq,
                   max_clusters=max_order,
                   em_runs=num_runs,
                   min_em_iter=min_iter,
                   max_em_iter=max_iter,
                   em_thresh=em_thresh,
                   checksums={ClusterBatchIO.btype(): checksums},
                   n_batches=len(checksums),
                   converged={runs.order: runs.converged for runs in orders},
                   log_likes={runs.order: runs.log_likes for runs in orders},
                   clusts_rmsds={runs.order: runs.rmsds for runs in orders},
                   clusts_meanr={runs.order: runs.meanr for runs in orders},
                   bic={runs.order: runs.bic for runs in orders},
                   best_order=find_best_order(orders),
                   began=began,
                   ended=ended)

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
