from datetime import datetime

from .compare import RunOrderResults, find_best_order
from .io import ClustIO, ClustBatchIO
from .uniq import UniqReads
from ..core import path
from ..core.report import (BatchedReport,
                           SampleF,
                           RefF,
                           SectF,
                           End5F,
                           End3F,
                           NumUniqReadKeptF,
                           MaxClustsF,
                           ClustNumRunsF,
                           MinIterClustF,
                           MaxIterClustF,
                           ClustConvThreshF,
                           ClustsConvF,
                           ClustsLogLikesF,
                           ClustsLikeMeanF,
                           ClustsLikeStdF,
                           ClustsVarInfoF,
                           ClustsBicF,
                           NumClustsF)


class ClustReport(BatchedReport, ClustIO):

    @classmethod
    def file_seg_type(cls):
        return path.ClustRepSeg

    @classmethod
    def _batch_types(cls):
        return ClustBatchIO,

    @classmethod
    def fields(cls):
        return [
            # Sample, reference, and section information.
            SampleF,
            RefF,
            SectF,
            End5F,
            End3F,
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
            ClustsLikeMeanF,
            ClustsLikeStdF,
            ClustsVarInfoF,
            ClustsBicF,
            NumClustsF,
        ] + super().fields()

    @classmethod
    def path_segs(cls):
        return (path.SampSeg,
                path.CmdSeg,
                path.RefSeg,
                path.SectSeg,
                path.ClustRepSeg)

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: path.CMD_CLS_DIR}

    @classmethod
    def get_batch_seg(cls):
        return path.ClustBatSeg

    @classmethod
    def from_clusters(cls,
                      ord_runs: dict[int, RunOrderResults],
                      uniq_reads: UniqReads,
                      max_order: int,
                      num_runs: int, *,
                      min_iter: int,
                      max_iter: int,
                      conv_thresh: float,
                      checksums: list[str],
                      began: datetime,
                      ended: datetime):
        """ Create a ClusterReport from EmClustering objects. """
        # Initialize a new ClusterReport.
        return cls(sample=uniq_reads.sample,
                   ref=uniq_reads.ref,
                   sect=uniq_reads.section.name,
                   end5=uniq_reads.section.end5,
                   end3=uniq_reads.section.end3,
                   n_uniq_reads=uniq_reads.num_uniq,
                   max_order=max_order,
                   num_runs=num_runs,
                   min_iter=min_iter,
                   max_iter=max_iter,
                   conv_thresh=conv_thresh,
                   checksums={ClustBatchIO.btype(): checksums},
                   n_batches=len(checksums),
                   converged={order: runs.converged
                              for order, runs in ord_runs.items()},
                   log_likes={order: runs.log_likes
                              for order, runs in ord_runs.items()},
                   log_like_mean={order: runs.log_like_mean
                                  for order, runs in ord_runs.items()},
                   log_like_std={order: runs.log_like_std
                                 for order, runs in ord_runs.items()},
                   var_info={order: runs.var_info
                             for order, runs in ord_runs.items()},
                   bic={order: runs.best.bic
                        for order, runs in ord_runs.items()},
                   best_order=find_best_order(ord_runs),
                   began=began,
                   ended=ended)

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
