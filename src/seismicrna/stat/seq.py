from .base import OneSeqGraph
from ..core.arg import CLUST_INDIV, CLUST_ORDER, CLUST_UNITE
from ..table.load import PosTableLoader, ClustPosTableLoader


def get_table_params(table: PosTableLoader | ClustPosTableLoader,
                     arrange: str,
                     es_type: type[OneSeqGraph] | None = None,
                     cs_type: type[OneSeqGraph] | None = None,
                     em_type: type[OneSeqGraph] | None = None,
                     cm_type: type[OneSeqGraph] | None = None):
    if isinstance(table, ClustPosTableLoader):
        single_type = cs_type
        multi_type = cm_type
        if arrange == CLUST_INDIV:
            # One file per cluster, with no subplots.
            clusters_params = [dict(order=order, clust=cluster)
                               for order, cluster in table.header.clusts]
        elif arrange == CLUST_ORDER:
            # One file per order, with one subplot per cluster.
            clusters_params = [dict(order=order) for order
                               in sorted(table.header.orders)]
        elif arrange == CLUST_UNITE:
            # One file, with one subplot per cluster for all orders.
            clusters_params = [dict()]
        else:
            raise ValueError(f"Invalid value for arrange: {repr(arrange)}")
    elif isinstance(table, PosTableLoader):
        single_type = es_type
        multi_type = em_type
        clusters_params = [dict()]
    else:
        single_type = None
        multi_type = None
        clusters_params = []
    return single_type, multi_type, clusters_params

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
