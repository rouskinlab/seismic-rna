"""
Cluster -- Names Module
========================================================================
Auth: Matty

Define names for the indexes of the cluster tables.
"""

# Bit vector counts
OBS_NAME = "Log Observed"
EXP_NAME = "Log Expected"

# Ensemble average name
AVERAGE_NAME = "Average"


def validate_order_cluster(order: int, cluster: int, allow_zero: bool = False):
    if not isinstance(order, int):
        raise TypeError(f"order must be an int, not {type(order).__name__}")
    if not isinstance(cluster, int):
        raise TypeError(f"cluster must be an int, not {type(cluster).__name__}")
    ineq = '≥' if allow_zero else '>'
    if order < 0 or (order == 0 and not allow_zero):
        raise ValueError(f"order must be {ineq} 0, but got {order}")
    if cluster < 0 or (cluster == 0 and not allow_zero):
        raise ValueError(f"cluster must be {ineq} 0, but got {cluster}")
    if cluster > order:
        raise ValueError(f"cluster ({cluster}) must be ≤ order ({order})")


def fmt_clust_name(order: int, cluster: int, allow_zero: bool = False):
    validate_order_cluster(order, cluster, allow_zero)
    return f"Cluster {order}-{cluster}" if order > 0 else AVERAGE_NAME

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
