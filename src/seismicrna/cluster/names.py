"""
Cluster -- Names Module
========================================================================
Auth: Matty

Define names for the indexes of the cluster tables.
"""


# Cluster indexes
ORD_NAME = "Order"
CLS_NAME = "Cluster"
ORD_CLS_NAME = ORD_NAME, CLS_NAME

# Cluster memberships
READ_NAME = "Read Name"

# Bit vector counts
OBS_NAME = "Log Observed"
EXP_NAME = "Log Expected"

# Ensemble name
ENSEMBLE_NAME = "Ensemble"


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
    return f"Cluster {order}-{cluster}" if order > 0 else ENSEMBLE_NAME
