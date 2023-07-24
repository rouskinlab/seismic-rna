from .base import OneSeqGraph
from ..core.cli import CLUST_INDIV, CLUST_ORDER, CLUST_UNITE
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
            clusters_params = [dict(order=order, cluster=cluster)
                               for order, cluster in table.ord_clust]
        elif arrange == CLUST_ORDER:
            # One file per order, with a subplot for each cluster.
            orders = sorted(table.orders)
            clusters_params = [dict(order=order) for order in orders]
        elif arrange == CLUST_UNITE:
            # One file, with subplots of all clusters of all orders.
            clusters_params = [dict()]
        else:
            raise ValueError(f"Invalid value for arrange: '{arrange}'")
    elif isinstance(table, PosTableLoader):
        single_type = es_type
        multi_type = em_type
        clusters_params = [dict()]
    else:
        single_type = None
        multi_type = None
        clusters_params = list()
    return single_type, multi_type, clusters_params
