from abc import ABC, abstractmethod
from functools import cached_property

from .base import BaseGraph, make_path_subject, make_title_action_sample
from .cgroup import ClusterGroupGraph
from ..core.header import K_CLUST_KEY


class OneSourceGraph(BaseGraph, ABC):
    """ Graph of data from one source of data (Dataset or Table). """

    @cached_property
    @abstractmethod
    def action(self) -> str:
        """ Action that generated the data. """

    @property
    def col_tracks(self):
        return None

    @cached_property
    def title_action_sample(self):
        return make_title_action_sample(self.action, self.sample)


class OneSourceClusterGroupGraph(OneSourceGraph, ClusterGroupGraph, ABC):

    def __init__(self, *,
                 k: int | None,
                 clust: int | None,
                 **kwargs):
        self.k_clust_list = kwargs.pop(K_CLUST_KEY, None)
        super().__init__(**kwargs)
        self.k = k
        self.clust = clust

    @cached_property
    def path_subject(self):
        return make_path_subject(self.action, self.k, self.clust)
