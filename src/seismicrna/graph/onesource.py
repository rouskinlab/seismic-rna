from abc import ABC, abstractmethod
from functools import cached_property

from .base import GraphBase, make_path_subject, make_title_action_sample


class OneSourceGraph(GraphBase, ABC):
    """ Graph of data from one source of data (Dataset or Table). """

    def __init__(self, *,
                 k: int | None,
                 clust: int | None,
                 **kwargs):
        super().__init__(**kwargs)
        self.k = k
        self.clust = clust

    @cached_property
    @abstractmethod
    def action(self) -> str:
        """ Action that generated the data. """

    @property
    def col_tracks(self):
        return None

    @cached_property
    def path_subject(self):
        return make_path_subject(self.action, self.k, self.clust)

    @cached_property
    def title_action_sample(self):
        return make_title_action_sample(self.action, self.sample)

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
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
