from abc import ABC
from functools import cached_property

from .base import GraphBase, GraphRunner
from ..core.arg import opt_window, opt_winmin
from ..core.seq import POS_NAME


class RollingGraph(GraphBase, ABC):

    def __init__(self, *, window: int, winmin: int, **kwargs):
        super().__init__(**kwargs)
        self._size = window
        self._min_count = winmin

    @property
    def x_title(self):
        return POS_NAME

    @cached_property
    def details(self):
        return super().details + [f"window = {self._size} nt",
                                  f"min = {self._min_count} nt"]

    @cached_property
    def predicate(self):
        return "_".join(
            [super().predicate,
             "-".join(map(str, [self._size, self._min_count]))]
        )


class RollingRunner(GraphRunner, ABC):

    @classmethod
    def var_params(cls):
        return super().var_params() + [opt_window, opt_winmin]

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
