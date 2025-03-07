from abc import ABC
from functools import cached_property

from .table import TableGraph, TableRunner
from ..core.arg import opt_window, opt_winmin
from ..core.seq import POS_NAME


class RollingGraph(TableGraph, ABC):

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
        return super().predicate + ["-".join(map(str, [self._size,
                                                       self._min_count]))]


class RollingRunner(TableRunner, ABC):

    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_window, opt_winmin]
