from abc import ABC
from functools import cached_property

from .base import GraphBase
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
