from abc import ABC

from .io import RelateFile
from .table import RelatePositionTableLoader
from ..core.lists import List, PositionList


class RelateList(List, RelateFile, ABC):
    pass


class RelatePositionList(PositionList, RelateList):

    @classmethod
    def get_table_type(cls):
        return RelatePositionTableLoader
