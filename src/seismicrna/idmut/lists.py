from abc import ABC

from .io import IDmutFile
from .table import IDmutPositionTableLoader
from ..core.lists import List, PositionList


class IDmutList(List, IDmutFile, ABC):
    pass


class IDmutPositionList(PositionList, IDmutList):

    @classmethod
    def get_table_type(cls):
        return IDmutPositionTableLoader
