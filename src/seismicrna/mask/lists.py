from abc import ABC

from .io import MaskFile
from .table import MaskPositionTableLoader
from ..core.lists import List, PositionList


class MaskList(List, MaskFile, ABC):
    pass


class MaskPositionList(PositionList, MaskList):

    @classmethod
    def get_table_type(cls):
        return MaskPositionTableLoader
