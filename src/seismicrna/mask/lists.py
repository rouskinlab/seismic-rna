from abc import ABC

from .io import MaskFile
from .table import MaskPositionTableLoader
from ..core import path
from ..core.lists import List, PositionList


class MaskList(List, MaskFile, ABC):

    def __init__(self, *, reg: str, **kwargs):
        super().__init__(**kwargs)
        self.reg = reg


class MaskPositionList(PositionList, MaskList):

    @classmethod
    def get_table_type(cls):
        return MaskPositionTableLoader

    @classmethod
    def list_init_table_attrs(cls):
        return [*super().list_init_table_attrs(), path.REG]
