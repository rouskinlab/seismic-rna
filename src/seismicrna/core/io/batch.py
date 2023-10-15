import re
from abc import ABC

from .file import BrickleIO
from ..batch import ReadBatch


class ReadBatchIO(ReadBatch, BrickleIO, ABC):
    """ Pickled file of a batch of data. """

    @classmethod
    def btype(cls):
        btype, = re.match("^([a-z]*)batchio", cls.__name__.lower()).groups()
        return btype
