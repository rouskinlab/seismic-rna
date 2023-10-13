import re
from abc import ABC

from .batch import Batch
from .iofile import SavedBrickle


class SavedBatch(SavedBrickle, Batch, ABC):
    """ Pickled file of a batch of data. """

    @classmethod
    def btype(cls):
        btype, = re.match("^saved([a-z]*)batch", cls.__name__.lower()).groups()
        return btype
