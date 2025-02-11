import re
from abc import ABC

from .brickle import BrickleIO
from ..batch import MutsBatch, ReadBatch
from ..seq import Region


class ReadBatchIO(ReadBatch, BrickleIO, ABC):
    """ Pickled file of a batch of data. """

    @classmethod
    def btype(cls):
        btype, = re.match("^([a-z_]*)batchio", cls.__name__.lower()).groups()
        return btype


class MutsBatchIO(MutsBatch, ReadBatchIO, ABC):
    """ Pickled file of a batch of mutational data. """

    def __init__(self, *args, region: Region, **kwargs):
        super().__init__(*args, **kwargs, region=region, ref=region.ref)
