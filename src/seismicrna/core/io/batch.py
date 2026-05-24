import re
from abc import ABC

from .brickle import BrickleIO
from ..batch import MutsBatch, ReadBatch
from ..seq import Region


class ReadBatchIO(ReadBatch, BrickleIO, ABC):
    """Pickled file of a batch of data."""

    @classmethod
    def btype(cls):
        (btype,) = re.match("^([a-z_]*)batchio", cls.__name__.lower()).groups()
        return btype

    @property
    def is_self_contained(self) -> bool:
        """Whether this batch file contains full mutation/coordinate data
        and does not require loading its predecessor batch to be useful."""
        return (
            getattr(self, "muts", None) is not None
            and getattr(self, "region", None) is not None
        )


class MutsBatchIO(MutsBatch, ReadBatchIO, ABC):
    """Pickled file of a batch of mutational data."""

    def __init__(self, *args, region: Region, **kwargs):
        super().__init__(*args, **kwargs, region=region, ref=region.ref)
