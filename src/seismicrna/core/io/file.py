from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Any

from .brickle import load_brickle, save_brickle
from .. import path

DEFAULT_BROTLI_LEVEL = 10
PICKLE_PROTOCOL = 5


class FileIO(path.HasFilePath, ABC):
    """ Object that can be saved to and loaded from a file path. """

    @classmethod
    @abstractmethod
    def load(cls, file: Path):
        """ Load an object from a file. """

    @abstractmethod
    def save(self, top: Path, **kwargs):
        """ Save the object to a file. """

    def __str__(self):
        return type(self).__name__


class SampleFileIO(FileIO, path.HasSampleFilePath, ABC):
    pass


class RefFileIO(SampleFileIO, path.HasRefFilePath, ABC):
    pass


class RegFileIO(RefFileIO, path.HasRegFilePath, ABC):
    pass


class BrickleIO(FileIO, ABC):
    """ Brotli-compressed file of a pickled object (brickle). """

    @classmethod
    def load(cls, file: Path, **kwargs):
        """ Load from a compressed pickle file. """
        return load_brickle(file, data_type=cls, **kwargs)

    def save(self, top: Path, *args, **kwargs):
        """ Save to a pickle file compressed with Brotli. """
        save_path = self.get_path(top)
        checksum = save_brickle(self, save_path, *args, **kwargs)
        return save_path, checksum

    def __getstate__(self):
        # Copy the __dict__ to avoid modifying this object's state.
        state = self.__dict__.copy()
        # Do not pickle cached properties.
        for name, value in vars(type(self)).items():
            if isinstance(value, cached_property):
                state.pop(name, None)
        return state

    def __setstate__(self, state: dict[str, Any]):
        self.__dict__.update(state)


class SampleBrickleIO(SampleFileIO, BrickleIO, ABC):

    def __init__(self,
                 *args,
                 sample: str,
                 branches: dict[str, str],
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.sample = sample
        self.branches = branches


class RefBrickleIO(SampleBrickleIO, RefFileIO, ABC):

    def __init__(self,
                 *args,
                 ref: str,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.ref = ref


class RegBrickleIO(RefBrickleIO, RegFileIO, ABC):

    def __init__(self,
                 *args,
                 reg: str,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.reg = reg
