from abc import ABC, abstractmethod
from pathlib import Path

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
