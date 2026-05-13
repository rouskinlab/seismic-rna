from abc import ABC

from ..core import path
from ..core.io import RegFileIO


class EnsemblesFile(path.HasRegFilePath, ABC):

    @classmethod
    def get_step(cls):
        return path.ENSEMBLES_STEP


class EnsemblesIO(EnsemblesFile, RegFileIO, ABC):
    pass
