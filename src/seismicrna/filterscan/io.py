from abc import ABC

from ..core import path
from ..core.io.file import RegFileIO


class FilterScanFile(path.HasRegFilePath, ABC):
    @classmethod
    def get_step(cls):
        return path.FILTERSCAN_STEP


class FilterScanIO(FilterScanFile, RegFileIO, ABC):
    pass
