from abc import ABC

from ..core import path
from ..core.io.file import RegFileIO


class ClusterScanFile(path.HasRegFilePath, ABC):
    @classmethod
    def get_step(cls):
        return path.CLUSTERSCAN_STEP


class ClusterScanIO(ClusterScanFile, RegFileIO, ABC):
    pass
