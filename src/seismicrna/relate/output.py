from abc import ABC

from ..core import path
from ..core.cmd import CMD_REL
from ..core.output import RefOutput


class RelateOutput(RefOutput, ABC):

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.CMD: CMD_REL}
