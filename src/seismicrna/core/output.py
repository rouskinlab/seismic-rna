from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Iterable

from . import path


class Output(ABC):
    """ Abstract base class for an output item. """

    @classmethod
    @abstractmethod
    def load(cls, file: Path) -> Output:
        """ Load an object from a file. """

    @classmethod
    @abstractmethod
    def dir_seg_types(cls) -> tuple[path.Segment, ...]:
        """ Type(s) of the directory segment(s) in the path. """
        return tuple()

    @classmethod
    @abstractmethod
    def file_seg_type(cls) -> path.Segment:
        """ Type of the last segment in the path. """

    @classmethod
    def seg_types(cls):
        return cls.dir_seg_types() + (cls.file_seg_type(),)

    @classmethod
    def auto_fields(cls) -> dict[str, Any]:
        """ Fields that are filled automatically. """
        if len(exts := cls.file_seg_type().exts) != 1:
            raise ValueError(f"Expected exactly one file extension, "
                             f"but got {cls.file_seg_type().exts}")
        return {path.EXT: exts[0]}

    @classmethod
    def build_path(cls, **path_fields):
        """ Build the file path from the given fields. """
        return path.buildpar(*cls.seg_types(),
                             **{**cls.auto_fields(),
                                **path_fields})

    def path_fields(self, top: Path, exclude: Iterable[str] = ()):
        """ Return the path fields as a dict. """
        return {"top": top} | {field: (getattr(self, field)
                                       if hasattr(self, field)
                                       else self.auto_fields()[field])
                               for segment in self.seg_types()
                               for field in segment.field_types
                               if field not in exclude}

    def get_path(self, top: Path):
        """ Return the file path. """
        return self.build_path(**self.path_fields(top))

    @abstractmethod
    def save(self, top: Path, **kwargs) -> Path:
        """ Save the object to a file. """


class RefOutput(Output, ABC):

    @classmethod
    def dir_seg_types(cls):
        return super().dir_seg_types() + (path.SampSeg,
                                          path.CmdSeg,
                                          path.RefSeg)


class SectOutput(RefOutput, ABC):

    @classmethod
    def dir_seg_types(cls):
        return super().dir_seg_types() + (path.SectSeg,)
