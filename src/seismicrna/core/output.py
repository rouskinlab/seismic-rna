from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Any, Iterable, TypeVar

from . import path
from .files import DEFAULT_BROTLI_LEVEL, load_pkl_br, save_pkl_br

logger = getLogger(__name__)


AnyOutput = TypeVar("AnyOutput")


class Output(ABC):
    """ Abstract base class for an output item. """

    @classmethod
    @abstractmethod
    def load(cls: type[AnyOutput], file: Path) -> AnyOutput:
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
                             **(cls.auto_fields() | path_fields))

    def path_fields(self, top: Path | None = None, exclude: Iterable[str] = ()):
        """ Return the path fields as a dict. """
        fields = {path.TOP: top} if top else dict()
        fields.update({field: (getattr(self, field) if hasattr(self, field)
                               else self.auto_fields()[field])
                       for segment in self.seg_types()
                       for field in segment.field_types})
        for field in exclude:
            fields.pop(field, None)
        return fields

    def get_path(self, top: Path):
        """ Return the file path. """
        return self.build_path(**self.path_fields(top))

    @abstractmethod
    def save(self, top: Path, **kwargs):
        """ Save the object to a file. """


class RefOutput(Output, ABC):

    @classmethod
    def dir_seg_types(cls):
        return super().dir_seg_types() + (path.SampSeg,
                                          path.CmdSeg,
                                          path.RefSeg)

    def __init__(self, sample: str, ref: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sample = sample
        self.ref = ref


class SectOutput(RefOutput, ABC):

    @classmethod
    def dir_seg_types(cls):
        return super().dir_seg_types() + (path.SectSeg,)

    def __init__(self, sect: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sect = sect


class PickleOutput(Output, ABC):

    @classmethod
    def load(cls, file: Path, checksum: str | None = None):
        """ Load from a compressed pickle file. """
        return load_pkl_br(file, check_type=cls, checksum=checksum)

    def save(self,
             top: Path,
             brotli_level: int = DEFAULT_BROTLI_LEVEL,
             overwrite: bool = False):
        """ Save to a pickle file compressed with Brotli. """
        checksum = save_pkl_br(self,
                               save_path := self.get_path(top),
                               brotli_level=brotli_level,
                               overwrite=overwrite)
        return save_path, checksum

    def __getstate__(self):
        # Copy the __dict__ to avoid modifying this object's state.
        state = self.__dict__.copy()
        # Do not pickle cached properties.
        for name, value in list(state.items()):
            if isinstance(value, cached_property):
                state.pop(name)
        return state

    def __setstate__(self, state: dict[str, Any]):
        self.__dict__.update(state)
