from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger
from pathlib import Path
from typing import Any, Iterable

from .util import load_pkl_br, save_pkl_br
from .. import path

logger = getLogger(__name__)

DEFAULT_BROTLI_LEVEL = 10
PICKLE_PROTOCOL = 5


class FileIO(ABC):
    """ Any file saved by SEISMIC-RNA, rather than by a dependency. """

    @classmethod
    @abstractmethod
    def load(cls, file: Path):
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
    def path_field_names(cls):
        """ Names of the fields of the path for the file type. """
        return tuple(field for segment in cls.seg_types()
                     for field in segment.field_types)

    @classmethod
    def auto_fields(cls) -> dict[str, Any]:
        """ Names and automatic values of selected fields. """
        try:
            return {path.EXT: cls.file_seg_type().exts[0]}
        except IndexError:
            raise ValueError(f"Got no file extensions for {cls.__name__}")

    @classmethod
    def parse_path(cls, file: Path):
        """ Parse a file path to determine the field values. """
        fields = path.parse(file, *cls.seg_types())
        return fields.pop(path.TOP), fields

    @classmethod
    def build_path(cls, **path_fields):
        """ Build the file path from the given field values. """
        return path.buildpar(*cls.seg_types(),
                             **(cls.auto_fields() | path_fields))

    @classmethod
    def normalize_fields(cls, **fields):
        """ Given arbitrary fields and values, process those for this
        class of file, giving preference to auto-fields. """
        return {field: cls.auto_fields().get(field, fields[field])
                for field in cls.path_field_names()}

    def path_fields(self, top: Path | None = None, exclude: Iterable[str] = ()):
        """ Return the path fields as a dict. """
        fields = {path.TOP: top} if top else dict()
        fields.update({field: (getattr(self, field) if hasattr(self, field)
                               else self.auto_fields()[field])
                       for field in self.path_field_names()})
        for field in exclude:
            fields.pop(field, None)
        return fields

    def get_path(self, top: Path):
        """ Return the file path. """
        return self.build_path(**self.path_fields(top))

    @abstractmethod
    def save(self, top: Path, **kwargs):
        """ Save the object to a file. """


class RefIO(FileIO, ABC):
    """ Saved file with a sample, command, and reference. """

    @classmethod
    def dir_seg_types(cls):
        return super().dir_seg_types() + (path.SampSeg,
                                          path.CmdSeg,
                                          path.RefSeg)

    def __init__(self, *, sample: str, ref: str, **kwargs):
        super().__init__(**kwargs)
        self.sample = sample
        self.ref = ref


class SectIO(RefIO, ABC):
    """ File with a section of a reference. """

    @classmethod
    def dir_seg_types(cls):
        return super().dir_seg_types() + (path.SectSeg,)

    def __init__(self, *, sect: str, **kwargs):
        super().__init__(**kwargs)
        self.sect = sect


class BrickleIO(FileIO, ABC):
    """ Brotli-compressed file of a Pickled object (Brickle). """

    @classmethod
    def load(cls, file: Path, checksum: str = ""):
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


def convert_path(file: Path, type1: type[FileIO], type2: type[FileIO]):
    """ Convert a path from that used by file type 1 to file type 2. """
    # Extract the fields from the path using file type 1.
    top, fields = type1.parse_path(file)
    # Normalize the fields to comply with file type 2.
    norm_fields = type2.normalize_fields(**fields)
    # Generate a new path for file type 2 from the normalized fields.
    return type2.build_path(top=top, **norm_fields)
