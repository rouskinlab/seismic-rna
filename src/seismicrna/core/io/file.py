from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Any, Iterable

from .brickle import load_brickle, save_brickle
from .. import path

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
    def dir_seg_types(cls) -> tuple[path.PathSegment, ...]:
        """ Types of the directory segments in the path. """
        return tuple()

    @classmethod
    @abstractmethod
    def file_seg_type(cls) -> path.PathSegment:
        """ Type of the last segment in the path. """

    @classmethod
    def seg_types(cls):
        """ Types of the segments in the path. """
        return cls.dir_seg_types() + (cls.file_seg_type(),)

    @classmethod
    def path_fields(cls):
        """ Path fields for the file type. """
        return path.get_fields_in_seg_types(cls.seg_types())

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
        return path.parse_top_separate(file, cls.seg_types())

    @classmethod
    def build_path(cls, path_fields: dict[str, Any]):
        """ Build the file path from the given field values. """
        return path.buildpar(cls.seg_types(), (cls.auto_fields() | path_fields))

    def path_field_values(self,
                          top: Path | None = None,
                          exclude: Iterable[str] = ()):
        """ Path field values as a dict. """
        fields = {path.TOP: top} if top else dict()
        fields.update({field: (getattr(self, field) if hasattr(self, field)
                               else self.auto_fields()[field])
                       for field in self.path_fields()})
        for field in exclude:
            fields.pop(field, None)
        return fields

    def get_path(self, top: Path):
        """ Return the file path. """
        return self.build_path(self.path_field_values(top))

    @abstractmethod
    def save(self, top: Path, **kwargs):
        """ Save the object to a file. """

    def __str__(self):
        return type(self).__name__


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


class RegIO(RefIO, ABC):
    """ File with a region of a reference. """

    @classmethod
    def dir_seg_types(cls):
        return super().dir_seg_types() + (path.RegSeg,)

    def __init__(self, *, reg: str, **kwargs):
        super().__init__(**kwargs)
        self.reg = reg


class BrickleIO(FileIO, ABC):
    """ Brotli-compressed file of a Pickled object (Brickle). """

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


def recast_file_path(input_path: Path,
                     input_type: type[FileIO],
                     output_type: type[FileIO]):
    """ Recast `input_path` from `input_type` to `output_type`.

    Parameters
    ----------
    input_path: Path
        Input path from which to take the path fields.
    input_type: type[FileIO]
        Type of file to use to determine the fields in `input_path`.
    output_type: type[FileIO]
        Type of file to use for fitting the path fields and building a
        new path.

    Returns
    -------
    Path
        Path for file of `output_type` made from fields in `input_path`
        (as determined by the file of `input_type`).
    """
    return path.cast_path(input_path,
                          input_type.seg_types(),
                          output_type.seg_types(),
                          output_type.auto_fields())
