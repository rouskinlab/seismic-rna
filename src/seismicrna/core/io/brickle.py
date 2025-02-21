import pickle
from abc import ABC
from functools import cached_property
from hashlib import sha512
from pathlib import Path
from typing import Any

import brotli

from .checksum import BadChecksumError, calc_sha512_digest
from .file import FileIO, SampleFileIO, RefFileIO, RegFileIO
from ..logs import logger
from ..write import write_mode

DEFAULT_BROTLI_LEVEL = 10
PICKLE_PROTOCOL = 5



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
        # Copy self.__dict__ to avoid modifying this object's state.
        state = self.__dict__.copy()
        # Do not pickle cached properties.
        for name, value in vars(type(self)).items():
            if isinstance(value, cached_property):
                state.pop(name, None)
        return state

    def __setstate__(self, state: dict[str, Any]):
        # All BrickleIO objects have a __dict__ rather than __slots__.
        # This method assumes that state has the correct attributes
        # because there is no easy, general way to verify that.
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


def save_brickle(item: BrickleIO,
                 file: str | Path,
                 brotli_level: int = DEFAULT_BROTLI_LEVEL,
                 force: bool = False):
    """ Pickle an object, compress with Brotli, and save to a file. """
    if not isinstance(item, BrickleIO):
        raise TypeError(
            f"item must be {BrickleIO.__name__}, but got {type(item).__name__}"
        )
    logger.routine(f"Began writing {item} to {file}")
    # Save the item's state rather than the item itself.
    state = item.__getstate__()
    logger.detail(f"Began compressing {item} with Brotli level {brotli_level}")
    data = brotli.compress(pickle.dumps(state, protocol=PICKLE_PROTOCOL),
                           quality=brotli_level)
    logger.detail(f"Ended compressing {item} with Brotli level {brotli_level}")
    with open(file, write_mode(force, binary=True)) as f:
        f.write(data)
    logger.action(f"Wrote {item} to {file}")
    checksum = calc_sha512_digest(data)
    logger.detail(f"Computed SHA-512 checksum of {file}: {checksum}")
    logger.routine(f"Ended writing {item} to {file}")
    return checksum


def load_brickle(file: str | Path,
                 data_type: type[BrickleIO],
                 checksum: str):
    """ Unpickle and return an object from a Brotli-compressed file. """
    logger.routine(f"Began loading {data_type} from {file}")
    if not issubclass(data_type, BrickleIO):
        raise TypeError(f"data_type must be subclass of {BrickleIO.__name__}, "
                        f"but got {type(data_type).__name__}")
    with open(file, "rb") as f:
        data = f.read()
    if checksum:
        sha512_digest = calc_sha512_digest(data)
        if sha512_digest != checksum:
            raise BadChecksumError(
                f"Expected SHA-512 digest of {file} to be {checksum}, "
                f"but got {sha512_digest}"
            )
    state = pickle.loads(brotli.decompress(data))
    logger.detail(f"{file} contains {type(state)}")
    if isinstance(state, data_type):
        item = state
    elif isinstance(state, dict):
        item = object.__new__(data_type)
        item.__setstate__(state)
    else:
        raise TypeError(f"Expected to unpickle {data_type}, "
                        f"but got {type(state)}")
    logger.routine(f"Ended loading {data_type} from {file}")
    return item
