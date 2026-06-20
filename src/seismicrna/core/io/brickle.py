import pickle
from abc import ABC
from functools import cached_property
from pathlib import Path
from typing import Any

from .checksum import BadChecksumError, calc_sha512_bytes
from .file import (
    FileIO,
    SampleFileIO,
    RefFileIO,
    RegFileIO,
    DEFAULT_BROTLI_LEVEL,
    PICKLE_PROTOCOL,
)
from ..logs import logger, format_sample_reference_region
from ..validate import require_isinstance, require_issubclass
from ..write import write_mode


class BrickleIO(FileIO, ABC):
    """Brotli-compressed file of a pickled object (brickle)."""

    @classmethod
    def load(cls, file: Path, **kwargs):
        """Load from a compressed pickle file."""
        return load_brickle(file, data_type=cls, **kwargs)

    def save(self, top: Path, *args, **kwargs):
        """Save to a pickle file compressed with Brotli."""
        save_path = self.get_path(top)
        checksum = save_brickle(self, save_path, *args, **kwargs)
        return save_path, checksum

    def __getstate__(self):
        # To reduce storage space, remove cached_property attributes.
        # Use getattr(cls) rather than getattr(self) because the latter
        # will return the value that is cached, while the former will
        # return the cached_property object that caches the value.
        cls = type(self)
        return {
            name: value
            for name, value in self.__dict__.items()
            if not isinstance(getattr(cls, name, None), cached_property)
        }

    def __setstate__(self, state: dict[str, Any]):
        # All BrickleIO objects have a __dict__ rather than __slots__.
        # This method assumes that state has the correct attributes
        # because there is no easy, general way to verify that.
        self.__dict__.update(state)


class SampleBrickleIO(SampleFileIO, BrickleIO, ABC):
    def __init__(self, *args, sample: str, branches: dict[str, str], **kwargs):
        super().__init__(*args, **kwargs)
        self.sample = sample
        self.branches = branches


class RefBrickleIO(SampleBrickleIO, RefFileIO, ABC):
    def __init__(self, *args, ref: str, **kwargs):
        super().__init__(*args, **kwargs)
        self.ref = ref

    def __str__(self):
        srr = format_sample_reference_region(self.sample, self.ref)
        return f"{type(self).__name__} of {srr}"


class RegBrickleIO(RefBrickleIO, RegFileIO, ABC):
    def __init__(self, *args, reg: str, **kwargs):
        super().__init__(*args, **kwargs)
        self.reg = reg

    def __str__(self):
        srr = format_sample_reference_region(self.sample, self.ref, self.reg)
        return f"{type(self).__name__} of {srr}"


def save_brickle(
    item: BrickleIO,
    file: str | Path,
    brotli_level: int = DEFAULT_BROTLI_LEVEL,
    force: bool = False,
):
    """Pickle an object, compress with Brotli, and save to a file.

    Parameters
    ----------
    item: BrickleIO
        Object to pickle and compress.
    file: str | Path
        Path at which to save the compressed pickle file.
    brotli_level: int
        Brotli compression quality (0–11); higher means better
        compression but slower.
    force: bool = False
        Overwrite the file if it already exists.

    Returns
    -------
    str
        SHA-512 checksum of the written data.
    """
    import brotli

    require_isinstance("item", item, BrickleIO)
    with logger.debug.single_context("Writing {} to {}", item, file):
        # Save the item's state rather than the item itself.
        state = item.__getstate__()
        logger.trace("State attributes of {}: {}", item, list(state))
        with logger.trace.single_context(
            "Compressing {} with Brotli level {}", item, brotli_level
        ):
            data = brotli.compress(
                pickle.dumps(state, protocol=PICKLE_PROTOCOL), quality=brotli_level
            )
        with open(file, write_mode(force, binary=True)) as f:
            f.write(data)
        logger.debug("Wrote {} to {}", item, file)
        with logger.trace.single_context("Computing SHA-512 digest of {}", file):
            checksum = calc_sha512_bytes(data)
            logger.trace("SHA-512 digest of {} is {}", file, checksum)
    return checksum


def load_brickle(file: str | Path, data_type: type[BrickleIO], checksum: str):
    """Unpickle and return an object from a Brotli-compressed file.

    Parameters
    ----------
    file: str | Path
        Path to the Brotli-compressed pickle file.
    data_type: type[BrickleIO]
        Expected type of the loaded object; must be a subclass of
        BrickleIO.
    checksum: str
        Expected SHA-512 checksum of the file's raw bytes; pass an
        empty string to skip checksum verification.

    Returns
    -------
    BrickleIO
        The loaded object.
    """
    import brotli

    require_issubclass("data_type", data_type, BrickleIO)
    with logger.debug.single_context("Loading {} from {}", data_type, file):
        with open(file, "rb") as f:
            data = f.read()
        if checksum:
            with logger.trace.single_context("Computing SHA-512 digest of {}", file):
                sha512_digest = calc_sha512_bytes(data)
                logger.trace("SHA-512 digest of {} is {}", file, sha512_digest)
            if sha512_digest != checksum:
                raise BadChecksumError(
                    f"Expected SHA-512 digest of {file} to be {checksum}, "
                    f"but got {sha512_digest}"
                )
        else:
            logger.trace(
                "Skipping the SHA-512 digest of {} because no checksum was given", file
            )
        state = pickle.loads(brotli.decompress(data))
        logger.trace(
            "{} contains an object of type {}", file, repr(type(state).__name__)
        )
        if isinstance(state, data_type):
            item = state
            state = item.__getstate__()
        elif isinstance(state, dict):
            item = object.__new__(data_type)
            item.__setstate__(state)
        else:
            raise TypeError(f"Expected to unpickle {data_type}, but got {type(state)}")
        logger.trace("State attributes of {}: {}", item, list(state))
    return item
