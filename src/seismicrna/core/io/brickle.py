import pickle
from hashlib import md5
from logging import getLogger
from pathlib import Path
from typing import Any

import brotli

from ..write import write_mode

logger = getLogger(__name__)

DEFAULT_BROTLI_LEVEL = 10
PICKLE_PROTOCOL = 5


def digest_data(data: bytes):
    """ Compute the MD5 digest of the data as a hexadecimal number. """
    return md5(data).hexdigest()


def save_brickle(item: Any,
                 file: Path,
                 brotli_level: int = DEFAULT_BROTLI_LEVEL,
                 force: bool = False):
    """ Pickle an object, compress with Brotli, and save to a file. """
    data = brotli.compress(pickle.dumps(item, protocol=PICKLE_PROTOCOL),
                           quality=brotli_level)
    with open(file, write_mode(force, binary=True)) as f:
        f.write(data)
    logger.info(f"Wrote {item} (brotli level {brotli_level}) to {file}")
    checksum = digest_data(data)
    logger.debug(f"Computed MD5 checksum of {file}: {checksum}")
    return checksum


def load_brickle(file: Path,
                 checksum: str,
                 check_type: None | type | tuple[type, ...] = None):
    """ Unpickle and return an object from a Brotli-compressed file. """
    with open(file, "rb") as f:
        data = f.read()
    if checksum != (digest := digest_data(data)):
        raise ValueError(
            f"Expected checksum of {file} to be {checksum}, but got {digest}")
    item = pickle.loads(brotli.decompress(data))
    if check_type is not None and not isinstance(item, check_type):
        raise TypeError(f"Expected to unpickle {check_type}, "
                        f"but got {type(item).__name__}")
    logger.debug(f"Loaded {item} from {file}")
    return item

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
