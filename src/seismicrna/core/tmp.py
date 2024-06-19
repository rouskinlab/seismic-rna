from functools import wraps
from logging import getLogger
from pathlib import Path
from shutil import move, rmtree
from typing import Callable

from .path import randdir, sanitize, transpath

logger = getLogger(__name__)

PENDING = "release"
WORKING = "working"


def release_to_out(out_dir: str | Path,
                   release_dir: str | Path,
                   initial_path: str | Path):
    """ Move temporary path(s) to the output directory. """
    # Determine the path in the output directory.
    out_path = transpath(out_dir, release_dir, initial_path)
    if initial_path.exists():
        # Delete the new path if it already exists.
        rmtree(out_path, ignore_errors=True)
        # Ensure the parent directory of the new path exists.
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Move the temporary path to the new location.
        move(initial_path, out_path.parent)
    if not out_path.exists():
        raise FileNotFoundError(out_path)
    return out_path


def get_release_working_dirs(tmp_dir: Path):
    release_dir = tmp_dir.joinpath(PENDING)
    working_dir = tmp_dir.joinpath(WORKING)
    release_dir.mkdir(parents=False, exist_ok=True)
    working_dir.mkdir(parents=False, exist_ok=True)
    return release_dir, working_dir


def with_tmp_dir(pass_keep_tmp: bool):
    """ Make a temporary directory, then delete it after the run. """

    def decorator(func: Callable):

        @wraps(func)
        def wrapper(*args,
                    tmp_pfx: str | Path,
                    keep_tmp: bool = False,
                    **kwargs):
            tmp_dir = None
            try:
                tmp_pfx = sanitize(tmp_pfx)
                tmp_dir = randdir(tmp_pfx.parent, prefix=tmp_pfx.name)
                logger.debug(f"Created temporary directory: {tmp_dir}")
                if pass_keep_tmp:
                    kwargs = dict(**kwargs, keep_tmp=keep_tmp)
                return func(*args, **kwargs, tmp_dir=tmp_dir)
            finally:
                if tmp_dir is not None and not keep_tmp:
                    try:
                        rmtree(tmp_dir)
                    except OSError as error:
                        logger.warning("Failed to delete temporary directory "
                                       f"{tmp_dir}:\n{error}")
                    else:
                        logger.debug(f"Deleted temporary directory {tmp_dir}")

        return wrapper

    return decorator

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
