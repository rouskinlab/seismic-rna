from functools import wraps
from inspect import Parameter, Signature
from pathlib import Path
from shutil import rmtree
from typing import Callable

from .logs import logger
from .path import randdir, sanitize, transpath

PENDING = "release"
WORKING = "working"


def release_to_out(out_dir: Path,
                   release_dir: Path,
                   initial_path: Path):
    """ Move temporary path(s) to the output directory. """
    # Determine the path in the output directory.
    out_path = transpath(out_dir, release_dir, initial_path)
    if initial_path.exists():
        # Ensure the parent directory of the new path exists.
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # If the output path already exists, then first rename it.
        delete_path = randdir(out_path.parent, f"{out_path.name}-")
        try:
            out_path.rename(delete_path)
        except FileNotFoundError:
            # The output path does not yet exist.
            deleted = False
            logger.detail("Output path {} does not yet exist", out_path)
        else:
            deleted = True
            logger.detail("Moved {} to {} (to be deleted)",
                          out_path, delete_path)
        try:
            # Move the initial path to the output location.
            initial_path.rename(out_path)
            logger.detail("Moved {} to {}", initial_path, out_path)
        except Exception:
            if deleted:
                # If an error occurred, then restore the original output
                # path before raising the exception.
                delete_path.rename(out_path)
                logger.detail("Restored {} from {}", out_path, delete_path)
            else:
                # No original files were moved to the delete directory,
                # which is therefore still empty. Delete it.
                delete_path.rmdir()
            raise
        # Once the initial path has been moved to its destination, the
        # original directory can be deleted safely.
        try:
            rmtree(delete_path)
        except Exception as error:
            logger.warning("Failed to delete {}; however, it is no longer "
                           "needed, and safe to delete yourself: {}",
                           delete_path, error)
        else:
            logger.detail("Deleted {}", delete_path)
    else:
        logger.detail("Skipped releasing {} (does not exist)", initial_path)
    if not out_path.exists():
        raise FileNotFoundError(out_path)
    logger.routine("Released {} to {}", initial_path, out_path)
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
                    keep_tmp: bool,
                    **kwargs):
            tmp_dir = None
            try:
                tmp_pfx = sanitize(tmp_pfx)
                tmp_dir = randdir(tmp_pfx.parent, prefix=tmp_pfx.name)
                logger.routine("Created temporary directory: {}", tmp_dir)
                if pass_keep_tmp:
                    kwargs = dict(**kwargs, keep_tmp=keep_tmp)
                return func(*args, **kwargs, tmp_dir=tmp_dir)
            finally:
                if tmp_dir is not None and not keep_tmp:
                    try:
                        rmtree(tmp_dir)
                    except OSError as error:
                        logger.warning(
                            "Failed to delete temporary directory {}:\n{}",
                            tmp_dir, error
                        )
                    else:
                        logger.routine("Deleted temporary directory: {}",
                                       tmp_dir)

        # Add tmp_pfx and keep_tmp to the signature of the wrapper, and
        # remove tmp_dir (functools.wraps does not do so automatically).
        params = dict(Signature.from_callable(func).parameters)
        params.pop("tmp_dir")
        for param in ["tmp_pfx", "keep_tmp"]:
            if param not in params:
                params[param] = Parameter(param, Parameter.KEYWORD_ONLY)
        wrapper.__signature__ = Signature(parameters=list(params.values()))
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
