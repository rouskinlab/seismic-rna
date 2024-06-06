import os
from functools import wraps
from logging import getLogger
from pathlib import Path
from shutil import rmtree
from typing import Callable

logger = getLogger(__name__)

# Lock directory (for function lock_output)
LOCK_DIR = ".seismic-rna_lock"


def lock_tmp_dir(run: Callable):
    @wraps(run)
    def wrapper(*args, tmp_dir: str | Path, keep_tmp: bool, **kwargs):
        lock_error = (f"The directory {tmp_dir} is currently being used by "
                      "another instance of SEISMIC-RNA. If SEISMIC-RNA is not "
                      "actually running, then please delete the directory with "
                      f"'rm -r {tmp_dir}'. Otherwise, please name another "
                      "directory with '-t /new/tmp/dir'.")
        tmp_error = (f"The directory {tmp_dir} exists. Please either delete "
                     f"it with 'rm -r {tmp_dir}' or name another directory "
                     f"with '-t /new/tmp/dir'.")
        # Determine whether the temporary directory and the lock exist.
        lock = os.path.join(tmp_dir, LOCK_DIR)
        try:
            os.mkdir(lock)
        except FileExistsError:
            # The lock already exists, which means another instance of
            # SEISMIC-RNA is using this temporary directory.
            logger.critical(lock_error)
            raise SystemExit()
        except FileNotFoundError:
            # The temporary directory does not exist yet, so create it
            # along with a lock.
            try:
                os.makedirs(lock, exist_ok=False)
            except FileExistsError:
                # If this error happens, it is due to a very unlikely
                # race condition wherein another instance of SEISMIC-RNA
                # raises a FileNotFoundError from os.mkdir(lock), then
                # this instance of SEISMIC-RNA does the same, then the
                # first run makes the directory with os.makedirs(lock),
                # and then this run tries to do the same thing but fails
                # because the directory was created moments before.
                logger.critical(lock_error)
                raise SystemExit()
            tmp_dir_existed_before = False
            logger.debug(f"Created and locked temporary directory: {tmp_dir}")
        else:
            # The temporary directory had existed, but the lock had not.
            tmp_dir_existed_before = True
            logger.debug(f"Locked temporary directory: {tmp_dir}")
        # The lock now exists, so any other instance of SEISMIC-RNA that
        # tries to use the same lock will exit before it can use the
        # temporary directory or delete the lock. Thus, this run must
        # delete the lock upon exiting, regardless of the circumstances.
        try:
            if tmp_dir_existed_before and not keep_tmp:
                logger.critical(tmp_error)
                raise SystemExit()
            try:
                # Run the wrapped function and return its result.
                return run(*args, **kwargs,
                           tmp_dir=tmp_dir,
                           keep_tmp=keep_tmp)
            finally:
                # Delete the temporary directory unless the option to
                # save it was enabled.
                if not keep_tmp:
                    rmtree(tmp_dir, ignore_errors=True)
                    logger.debug(f"Deleted temporary directory: {tmp_dir}")
        finally:
            # Always ensure that the temporary directory is unlocked
            # upon exiting.
            try:
                os.rmdir(lock)
                logger.debug(f"Unlocked temporary directory: {tmp_dir}")
            except FileNotFoundError:
                pass

    # Return the decorator.
    return wrapper

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
