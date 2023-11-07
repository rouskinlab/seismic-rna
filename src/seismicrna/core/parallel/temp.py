import os
from functools import wraps
from logging import getLogger
from pathlib import Path
from shutil import rmtree
from typing import Callable

logger = getLogger(__name__)

# Lock directory (for function lock_output)
LOCK_DIR = ".seismic-rna_lock"


def lock_temp_dir(run: Callable):
    @wraps(run)
    def wrapper(*args, temp_dir: str | Path, keep_temp: bool, **kwargs):
        lock_error = (f"The directory {temp_dir} is currently being used by "
                      f"another instance of SEISMIC-RNA. If possible, please "
                      f"name a temporary directory that does not yet exist "
                      f"with '--temp-dir /path/to/new/temp-dir/'. If a former "
                      f"run of SEISMIC-RNA failed to unlock this directory, "
                      f"then please delete it with 'rm -r {temp_dir}'.")
        temp_error = (f"The temporary directory {temp_dir} exists. If any "
                      f"needed files reside in {temp_dir}, then please "
                      f"specify a nonexistent temporary directory with "
                      f"'--temp-dir /new/temp/dir'. Otherwise, please delete "
                      f"the directory with 'rm -r {temp_dir}' and then force.")
        # Determine whether the temporary directory and the lock exist.
        lock = os.path.join(temp_dir, LOCK_DIR)
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
            temp_dir_existed_before = False
            logger.debug(f"Created and locked temporary directory: {temp_dir}")
        else:
            # The temporary directory had existed, but the lock had not.
            temp_dir_existed_before = True
            logger.debug(f"Locked temporary directory: {temp_dir}")
        # The lock now exists, so any other instance of SEISMIC-RNA that
        # tries to use the same lock will exit before it can use the
        # temporary directory or delete the lock. Thus, this run must
        # delete the lock upon exiting, regardless of the circumstances.
        try:
            if temp_dir_existed_before and not keep_temp:
                logger.critical(temp_error)
                raise SystemExit()
            try:
                # Run the wrapped function and return its result.
                return run(*args, **kwargs,
                           temp_dir=temp_dir,
                           keep_temp=keep_temp)
            finally:
                # Delete the temporary directory unless the option to
                # save it was enabled.
                if not keep_temp:
                    rmtree(temp_dir, ignore_errors=True)
                    logger.debug(f"Deleted temporary directory: {temp_dir}")
        finally:
            # Always ensure that the temporary directory is unlocked
            # upon exiting.
            try:
                os.rmdir(lock)
                logger.debug(f"Unlocked temporary directory: {temp_dir}")
            except FileNotFoundError:
                pass

    # Return the decorator.
    return wrapper

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
