from functools import wraps
from typing import Callable, Optional

from .arg import docdef
from .logs import logger, log_exceptions
from .tmp import with_tmp_dir


def log_command(command: str):
    """ Log the name of the command. """

    def decorator(func: Callable):

        @wraps(func)
        def wrapper(*args, **kwargs):
            logger.command(f"Began {command}")
            result = func(*args, **kwargs)
            logger.command(f"Ended {command}")
            return result

        return wrapper

    return decorator


def run_func(command: str,
             default: Optional[Callable] = list,
             with_tmp: bool = False,
             pass_keep_tmp: bool = False,
             *args,
             **kwargs):
    """ Decorator for a run function. """

    def decorator(func: Callable):
        # Apply each decorator to the run function.
        if with_tmp:
            func = with_tmp_dir(pass_keep_tmp)(func)
        func = docdef.auto(*args, **kwargs)(func)
        func = log_exceptions(default)(func)
        func = log_command(command)(func)
        return func

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
