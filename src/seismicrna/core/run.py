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
            logger.status(f"Began {command}")
            result = func(*args, **kwargs)
            logger.status(f"Ended {command}")
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
