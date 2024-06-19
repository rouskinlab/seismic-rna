from typing import Callable

from .arg import docdef
from .logs import log_exceptions
from .tmp import with_tmp_dir


def run_func(logging_method: Callable,
             with_tmp: bool = False,
             pass_keep_tmp: bool = False,
             *args,
             **kwargs):
    """ Decorator for a run function. """

    docdef_decorator = docdef.auto(*args, **kwargs)
    tmp_decorator = with_tmp_dir(pass_keep_tmp) if with_tmp else None
    exceptions_decorator = log_exceptions(logging_method)

    def decorator(func: Callable):
        # Apply each decorator to the run function.
        func = docdef_decorator(func)
        if tmp_decorator is not None:
            func = tmp_decorator(func)
        func = exceptions_decorator(func)
        return func

    return decorator
