from typing import Callable, Optional

from .arg import docdef
from .logs import log_exceptions
from .tmp import with_tmp_dir


def run_func(logging_method: Callable,
             default: Optional[Callable] = list,
             with_tmp: bool = False,
             pass_keep_tmp: bool = False,
             *args,
             **kwargs):
    """ Decorator for a run function. """

    tmp_decorator = with_tmp_dir(pass_keep_tmp) if with_tmp else None
    docdef_decorator = docdef.auto(*args, **kwargs)
    exceptions_decorator = log_exceptions(logging_method, default)

    def decorator(func: Callable):
        # Apply each decorator to the run function.
        if tmp_decorator is not None:
            func = tmp_decorator(func)
        func = docdef_decorator(func)
        func = exceptions_decorator(func)
        return func

    return decorator
