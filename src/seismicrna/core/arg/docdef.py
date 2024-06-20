from functools import wraps
from inspect import Parameter, Signature
from textwrap import dedent
from typing import Any, Callable

from . import cli
from .default import cli_defaults, cli_opts

# Ignore special parameters with reserved names.
reserved_params = ["self", "cls"]

# Get the default value for every parameter.
api_defaults = dict(n_procs=cli.NUM_CPUS)
all_defaults = cli_defaults | api_defaults

# Get the documentation for every CLI option.
cli_docstrs = {option.name: option.help for option in cli_opts.values()}


def get_param_default(param: Parameter,
                      defaults: dict[str, Any],
                      exclude_defaults: tuple[str, ...]):
    """ Return the parameter, possibly with a new default value. """
    if param.name in exclude_defaults:
        return param
    if param.name in reserved_params:
        return param
    try:
        # Return a copy of the parameter with a new default value.
        return param.replace(default=defaults[param.name])
    except KeyError:
        return param


def param_defaults(defaults: dict[str, Any],
                   exclude_defaults: tuple[str, ...]):
    """ Give the keyword argments of a function default values. """

    def decorator(func: Callable):
        # List all the parameters of the function, replacing the default
        # value of those parameters with defaults given in defaults.
        sig = Signature.from_callable(func)
        new_params = [get_param_default(param, defaults, exclude_defaults)
                      for param in sig.parameters.values()]
        # Update the help text (does not affect actual default values).
        try:
            func.__signature__ = Signature(parameters=new_params)
        except ValueError as error:
            raise ValueError(f"Failed to set signature of {func.__name__} to "
                             f"({', '.join(map(str, new_params))}): {error}")

        # Update the actual default values of keyword-only arguments
        # (does not affect help text).
        default_kwargs = {param.name: param.default for param in new_params
                          if param.kind == Parameter.KEYWORD_ONLY
                          and param.default is not Parameter.empty}

        @wraps(func)
        def new_func(*args, **kwargs):
            return func(*args, **(default_kwargs | kwargs))

        return new_func

    return decorator


def auto_defaults(extra_defaults: dict[str, Any] | None = None,
                  exclude_defaults: tuple[str, ...] = ()):
    """ Call `paramdef` and automatically infer default values from
    the CLI and API. Extra defaults (if needed) may be given as keyword
    arguments. """
    return param_defaults((all_defaults if extra_defaults is None
                           else all_defaults | extra_defaults),
                          exclude_defaults)


def get_param_lines(func: Callable, docstrs: dict[str, str]):
    sig = Signature.from_callable(func)
    # Add information about every parameter to the docstring.
    param_lines = list()
    for name, param in sig.parameters.items():
        if name in reserved_params:
            # Ignore reserved parameters (if any).
            continue
        if docstr := docstrs.get(name):
            # Add the type annotation (if any) after the name of the
            # parameter.
            if param.annotation is param.empty:
                name_type = name
            else:
                try:
                    name_type = f"{name}: {param.annotation.__name__}"
                except AttributeError:
                    # Some types (e.g. UnionType) have no name.
                    name_type = f"{name}: {param.annotation}"
            # Add the kind of parameter and its default value (if any)
            # in brackets after the documentation of the parameter.
            docstr = (f"{docstr} [{param.kind.description}]"
                      if param.default is param.empty
                      else (f"{docstr} [{param.kind.description}, "
                            f"default: {repr(param.default)}]"))
            # Add the parameter's name, type, kind, and documentation to
            # the docstring.
            param_lines.extend([f"{name_type}",
                                f"    {docstr}"])
    return param_lines


def get_docstr_lines(func: Callable,
                     param_lines: list[str],
                     return_docstr: str):
    sig = Signature.from_callable(func)
    docstr_lines = list()
    if func.__doc__:
        # Use the existing docstring to start the new docstring.
        docstr_lines.append(dedent(func.__doc__))
    if param_lines:
        if docstr_lines:
            docstr_lines.append("")
        docstr_lines.extend(["Parameters",
                             "----------"])
        docstr_lines.extend(param_lines)
    if sig.return_annotation is not sig.empty:
        if return_docstr:
            if docstr_lines:
                docstr_lines.append("")
            docstr_lines.extend(["Return",
                                 "------",
                                 f"{sig.return_annotation}",
                                 f"    {return_docstr}"])
    return docstr_lines


def param_docstrs(docstrs: dict[str, str], return_docstr: str):
    """
    Give a function a new docstring where each parameter gets annotated
    with the text given in the keyword arguments (`param_docs`).

    Parameters
    ----------
    docstrs: dict[str, str]
        Description of each parameter, keyed by name
    return_docstr: str
        Description of the return value; will occur at end of docstring

    Return
    ------
    callable
        Function with new docstring
    """

    def decorator(func: Callable):
        param_lines = get_param_lines(func, docstrs)
        docstr_lines = get_docstr_lines(func, param_lines, return_docstr)
        func.__doc__ = "\n".join(docstr_lines)
        return func

    return decorator


def auto_docstrs(extra_docstrs: dict[str, str] | None = None,
                 return_docstr: str = ""):
    """ Call `param_docstrs` and automatically infer descriptions and
    type annotations about all parameters from the CLI and API.
    Documentation of any extra parameters may also be given. """
    return param_docstrs((cli_docstrs if extra_docstrs is None
                          else cli_docstrs | extra_docstrs),
                         return_docstr)


def auto(*,
         extra_defaults: dict[str, Any] | None = None,
         exclude_defaults: tuple[str, ...] = (),
         extra_docstrs: dict[str, str] | None = None,
         return_docstr: str = ""):
    """ Combine `auto_defaults` and `auto_docstrs`, in that order. """

    def decorator(func: Callable):
        func = auto_defaults(extra_defaults, exclude_defaults)(func)
        func = auto_docstrs(extra_docstrs, return_docstr)(func)
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
