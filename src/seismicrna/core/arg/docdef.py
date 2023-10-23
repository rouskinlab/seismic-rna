from functools import wraps
from inspect import getmembers, Parameter, Signature
from textwrap import dedent
from typing import Any, Callable

from click import Argument, Option

from . import cli as cpar


# Ignore special parameters with reserved names.
reserved_params = "self", "cls"

# Get every parameter defined for the command line interface.
cli_args = dict(getmembers(cpar, lambda member: isinstance(member, Argument)))
cli_opts = dict(getmembers(cpar, lambda member: isinstance(member, Option)))

# Get the default value for every parameter.
api_defs = {"n_procs": cpar.NUM_CPUS}
cli_defs = {param.name: param.default
            for param in (cli_args | cli_opts).values()
            if param.default is not None}
all_defs = cli_defs | api_defs

# Get the documentation for every CLI option.
cli_docs = {option.name: option.help for option in cli_opts.values()}


def get_param_default(param: Parameter,
                      defaults: dict[str, Any],
                      exclude_defs: tuple[str, ...]) -> Parameter:
    """ Return the parameter, possibly with a new default value. """
    if param.name in exclude_defs:
        return param
    if param.name in reserved_params:
        return param
    try:
        default = defaults[param.name]
    except KeyError:
        return param
    # Return a copy of the parameter with a new default value.
    return param.replace(default=default)


def paramdef(defaults: dict[str, Any], exclude_defs: tuple[str, ...]):
    """ Give the keyword argments of a function default values. """
    if defaults is None:
        defaults = dict()

    def decorator(func: Callable):
        # List all the parameters of the function, replacing the default
        # value of those parameters with defaults given in defaults.
        sig = Signature.from_callable(func)
        new_params = [get_param_default(param, defaults, exclude_defs)
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
                          if param.kind == Parameter.KEYWORD_ONLY}

        @wraps(func)
        def new_func(*args, **kwargs):
            return func(*args, **{**default_kwargs, **kwargs})
        return new_func

    return decorator


def autodef(extra_defs: dict[str, Any] | None = None,
            exclude_defs: tuple[str, ...] = ()):
    """ Call `paramdef` and automatically infer default values from
    the CLI and API. Extra defaults (if needed) may be given as keyword
    arguments. """
    return paramdef(all_defs if extra_defs is None else all_defs | extra_defs,
                    exclude_defs)


def get_param_lines(func: Callable, param_docs: dict[str, str]):
    sig = Signature.from_callable(func)
    # Add information about every parameter to the docstring.
    param_lines = list()
    for name, param in sig.parameters.items():
        if name in reserved_params:
            # Ignore reserved parameters (if any).
            continue
        if doc := param_docs.get(name):
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
            doc = (f"{doc} [{param.kind.description}]"
                   if param.default is param.empty
                   else f"{doc} [{param.kind.description}, "
                        f"default: {repr(param.default)}]")
            # Add the parameter's name, type, kind, and documentation to
            # the docstring.
            param_lines.extend([f"{name_type}",
                                f"    {doc}"])
    return param_lines


def get_doc_lines(func: Callable, param_lines: list[str], return_doc: str):
    sig = Signature.from_callable(func)
    doc_lines = list()
    if func.__doc__:
        # Use the existing docstring to start the new docstring.
        doc_lines.append(dedent(func.__doc__))
    if param_lines:
        if doc_lines:
            doc_lines.append("")
        doc_lines.extend(["Parameters",
                          "----------"])
        doc_lines.extend(param_lines)
    if sig.return_annotation is not sig.empty:
        if return_doc:
            if doc_lines:
                doc_lines.append("")
            doc_lines.extend(["Return",
                              "------",
                              f"{sig.return_annotation}",
                              f"    {return_doc}"])
    return doc_lines


def paramdoc(param_docs: dict[str, str], return_doc: str):
    """
    Give a function a new docstring where each parameter gets annotated
    with the text given in the keyword arguments (`param_docs`).

    Parameters
    ----------
    param_docs: dict[str, str]
        Description of each parameter, keyed by name
    return_doc: str
        Description of the return value; will occur at end of docstring

    Return
    ------
    callable
        Function with new docstring
    """

    def decorator(func: Callable):
        param_lines = get_param_lines(func, param_docs)
        doc_lines = get_doc_lines(func, param_lines, return_doc)
        func.__doc__ = "\n".join(doc_lines)
        return func

    return decorator


def autodoc(extra_docs: dict[str, str] | None = None, return_doc: str = ""):
    """ Call `paramdoc` and automatically infer descriptions and
    type annotations about all parameters from the CLI and API.
    Documentation of any extra parameters may also be given. """
    return paramdoc(cli_docs if extra_docs is None else cli_docs | extra_docs,
                    return_doc)


def auto(*,
         extra_defs: dict[str, Any] | None = None,
         exclude_defs: tuple[str, ...] = (),
         extra_docs: dict[str, str] | None = None,
         return_doc: str = ""):
    """ Combine `autodef` and `autodoc`, in that order. """

    def decorator(func: Callable):
        func = autodef(extra_defs, exclude_defs)(func)
        func = autodoc(extra_docs, return_doc)(func)
        return func

    return decorator

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
