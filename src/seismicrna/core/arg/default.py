from inspect import getmembers

from click import Argument, Option

from . import cli

# Get every parameter defined for the command line interface.
cli_args = dict(getmembers(cli, lambda member: isinstance(member, Argument)))
cli_opts = dict(getmembers(cli, lambda member: isinstance(member, Option)))

# Get the default value for every parameter.
cli_defaults = {param.name: param.default
                for param in (cli_args | cli_opts).values()
                if param.default is not None}
