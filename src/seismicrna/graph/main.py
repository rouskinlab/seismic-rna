from click import group

from . import (abundance,
               aucroll,
               corroll,
               delprof,
               giniroll,
               histpos,
               histread,
               mutdist,
               poscorr,
               profile,
               roc,
               scatter,
               snrroll)
from ..core.arg import CMD_GRAPH


# Group for all graph commands
@group(CMD_GRAPH)
def cli():
    """ Graph and compare data from tables and/or structures. """


# Add graph commands to the CLI.
for module in (abundance,
               aucroll,
               corroll,
               delprof,
               giniroll,
               histpos,
               histread,
               mutdist,
               poscorr,
               profile,
               roc,
               scatter,
               snrroll):
    cli.add_command(module.cli)
