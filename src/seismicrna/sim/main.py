from click import group

from . import (abstract,
               clusts,
               ends,
               fastq,
               fold,
               muts,
               params,
               ref,
               relate,
               total)
from ..core.arg import CMD_SIM


# Group for all sim commands
@group(CMD_SIM)
def cli():
    """ Simulate sequences, structures, and samples. """


# Add simulation commands to the CLI.
for module in [abstract,
               clusts,
               ends,
               fastq,
               fold,
               muts,
               params,
               ref,
               relate,
               total]:
    cli.add_command(module.cli)
