from click import group

from . import seqbar, seqcorr, seqdiff, scatter
from ..core.cmd import CMD_GRAPH


# Group for all graph commands
@group(CMD_GRAPH)
def cli():
    """ Graphing command line interface """


# Add graph commands to the CLI.
cli.add_command(seqbar.cli)
cli.add_command(seqcorr.cli)
cli.add_command(seqdiff.cli)
cli.add_command(scatter.cli)
