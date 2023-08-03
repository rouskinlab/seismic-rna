from click import group

from . import seqbar as seqbar_mod
from . import seqcorr as seqcorr_mod
from . import seqdiff as seqdiff_mod
from . import scatter as scatter_mod
from ..core import path


# Group for all graph commands
@group(path.MOD_GRAPH)
def cli():
    """ Graphing command line interface """


# Add graph commands to the CLI.
cli.add_command(seqbar_mod.cli)
cli.add_command(seqcorr_mod.cli)
cli.add_command(seqdiff_mod.cli)
cli.add_command(scatter_mod.cli)
