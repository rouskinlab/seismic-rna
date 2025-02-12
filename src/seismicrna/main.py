from pathlib import Path

from click import group, version_option

from . import (wf,
               demult,
               align,
               relate,
               pool,
               mask,
               cluster,
               join,
               ensembles,
               table,
               lists,
               fold,
               graph,
               export,
               draw,
               migrate,
               test,
               sim,
               cleanfa,
               renumct,
               __version__)
from .align import split
from .core import path, rna
from .core.arg import (opt_exit_on_error,
                       opt_log,
                       opt_log_color,
                       opt_quiet,
                       opt_verbose)
from .core.logs import logger, set_config
from .urls import (cli_docs,
                   cli_github,
                   cli_pypi,
                   cli_conda,
                   cli_biorxiv)

params = [
    opt_verbose,
    opt_quiet,
    opt_log,
    opt_log_color,
    opt_exit_on_error,
]


# Group for main commands
@group(params=params, context_settings={"show_default": True})
@version_option(__version__)
def cli(verbose: int,
        quiet: int,
        log: str | Path,
        log_color: bool,
        exit_on_error: bool):
    """ Command line interface of SEISMIC-RNA. """
    if log:
        log_file_path = path.sanitize(log)
        log_file_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        log_file_path = None
    set_config(verbose - quiet, log_file_path, log_color, exit_on_error)
    logger.detail(f"This is SEISMIC-RNA version {__version__}")


# Add all commands to the main CLI command group.

for module in [wf,
               demult,
               split,
               align,
               relate,
               pool,
               mask,
               cluster,
               join,
               ensembles,
               table,
               lists,
               fold,
               graph,
               export,
               draw,
               migrate,
               test,
               sim,
               cleanfa,
               renumct]:
    cli.add_command(module.cli)

cli.add_command(rna.convert.cli_ct2db)
cli.add_command(rna.convert.cli_db2ct)
cli.add_command(fold.cli_datapath)
cli.add_command(cli_docs)
cli.add_command(cli_github)
cli.add_command(cli_pypi)
cli.add_command(cli_conda)
cli.add_command(cli_biorxiv)
