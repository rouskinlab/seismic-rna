import webbrowser

from click import command

from .core.logs import logger


def open_url(url: str):
    try:
        webbrowser.open(url)
    except Exception as error:
        logger.error(error)


@command("docs")
def cli_docs():
    """ Open the documentation for SEISMIC-RNA. """
    open_url("https://rouskinlab.github.io/seismic-rna/index.html")


@command("github")
def cli_github():
    """ Open SEISMIC-RNA's source code on GitHub. """
    open_url("https://github.com/rouskinlab/seismic-rna")


@command("pypi")
def cli_pypi():
    """ Open SEISMIC-RNA on the Python Packaging Index. """
    open_url("https://pypi.org/project/seismic-rna")


@command("conda")
def cli_conda():
    """ Open SEISMIC-RNA on Anaconda's website. """
    open_url("https://anaconda.org/bioconda/seismic-rna")


@command("biorxiv")
def cli_biorxiv():
    """ Open the preprint for SEISMIC-RNA in bioRxiv. """
    open_url("https://www.biorxiv.org/content/10.1101/2024.04.29.591762v2")

########################################################################
#                                                                      #
# © Copyright 2022-2025, the Rouskin Lab.                              #
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
