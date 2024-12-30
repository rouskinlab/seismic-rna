import webbrowser

from click import command


@command("docs")
def cli_docs():
    """ Open the documentation for SEISMIC-RNA. """
    webbrowser.open("https://rouskinlab.github.io/seismic-rna/index.html")


@command("github")
def cli_github():
    """ Open SEISMIC-RNA's source code on GitHub. """
    webbrowser.open("https://github.com/rouskinlab/seismic-rna")


@command("pypi")
def cli_pypi():
    """ Open SEISMIC-RNA on the Python Packaging Index. """
    webbrowser.open("https://pypi.org/project/seismic-rna")


@command("conda")
def cli_conda():
    """ Open SEISMIC-RNA on Anaconda's website. """
    webbrowser.open("https://anaconda.org/bioconda/seismic-rna")


@command("biorxiv")
def cli_biorxiv():
    """ Open the preprint for SEISMIC-RNA in bioRxiv. """
    webbrowser.open(
        "https://www.biorxiv.org/content/10.1101/2024.04.29.591762v2"
    )

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
