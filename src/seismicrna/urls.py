import webbrowser

from click import command

from .core.logs import logger


def open_url(url: str):
    try:
        webbrowser.open(url)
    except Exception as error:
        logger.error(error)


@command("biorxiv")
def cli_biorxiv():
    """ Open the preprint for SEISMIC-RNA in bioRxiv. """
    open_url("https://www.biorxiv.org/content/10.1101/2024.04.29.591762v2")


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
