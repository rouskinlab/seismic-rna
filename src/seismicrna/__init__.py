import warnings
from importlib import metadata

from . import demult, align, relate, cluster, table, draw
from .main import run, main_cli

warnings.simplefilter(action='ignore', category=FutureWarning)


try:
    __version__ = metadata.version(__package__ or __name__)
except metadata.PackageNotFoundError:
    __version__ = None
