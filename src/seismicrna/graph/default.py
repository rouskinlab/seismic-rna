import os
from pathlib import Path

from click import command

from .base import find_tables, GraphWriter
from .seq import run as run_pos
from ..core import docdef, path
from ..core.cli import (opt_table, opt_csv, opt_html, opt_pdf,
                        opt_max_procs, opt_parallel)


params = [
    opt_table,
    opt_csv,
    opt_html,
    opt_pdf,
    opt_max_procs,
    opt_parallel,
]


@command(__name__.split(os.path.extsep)[-1], params=params)
def cli(*args, **kwargs):
    """ Create the canonical set of graphs. """
    return run(*args, **kwargs)


@docdef.auto()
def run(table: tuple[str, ...],
        csv: bool,
        html: bool,
        pdf: bool,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Create the canonical set of graphs. """
    # Positional graphs.
    run_pos()


class CanonicalGraphWriter(GraphWriter):

    def iter(self, fields: str, count: bool, stack: bool):
        if isinstance(self.table, RelPosTableLoader):
            if stack:
                if count:
                    yield RelCountStackedPosBarGraph(table=self.table,
                                                     codes=fields)
                else:
                    yield RelFracStackedPosBarGraph(table=self.table,
                                                    codes=fields)
            else:
                for field in fields:
                    if count:
                        yield RelCountSerialPosBarGraph(table=self.table,
                                                        codes=field)
                    else:
                        yield RelFracSerialPosBarGraph(table=self.table,
                                                       codes=field)
        elif isinstance(self.table, MaskPosTableLoader):
            if stack:
                if count:
                    yield MaskCountStackedPosBarGraph(table=self.table,
                                                      codes=fields)
                else:
                    yield MaskFracStackedPosBarGraph(table=self.table,
                                                     codes=fields)
            else:
                for field in fields:
                    if count:
                        yield MaskCountSerialPosBarGraph(table=self.table,
                                                         codes=field)
                    else:
                        yield MaskFracSerialPosBarGraph(table=self.table,
                                                        codes=field)
        elif isinstance(self.table, ClusterPosTableLoader):
            for cluster in self.table.clusters:
                yield ClustPosBarGraph(table=self.table, cluster=cluster)
