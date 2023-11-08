from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain, product
from logging import getLogger
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .base import PRECISION, CartesianGraph, TwoTableSeqGraph, make_subject
from .color import SeqColorMap
from .seq import get_table_params
from .write import TwoTableGraphWriter
from ..core.arg import (docdef,
                        arg_input_path,
                        opt_rels,
                        opt_y_ratio,
                        opt_quantile,
                        opt_arrange,
                        opt_csv,
                        opt_html,
                        opt_pdf,
                        opt_force,
                        opt_max_procs,
                        opt_parallel)
from ..core.header import make_header
from ..core.parallel import dispatch
from ..core.seq import POS_NAME
from ..table.base import get_rel_name
from ..table.load import (PosTableLoader,
                          RelPosTableLoader,
                          MaskPosTableLoader,
                          ClustPosTableLoader,
                          find_tables)

logger = getLogger(__name__)

# Sample index level name
SAMPLE_NAME = "Sample"
ROW_NAME = "Row"
COL_NAME = "Column"


def _get_clusts(max_order: int,
                min_order: int,
                order: int | None,
                clust: int | None):
    return (make_header(max_order=max_order,
                        min_order=min_order).select(order=order,
                                                    clust=clust)
            if max_order > 0
            else None)


# Base Sequence Pair Graph #############################################

class SeqPairGraph(CartesianGraph, TwoTableSeqGraph, ABC):
    """ Graph that compares a pair of samples of the same sequence. """

    @classmethod
    def get_cmap_type(cls):
        return SeqColorMap

    @classmethod
    def sources(cls) -> dict[type[PosTableLoader], str]:
        """ Names of the sources of data. """
        return {RelPosTableLoader: "Related",
                MaskPosTableLoader: "Masked",
                ClustPosTableLoader: "Clustered"}

    def __init__(self,
                 *args,
                 rel: str,
                 y_ratio: bool,
                 quantile: float,
                 order1: int | None = None,
                 clust1: int | None = None,
                 order2: int | None = None,
                 clust2: int | None = None,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.rel_code = rel
        self.y_ratio = y_ratio
        self.quantile = quantile
        self._order1 = order1
        self._clust1 = clust1
        self._order2 = order2
        self._clust2 = clust2

    @property
    def rel(self):
        """ Relationship to graph. """
        return get_rel_name(self.rel_code)

    def _get_source(self, table: PosTableLoader):
        table_type = type(table)
        try:
            return self.sources()[table_type]
        except KeyError:
            raise TypeError(f"Invalid table for {self}: {table_type.__name__}")

    @property
    def source1(self):
        """ Source of the first data set. """
        return self._get_source(self.table1)

    @property
    def source2(self):
        """ Source of the second data set. """
        return self._get_source(self.table2)

    @property
    def title(self):
        return (f"{self.quantity} of {self.rel} bases "
                f"in {self.source1} reads from {self.sample1} "
                f"vs. {self.source2} reads from {self.sample2} "
                f"per position in {self.ref}:{self.sect}")

    @property
    def subject1(self):
        return (make_subject(self.source1, self._order1, self._clust1)
                if isinstance(self.table1, ClustPosTableLoader)
                else self.source1)

    @property
    def subject2(self):
        return (make_subject(self.source2, self._order2, self._clust2)
                if isinstance(self.table2, ClustPosTableLoader)
                else self.source2)

    @property
    def subject(self):
        return f"{self.subject1}__and__{self.subject2}"

    @property
    def quantity(self):
        return "Ratio" if self.y_ratio else "Count"

    @property
    def predicate(self):
        return f"{self.rel_code}-{self.quantity}"

    @classmethod
    @abstractmethod
    def graph_type(cls):
        """ Type of the graph. """

    @property
    def graph_filename(self):
        return f"{self.subject}_{self.predicate}_{self.graph_type()}".lower()

    def _fetch_data(self,
                    table: PosTableLoader,
                    order: int | None,
                    clust: int | None):
        # Determine which relationships, orders, and clusters to use.
        all_kwargs = dict(rel=self.rel, order=order, clust=clust)
        return (table.fetch_ratio(quantile=self.quantile,
                                  precision=PRECISION,
                                  **all_kwargs)
                if self.y_ratio
                else table.fetch_count(**all_kwargs))

    @cached_property
    def data1(self):
        return self._fetch_data(self.table1, self._order1, self._clust1)

    @cached_property
    def data2(self):
        return self._fetch_data(self.table2, self._order2, self._clust2)

    @cached_property
    def row_index(self):
        return _get_clusts(self.table2.header.max_order,
                           self.table2.header.min_order,
                           self._order2,
                           self._clust2)

    @cached_property
    def col_index(self):
        return _get_clusts(self.table1.header.max_order,
                           self.table1.header.min_order,
                           self._order1,
                           self._clust1)


class SeqPairOneAxisGraph(SeqPairGraph, ABC):
    """ Graph of a pair of datasets over the same sequence in which the
    data series are merged in some fashion into another series, and the
    original data are not graphed directly. """

    @property
    def x_title(self) -> str:
        return POS_NAME

    @classmethod
    @abstractmethod
    def _trace_function(cls) -> Callable:
        """ Function to generate the graph's traces. """

    @property
    @abstractmethod
    def _merge_data(self) -> Callable:
        """ Function to merge the two datasets into one. """

    @cached_property
    def data(self):
        # Merge each pair in the Cartesian product of datasets 1 and 2.
        data = pd.DataFrame.from_dict(
            {(row, col): self._merge_data(vals1, vals2)
             for (col, (key1, vals1)), (row, (key2, vals2))
             in product(enumerate(self.data1.items(), start=1),
                        enumerate(self.data2.items(), start=1))}
        )
        # Indicate that the column index levels now correspond to the
        # rows and columns of the graph.
        data.columns.rename([ROW_NAME, COL_NAME], inplace=True)
        return data

    def get_traces(self):
        for (row, col), series in self.data.items():
            for trace in self._trace_function()(series, self.cmap):
                yield (row, col), trace


class SeqPairTwoAxisGraph(SeqPairGraph, ABC):
    """ Graph of a pair of datasets over the same sequence in which one
    dataset is graphed on each axis. """

    @property
    def x_title(self):
        return self.sample1

    @property
    def y_title(self):
        return self.sample2

    @cached_property
    def data(self):
        # Join data tables 1 and 2 horizontally.
        data = pd.concat([self.data1, self.data2], axis=1, join="inner")
        # Add the sample names as the first level of the columns.
        samples = np.hstack([np.repeat(self.sample1, self.data1.columns.size),
                             np.repeat(self.sample2, self.data2.columns.size)])
        names = [SAMPLE_NAME] + list(data.columns.names)
        data.columns = pd.MultiIndex.from_arrays(
            [(samples if name == SAMPLE_NAME
              else data.columns.get_level_values(name).values)
             for name in names],
            names=names
        )
        return data


# Sequence Pair Graph Writer ###########################################

class SeqPairGraphWriter(TwoTableGraphWriter, ABC):

    @property
    @abstractmethod
    def graph_type(self) -> type[SeqPairGraph]:
        """ Type of the graph to iterate. """

    def iter(self, rels_sets: tuple[str, ...],
             arrange: str, y_ratio: bool, quantile: float):
        _, _, csparams1 = get_table_params(self.table1, arrange)
        _, _, csparams2 = get_table_params(self.table2, arrange)
        for cparams1, cparams2 in product(csparams1, csparams2):
            for rels in rels_sets:
                yield self.graph_type(table1=self.table1,
                                      table2=self.table2,
                                      rel=rels,
                                      y_ratio=y_ratio,
                                      quantile=quantile,
                                      order1=cparams1.get("order"),
                                      clust1=cparams1.get("clust"),
                                      order2=cparams2.get("order"),
                                      clust2=cparams2.get("clust"))


# Helper functions #####################################################


class SeqPairGraphRunner(object):
    params = [
        arg_input_path,
        opt_rels,
        opt_y_ratio,
        opt_quantile,
        opt_arrange,
        opt_csv,
        opt_html,
        opt_pdf,
        opt_force,
        opt_max_procs,
        opt_parallel,
    ]

    @classmethod
    @abstractmethod
    def writer_type(cls) -> type[SeqPairGraphWriter]:
        """ Type of SeqPairGraphWriter. """

    @classmethod
    @docdef.auto()
    def run(cls,
            input_path: tuple[str, ...],
            rels: tuple[str, ...], *,
            y_ratio: bool,
            quantile: float,
            arrange: str,
            csv: bool,
            html: bool,
            pdf: bool,
            force: bool,
            max_procs: int,
            parallel: bool) -> list[Path]:
        """ Run the graph seqpair module. """
        tables = list(find_tables(input_path))
        if len(tables) % 2 != 0:
            raise ValueError(f"Number of files must be even, but got {tables}")
        writers = [cls.writer_type()(table1_file=t1, table2_file=t2)
                   for t1, t2 in zip(tables[0::2], tables[1::2], strict=True)]
        return list(chain(*dispatch([writer.write for writer in writers],
                                    max_procs,
                                    parallel,
                                    pass_n_procs=False,
                                    kwargs=dict(rels_sets=rels,
                                                y_ratio=y_ratio,
                                                quantile=quantile,
                                                arrange=arrange,
                                                csv=csv,
                                                html=html,
                                                pdf=pdf,
                                                force=force))))

########################################################################
#                                                                      #
# Copyright ©2023, the Rouskin Lab.                                    #
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
