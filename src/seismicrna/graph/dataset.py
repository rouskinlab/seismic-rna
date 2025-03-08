from abc import ABC, abstractmethod
from functools import cached_property
from itertools import chain
from pathlib import Path
from typing import Iterable

from .base import (BaseWriter,
                   get_action_name,
                   make_path_subject,
                   make_title_action_sample)
from .cgroup import ClusterGroupRunner, cgroup_table
from .onesource import OneSourceClusterGroupGraph
from .rel import OneRelGraph, RelRunner
from ..core.arg import opt_verify_times
from ..core.dataset import MutsDataset
from ..core.table import get_subpattern
from ..core.task import dispatch
from ..table import load_all_datasets


class DatasetGraph(OneRelGraph, OneSourceClusterGroupGraph, ABC):
    """ Graph based on one Dataset. """

    def __init__(self, *, dataset: MutsDataset, num_cpus: int, **kwargs):
        super().__init__(**kwargs)
        self.dataset = dataset
        self.num_cpus = num_cpus

    @property
    def top(self):
        return self.dataset.top

    @property
    def branches(self):
        return self.dataset.branches

    @property
    def sample(self):
        return self.dataset.sample

    @property
    def ref(self):
        return self.dataset.ref

    @property
    def reg(self):
        return self.dataset.region.name

    @property
    def seq(self):
        return self.dataset.region.seq

    @cached_property
    def action(self):
        return get_action_name(self.dataset)

    @cached_property
    def path_subject(self):
        return make_path_subject(self.action, self.k, self.clust)

    @cached_property
    def title_action_sample(self):
        return make_title_action_sample(self.action, self.sample)

    @cached_property
    def _title_main(self):
        return [f"{self.what()} "
                f"{self.relationships} bases "
                f"in {self.title_action_sample} "
                f"over reference {repr(self.ref)} "
                f"region {repr(self.reg)}"]

    @property
    def details(self):
        return []

    @cached_property
    def predicate(self):
        return super().predicate + [self.codestring]

    @cached_property
    def pattern(self):
        """ Relationship pattern for the graph. """
        return get_subpattern(self.rel_name, self.dataset.pattern)


class DatasetWriter(BaseWriter, ABC):

    def __init__(self, *, dataset: MutsDataset, **kwargs):
        super().__init__(**kwargs)
        self.dataset = dataset

    @abstractmethod
    def get_graph(self, *args, **kwargs) -> DatasetGraph:
        """ Return a graph instance. """

    def iter_graphs(self, *, rels: list[str], cgroup: str, **kwargs):
        for cparams in cgroup_table(self.dataset, cgroup):
            for rel in rels:
                yield self.get_graph(rel, **kwargs | cparams)


class DatasetRunner(RelRunner, ClusterGroupRunner, ABC):

    @classmethod
    @abstractmethod
    def get_writer_type(cls) -> type[DatasetWriter]:
        pass

    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_verify_times]

    @classmethod
    def get_input_loader(cls):
        return load_all_datasets

    @classmethod
    def run(cls,
            input_path: Iterable[str | Path], *,
            verify_times: bool,
            num_cpus: int,
            **kwargs):
        # Generate a table writer for each table.
        writer_type = cls.get_writer_type()
        writers = [writer_type(dataset=dataset)
                   for dataset
                   in cls.load_input_files(input_path,
                                           verify_times=verify_times)]
        return list(chain(*dispatch([writer.write for writer in writers],
                                    num_cpus=num_cpus,
                                    pass_num_cpus=True,
                                    as_list=False,
                                    ordered=False,
                                    raise_on_error=False,
                                    kwargs=kwargs)))
