from __future__ import annotations
from functools import cached_property
from pathlib import Path

from .io import DuplexRefseqIO
from .report import DuplexReport
from ..core import path
from ..core.dataset import LoadFunction, RegionDataset
from ..core.header import NO_K, NO_KS
from ..core.report import DuplexCutF, DuplexKsF, DuplexSourcesF, RefseqChecksumF, RegF
from ..core.seq.region import Region

# Token joining the identifiers of the two combined sources.
FUSE_SEP = "__"

# Number of unpaired nucleotides inserted between the two strands so the
# duplex is a single drawable/graphable molecule; matches RNAstructure's
# intermolecular linker length.
LINKER_LENGTH = 3
LINKER_BASE = "N"
LINKER = LINKER_BASE * LINKER_LENGTH


def fuse_names(name_a: str, name_b: str, collapse_identical: bool = True):
    """Join two identifiers. An empty name is dropped (returns the other).
    Two identical names collapse to one when ``collapse_identical`` is
    True (e.g. the shared sample of a homodimer); pass False to keep them
    joined (e.g. a homodimer reference ``X__X``, so its fold is not
    mistaken for a monomer fold of ``X``)."""
    if not name_a:
        return name_b
    if not name_b:
        return name_a
    if collapse_identical and name_a == name_b:
        return name_a
    return f"{name_a}{FUSE_SEP}{name_b}"


def load_source_table(table_file: str | Path):
    """Load a source position table (filter or cluster) from its CSV."""
    from ..cluster.data import ClusterPositionTableLoader
    from ..filter.table import FilterPositionTableLoader

    errors = dict()
    for loader in [FilterPositionTableLoader, ClusterPositionTableLoader]:
        try:
            return loader(table_file)
        except Exception as error:
            errors[loader.__name__] = error
    raise RuntimeError(
        f"Could not load source table {table_file}:\n"
        + "\n".join(f"  {name}: {error}" for name, error in errors.items())
    )


class DuplexDataset(RegionDataset):
    """Dataset made by combining two source datasets from different
    references into one duplex, without duplicating their read-level
    data (the sources are reached through table1/table2)."""

    @classmethod
    def get_report_type(cls):
        return DuplexReport

    @cached_property
    def refseq(self):
        return DuplexRefseqIO.load(
            DuplexRefseqIO.build_path(
                {
                    path.TOP: self.top,
                    path.SAMPLE: self.sample,
                    path.BRANCHES: self.branches,
                    path.REF: self.ref,
                    path.REG: self.report.get_field(RegF),
                }
            ),
            checksum=self.report.get_field(RefseqChecksumF),
        ).refseq

    @cached_property
    def region(self):
        # A duplex spans its whole fused reference, so build a full region
        # (named "full") from the fused reference sequence.
        return Region(self.ref, self.refseq)

    @property
    def cut(self) -> int:
        """Length of the 5' strand (position of the strand break)."""
        return self.report.get_field(DuplexCutF)

    @property
    def duplex_ks(self) -> list[int]:
        """Numbers of clusters in the duplex (empty if unclustered)."""
        return [int(k) for k in self.report.get_field(DuplexKsF)]

    @property
    def header_depth(self) -> int:
        """Number of column-header rows in the duplex CSV: 3 when the duplex is
        clustered (relationship, K, cluster), else 1 (relationship)."""
        return 3 if self.duplex_ks else 1

    @property
    def probe(self) -> str:
        # The chemical probe is set upstream (at the filter step); like
        # cluster, delegate to the 5' strand's source rather than storing
        # a copy in the duplex report.
        return self.dataset1.probe

    @cached_property
    def source_tables(self):
        """The two source position tables (loaded lazily); an entry is
        None for a data-less strand (e.g. a bare sequence)."""
        return [
            load_source_table(self.top.joinpath(rel)) if rel else None
            for rel in self.report.get_field(DuplexSourcesF)
        ]

    @property
    def table1(self):
        """Position table of the 5' strand's source (always present)."""
        return self.source_tables[0]

    @property
    def table2(self):
        """Position table of the 3' strand's source, or None if it has
        no data (a bare sequence)."""
        return self.source_tables[1]

    @property
    def dataset1(self):
        return self.table1._dataset

    @property
    def dataset2(self):
        return self.table2._dataset if self.table2 is not None else None

    @property
    def pattern(self):
        return self.dataset1.pattern

    @cached_property
    def ks(self):
        return self.duplex_ks or NO_KS

    @cached_property
    def best_k(self):
        ks = self.duplex_ks
        return ks[-1] if ks else NO_K

    @property
    def data_dirs(self):
        return [self.dir]

    @property
    def num_batches(self):
        # A duplex has no read-level batches of its own; its read-level data
        # live in the source datasets (self.dataset1/self.dataset2).
        return 0

    def get_batch(self, batch_num: int):
        raise NotImplementedError(
            "A duplex has no read-level batches; use table1/table2 to reach "
            "the source datasets"
        )


load_duplex_dataset = LoadFunction(DuplexDataset)
