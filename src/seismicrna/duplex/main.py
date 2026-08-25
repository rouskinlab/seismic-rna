from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from itertools import combinations, product
from pathlib import Path
from typing import Iterable

from click import command

from .dataset import LINKER, LINKER_LENGTH, fuse_names
from .io import DuplexRefseqIO
from .report import DuplexReport
from .table import DuplexPositionTable, DuplexPositionTableLoader
from ..core import path
from ..core.arg.cmd import CMD_DUPLEX
from ..core.arg.cli import (
    arg_input_path,
    opt_branch,
    opt_duplex_pair,
    opt_dimer,
    opt_duplex_file,
    opt_duplex_sequence,
    opt_duplex_table_region,
    opt_duplex_full_ref,
    opt_duplex_region_ref,
    opt_verify_times,
    opt_num_cpus,
    opt_force,
)
from ..core.error import IncompatibleOptionsError, IncompatibleValuesError
from ..core.header import CLUST_NAME, NUM_CLUSTS_NAME, REL_NAME
from ..core.logs import logger
from ..core.run import run_func
from ..core.task import dispatch
from ..core.seq.fasta import parse_fasta
from ..core.seq.region import FULL_NAME, SEQ_INDEX_NAMES
from ..core.seq.xna import DNA
from ..core.table.write import PRECISION
from ..core.write import need_write


@dataclass
class Strand:
    """One strand of a duplex to be combined."""

    seq: DNA
    ref: str
    reg: str
    sample: str
    branches: dict
    # Source position table if the strand has data; None for a data-less
    # partner (e.g. a bare sequence). Its source table path, relative to
    # the duplex top, or "" if it has no source table.
    table: object | None
    source: str
    # Whether the strand spans only the source table's region (True) or
    # its full reference sequence with the region's data mapped onto it
    # (False, the default); irrelevant for a data-less strand.
    table_region: bool = False

    @classmethod
    def from_table(cls, table, top: Path, table_region: bool = False):
        return cls(
            seq=table.region.seq if table_region else table.refseq,
            ref=table.ref,
            reg=table.reg,
            sample=table.sample,
            branches=table.branches,
            table=table,
            source=str(table.path.relative_to(top)),
            table_region=table_region,
        )

    @classmethod
    def from_sequence(cls, seq: DNA, name: str):
        return cls(
            seq=seq,
            ref=name,
            reg=name,
            sample="",
            branches=dict(),
            table=None,
            source="",
        )


def strand_clusters(strand: Strand):
    """List the per-cluster data of one strand as (label, DataFrame)
    pairs, one per cluster at the table's best number of clusters.

    Each DataFrame is indexed by position (spanning either the table's
    region or the full reference, per ``strand.table_region``) and has
    one column per relationship. An unclustered (average) table yields a
    single ``(None, data)`` pair, and a data-less strand yields a single
    ``(None, None)`` placeholder (its positions get no data).
    """
    table = strand.table
    if table is None:
        return [(None, None)]
    header = table.header
    if not header.get_is_clustered():
        clusters = [(None, table.data)]
    else:
        # Prefer the number of clusters that the cluster step chose as
        # best. That is 0 if no number of clusters passed its filters, in
        # which case it wrote a fallback (the ensemble average) and left
        # best_k at 0 to record that none was actually best; use the one
        # number of clusters that such a table contains.
        best_k = table._dataset.best_k
        if best_k < 1:
            ks = list(header.ks)
            if len(ks) != 1:
                raise IncompatibleValuesError(
                    f"{table} has no best number of clusters, and "
                    f"{len(ks)} to choose from ({ks}), so which to "
                    f"duplex is ambiguous"
                )
            best_k = ks[0]
            logger.warning(
                "{} has no best number of clusters, since none passed the "
                "filters of the cluster step; using its only number of "
                "clusters (K = {})",
                table,
                best_k,
            )
        clusters = list()
        for clust in range(1, best_k + 1):
            cols = header.select(k=best_k, clust=clust)
            sub = table.data.loc[:, cols]
            # Flatten the (relationship, K, cluster) columns to relationship.
            sub.columns = sub.columns.get_level_values(REL_NAME)
            clusters.append((clust, sub))
    if not strand.table_region:
        # Map the region's data onto the full reference sequence, leaving
        # positions outside the region empty (so cofolding uses the whole
        # reference while only the table's region carries reactivities).
        import pandas as pd

        full = table.refseq
        target = pd.MultiIndex.from_arrays(
            [list(range(1, len(full) + 1)), list(str(full))], names=SEQ_INDEX_NAMES
        )
        clusters = [(label, sub.reindex(index=target)) for label, sub in clusters]
    return clusters


def make_duplex(
    strand1: Strand, strand2: Strand, top: Path, *, branch: str, force: bool
):
    """Fuse two strands into one duplex position table.

    The 5' strand (strand1) must carry data (its columns define the
    table); the 3' strand may be data-less (e.g. a bare sequence). When
    either strand is clustered, the duplex is the cross-product of the two
    strands' clusters (at each strand's best K): one fused profile per
    combination of a 5'-strand cluster and a 3'-strand cluster.
    """
    if strand1.table is None:
        raise ValueError("The 5' strand of a duplex must have data")
    import numpy as np
    import pandas as pd

    # Propagate the parent branches. Like pool/join, this multi-parent
    # step requires its data-bearing sources to share one branch chain
    # (the linear branch model has no room for two different branches at
    # the same step); a data-less partner contributes no branches.
    parent_branches = strand1.branches
    if strand2.table is not None and strand2.branches != parent_branches:
        raise IncompatibleOptionsError(
            "Cannot duplex references from different branches: "
            f"{parent_branches} (5') vs {strand2.branches} (3')"
        )
    # The duplex is folded with a single energy method, so both
    # data-bearing strands must use the same chemical probe.
    if strand2.table is not None:
        probe1 = strand1.table._dataset.probe
        probe2 = strand2.table._dataset.probe
        if probe1 != probe2:
            raise IncompatibleOptionsError(
                "Cannot duplex references probed with different chemicals: "
                f"{probe1} (5') vs {probe2} (3')"
            )
    fused_seq = strand1.seq + DNA(LINKER) + strand2.seq
    cut = len(strand1.seq)
    length = len(fused_seq)
    fields = {
        path.TOP: top,
        path.SAMPLE: fuse_names(strand1.sample, strand2.sample),
        path.BRANCHES: path.add_branch(path.DUPLEX_STEP, branch, parent_branches),
        # Keep identical reference names joined (a homodimer is "X__X",
        # not "X") so its duplex/fold output is never confused with a
        # monomer's.
        path.REF: fuse_names(strand1.ref, strand2.ref, collapse_identical=False),
        # A duplex spans its whole fused reference, so its region is
        # "full" (whether each strand came from a full reference or a
        # sub-region), keeping fold/graph/draw region names consistent.
        path.REG: FULL_NAME,
    }
    report_file = DuplexReport.build_path(fields)
    if need_write(report_file, force):
        began = datetime.now()
        index = pd.MultiIndex.from_arrays(
            [list(range(1, length + 1)), list(str(fused_seq))], names=SEQ_INDEX_NAMES
        )
        # Cross the two strands' clusters: each combination is one column
        # block whose 5' rows come from a strand-1 cluster and whose 3'
        # rows come from a strand-2 cluster (empty for a data-less strand).
        clusters1 = strand_clusters(strand1)
        clusters2 = strand_clusters(strand2)
        rels = list(clusters1[0][1].columns)
        combos = [(d1, d2) for (_, d1), (_, d2) in product(clusters1, clusters2)]
        k_duplex = len(combos)
        # Fill a (position, relationship x cluster) array in the column
        # order that a relationship/cluster header expects: relationship
        # varies slowest, cluster fastest.
        values = np.full((length, len(rels) * k_duplex), np.nan)
        for n, (d1, d2) in enumerate(combos):
            top5 = d1.reindex(columns=rels).values
            bot3 = d2.reindex(columns=rels).values if d2 is not None else None
            for c, _ in enumerate(rels):
                col = c * k_duplex + n
                values[:cut, col] = top5[:, c]
                if bot3 is not None:
                    values[cut + LINKER_LENGTH :, col] = bot3[:, c]
        if k_duplex > 1:
            # A cluster product: label the combinations 1..k_duplex under a
            # single number-of-clusters K = k_duplex.
            columns = pd.MultiIndex.from_tuples(
                [(rel, k_duplex, n + 1) for rel in rels for n in range(k_duplex)],
                names=[REL_NAME, NUM_CLUSTS_NAME, CLUST_NAME],
            )
            duplex_ks = [k_duplex]
        else:
            columns = pd.Index(rels, name=REL_NAME)
            duplex_ks = list()
        data = pd.DataFrame(values, index=index, columns=columns)
        csv_file = DuplexPositionTableLoader.build_path(fields)
        csv_file.parent.mkdir(parents=True, exist_ok=True)
        data.round(PRECISION).to_csv(csv_file)
        # Store the fused reference sequence like any reference sequence,
        # in the same (region) directory as the table and report.
        _, refseq_checksum = DuplexRefseqIO(
            sample=fields[path.SAMPLE],
            branches=fields[path.BRANCHES],
            ref=fields[path.REF],
            reg=fields[path.REG],
            refseq=fused_seq,
        ).save(top, force=force)
        report = DuplexReport(
            sample=fields[path.SAMPLE],
            branches=fields[path.BRANCHES],
            ref=fields[path.REF],
            reg=fields[path.REG],
            end5=1,
            end3=length,
            refseq_checksum=refseq_checksum,
            duplex_cut=cut,
            duplex_ks=duplex_ks,
            duplex_sources=[strand1.source, strand2.source],
            began=began,
            ended=datetime.now(),
        )
        report.save(top, force=force)
    return report_file


def iter_sequence_strands(fasta: Path):
    """Yield a data-less Strand for every sequence in a FASTA file."""
    for ref, seq in parse_fasta(fasta, DNA):
        yield Strand.from_sequence(seq, ref)


def iter_raw_strands(sequences: Iterable[tuple[str, DNA]]):
    """Yield a data-less Strand for every (name, sequence) pair (each
    sequence given 5' to 3')."""
    for name, seq in sequences:
        yield Strand.from_sequence(DNA.from_any_seq(seq), name)


def strand_uses_table_region(
    ref: str, *, default: bool, region_refs: set[str], full_refs: set[str]
) -> bool:
    """Decide whether the strand for one reference spans only its table
    region (True) or its full reference (False), letting per-reference
    overrides take precedence over the global default."""
    if ref in full_refs and ref in region_refs:
        raise IncompatibleOptionsError(
            f"Reference {repr(ref)} was given to both --duplex-full-ref "
            "and --duplex-region-ref"
        )
    if ref in full_refs:
        return False
    if ref in region_refs:
        return True
    return default


def _table_order(table):
    return table.ref, table.reg, table.sample


def iter_duplex_pairs(tables: Iterable):
    """Yield every pair of tables whose branches match.

    Unlike a comparison graph, whose two tables must share a reference
    and region, a duplex is two different references, so its two strands
    are grouped by their branches instead: ``make_duplex`` refuses to
    duplex references from different branches, so pairing only within a
    branch avoids building pairs that would just be rejected.
    """
    tables = list(tables)
    # Group the tables by their branches. The key must be the branches
    # themselves, not their flattened form: flattening drops the step of
    # every branch and every empty branch name, so two tables from
    # different steps (e.g. a filter table and a duplex table, both on
    # default branches) would flatten alike but still be rejected.
    table_groups = defaultdict(list)
    for table in tables:
        key = tuple(sorted(table.branches.items()))
        if any(
            _table_order(other) == _table_order(table) for other in table_groups[key]
        ):
            logger.warning(
                "Duplicate reference, region, and sample: {}", _table_order(table)
            )
        else:
            table_groups[key].append(table)
    # Yield every pair of tables from each table group.
    for branches, table_group in table_groups.items():
        n_files = len(table_group)
        n_pairs = n_files * (n_files - 1) // 2
        logger.trace(
            "Found {} table(s) and {} pair(s) with branches {}",
            n_files,
            n_pairs,
            dict(branches),
        )
        # Sort the tables by reference to ensure the order of combinations
        # is consistent no matter the order of the "tables" argument.
        yield from combinations(sorted(table_group, key=_table_order), 2)


@run_func(CMD_DUPLEX)
def run(
    input_path: Iterable[str | Path],
    *,
    duplex_pair: bool,
    dimer: bool,
    duplex_file: str | None,
    duplex_sequence: Iterable[str],
    duplex_table_region: bool,
    duplex_full_ref: Iterable[str],
    duplex_region_ref: Iterable[str],
    branch: str,
    verify_times: bool,
    num_cpus: int,
    force: bool,
):
    """Combine two references into a duplex for cofolding.

    Each of these ways of choosing the 3' strand contributes its own
    duplexes, so any combination of them can be used at once:

    - ``--duplex-pair`` (the default): every pairwise combination of the
      input tables that share the same branches.
    - ``--dimer``: the input table itself (a homodimer).
    - ``--duplex-file`` / ``--duplex-sequence``: a partner without data
      (a FASTA file or a raw 5'-to-3' sequence).

    By default each strand spans its full reference (with the table's
    data mapped onto its region); ``--duplex-table-region`` restricts
    every strand to its table region, and ``--duplex-full-ref`` /
    ``--duplex-region-ref`` override that choice per reference.
    """
    from ..fold.main import load_foldable_tables

    # A duplex is made from single-molecule tables; drop any duplex table
    # among the inputs, so that re-running duplex over a directory that
    # already holds its own output does not fuse duplexes into chimeras.
    tables = [
        table
        for table in load_foldable_tables(input_path, verify_times=verify_times)
        if not isinstance(table, DuplexPositionTable)
    ]
    full_refs = set(duplex_full_ref)
    region_refs = set(duplex_region_ref)
    # Build each table's strand once, so that every duplex containing a
    # table applies the same choice of region to it.
    strands = {
        table: Strand.from_table(
            table,
            table.top,
            strand_uses_table_region(
                table.ref,
                default=duplex_table_region,
                region_refs=region_refs,
                full_refs=full_refs,
            ),
        )
        for table in tables
    }
    partners = list()
    if duplex_file:
        partners.extend(iter_sequence_strands(Path(duplex_file)))
    if duplex_sequence:
        partners.extend(iter_raw_strands(duplex_sequence))
    # Determine all pairs of strands to duplex. Each source of pairs adds
    # its own duplexes, so they can all be used at once.
    args = list()
    if duplex_pair:
        # Duplex every pair of two different tables.
        args.extend(
            (strands[table1], strands[table2], table1.top)
            for table1, table2 in iter_duplex_pairs(tables)
        )
    if dimer:
        # Duplex every table with itself.
        args.extend((strands[table], strands[table], table.top) for table in tables)
    # Duplex every table with every data-less partner.
    args.extend(
        (strands[table], partner, table.top) for table in tables for partner in partners
    )
    if not args:
        # Dispatching no tasks would warn that it was given none, so say
        # what would have made a duplex instead.
        logger.warning(
            "Nothing to duplex: give at least two tables, or use --dimer, "
            "--duplex-file, or --duplex-sequence for the second strand"
        )
        return list()
    return dispatch(
        make_duplex,
        num_cpus=num_cpus,
        pass_num_cpus=False,
        as_list=True,
        ordered=False,
        raise_on_error=False,
        label="making duplexes",
        unit="duplex",
        args=args,
        kwargs=dict(branch=branch, force=force),
    )


params = [
    arg_input_path,
    opt_duplex_pair,
    opt_dimer,
    opt_duplex_file,
    opt_duplex_sequence,
    opt_duplex_table_region,
    opt_duplex_full_ref,
    opt_duplex_region_ref,
    opt_branch,
    opt_verify_times,
    opt_num_cpus,
    opt_force,
]


@command(CMD_DUPLEX, params=params)
def cli(*args, **kwargs):
    """Combine two references into a duplex for cofolding."""
    return run(*args, **kwargs)
