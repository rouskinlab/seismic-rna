from abc import ABC, abstractmethod
from functools import cached_property

from .base import BaseGraph, BaseRunner
from ..core.arg.cli import NO_GROUP, GROUP_BY_K, GROUP_ALL, opt_cgroup
from ..core.dataset import Dataset
from ..core.header import NO_KS, NO_CLUSTS, format_clust_names, list_ks_clusts
from ..core.table.base import Table


def get_ks(source: Dataset | Table):
    """List the numbers of clusters for a source of data."""
    if isinstance(source, Dataset):
        return source.ks
    if isinstance(source, Table):
        return source.header.ks
    raise TypeError(source)


def get_ks_clusts(source: Dataset | Table):
    """List the clusters for a source of data."""
    ks = get_ks(source)
    if ks == NO_KS:
        return NO_CLUSTS
    return list_ks_clusts(ks)


def make_tracks(
    source: Dataset | Table,
    k: int | None,
    clust: int | None,
    *,
    k_clust_list: list[tuple[int, int]] | None = None,
):
    """Make an index for the rows or columns of a graph."""
    clusts = get_ks_clusts(source)

    if k is None and clust is None:
        tracks = clusts
    else:
        tracks = [
            (k_, clust_)
            for k_, clust_ in clusts
            if ((k is None or k_ == k) and (clust is None or clust_ == clust))
        ]
    if k_clust_list:
        assert isinstance(k_clust_list, list), "k_clust_list must be a list of tuples."
        selected = set(tracks)
        combo_selected = set(k_clust_list)
        if tracks == clusts:
            selected &= combo_selected
        else:
            selected |= combo_selected
        tracks = sorted(list(selected))
    return tracks


def _track_count(tracks: list[tuple[int, int]] | None):
    return len(tracks) if tracks is not None else 1


def get_clust_proportions(top, sample, branches, ref, reg):
    """Load the proportion of each cluster, keyed by (k, cluster), from
    the cluster abundance table identified by these path fields. Returns
    None if no such abundance table exists (e.g. the source is not
    clustered)."""
    from ..core import path
    from ..cluster.data import ClusterAbundanceTableLoader

    path_fields = {
        path.TOP: top,
        path.SAMPLE: sample,
        path.BRANCHES: branches,
        path.REF: ref,
        path.REG: reg,
    }
    abundance_file = ClusterAbundanceTableLoader.build_path(path_fields)
    if not abundance_file.is_file():
        return None
    return ClusterAbundanceTableLoader(abundance_file, verify_times=False).proportions


def _format_proportion(proportion: float):
    """Format a cluster proportion as a percentage."""
    return f"{round(proportion * 100.0, 1)}%"


def _track_titles(tracks: list[tuple[int, int]] | None, proportions=None):
    """Format the title of each track, optionally annotating each
    cluster with its proportion (as a percentage)."""
    if tracks is None:
        return None
    names = format_clust_names(tracks, allow_duplicates=False)
    if proportions is None:
        return names
    titles = list()
    for name, k_clust in zip(names, tracks, strict=True):
        proportion = proportions.get(k_clust)
        titles.append(
            f"{name} ({_format_proportion(proportion)})"
            if proportion is not None
            else name
        )
    return titles


def cgroup_table(source: Dataset | Table, cgroup: str):
    if cgroup == NO_GROUP:
        # One file per cluster, with no subplots.
        return [dict(k=k, clust=clust) for k, clust in get_ks_clusts(source)]
    elif cgroup == GROUP_BY_K:
        # One file per k, with one subplot per cluster.
        return [dict(k=k, clust=None) for k in sorted(get_ks(source))]
    elif cgroup == GROUP_ALL:
        # One file, with one subplot per cluster.
        return [dict(k=None, clust=None)]
    raise ValueError(f"Invalid value for cgroup: {repr(cgroup)}")


class ClusterGroupGraph(BaseGraph, ABC):
    """Graph in which clusters can be placed in subplots."""

    @property
    @abstractmethod
    def row_tracks(self) -> list[tuple[int, int]] | None:
        """Track for each row of subplots."""

    @property
    @abstractmethod
    def col_tracks(self) -> list[tuple[int, int]] | None:
        """Track for each column of subplots."""

    @property
    def nrows(self):
        """Number of rows of subplots."""
        return _track_count(self.row_tracks)

    @property
    def ncols(self):
        """Number of columns of subplots."""
        return _track_count(self.col_tracks)

    @property
    def clust_proportions(self):
        """Mapping from each (k, cluster) to its proportion, or None if
        cluster proportions are not available for this graph."""
        return None

    @property
    def row_proportions(self):
        """Cluster proportions for the row tracks."""
        return self.clust_proportions

    @property
    def col_proportions(self):
        """Cluster proportions for the column tracks."""
        return self.clust_proportions

    @cached_property
    def row_titles(self):
        """Titles of the rows."""
        return _track_titles(self.row_tracks, self.row_proportions)

    @cached_property
    def col_titles(self):
        """Titles of the columns."""
        return _track_titles(self.col_tracks, self.col_proportions)

    @cached_property
    def _subplots_params(self):
        return super()._subplots_params | dict(
            rows=self.nrows,
            cols=self.ncols,
            row_titles=self.row_titles,
            column_titles=self.col_titles,
            shared_xaxes="all",
            shared_yaxes="all",
        )


class ClusterGroupRunner(BaseRunner, ABC):
    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_cgroup]
