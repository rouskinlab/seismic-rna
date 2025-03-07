from abc import ABC
from functools import cached_property
from pathlib import Path
from typing import Iterable

from .cgroup import cgroup_table
from .onetable import (OneTableRelClusterGroupGraph,
                       OneTableRelClusterGroupWriter,
                       OneTableRelClusterGroupRunner)
from .rel import OneRelGraph
from ..core import path
from ..core.arg import (opt_struct_file,
                        opt_branch,
                        opt_fold_regions_file,
                        opt_fold_coords,
                        opt_fold_primers,
                        opt_fold_full,
                        optional_path)
from ..core.logs import logger
from ..core.rna import RNAState, from_ct
from ..core.seq import DNA, RefRegions
from ..fold.main import load_foldable_tables


class StructOneTableGraph(OneTableRelClusterGroupGraph, OneRelGraph, ABC):
    """ Graph of data from one Table applied to RNA structure(s). """

    def __init__(self, *,
                 struct_file: Path | None,
                 struct_reg: str | None,
                 branch: str,
                 **kwargs):
        super().__init__(**kwargs)
        self._struct_file = struct_file
        self._struct_reg = struct_reg
        self.branch = branch

    @property
    def branches(self):
        return path.add_branch(path.FOLD_STEP, self.branch, super().branches)

    @property
    def struct_reg(self):
        """ Name of the region from which the structure comes. """
        if self._struct_file is not None:
            # Use the region from the given structure file.
            fields = path.parse(self._struct_file, path.CT_FILE_LAST_SEGS)
            return fields[path.REG]
        if self._struct_reg is None:
            raise ValueError("A structure region is required if no structure "
                             "file is given")
        return self._struct_reg

    def get_path_fields(self):
        fields = super().get_path_fields()
        # Replace the region with the structure region.
        fields[path.REG] = self.struct_reg
        return fields

    def _get_struct_file(self, profile: str):
        """ Get the path of the CT file of RNA structures for a specific
        RNA profile. """
        if self._struct_file is not None:
            # Use the given structure file for every profile.
            return self._struct_file
        # Determine the path of the structure file from the profile.
        return path.build(path.CT_FILE_ALL_SEGS,
                          {path.TOP: self.top,
                           path.SAMPLE: self.sample,
                           path.STEP: path.FOLD_STEP,
                           path.BRANCHES: self.branches,
                           path.REF: self.ref,
                           path.REG: self.struct_reg,
                           path.PROFILE: profile,
                           path.EXT: path.CT_EXT})

    @cached_property
    def path_subject(self):
        return f"{self.table.reg}__{super().path_subject}"

    @cached_property
    def _title_main(self):
        return [f"{self.what()} of structures "
                f"for reference {repr(self.ref)} "
                f"region {repr(self._struct_reg)} "
                f"vs. {self.data_kind}s "
                f"of {self.relationships} bases "
                f"in {self.title_action_sample} "
                f"over region {repr(self.reg)}"]

    def iter_profiles(self):
        """ Yield each RNAProfile from the table. """
        yield from self.table.iter_profiles(quantile=self.quantile,
                                            rel=self.rel_name,
                                            k=self.k,
                                            clust=self.clust)

    def iter_states(self):
        """ Yield each RNAState. """
        for profile in self.iter_profiles():
            ct_file = self._get_struct_file(profile.profile)
            try:
                for struct in from_ct(ct_file):
                    yield RNAState.from_struct_profile(struct, profile)
            except FileNotFoundError:
                logger.error(f"Structure file {ct_file} does not exist; "
                             "please obtain the file, e.g. using seismic fold")


class StructOneTableWriter(OneTableRelClusterGroupWriter, ABC):

    def iter_graphs(self, *,
                    rels: list[str],
                    cgroup: str,
                    struct_file: Iterable[str | Path] = (),
                    branch: str = "",
                    fold_coords: Iterable[tuple[str, int, int]] = (),
                    fold_primers: Iterable[tuple[str, DNA, DNA]] = (),
                    fold_regions_file: str | None = None,
                    fold_full: bool = opt_fold_full.default,
                    **kwargs):
        struct_files = list()
        for file in path.find_files_chain(struct_file, path.CT_FILE_LAST_SEGS):
            # Use a given CT file of an RNA structure, and determine the
            # structure region name from the file path.
            fields = path.parse(file, path.CT_FILE_LAST_SEGS)
            ref = fields[path.REF]
            if ref == self.table.ref:
                struct_files.append(file)
            else:
                logger.warning(f"Skipped CT file {file} in region directory "
                               f"{repr(path.REG)} in reference directory "
                               f"{repr(ref)}, which differs from the reference "
                               f"name of the table file {repr(self.table.ref)}")
        if struct_files:
            fold_regs = list()
        else:
            # Add the regions from the given coordinates/primers.
            ref_regions = RefRegions([(self.table.ref, self.table.refseq)],
                                     regs_file=optional_path(fold_regions_file),
                                     coords=fold_coords,
                                     primers=fold_primers,
                                     default_full=fold_full)
            fold_regs = [region.name
                         for region in ref_regions.list(self.table.ref)]
            if not fold_regs:
                # Add the table's region if no other regions were defined.
                fold_regs.append(self.table.reg)
        # Generate a graph for each cluster, relationship, and region.
        for cparams in cgroup_table(self.table, cgroup):
            kwparams = kwargs | cparams
            for rels_group in rels:
                for file in struct_files:
                    yield self.get_graph(rels_group,
                                         struct_file=file,
                                         struct_reg=None,
                                         branch=branch,
                                         **kwparams)
                for reg in fold_regs:
                    yield self.get_graph(rels_group,
                                         struct_file=None,
                                         struct_reg=reg,
                                         branch=branch,
                                         **kwparams)


class StructOneTableRunner(OneTableRelClusterGroupRunner, ABC):

    @classmethod
    def get_input_loader(cls):
        return load_foldable_tables

    @classmethod
    def get_var_params(cls):
        return super().get_var_params() + [opt_struct_file,
                                           opt_branch,
                                           opt_fold_regions_file,
                                           opt_fold_coords,
                                           opt_fold_primers,
                                           opt_fold_full]
