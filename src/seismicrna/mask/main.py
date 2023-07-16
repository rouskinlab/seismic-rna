from logging import getLogger
from pathlib import Path
from typing import Iterable

from click import command

from .write import mask_section
from ..relate.load import open_reports
from ..core import docdef, path
from ..core.cli import (arg_input_file,
                        opt_coords, opt_primers, opt_primer_gap, opt_sections_file,
                        opt_count_del, opt_count_ins, opt_discount_mut,
                        opt_exclude_polya, opt_exclude_gu, opt_exclude_pos,
                        opt_min_finfo_read, opt_max_fmut_read, opt_min_mut_gap,
                        opt_min_ninfo_pos, opt_max_fmut_pos,
                        opt_max_procs, opt_parallel, opt_rerun)
from ..core.parallel import dispatch
from ..core.sect import encode_primers, RefSections
from ..core.seq import DNA

logger = getLogger(__name__)

params = [
    # Input/output paths
    arg_input_file,
    # Sections
    opt_coords, opt_primers, opt_primer_gap, opt_sections_file,
    # Mutation counting
    opt_count_del, opt_count_ins, opt_discount_mut,
    # Filtering
    opt_exclude_polya, opt_exclude_gu, opt_exclude_pos,
    opt_min_finfo_read, opt_max_fmut_read, opt_min_mut_gap,
    opt_min_ninfo_pos, opt_max_fmut_pos,
    # Parallelization
    opt_max_procs, opt_parallel,
    # Effort
    opt_rerun,
]


@command(path.MOD_MASK, params=params)
def cli(*args, **kwargs):
    """ Select a section of the reference, define which relationships
    count as mutations, and filter out unusable positions and reads. """
    return run(*args, **kwargs)


@docdef.auto()
def run(input_file: tuple[str, ...], *,
        # Sections
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        sections_file: str,
        # Mutation counting
        count_del: bool,
        count_ins: bool,
        discount_mut: tuple[str, ...],
        # Filtering
        exclude_polya: int,
        exclude_gu: bool,
        exclude_pos: tuple[tuple[str, int], ...],
        min_finfo_read: float,
        max_fmut_read: int,
        min_mut_gap: int,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        # Parallelization
        max_procs: int,
        parallel: bool,
        # Effort
        rerun: bool) -> list[Path]:
    """ Run the filtering module. """
    # Open all relation vector loaders and get the sections for each.
    loaders, sections = open_sections(map(Path, input_file),
                                      coords=coords,
                                      primers=encode_primers(primers),
                                      primer_gap=primer_gap,
                                      sections_file=(Path(sections_file)
                                                     if sections_file
                                                     else None))
    # List the relation loaders and their sections.
    args = [(loader, section) for loader in loaders
            for section in sections.list(loader.ref)]
    # Define the keyword arguments.
    kwargs = dict(count_del=count_del,
                  count_ins=count_ins,
                  discount=discount_mut,
                  exclude_polya=exclude_polya,
                  exclude_gu=exclude_gu,
                  exclude_pos=exclude_pos,
                  min_finfo_read=min_finfo_read,
                  max_fmut_read=max_fmut_read,
                  min_mut_gap=min_mut_gap,
                  min_ninfo_pos=min_ninfo_pos,
                  max_fmut_pos=max_fmut_pos,
                  rerun=rerun)
    # Call the mutations and filter the mutation vectors.
    reports = dispatch(mask_section, max_procs=max_procs, parallel=parallel,
                       pass_n_procs=False, args=args, kwargs=kwargs)
    return list(map(Path, reports))


def open_sections(report_paths: Iterable[Path],
                  coords: Iterable[tuple[str, int, int]],
                  primers: Iterable[tuple[str, DNA, DNA]],
                  primer_gap: int,
                  sections_file: Path | None = None):
    """ Open sections of relate reports. """
    report_files = path.find_files_chain(report_paths, [path.RelateRepSeg])
    loaders = open_reports(report_files)
    sections = RefSections({(rep.ref, rep.seq) for rep in loaders},
                           coords=coords, primers=primers, primer_gap=primer_gap,
                           sects_file=sections_file)
    return loaders, sections
