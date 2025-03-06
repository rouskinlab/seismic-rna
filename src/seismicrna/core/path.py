from __future__ import annotations

import os
import pathlib
import re
import shutil
from abc import ABC, abstractmethod
from collections import Counter
from functools import cache, cached_property, partial, wraps
from itertools import chain, islice, product
from string import ascii_letters, digits, printable
from tempfile import mkdtemp
from typing import Any, Callable, Iterable, Sequence

from .logs import logger
from .validate import (require_isinstance,
                       require_issubclass,
                       require_isin,
                       require_equal,
                       require_atleast)

# Constants ############################################################

BRANCH_SEP = "_"
VERSUS_BRANCH = "VS"

# Valid/invalid characters in fields

PATH_CHARS = printable
ALPHANUM_CHARS = ascii_letters + digits
STR_CHARS = ALPHANUM_CHARS + "_.=+-"
STR_CHARS_SET = frozenset(STR_CHARS)
BRANCH_CHARS_SET = frozenset(STR_CHARS_SET - {BRANCH_SEP})
INT_CHARS = digits
PATH_PATTERN = f"([{PATH_CHARS}]+)"
STR_PATTERN = f"([{STR_CHARS}]+)"
BRANCHES_PATTERN = f"([{STR_CHARS}]*)"
INT_PATTERN = f"([{INT_CHARS}]+)"
STEP_PATTERN = f"([{ALPHANUM_CHARS}]+)"
RE_PATTERNS = {str: STR_PATTERN,
               int: INT_PATTERN,
               pathlib.Path: PATH_PATTERN}

# Names of steps
ALIGN_STEP = "align"
RELATE_STEP = "relate"
NAMES_BATCH = "names"
MASK_STEP = "mask"
CLUSTER_STEP = "cluster"
FOLD_STEP = "fold"
GRAPH_STEP = "graph"
LIST_STEP = "list"

# Directories for simulation

SIM_REFS_DIR = "refs"
SIM_PARAM_DIR = "params"
SIM_SAMPLES_DIR = "samples"

# Directories for stages

STAGE_ALIGN_INDEX = "index"
STAGE_ALIGN_INDEX_DEMULT = "index-demult"
STAGE_ALIGN_TRIM = "trim"
STAGE_ALIGN_MAP = "map"
STAGE_ALIGN_SORT = "sort"

STAGE_REL_SAMS = "sams"

STAGES = (STAGE_ALIGN_INDEX,
          STAGE_ALIGN_INDEX_DEMULT,
          STAGE_ALIGN_TRIM,
          STAGE_ALIGN_MAP,
          STAGE_ALIGN_SORT,
          STAGE_REL_SAMS)

# Cluster information
CLUST_PARAM_PIS = "pis"
CLUST_PARAM_MUS = "mus"
CLUST_PARAMS = (CLUST_PARAM_PIS,
                CLUST_PARAM_MUS)
CLUST_PARAMS_DIR = "parameters"
CLUST_STATS_DIR = "statistics"
CLUST_COUNTS_DIR = "read-counts"

TABLES = (RELATE_STEP, MASK_STEP, CLUSTER_STEP)

# File extensions

GZIP_EXT = ".gz"
TXT_EXT = ".txt"
CSV_EXT = ".csv"
CSVZIP_EXT = f"{CSV_EXT}{GZIP_EXT}"
CSV_EXTS = CSV_EXT, CSVZIP_EXT
BRICKLE_EXT = ".brickle"
JSON_EXT = ".json"
FASTA_EXTS = ".fa", ".fna", ".fasta"
BOWTIE2_INDEX_EXTS = (".1.bt2",
                      ".2.bt2",
                      ".3.bt2",
                      ".4.bt2",
                      ".rev.1.bt2",
                      ".rev.2.bt2")
FQ_EXTS = (".fq.gz",
           ".fastq.gz",
           ".fq",
           ".fastq",
           "_001.fq.gz",
           "_001.fastq.gz",
           "_001.fq",
           "_001.fastq")
FQ_PAIRED_EXTS_TEMPLATES = "_R{}{}", "_{}{}", "_mate{}{}", "_{}_sequence{}"
FQ1_EXTS = tuple(template.format(1, ext) for template, ext in
                 product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
FQ2_EXTS = tuple(template.format(2, ext) for template, ext in
                 product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
SAM_EXT = ".sam"
BAM_EXT = ".bam"
CRAM_EXT = ".cram"
XAM_EXTS = SAM_EXT, BAM_EXT, CRAM_EXT
FAI_EXT = ".fai"
CT_EXT = ".ct"
DB_EXT = ".db"
DBN_EXT = ".dbn"
DOT_EXT = ".dot"
DOT_EXTS = DB_EXT, DBN_EXT, DOT_EXT
DMS_EXT = ".dms"
KTS_EXT = ".kts"
HTML_EXT = ".html"
SVG_EXT = ".svg"
PDF_EXT = ".pdf"
PNG_EXT = ".png"
GRAPH_EXTS = CSV_EXT, HTML_EXT, SVG_EXT, PDF_EXT, PNG_EXT
PARAM_MUTS_EXT = f".muts{CSV_EXT}"
PARAM_ENDS_EXT = f".ends{CSV_EXT}"
PARAM_CLUSTS_EXT = f".clusts{CSV_EXT}"


# Path Exceptions ######################################################

class PathError(Exception):
    """ Any error involving a path. """


class PathTypeError(PathError, TypeError):
    """ Use of the wrong type of path or segment. """


class PathValueError(PathError, ValueError):
    """ Invalid value of a path segment field. """


class PathNotFoundError(PathError, FileNotFoundError):
    """ Path does not exist but should. """


class PathExistsError(PathError, FileExistsError):
    """ Path exists but should not. """


class WrongFileExtensionError(PathValueError):
    """ A file has the wrong extension. """


# Path Functions #######################################################

def fill_whitespace(path: str | pathlib.Path,
                    fill: str = "_") -> str | pathlib.Path:
    """ Replace all whitespace in `path` with `fill`. """
    return path.__class__(fill.join(str(path).split()))


def sanitize(path: str | pathlib.Path, strict: bool = False):
    """ Sanitize a path-like object by ensuring it is an absolute path,
    eliminating symbolic links and redundant path separators/references,
    and returning a Path object.

    Parameters
    ----------
    path: str | pathlib.Path
        Path to sanitize.
    strict: bool = False
        Require the path to exist and contain no symbolic link loops.

    Returns
    -------
    pathlib.Path
        Absolute, normalized, symlink-free path.
    """
    return pathlib.Path(path).resolve(strict=strict)


@cache
def get_seismicrna_source_dir():
    """ SEISMIC-RNA source directory, named seismicrna, containing
    __init__.py and the top-level modules and subpackages. """
    seismicrna_src_dir = sanitize(__file__, strict=True).parent.parent
    try:
        from seismicrna import __file__ as seismicrna_file
    except ImportError:
        seismicrna_file = None
    if seismicrna_file:
        require_equal("source directory",
                      sanitize(seismicrna_file).parent,
                      seismicrna_src_dir,
                      classes=pathlib.Path,
                      error_type=PathValueError)
    else:
        logger.warning("seismicrna is not installed: skipped verifying path")
    require_equal("source directory name",
                  seismicrna_src_dir.name,
                  "seismicrna",
                  error_type=PathValueError)
    return seismicrna_src_dir


@cache
def get_seismicrna_project_dir():
    """ SEISMIC-RNA project directory, named seismic-rna, containing
    src, pyproject.toml, and all other project files. Will exist if the
    entire SEISMIC-RNA project has been downloaded, e.g. from GitHub,
    but not if SEISMIC-RNA was only installed using pip or conda. """
    seismicrna_prj_dir = get_seismicrna_source_dir().parent.parent
    name = "seismic-rna"
    if seismicrna_prj_dir.name != name:
        # It is fine if the project directory does not exist because
        # installing SEISMIC-RNA using pip or conda installs only the
        # source directory, but not the project directory.
        return None
    return seismicrna_prj_dir


# Path Fields ##########################################################

# Field validation functions

def validate_str(txt: str):
    require_isinstance("txt", txt, str, error_type=PathTypeError)
    if not txt:
        raise PathValueError("txt cannot be empty string")
    if illegal := sorted(set(txt) - STR_CHARS_SET):
        raise PathValueError(
            f"txt {repr(txt)} has illegal characters: {illegal}"
        )


def validate_top(top: pathlib.Path):
    require_isinstance("top", top, pathlib.Path, error_type=PathTypeError)
    if not top.parent.is_dir():
        raise PathNotFoundError(top.parent)
    if top.is_file():
        raise PathExistsError(top)


def validate_int(num: int):
    require_isinstance("num", num, int, error_type=PathTypeError)
    require_atleast("num", num, 0, error_type=PathValueError)


def validate_branch(branch: str):
    require_isinstance("branch", branch, str, error_type=PathTypeError)
    if illegal := sorted(set(branch) - BRANCH_CHARS_SET):
        raise PathValueError(
            f"branch {repr(branch)} has illegal characters: {illegal}"
        )


def validate_branches(branches: dict[str, str]):
    require_isinstance("branches", branches, dict, error_type=PathTypeError)
    for step, branch in branches.items():
        validate_str(step)
        validate_branch(branch)


def get_ancestors(branches: dict[str, str]):
    """ Get all but the last branch in a dict of branches. """
    validate_branches(branches)
    if not branches:
        return dict()
    return dict(islice(branches.items(), len(branches) - 1))


def add_branch(step: str, branch: str, ancestors: dict[str, str]):
    """ Add a new branch to a dict of branches. """
    validate_str(step)
    validate_branch(branch)
    validate_branches(ancestors)
    if step in ancestors:
        raise PathValueError(
            f"A step named {repr(step)} already exists in {ancestors}"
        )
    return {**ancestors, step: branch}


def flatten_branches(branches: dict[str, str]):
    validate_branches(branches)
    return list(filter(None, branches.values()))


def validate_branches_flat(branches_flat: list[str]):
    require_isinstance("branches_flat",
                       branches_flat,
                       list,
                       error_type=PathTypeError)
    for branch in branches_flat:
        validate_branch(branch)
        if not branch:
            raise PathValueError("branch cannot be the empty string")


VALIDATE = {int: validate_int,
            str: validate_str,
            pathlib.Path: validate_top}


# Field class

class PathField(object):
    __slots__ = ["dtype", "options", "is_ext", "pattern"]

    def __init__(self,
                 dtype: type[str | int | pathlib.Path | list],
                 options: Iterable = (),
                 is_ext: bool = False,
                 pattern: str = ""):
        self.dtype = dtype
        self.options = list(options)
        if self.dtype is list and self.options:
            raise PathValueError("Cannot take options if the data type is list")
        for i, option in enumerate(self.options):
            require_isinstance(f"option {i}", option, self.dtype,
                               error_type=PathTypeError)
        self.is_ext = is_ext
        if self.is_ext:
            require_issubclass("Extension field dtype", self.dtype, str,
                               error_type=PathValueError)
            if not self.options:
                raise PathValueError("Extension field must have options")
            if pattern:
                raise PathValueError("Extension field cannot have a pattern")
            self.pattern = ""
        else:
            require_isinstance("pattern", pattern, str,
                               error_type=PathTypeError)
            if pattern:
                self.pattern = pattern
            else:
                try:
                    self.pattern = RE_PATTERNS[self.dtype]
                except KeyError:
                    raise PathTypeError(
                        f"No default pattern for {self.dtype}"
                    ) from None

    def validate(self, val: Any):
        """ Validate a value before turning it into a string. """
        validate_func = VALIDATE[self.dtype]
        validate_func(val)
        if self.options:
            require_isin("option", val, self.options,
                         error_type=PathValueError)

    def build(self, val: Any):
        """ Validate a value and return it as a string. """
        self.validate(val)
        return str(val)

    def parse(self, text: str) -> Any:
        """ Parse a value from a string, validate it, and return it. """
        require_isinstance("text", text, str, error_type=PathTypeError)
        try:
            val = self.dtype(text)
        except Exception as error:
            raise PathValueError(f"Could not interpret {repr(text)} as type "
                                 f"{self.dtype}: {error}") from None
        self.validate(val)
        return val

    def __str__(self):
        return f"{type(self).__name__} {self.dtype.__name__}"


class BranchesPathField(PathField):
    """ The field for branches requires special functions. """

    def __init__(self):
        super().__init__(list, (), False, BRANCHES_PATTERN)

    def validate(self, val: Any):
        validate_branches_flat(val)

    def build(self, val: Any):
        if isinstance(val, dict):
            val = flatten_branches(val)
        else:
            self.validate(val)
        return "".join(f"{BRANCH_SEP}{branch}" for branch in val)

    def parse(self, text: str):
        require_isinstance("text", text, str, error_type=PathTypeError)
        val = list(filter(None, text.split(BRANCH_SEP)))
        self.validate(val)
        return val


# Fields
TopField = PathField(pathlib.Path)
NameField = PathField(str)
StepField = PathField(str,
                      [ALIGN_STEP,
                       RELATE_STEP,
                       MASK_STEP,
                       CLUSTER_STEP,
                       FOLD_STEP,
                       GRAPH_STEP],
                      pattern=STEP_PATTERN)
BranchesField = BranchesPathField()
StageField = PathField(str, STAGES)
IntField = PathField(int)
ClustRunResultsField = PathField(str, CLUST_PARAMS)
PosTableField = PathField(str, TABLES)
ReadTableField = PathField(str, TABLES)
AbundanceField = PathField(str, [CLUSTER_STEP])

# File extensions
TextExt = PathField(str, [TXT_EXT], is_ext=True)
ReportExt = PathField(str, [JSON_EXT], is_ext=True)
RefseqFileExt = PathField(str, [BRICKLE_EXT], is_ext=True)
BatchExt = PathField(str, [BRICKLE_EXT], is_ext=True)
ClustTabExt = PathField(str, CSV_EXTS, is_ext=True)
PosTableExt = PathField(str, [CSV_EXT], is_ext=True)
ReadTableExt = PathField(str, [CSVZIP_EXT], is_ext=True)
AbundanceExt = PathField(str, [CSV_EXT], is_ext=True)
FastaExt = PathField(str, FASTA_EXTS, is_ext=True)
FastaIndexExt = PathField(str, BOWTIE2_INDEX_EXTS, is_ext=True)
FastqExt = PathField(str, FQ_EXTS, is_ext=True)
Fastq1Ext = PathField(str, FQ1_EXTS, is_ext=True)
Fastq2Ext = PathField(str, FQ2_EXTS, is_ext=True)
XamExt = PathField(str, XAM_EXTS, is_ext=True)
ConnectTableExt = PathField(str, [CT_EXT], is_ext=True)
DotBracketExt = PathField(str, DOT_EXTS, is_ext=True)
DmsReactsExt = PathField(str, [DMS_EXT], is_ext=True)
GraphExt = PathField(str, GRAPH_EXTS, is_ext=True)
WebAppFileExt = PathField(str, [JSON_EXT], is_ext=True)
SvgExt = PathField(str, [SVG_EXT], is_ext=True)
PngExt = PathField(str, [PNG_EXT], is_ext=True)
KtsExt = PathField(str, [KTS_EXT], is_ext=True)


def check_file_extension(file: str | pathlib.Path,
                         extensions: Iterable[str] | PathField):
    if isinstance(extensions, PathField):
        if not extensions.is_ext:
            raise PathValueError(f"{extensions} is not an extension field")
        extensions = extensions.options
    elif not isinstance(extensions, (tuple, list, set, dict)):
        extensions = set(extensions)
    if not isinstance(file, pathlib.Path):
        file = pathlib.Path(file)
    require_isin("file extension", file.suffix, extensions,
                 error_type=WrongFileExtensionError)


# Path Segments ########################################################

# Segment class

class PathSegment(object):

    def __init__(self,
                 segment_name: str,
                 field_types: dict[str, PathField], *,
                 order: int = 0,
                 frmt: str | None = None):
        self.name = segment_name
        self.field_types = field_types
        # Verify that a field has the key EXT if and only if it is an
        # extension and is the last field in the segment.
        for i, (name, field) in enumerate(self.field_types.items(), start=1):
            if name == EXT:
                if not field.is_ext:
                    raise PathValueError(
                        f"Field {repr(EXT)} is not an extension"
                    )
                if i != len(self.field_types):
                    raise PathValueError(
                        f"Extension of {self} is not the last field"
                    )
                require_equal("order of segment with an extension", order, 0,
                              error_type=PathValueError)
            elif field.is_ext:
                raise PathValueError(f"{self} extension has name {repr(name)}")
        if not any(ft in self.field_types for ft in [EXT, TOP]):
            require_atleast("order of segment without an extension", order, 1,
                            error_type=PathValueError)
        self.order = order
        # Determine the format string.
        if frmt is None:
            # Default format is to concatenate all the fields.
            frmt = "".join("{" + name + "}" for name, field
                           in self.field_types.items())
        self.frmt = frmt
        # Generate the pattern string using the format string, excluding
        # the extension (if any) because before parsing, it is removed
        # from the end of the string.
        patterns = {name: field.pattern
                    for name, field in self.field_types.items()}
        self.ptrn = re.compile(self.frmt.format(**patterns))

    @property
    def ext_type(self):
        """ Type of the segment's file extension, or None if it has no
        file extension. """
        return self.field_types.get(EXT)

    @cached_property
    def exts(self) -> list[str]:
        """ Valid file extensions of the segment. """
        if self.ext_type is None:
            return list()
        if not self.ext_type.options:
            raise PathValueError(
                f"{self} extension {self.ext_type} has no options"
            )
        return self.ext_type.options

    def match_longest_ext(self, text: str):
        """ Find the longest extension of the given text that matches a
        valid file extension. If none match, return None. """
        # Iterate over the file extensions from longest to shortest.
        for ext in sorted(self.exts, key=len, reverse=True):
            if text.endswith(ext):
                # The text ends with this extension, so it must be the
                # longest valid extension in which the text ends.
                return ext
        return

    def build(self, vals: dict[str, Any]):
        require_isinstance("vals", vals, dict, error_type=PathTypeError)
        # Verify that a value is given for every field, with no extras.
        require_equal("segment fields",
                      sorted(vals.keys()),
                      sorted(self.field_types.keys()),
                      error_type=PathValueError)
        # Validate the value passed to every field.
        fields = {name: field.build(vals[name])
                  for name, field in self.field_types.items()}
        # Return the formatted segment.
        return self.frmt.format(**fields)

    def parse(self, text: str):
        ext = None
        if self.ext_type is not None:
            # If the segment has a file extension, then determine the
            # longest valid file extension that matches the text.
            ext = self.match_longest_ext(text)
            if ext is None:
                raise PathValueError(f"text {repr(text)} has no file "
                                     f"extension; expected one of {self.exts}")
            # Remove the file extension from the end of the text.
            text = text[: -len(ext)]
        # Try to parse the text (with the extension, if any, removed).
        match = self.ptrn.match(text)
        if not match:
            raise PathValueError(f"Could not parse fields in text {repr(text)} "
                                 f"using pattern {repr(self.ptrn)}")
        vals = list(match.groups())
        # If there is an extension field, add its value back to the end
        # of the parsed values.
        if ext is not None:
            vals.append(ext)
        # Return a dict of the names of the fields in the segment and
        # their parsed values.
        return {name: field.parse(group) for (name, field), group
                in zip(self.field_types.items(), vals, strict=True)}

    def __str__(self):
        return f"{type(self).__name__} {repr(self.name)}"


# Field names

TOP = "top"
STAGE = "stage"
STEP = "step"
BRANCHES = "branches"
SAMPLE = "sample"
REF = "ref"
REG = "reg"
BATCH = "batch"
TABLE = "table"
LIST = "list"
NCLUST = "k"
RUN = "run"
PROFILE = "profile"
GRAPH = "graph"
EXT = "ext"
STRUCT = "struct"

# Directory segments

TopSeg = PathSegment("top-dir", {TOP: TopField}, order=-1)
StageSeg = PathSegment("stage-dir", {STAGE: StageField}, order=70)
SampSeg = PathSegment("sample-dir", {SAMPLE: NameField}, order=60)
StepSeg = PathSegment("command-dir",
                      {STEP: StepField, BRANCHES: BranchesField},
                      order=50)
RefSeg = PathSegment("ref-dir", {REF: NameField}, order=30)
RegSeg = PathSegment("reg-dir", {REG: NameField}, order=20)

# File segments

# FASTA
FastaSeg = PathSegment("fasta", {REF: NameField, EXT: FastaExt})
FastaIndexSeg = PathSegment("fasta-index", {REF: NameField, EXT: FastaIndexExt})

# FASTQ
FastqSeg = PathSegment("fastq", {SAMPLE: NameField, EXT: FastqExt})
Fastq1Seg = PathSegment("fastq1", {SAMPLE: NameField, EXT: Fastq1Ext})
Fastq2Seg = PathSegment("fastq2", {SAMPLE: NameField, EXT: Fastq2Ext})

# Demultiplexed FASTQ
DmFastqSeg = PathSegment("dm-fastq", {REF: NameField, EXT: FastqExt})
DmFastq1Seg = PathSegment("dm-fastq1", {REF: NameField, EXT: Fastq1Ext})
DmFastq2Seg = PathSegment("dm-fastq2", {REF: NameField, EXT: Fastq2Ext})

# Align
XamSeg = PathSegment("xam", {REF: NameField, EXT: XamExt})
AlignSampleRepSeg = PathSegment("align-samp-rep",
                                {EXT: ReportExt},
                                frmt="align-report{ext}")
AlignRefRepSeg = PathSegment("align-ref-rep",
                             {REF: NameField, EXT: ReportExt},
                             frmt="{ref}__align-report{ext}")
SplitRepSeg = PathSegment("split-rep",
                          {EXT: ReportExt},
                          frmt="split-report{ext}")

# Relate
RefseqFileSeg = PathSegment("refseq-file",
                            {EXT: RefseqFileExt},
                            frmt="refseq{ext}")
ReadNamesBatSeg = PathSegment("names-bat",
                              {BATCH: IntField, EXT: BatchExt},
                              frmt=NAMES_BATCH + "-batch-{batch}{ext}")
RelateBatSeg = PathSegment("relate-bat",
                           {BATCH: IntField, EXT: BatchExt},
                           frmt=RELATE_STEP + "-batch-{batch}{ext}")
RelateRepSeg = PathSegment("relate-rep",
                           {EXT: ReportExt},
                           frmt=RELATE_STEP + "-report{ext}")

# Mask
MaskBatSeg = PathSegment(f"{MASK_STEP}-bat",
                         {BATCH: IntField, EXT: BatchExt},
                         frmt=MASK_STEP + "-batch-{batch}{ext}")
MaskRepSeg = PathSegment("mask-rep",
                         {EXT: ReportExt},
                         frmt=MASK_STEP + "-report{ext}")

# Cluster
ClustParamsDirSeg = PathSegment("cluster-run-res-dir",
                                {},
                                frmt=CLUST_PARAMS_DIR,
                                order=10)
ClustParamsFileSeg = PathSegment("cluster-run-res",
                                 {TABLE: ClustRunResultsField,
                                  NCLUST: IntField,
                                  RUN: IntField,
                                  EXT: ClustTabExt},
                                 frmt="k{k}-r{run}_{table}{ext}")
ClustBatSeg = PathSegment("cluster-bat",
                          {BATCH: IntField, EXT: BatchExt},
                          frmt=CLUSTER_STEP + "-batch-{batch}{ext}")
ClustRepSeg = PathSegment("cluster-rep",
                          {EXT: ReportExt},
                          frmt=CLUSTER_STEP + "-report{ext}")

# Table
PositionTableSeg = PathSegment("position-table",
                               {TABLE: PosTableField, EXT: PosTableExt},
                               frmt="{table}-position-table{ext}")
ReadTableSeg = PathSegment("read-table",
                           {TABLE: ReadTableField, EXT: ReadTableExt},
                           frmt="{table}-read-table{ext}")
AbundanceTableSeg = PathSegment("abundance-table",
                                {TABLE: AbundanceField, EXT: AbundanceExt},
                                frmt="{table}-abundance-table{ext}")

# List
PositionListSeg = PathSegment("position-list",
                              {LIST: PosTableField, EXT: PosTableExt},
                              frmt="{list}-position-list{ext}")
ReadListSeg = PathSegment("read-list",
                          {LIST: ReadTableField, EXT: ReadTableExt},
                          frmt="{list}-read-list{ext}")

# Fold
FoldRepSeg = PathSegment("fold-rep",
                         {PROFILE: NameField, EXT: ReportExt},
                         frmt="{profile}__fold-report{ext}")
ConnectTableSeg = PathSegment("rna-ct",
                              {PROFILE: NameField, EXT: ConnectTableExt})
DotBracketSeg = PathSegment("rna-dot",
                            {PROFILE: NameField, EXT: DotBracketExt})
DmsReactsSeg = PathSegment("dms-reacts",
                           {PROFILE: NameField, EXT: DmsReactsExt})
VarnaColorSeg = PathSegment("varna-color",
                            {PROFILE: NameField, EXT: TextExt},
                            frmt="{profile}__varna-color{ext}")

# Draw
SvgSeg = PathSegment("svg",
                     {PROFILE: NameField, STRUCT: IntField, EXT: SvgExt},
                     frmt="{profile}-{struct}{ext}")
PngSeg = PathSegment("png",
                     {PROFILE: NameField, STRUCT: IntField, EXT: PngExt},
                     frmt="{profile}-{struct}{ext}")
KtsSeg = PathSegment("kts",
                     {PROFILE: NameField, STRUCT: IntField, EXT: KtsExt},
                     frmt="{profile}-{struct}{ext}")

# Graphs
GraphSeg = PathSegment("graph", {GRAPH: NameField, EXT: GraphExt})

# Web App Export
WebAppFileSeg = PathSegment("webapp",
                            {SAMPLE: NameField, EXT: WebAppFileExt},
                            frmt="{sample}__webapp{ext}")

# Path segment patterns
STEP_DIR_SEGS = SampSeg, StepSeg
REF_DIR_SEGS = STEP_DIR_SEGS + (RefSeg,)
REG_DIR_SEGS = REF_DIR_SEGS + (RegSeg,)
STAGE_DIR_SEGS = SampSeg, StepSeg, StageSeg
FASTA_STAGE_SEGS = StageSeg, FastaSeg
FASTA_INDEX_DIR_STAGE_SEGS = StageSeg, RefSeg
FASTQ_SEGS = FastqSeg,
FASTQ1_SEGS = Fastq1Seg,
FASTQ2_SEGS = Fastq2Seg,
DMFASTQ_SEGS = SampSeg, DmFastqSeg
DMFASTQ1_SEGS = SampSeg, DmFastq1Seg
DMFASTQ2_SEGS = SampSeg, DmFastq2Seg
XAM_SEGS = STEP_DIR_SEGS + (XamSeg,)
XAM_STAGE_SEGS = STAGE_DIR_SEGS + (XamSeg,)
CLUST_TAB_SEGS = REG_DIR_SEGS + (ClustParamsDirSeg, ClustParamsFileSeg)
CT_FILE_ALL_SEGS = REG_DIR_SEGS + (ConnectTableSeg,)
CT_FILE_LAST_SEGS = CT_FILE_ALL_SEGS[-3:]
DB_FILE_ALL_SEGS = REG_DIR_SEGS + (DotBracketSeg,)
DB_FILE_LAST_SEGS = DB_FILE_ALL_SEGS[-3:]


# Paths ################################################################


# Path class

class Path(object):

    def __init__(self, seg_types: Iterable[PathSegment]):
        # Sort the non-redundant segment types in the path from largest
        # to smallest value of their order attribute.
        self.seg_types = sorted(set(seg_types),
                                key=lambda segt: segt.order,
                                reverse=True)
        # Check for TopSeg.
        if TopSeg in self.seg_types:
            raise PathValueError(f"{TopSeg} may not be given in seg_types")
        self.seg_types.insert(0, TopSeg)
        # Check for duplicate orders.
        if max(Counter(segt.order for segt in self.seg_types).values()) > 1:
            raise PathValueError(
                f"Got duplicate order values in {self.seg_types}"
            )

    def build(self, fields: dict[str, Any]):
        """ Return a `pathlib.Path` instance by assembling the given
        `fields` into a full path. """
        # Build the new path one segment at a time.
        segments = list()
        used_fields = set()
        for seg_type in self.seg_types:
            # For each type of segment in the path, try to get the names
            # and values of all fields of the segment.
            try:
                seg_fields = {name: fields[name]
                              for name in seg_type.field_types}
            except KeyError as error:
                raise PathValueError(f"Missing field for {seg_type}: {error}")
            used_fields.update(seg_fields)
            # Generate a string representation of the segment using the
            # values of its fields, and add it to the growing path.
            segments.append(seg_type.build(seg_fields))
        # Check whether any fields were given but not used by the path.
        extras = fields.keys() - used_fields
        if extras:
            raise PathValueError(f"Extra fields for {type(self)}: {extras}")
        # Assemble the segment strings into a path, and return it.
        path = pathlib.Path(*segments)
        return path

    def parse(self, path: str | pathlib.Path):
        """ Return the field names and values from a given path. """
        # Convert the given path into a canonical, absolute path.
        path = str(sanitize(path))
        # Get the field names and values one segment at a time.
        fields = dict()
        # Iterate from the deepest (last) to shallowest (first) segment.
        for seg_type in reversed(self.seg_types):
            if seg_type is TopSeg:
                # The top-most segment has been reached and must be used
                # to parse the entire remaining path.
                tail = path
            else:
                # The top-most segment of the path has not been reached.
                # Split off the deepest part of the path (tail), and
                # parse it using the current segment type.
                path, tail = os.path.split(path)
            # Verify that the entire path has not been consumed.
            if not tail:
                raise PathValueError(f"No path remains to parse {seg_type}")
            # Parse the deepest part of the path to obtain the fields,
            # and use them to update the field names and values.
            fields.update(seg_type.parse(tail))
        return fields

    def __str__(self):
        return f"{type(self).__name__}: {list(map(str, self.seg_types))}"


# mkdir/symlink/rmdir.


def mkdir_if_needed(path: pathlib.Path | str):
    """ Create a directory and log that event if it does not exist. """
    path = sanitize(path, strict=False)
    try:
        path.mkdir(parents=True)
    except FileExistsError:
        if not path.is_dir():
            # Raise an error if the existing path is not a directory,
            # e.g. if it is a file.
            raise NotADirectoryError(path) from None
        return path
    logger.action(f"Created directory {path}")
    return path


def symlink_if_needed(link_path: pathlib.Path | str,
                      target_path: pathlib.Path | str):
    """ Make link_path a link pointing to target_path and log that event
    if it does not exist. """
    link_path = pathlib.Path(link_path)
    target_path = sanitize(target_path, strict=True)
    try:
        link_path.symlink_to(target_path)
    except FileExistsError:
        # link_path already exists, so make sure it is a symbolic link
        # that points to target_path.
        try:
            readlink = link_path.readlink()
        except OSError:
            raise OSError(f"{link_path} is not a symbolic link") from None
        if readlink != target_path:
            raise OSError(f"{link_path} is a symbolic link to {readlink}, "
                          f"not to {target_path}") from None
        return link_path
    logger.action(f"Made {link_path} a symbolic link to {target_path}")
    return link_path


def rmdir_if_needed(path: pathlib.Path | str,
                    rmtree: bool = False,
                    rmtree_ignore_errors: bool = False,
                    raise_on_rmtree_error: bool = True):
    """ Remove a directory and log that event if it exists. """
    path = sanitize(path, strict=False)
    try:
        path.rmdir()
    except FileNotFoundError:
        # The path does not exist, so there is no need to delete it.
        # FileNotFoundError is a subclass of OSError, so need to handle
        # this exception before OSError.
        logger.detail(f"Skipped removing directory {path}: does not exist")
        return path
    except NotADirectoryError:
        # Trying to rmdir() something that is not a directory should
        # always raise an error. NotADirectoryError is a subclass of
        # OSError, so need to handle this exception before OSError.
        raise
    except OSError:
        # The directory exists but could not be removed for some reason,
        # probably that it is not empty.
        if not rmtree:
            # For safety, remove directories recursively only if given
            # explicit permission to do so; if not, re-raise the error.
            raise
        try:
            shutil.rmtree(path, ignore_errors=rmtree_ignore_errors)
        except Exception as error:
            if raise_on_rmtree_error:
                raise
            # If not raising errors, then log a warning but return now
            # to avoid logging that the directory was removed.
            logger.warning(error)
            return path
    logger.action(f"Deleted directory {path}")


# Path creation routines


# The purpose of this function (which just wraps Path(segment_types)
# is to cache every type of Path; thus, segment_types must be a hashable
# sequence, i.e. a tuple.
@cache
def create_path_type(segment_types: tuple[PathSegment, ...]):
    """ Create and cache a Path instance from the segment types. """
    return Path(segment_types)


def build(segment_types: Iterable[PathSegment], field_values: dict[str, Any]):
    """ Return a `pathlib.Path` from the segment types and field values.
    """
    return create_path_type(tuple(segment_types)).build(field_values)


def builddir(segment_types: Iterable[PathSegment], field_values: dict[str, Any]):
    """ Build the path and create it on the file system as a directory
    if it does not already exist. """
    return mkdir_if_needed(build(segment_types, field_values))


def buildpar(segment_types: Iterable[PathSegment], field_values: dict[str, Any]):
    """ Build a path and create its parent directory if it does not
    already exist. """
    path = build(segment_types, field_values)
    mkdir_if_needed(path.parent)
    return path


def randdir(parent: str | pathlib.Path | None = None,
            prefix: str = "",
            suffix: str = ""):
    """ Build a path of a new directory that does not exist and create
    it on the file system. """
    parent = sanitize(parent) if parent is not None else pathlib.Path.cwd()
    path = pathlib.Path(mkdtemp(dir=parent, prefix=prefix, suffix=suffix))
    logger.action(f"Created directory {path}")
    return path


# Path parsing routines

def get_fields_in_seg_types(segment_types: Iterable[PathSegment],
                            include_top: bool = False) -> dict[str, PathField]:
    """ Get all fields among the given segment types. """
    fields_no_top = {field_name: field
                     for segment_type in segment_types
                     for field_name, field in segment_type.field_types.items()}
    return {TOP: TopField, **fields_no_top} if include_top else fields_no_top


def deduplicate(paths: Iterable[str | pathlib.Path], warn: bool = True):
    """ Yield the non-redundant paths. """
    total = 0
    seen = set()
    for path in map(sanitize, paths):
        total += 1
        if path in seen:
            if warn:
                logger.warning(f"Duplicate path: {path}")
        else:
            seen.add(path)
            yield path


def deduplicated(func: Callable):
    """ Decorate a Path generator to yield non-redundant paths. """

    @wraps(func)
    def wrapper(*args, **kwargs):
        yield from deduplicate(func(*args, **kwargs))

    return wrapper


def parse(path: str | pathlib.Path, segment_types: Iterable[PathSegment]):
    """ Return the fields of a path based on the segment types. """
    return create_path_type(tuple(segment_types)).parse(path)


def parse_top_separate(path: str | pathlib.Path,
                       segment_types: Iterable[PathSegment]):
    """ Return the fields of a path, and the `top` field separately. """
    field_values = parse(path, segment_types)
    return field_values.pop(TOP), field_values


def path_matches(path: str | pathlib.Path, segments: Iterable[PathSegment]):
    """ Check if a path matches a sequence of path segments.

    Parameters
    ----------
    path: str | pathlib.Path
        Path of the file/directory.
    segments: Iterable[PathSegment]
        Sequence of path segments to check if the file matches.

    Returns
    -------
    bool
        Whether the path matches any given sequence of path segments.
    """
    # Parsing the path will succeed if and only if it matches the
    # sequence of path segments.
    try:
        parse(path, segments)
    except PathError:
        # The path does not match this sequence of path segments.
        return False
    else:
        # The path matches this sequence of path segments.
        return True


@deduplicated
def find_files(path: str | pathlib.Path,
               segments: Sequence[PathSegment],
               pre_sanitize: bool = True):
    """ Yield all files that match a sequence of path segments.
    The behavior depends on what `path` is:

    - If it is a file, then yield `path` if it matches the segments;
      otherwise, yield nothing.
    - If it is a directory, then search it recursively and yield every
      matching file in the directory and its subdirectories.

    Parameters
    ----------
    path: str | pathlib.Path
        Path of a file to check or a directory to search recursively.
    segments: Sequence[PathSegment]
        Sequence(s) of Path segments to check if each file matches.
    pre_sanitize: bool
        Whether to sanitize the path before searching it.

    Returns
    -------
    Generator[Path, Any, None]
        Paths of files matching the segments.
    """
    if pre_sanitize:
        path = sanitize(path, strict=True)
    if path.is_file():
        # Check if the file matches the segments.
        if path_matches(path, segments):
            # If so, then yield it.
            logger.detail(f"Found file {path}")
            yield path
    else:
        # Search the directory for files matching the segments.
        logger.routine(f"Began recursively searching directory {path}")
        yield from chain(*map(partial(find_files,
                                      segments=segments,
                                      pre_sanitize=False),
                              path.iterdir()))
        logger.routine(f"Ended recursively searching directory {path}")


@deduplicated
def find_files_chain(paths: Iterable[str | pathlib.Path],
                     segments: Sequence[PathSegment]):
    """ Yield from `find_files` called on every path in `paths`. """
    for path in deduplicate(paths):
        try:
            yield from find_files(path, segments)
        except Exception as error:
            logger.error(error)


# Path transformation routines


def cast_path(input_path: str | pathlib.Path,
              input_segments: Sequence[PathSegment],
              output_segments: Sequence[PathSegment],
              override: dict[str, Any] | None = None):
    """ Cast `input_path` made of `input_segments` to a new path made of
    `output_segments`.

    Parameters
    ----------
    input_path: str | pathlib.Path
        Input path from which to take the path fields.
    input_segments: Sequence[PathSegment]
        Path segments to use to determine the fields in `input_path`.
    output_segments: Sequence[PathSegment]
        Path segments to use to determine the fields in `output_path`.
    override: dict[str, Any] | None
        Override and supplement the fields in `input_path`.

    Returns
    -------
    pathlib.Path
        Path comprising `output_segments` made of fields in `input_path`
        (as determined by `input_segments`).
    """
    output_segments = tuple(output_segments)
    # Extract the fields from the input path using the input segments.
    top, fields = parse_top_separate(input_path, input_segments)
    if override:
        # Override and supplement the fields in the input path.
        fields |= override
    # Normalize the fields to comply with the output segments.
    fields = {field_name: fields[field_name]
              for field_name in get_fields_in_seg_types(output_segments)}
    # Generate a new output path from the normalized fields.
    fields[TOP] = top
    output_path = build(output_segments, fields)
    return output_path


def transpath(to_dir: str | pathlib.Path,
              from_dir: str | pathlib.Path,
              path: str | pathlib.Path,
              strict: bool = False):
    """ Return the path that would be produced by moving `path` from
    `from_dir` to `to_dir` (but do not actually move the path on the
    file system). This function does not require that any of the given
    paths exist, unless `strict` is True.

    Parameters
    ----------
    to_dir: str | pathlib.Path
        Directory to which to move `path`.
    from_dir: str | pathlib.Path
        Directory from which to move `path`; must contain `path` but not
        necessarily be the direct parent directory of `path`.
    path: str | pathlib.Path
        Path to move; can be a file or directory.
    strict: bool = False
        Require that all paths exist and contain no symbolic link loops.

    Returns
    -------
    pathlib.Path
        Hypothetical path after moving `path` from `indir` to `outdir`.
    """
    # Ensure from_dir is sanitized.
    from_dir = sanitize(from_dir, strict)
    # Find the part of the given path relative to from_dir.
    relpath = sanitize(path, strict).relative_to(from_dir)
    if relpath == pathlib.Path():
        # If the relative path is empty, then use the parent directory
        # of from_dir instead.
        return transpath(to_dir, from_dir.parent, path, strict)
    # Append the relative part of the path to to_dir.
    return sanitize(to_dir, strict).joinpath(relpath)


def transpaths(to_dir: str | pathlib.Path,
               paths: Iterable[str | pathlib.Path],
               strict: bool = False):
    """ Return all paths that would be produced by moving all paths in
    `paths` from their longest common sub-path to `to_dir` (but do not
    actually move the paths on the file system). This function does not
    require that any of the given paths exist, unless `strict` is True.

    Parameters
    ----------
    to_dir: str | pathlib.Path
        Directory to which to move every path in `path`.
    paths: Iterable[str | pathlib.Path]
        Paths to move; can be files or directories. A common sub-path
        must exist among all of these paths.
    strict: bool = False
        Require that all paths exist and contain no symbolic link loops.

    Returns
    -------
    tuple[pathlib.Path, ...]
        Hypothetical paths after moving all paths in `path` to `outdir`.
    """
    if not paths:
        # There are no paths to transplant.
        return tuple()
    # Determine the longest common sub-path of all given paths.
    common_path = os.path.commonpath([sanitize(p, strict) for p in paths])
    # Move each path from that common path to the given directory.
    return tuple(transpath(to_dir, common_path, p, strict) for p in paths)


# Classes each of whose instances corresponds to a path

class HasFilePath(ABC):
    """ Object that corresponds to the path of a file (which may or may
    not actually exist on the file system). """

    @classmethod
    @abstractmethod
    def get_dir_seg_types(cls) -> list[PathSegment]:
        """ Types of the directory segments in the path. """
        return list()

    @classmethod
    @abstractmethod
    def get_file_seg_type(cls) -> PathSegment:
        """ Type of the last segment in the path. """

    @classmethod
    def get_path_seg_types(cls):
        """ Types of the segments in the path. """
        return [*cls.get_dir_seg_types(), cls.get_file_seg_type()]

    @classmethod
    def get_path_fields(cls):
        """ Path fields for the file type. """
        return get_fields_in_seg_types(cls.get_path_seg_types())

    @classmethod
    def get_ext(cls):
        """ File extension. """
        try:
            return cls.get_file_seg_type().exts[0]
        except IndexError:
            raise ValueError(f"{cls} has no file extension") from None

    @classmethod
    def get_auto_path_fields(cls) -> dict[str, Any]:
        """ Names and path fields that have automatic values. """
        return {EXT: cls.get_ext()}

    @classmethod
    def parse_path(cls, file: str | Path, exclude_auto: bool = False):
        """ Parse a file path to determine the field values. """
        top, field_values = parse_top_separate(file, cls.get_path_seg_types())
        if exclude_auto:
            # Taking the union creates a new set, which avoids iterating
            # over field_values itself as it is being modified.
            exclude = field_values.keys() | cls.get_auto_path_fields().keys()
            for field in exclude:
                field_values.pop(field)
        return top, field_values

    @classmethod
    def build_path(cls, path_fields: dict[str, Any]):
        """ Build the file path from the given field values. """
        return buildpar(cls.get_path_seg_types(),
                        {**cls.get_auto_path_fields(), **path_fields})

    def get_path_field_values(self,
                              top: str | Path | None = None,
                              exclude_auto: bool = False,
                              exclude: Iterable[str] = ()):
        """ Path field values as a dict. """
        auto_path_fields = self.get_auto_path_fields()
        # Make exclude a set for fast membership checking, and to ensure
        # it is not an exhaustible generator.
        if not isinstance(exclude, set):
            exclude = set(exclude)
        if exclude_auto:
            exclude.update(auto_path_fields)
        fields = {TOP: pathlib.Path(top)} if top else dict()
        for field in self.get_path_fields():
            if field not in exclude:
                try:
                    fields[field] = getattr(self, field)
                except AttributeError:
                    # If the field is in neither self.path_fields() nor
                    # self.get_auto_path_fields(), then there is a bug
                    # in how the class itself is designed.
                    assert field in auto_path_fields
                    fields[field] = auto_path_fields[field]
        return fields

    def get_path(self, top: str | Path):
        """ Return the file path. """
        return self.build_path(self.get_path_field_values(top))

    def __str__(self):
        return f"{type(self).__name__}({self.get_path_field_values()})"

    def __repr__(self):
        return str(self)


class HasSampleFilePath(HasFilePath, ABC):
    """ Object that has a path with a sample, step, and branches. """

    @classmethod
    def get_dir_seg_types(cls):
        return [*super().get_dir_seg_types(),
                SampSeg,
                StepSeg]

    @classmethod
    @abstractmethod
    def get_step(cls) -> str:
        """ Step of the workflow. """

    @classmethod
    def get_auto_path_fields(cls):
        return {**super().get_auto_path_fields(),
                STEP: cls.get_step()}


class HasRefFilePath(HasSampleFilePath, ABC):
    """ Object that has a path with a reference. """

    @classmethod
    def get_dir_seg_types(cls):
        return [*super().get_dir_seg_types(),
                RefSeg]


class HasRegFilePath(HasRefFilePath, ABC):
    """ Object that has a path with a region. """

    @classmethod
    def get_dir_seg_types(cls):
        return [*super().get_dir_seg_types(),
                RegSeg]
