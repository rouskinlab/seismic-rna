"""
Core -- Path Module
========================================================================
Auth: Matty

Purpose
========================================================================
Several steps in SEISMIC-RNA produce files that other steps use. For
example, the 'align' step writes alignment map (BAM) files, from which
the 'relate' step writes mutation vector files, which both the 'mask'
and 'table' steps use.

Steps that pass files to each other must agree on
- the path to the file, so that the second step can find the file
- the meaning of each part of the path, so that the second step can
  parse information contained in the path

Although these path conventions could be written separately in each
subpackage or module, this strategy is not ideal for several reasons:
- It would risk inconsistencies among the modules, causing bugs.
- Changing the conventions would require modifying every module, which
  would be not only tedious but also risky for the first reason.
- Defining all the conventions in one place would reduce the size of the
  code base, improving readability, maintainability, and distribution.  

This module defines all file path conventions for all other modules.
"""

from __future__ import annotations

from collections import Counter
from functools import cache, cached_property, partial
from itertools import chain, product
from logging import getLogger
import os
import pathlib as pl
import re
from string import ascii_letters, digits, printable
from typing import Any, Iterable, Sequence

# Constants ############################################################

logger = getLogger(__name__)

# Valid/invalid characters in fields

PATH_CHARS = printable
STR_CHARS = ascii_letters + digits + "_.=+-"
STR_CHARS_SET = frozenset(STR_CHARS)
INT_CHARS = digits
PATH_PATTERN = f"([{PATH_CHARS}]+)"
STR_PATTERN = f"([{STR_CHARS}]+)"
INT_PATTERN = f"([{INT_CHARS}]+)"
RE_PATTERNS = {str: STR_PATTERN, int: INT_PATTERN, pl.Path: PATH_PATTERN}

MOD_DEMULT = "demult"
MOD_QC = "qc"
MOD_ALIGN = "align"
MOD_REL = "relate"
MOD_MASK = "mask"
MOD_CLUST = "cluster"
MOD_TABLE = "table"
MOD_STRUCT = "struct"
MOD_GRAPH = "graph"
MOD_TEST = "test"
MODULES = (MOD_DEMULT, MOD_QC, MOD_ALIGN, MOD_REL, MOD_MASK, MOD_CLUST,
           MOD_TABLE, MOD_STRUCT, MOD_GRAPH, MOD_TEST)

STEPS_QC = "input", "trimmed"
STEPS_ALIGN = ("align-0_refs", "align-1_trim", "align-2_align",
               "align-3_dedup", "align-4_sort", "align-5_split")
STEPS_VECT = "vector-0_bams",
STEPS = STEPS_QC + STEPS_ALIGN + STEPS_VECT

CLUST_PROP_RUN_TABLE = "props"
CLUST_MUS_RUN_TABLE = "mus"
CLUST_RESP_RUN_TABLE = "resps"
CLUST_COUNT_RUN_TABLE = "counts"
CLUST_TABLES = (CLUST_PROP_RUN_TABLE, CLUST_MUS_RUN_TABLE,
                CLUST_RESP_RUN_TABLE, CLUST_COUNT_RUN_TABLE)

RELATE_POS_TAB = "relate-per-pos"
RELATE_READ_TAB = "relate-per-read"
MASKED_POS_TAB = "mask-per-pos"
MASKED_READ_TAB = "mask-per-read"
CLUST_POS_TAB = "clust-per-pos"
CLUST_FREQ_TAB = "clust-freq"
CLUST_READ_TAB = "clust-per-read"
COUNT_TABLES = (RELATE_POS_TAB, RELATE_READ_TAB,
                MASKED_POS_TAB, MASKED_READ_TAB,
                CLUST_POS_TAB, CLUST_FREQ_TAB, CLUST_READ_TAB)

# File extensions

TXT_EXT = ".txt"
CSV_EXT = ".csv"
CSVZIP_EXT = ".csv.gz"
CSV_EXTS = CSV_EXT, CSVZIP_EXT
PARQ_EXTS = ".parquet", ".parq"
JSON_EXT = ".json"
FASTA_EXTS = ".fasta", ".fna", ".fa"
BOWTIE2_INDEX_EXTS = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                      ".rev.1.bt2", ".rev.2.bt2")
FQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz",
           "_001.fastq", "_001.fq", "_001.fastq.gz", "_001.fq.gz")
FQ_PAIRED_EXTS_TEMPLATES = ("_R{}{}", "_mate{}{}", "_{}_sequence{}")
FQ1_EXTS = tuple(template.format(1, ext) for template, ext in
                 product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
FQ2_EXTS = tuple(template.format(2, ext) for template, ext in
                 product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
SAM_EXT = ".sam"
BAM_EXT = ".bam"
XAM_EXTS = SAM_EXT, BAM_EXT
BAI_EXT = f"{BAM_EXT}.bai"
CT_EXT = ".ct"
DOT_EXT = ".dot"
DBN_EXT = ".dbn"
DOT_EXTS = DOT_EXT, DBN_EXT
DMS_EXT = ".dms"
HTML_EXT = ".html"
PDF_EXT = ".pdf"
PNG_EXT = ".png"
GRAPH_EXTS = CSV_EXT, HTML_EXT, PDF_EXT, PNG_EXT
IMAGE_EXTS = PDF_EXT, PNG_EXT


# Path Exceptions ######################################################

class PathError(Exception):
    """ Any error involving a path """


class PathTypeError(PathError, TypeError):
    """ Use of the wrong type of path or segment """


class PathValueError(PathError, ValueError):
    """ Invalid value of a path segment field """


# Path Functions #######################################################

def fill_whitespace(path: str | Path, fill: str = "_"):
    """ Replace all whitespace in `path` with `fill`. """
    return type(path)(fill.join(str(path).split()))


def sanitize(path: str | pl.Path):
    return pl.Path(os.path.realpath(os.path.normpath(os.path.abspath(path))))


# Path Fields ##########################################################

# Field validation functions

def validate_str(txt: str):
    if not isinstance(txt, str):
        raise PathTypeError(f"Expected 'str', got '{type(txt).__name__}'")
    if not txt:
        raise PathValueError(f"Text is empty")
    if illegal := "".join(sorted(set(txt) - STR_CHARS_SET)):
        raise PathValueError(f"Text '{txt}' has illegal characters: {illegal}")


def validate_top(top: pl.Path):
    if not isinstance(top, pl.Path):
        raise PathTypeError(f"Expected 'Path', got '{type(top).__name__}'")
    if not top.parent.is_dir():
        raise PathValueError(f"Not a directory: {top.parent}")
    if top.is_file():
        raise PathValueError(f"File exists: {top}")


def validate_int(num: int):
    if num < 0:
        raise PathValueError(f"Expected integer ≥ 0, got {num}")


VALIDATE = {int: validate_int,
            str: validate_str,
            pl.Path: validate_top}


# Field class

class Field(object):
    def __init__(self,
                 dtype: type[str | int | pl.Path],
                 options: Iterable = (),
                 is_ext: bool = False):
        self.dtype = dtype
        self.options = tuple(options)
        if not all(isinstance(option, self.dtype) for option in self.options):
            raise PathTypeError("All options of a field must be of its type")
        self.is_ext = is_ext
        if self.is_ext:
            if self.dtype is not str:
                raise PathTypeError("Extension field must be type 'str', "
                                    f"but got type '{self.dtype.__name__}'")
            if not self.options:
                raise PathValueError("Extension field must have options")

    def validate(self, val: Any):
        if not isinstance(val, self.dtype):
            raise PathTypeError(f"Expected a '{self.dtype.__name__}', but got "
                                f"{repr(val)} of type '{type(val).__name__}'")
        if self.options and val not in self.options:
            raise PathValueError(
                f"Invalid option {repr(val)}; expected one of {self.options}")
        VALIDATE[self.dtype](val)

    def build(self, val: Any):
        """ Validate a value and return it as a string. """
        self.validate(val)
        return str(val)

    def parse(self, text: str) -> Any:
        """ Parse a value from a string, validate it, and return it. """
        try:
            val = self.dtype(text)
        except Exception as error:
            raise PathValueError(f"Failed to interpret '{text}' as type "
                                 f"'{self.dtype.__name__}': {error}")
        self.validate(val)
        return val

    def __str__(self):
        return f"Path Field '{self.dtype.__name__}'"


# Fields
TopField = Field(pl.Path)
NameField = Field(str)
ModField = Field(str, MODULES)
StepField = Field(str, STEPS)
IntField = Field(int)
CountTabField = Field(str, COUNT_TABLES)
ClustTabField = Field(str, CLUST_TABLES)

# File extensions
TextExt = Field(str, [TXT_EXT], is_ext=True)
ReportExt = Field(str, [JSON_EXT], is_ext=True)
RelVecBatExt = Field(str, PARQ_EXTS, is_ext=True)
MaskRepExt = Field(str, [JSON_EXT], is_ext=True)
MaskBatExt = Field(str, CSV_EXTS, is_ext=True)
ClustTabExt = Field(str, CSV_EXTS, is_ext=True)
ClustCountExt = Field(str, CSV_EXTS, is_ext=True)
ClustBatExt = Field(str, CSV_EXTS, is_ext=True)
MutTabExt = Field(str, CSV_EXTS, is_ext=True)
FastaExt = Field(str, FASTA_EXTS, is_ext=True)
FastaIndexExt = Field(str, BOWTIE2_INDEX_EXTS, is_ext=True)
FastqExt = Field(str, FQ_EXTS, is_ext=True)
Fastq1Ext = Field(str, FQ1_EXTS, is_ext=True)
Fastq2Ext = Field(str, FQ2_EXTS, is_ext=True)
XamExt = Field(str, XAM_EXTS, is_ext=True)
BamIndexExt = Field(str, [BAI_EXT], is_ext=True)
ConnectTableExt = Field(str, [CT_EXT], is_ext=True)
DotBracketExt = Field(str, DOT_EXTS, is_ext=True)
DmsReactsExt = Field(str, [DMS_EXT], is_ext=True)
GraphExt = Field(str, GRAPH_EXTS, is_ext=True)


# Path Segments ########################################################

# Segment class

class Segment(object):
    def __init__(self, segment_name: str,
                 field_types: dict[str, Field], *,
                 order: int = 0,
                 frmt: str | None = None):
        self.name = segment_name
        self.field_types = field_types
        if not self.field_types:
            raise PathValueError(f"Segment got no fields")
        # Verify that a field has the key EXT if and only if it is an
        # extension and is the last field in the segment.
        for i, (name, field) in enumerate(self.field_types.items(), start=1):
            if name == EXT:
                if not field.is_ext:
                    raise PathValueError(f"Field '{EXT}' is not an extension")
                if i != len(self.field_types):
                    raise PathValueError(
                        f"Extension of {self} is not the last field")
                if order != 0:
                    raise ValueError("Segments with extensions must have order "
                                     f"= 0, but {self.name} has order {order}")
            elif field.is_ext:
                raise PathValueError(f"{self} extension has name '{name}'")
        if order <= 0 and not any(ft in self.field_types for ft in [EXT, TOP]):
            raise ValueError("Segments without extensions must have order > 0, "
                             f"but {self.name} has order {order}")
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
        patterns = {name: "" if field.is_ext else RE_PATTERNS[field.dtype]
                    for name, field in self.field_types.items()}
        self.ptrn = re.compile(self.frmt.format(**patterns))

    @property
    def ext_type(self):
        """ Type of the segment's file extension, or None if it has no
        file extension. """
        return self.field_types.get(EXT)

    @cached_property
    def exts(self) -> tuple[str, ...]:
        """ Valid file extensions of the segment. """
        if self.ext_type is None:
            return tuple()
        if not self.ext_type.options:
            raise ValueError(f"{self} extension {self.ext_type} has no options")
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

    def build(self, **vals: Any):
        # Verify that a value is given for every field, with no extras.
        if (v := sorted(vals.keys())) != (f := sorted(self.field_types.keys())):
            raise PathValueError(f"{self} expected fields {f}, but got {v}")
        # Validate the value passed to every field.
        fields = {name: field.build(vals[name])
                  for name, field in self.field_types.items()}
        # Return the formatted segment.
        segment = self.frmt.format(**fields)
        return segment

    def parse(self, text: str):
        ext = None
        if self.ext_type is not None:
            # If the segment has a file extension, then determine the
            # longest valid file extension that matches the text.
            if (ext := self.match_longest_ext(text)) is None:
                raise PathValueError(f"Segment '{text}' is missing a file "
                                     f"extension; expected one of {self.exts}")
            # Remove the file extension from the end of the text.
            text = text[: -len(ext)]
        # Try to parse the text (with the extension, if any, removed).
        if not (match := self.ptrn.match(text)):
            raise PathValueError(f"Could not parse fields in text '{text}' "
                                 f"using pattern '{self.ptrn}'")
        vals = list(match.groups())
        # If there is an extension field, add its value back to the end
        # of the parsed values.
        if ext is not None:
            vals.append(ext)
        # Return a dict of the names of the fields in the segment and
        # their parsed values.
        fields = {name: field.parse(group) for (name, field), group
                  in zip(self.field_types.items(), vals, strict=True)}
        return fields

    def __str__(self):
        return f"Path Segment '{self.name}'"


# Field names
TOP = "top"
STEP = "step"
MOD = "module"
SAMP = "sample"
REF = "ref"
SECT = "sect"
FOLD_SECT = "fold_sect"
BATCH = "batch"
TABLE = "table"
NCLUST = "k"
RUN = "run"
STRUCT = "struct"
REACTS = "reacts"
GRAPH = "graph"
EXT = "ext"

# Directory segments
TopSeg = Segment("top-dir", {TOP: TopField}, order=-1)
ModSeg = Segment("module-dir", {MOD: ModField}, order=60)
StepSeg = Segment("step-dir", {STEP: StepField}, order=50)
SampSeg = Segment("sample-dir", {SAMP: NameField}, order=40)
RefSeg = Segment("ref-dir", {REF: NameField}, order=30)
SectSeg = Segment("section-dir", {SECT: NameField}, order=20)
FoldSectSeg = Segment("fold-section-dir", {FOLD_SECT: NameField}, order=10)

# File segments
# FASTA
FastaSeg = Segment("fasta", {REF: NameField, EXT: FastaExt})
FastaIndexSeg = Segment("fasta-index", {REF: NameField, EXT: FastaIndexExt})
# FASTQ
FastqSeg = Segment("fastq", {SAMP: NameField, EXT: FastqExt})
Fastq1Seg = Segment("fastq1", {SAMP: NameField, EXT: Fastq1Ext})
Fastq2Seg = Segment("fastq2", {SAMP: NameField, EXT: Fastq2Ext})
# Demultiplexed FASTQ
DmFastqSeg = Segment("dm-fastq", {REF: NameField, EXT: FastqExt})
DmFastq1Seg = Segment("dm-fastq1", {REF: NameField, EXT: Fastq1Ext})
DmFastq2Seg = Segment("dm-fastq2", {REF: NameField, EXT: Fastq2Ext})
# Alignment
XamSeg = Segment("xam", {REF: NameField, EXT: XamExt})
BamIndexSeg = Segment("bai", {REF: NameField, EXT: BamIndexExt})
AlignRepSeg = Segment("align-rep", {EXT: ReportExt}, frmt="report-align{ext}")
# Relation Vectors
RelateBatSeg = Segment("rel-bat", {BATCH: IntField, EXT: RelVecBatExt},
                       frmt="batch-relate-{batch}{ext}")
RelateRepSeg = Segment("rel-rep", {EXT: ReportExt}, frmt="report-relate{ext}")
# Masking
MaskBatSeg = Segment("mask-bat", {BATCH: IntField, EXT: MaskBatExt},
                     frmt="batch-mask-{batch}{ext}")
MaskRepSeg = Segment("mask-rep", {EXT: ReportExt}, frmt="report-mask{ext}")
# Clustering
ClustTabSeg = Segment("clust-tab", {TABLE: ClustTabField,
                                    NCLUST: IntField,
                                    RUN: IntField,
                                    EXT: ClustTabExt},
                      frmt="{table}-k{k}-r{run}{ext}")
ClustCountSeg = Segment("clust-count", {EXT: ClustCountExt}, frmt="counts{ext}")
ClustBatSeg = Segment("clust-bat", {BATCH: IntField, EXT: ClustBatExt},
                      frmt="batch-cluster-{batch}{ext}")
ClustRepSeg = Segment("clust-rep", {EXT: ReportExt},
                      frmt="report-cluster{ext}")
# Mutation Tables
MutTabSeg = Segment("mut-tab", {TABLE: CountTabField, EXT: MutTabExt})
# RNA Structure Formats
ConnectTableSeg = Segment("rna-ct", {STRUCT: NameField, EXT: ConnectTableExt})
DotBracketSeg = Segment("rna-dot", {STRUCT: NameField, EXT: DotBracketExt})
DmsReactsSeg = Segment("dms-reacts", {REACTS: NameField, EXT: DmsReactsExt})
VarnaColorSeg = Segment("varna-color", {REACTS: NameField, EXT: TextExt},
                        frmt="{reacts}_varna-color{ext}")
# Graphs
GraphSeg = Segment("graph", {GRAPH: NameField, EXT: GraphExt})


class Path(object):
    def __init__(self, *seg_types: Segment):
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
            raise ValueError(f"Got duplicate order values in {self.seg_types}")

    def build(self, **fields: Any):
        """ Return a `pathlib.Path` instance by assembling the given
        `fields` into a full path. """
        fstr = str(fields)
        # Build the new path one segment at a time.
        segments = list()
        for seg_type in self.seg_types:
            # For each type of segment in the path, try to get the names
            # and values of all fields of the segment.
            try:
                seg_fields = {name: fields.pop(name)
                              for name in seg_type.field_types}
            except KeyError as error:
                raise PathValueError(f"Missing field for {seg_type}: {error}")
            # Generate a string representation of the segment using the
            # values of its fields, and add it to the growing path.
            segments.append(seg_type.build(**seg_fields))
        # Check whether any fields were given but not used by the path.
        if fields:
            exp = [ft for seg in self.seg_types for ft in seg.field_types]
            segs = [str(seg) for seg in self.seg_types]
            raise PathValueError(f"Unexpected fields: {fields}; expected "
                                 f"fields {exp} for segment types {segs}")
        # Assemble the segment strings into a path, and return it.
        path = pl.Path(*segments)
        logger.debug(f"Built path: {fstr} + {tuple(map(str, self.seg_types))} "
                     f"-> {path}")
        return path

    def parse(self, path: str | pl.Path):
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
        logger.debug(f"Parsed path: {path}, {tuple(map(str, self.seg_types))} "
                     f"-> {fields}")
        return fields


@cache
def create_path_type(*segment_types: Segment):
    """ Create and cache a Path instance from the segment types. """
    return Path(*segment_types)


def build(*segment_types: Segment, **field_values: Any):
    """ Return a `pathlib.Path` from the given segment types and
    field values. """
    return create_path_type(*segment_types).build(**field_values)


def builddir(*segment_types: Segment, **field_values: Any):
    """ Build the path and create it on the file system as a directory
    if it does not already exist. """
    path = build(*segment_types, **field_values)
    path.mkdir(parents=True, exist_ok=True)
    return path


def buildpar(*segment_types: Segment, **field_values: Any):
    """ Build a path and create its parent directory if it does not
    already exist. """
    path = build(*segment_types, **field_values)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def parse(path: str | pl.Path, /, *segment_types: Segment):
    """ Return the fields of a path given as a `str` based on the
    segment types. """
    return create_path_type(*segment_types).parse(path)


def find_files(path: pl.Path, segments: Sequence[Segment]):
    """ Yield all files that match a given sequence of path segments.
    The behavior depends on what `path` is:

    - If it is a file, then yield `path` if it matches the segments.
      Otherwise, yield nothing.
    - If it is a directory, then search it recursively and yield every
      matching file in the directory and its subdirectories.
    - If it does not exist, then raise `FileNotFoundError` via calling
      `path.iterdir()`.
    """
    if path.is_file():
        # Check if the given path is a file.
        try:
            # Determine if the file matches the segments.
            parse(path, *segments)
        except PathError:
            # If not, skip it.
            logger.debug(f"File {path} does not match {segments}")
        else:
            # If so, return it.
            logger.debug(f"File {path} matches {segments}")
            yield path
    else:
        # Otherwise, assume it is a directory and search it for reports.
        logger.debug(f"Searching {path} for files matching {segments}")
        yield from chain(*map(partial(find_files, segments=segments),
                              path.iterdir()))


def find_files_chain(paths: Iterable[pl.Path], segments: Sequence[Segment]):
    """ Yield from `find_files` called on every path in `paths`. """
    for path, count in Counter(paths).items():
        if count > 1:
            logger.warning(f"Path {path} was given {count} times")
        try:
            yield from find_files(path, segments)
        except FileNotFoundError:
            logger.error(f"Path does not exist: {path}")
        except Exception as error:
            logger.error(f"Failed search for {path}: {error}")
