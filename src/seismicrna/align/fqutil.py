from functools import cached_property
from itertools import chain
from logging import getLogger
from pathlib import Path

from ..core import path
from ..core.cli import BOWTIE2_ORIENT
from ..core.shell import run_cmd, BOWTIE2_CMD, CUTADAPT_CMD, FASTQC_CMD

logger = getLogger(__name__)

# Bowtie2 parameters
# Consider this example: Ref = ACGT, Read = AG
# Assume that we want to minimize the number of edits needed to convert
# the reference into the read sequence. The smallest number of edits is
# two, specifically these two deletions (/) from the reference: [A/G/]
# which gets a score of (2 * match - 2 * gap_open - 2 * gap_extend).
# But there are two alternative alignment, each with 3 edits:
# [Ag//] and [A//g] (substitutions marked in lowercase). Each gets the
# score (match - substitution - gap_open - 2 * gap_extend).
# In order to favor the simpler alignment with two edits,
# (2 * match - 2 * gap_open - 2 * gap_extend) must be greater than
# (match - substitution - gap_open - 2 * gap_extend); this simplifies to
# (substitution > gap_open - match). Thus, the substitution penalty and
# match bonus must be relatively large, and the gap open penalty small.
# We want to avoid introducing too many gaps, especially to prevent the
# introduction of an insertion and a deletion from scoring better than
# one substitution. Consider this example: Ref = ATAT, Read = ACTT
# The simplest alignment (the smallest number of mutations) is ActT,
# which gets a score of (2 * match - 2 * substitution).
# Another alignment with indels is A{C}T/T, where {C} means a C was
# inserted into the read and the / denotes an A deleted from the read.
# This alignment scores (3 * match - 2 * gap_open - 2 * gap_extend).
# Thus, (2 * match - 2 * substitution) must be greater than
# (3 * match - 2 * gap_open - 2 * gap_extend), which simplifies to
# (2 * gap_open + 2 * gap_extend > match + 2 * substitution).
# There are two easy solutions to these inequalities:
# - Bowtie v2.5 defaults: 6 > 5 - 2 and 2*5 + 2*3 > 2 + 2*6
# - Set every value to 1: 1 > 1 - 1 and 2*1 + 2*1 > 1 + 2*1
MATCH_BONUS = "1"
MISMATCH_PENALTY = "1,1"
N_PENALTY = "0"
REF_GAP_PENALTY = "1,1"
READ_GAP_PENALTY = "1,1"
METRICS_INTERVAL = 60  # Write metrics once every 60 seconds.


class FastqUnit(object):
    """
    Unified interface for the following sets of sequencing reads:

    - One FASTQ file of single-end reads from one sample
    - One FASTQ file of interleaved, paired-end reads from one sample
    - Two FASTQ files of mate 1 and 2 paired-end reads from one sample
    - One FASTQ file of single-end reads originating from one reference
      sequence in one sample
    - One FASTQ file of interleaved, paired-end reads originating from
      one reference sequence in one sample
    - Two FASTQ files of mate 1 and mate 2 paired-end reads originating
      from one reference sequence in one sample
    """

    MAX_PHRED_ENC = 127  # 2^7 - 1

    KEY_SINGLE = "fastqs"
    KEY_INTER = "fastqi"
    KEY_MATED = "fastqp"
    KEY_DSINGLE = "dmfastqs"
    KEY_DINTER = "dmfastqi"
    KEY_DMATED = "dmfastqp"
    KEY_MATE1 = "fastq1"
    KEY_MATE2 = "fastq2"

    BOWTIE2_FLAGS = {KEY_SINGLE: "-U",
                     KEY_INTER: "--interleaved",
                     KEY_MATE1: "-1",
                     KEY_MATE2: "-2"}

    def __init__(self, *,
                 fastqs: Path | None = None,
                 fastqi: Path | None = None,
                 fastq1: Path | None = None,
                 fastq2: Path | None = None,
                 phred_enc: int,
                 one_ref: bool):
        if fastqs:
            if fastqi or fastq1 or fastq2:
                raise TypeError("Got too many FASTQ files")
            self.paths: dict[str, Path] = {self.KEY_SINGLE: fastqs}
            self.paired = False
            self.interleaved = False
        elif fastqi:
            if fastq1 or fastq2:
                raise TypeError("Got too many FASTQ files")
            self.paths: dict[str, Path] = {self.KEY_INTER: fastqi}
            self.paired = True
            self.interleaved = True
        elif fastq1:
            if not fastq2:
                raise TypeError("Got fastq1 but not fastq2")
            self.paths: dict[str, Path] = {self.KEY_MATE1: fastq1,
                                           self.KEY_MATE2: fastq2}
            self.paired = True
            self.interleaved = False
        elif fastq2:
            raise TypeError("Got fastq2 but not fastq1")
        if phred_enc < 0 or phred_enc > self.MAX_PHRED_ENC:
            raise ValueError(f"Invalid Phred encoding: {phred_enc}")
        self.phred_enc = phred_enc
        self.one_ref = one_ref
        self.sample, self.ref, self.exts = self.get_sample_ref_exts()
        logger.debug(f"Instantiated a {self.__class__.__name__} with "
                     + ", ".join(f"{k} = {v} (type '{type(v).__name__}')"
                                 for k, v in self.paths.items())
                     + f", phred_enc = {phred_enc}, one_ref = {one_ref}")

    @property
    def phred_arg(self):
        return f"--phred{self.phred_enc}"

    @property
    def kind(self):
        if self.paired:
            if self.interleaved:
                return "interleaved paired-end FASTQ file"
            return "separate paired-end FASTQ files"
        return "single-end FASTQ file"

    @cached_property
    def parent(self):
        """ Return the parent directory of the FASTQ file(s). """
        parents = [inp.parent for inp in self.paths.values()]
        if not parents:
            raise TypeError("Not parent directory")
        if any(parent != parents[0] for parent in parents[1:]):
            raise ValueError("More than one parent directory")
        return parents[0]

    @cached_property
    def seg_types(self) -> dict[str, tuple[path.Segment, ...]]:
        if self.one_ref:
            seg_types = {self.KEY_SINGLE: (path.SampSeg, path.DmFastqSeg),
                         self.KEY_INTER: (path.SampSeg, path.DmFastqSeg),
                         self.KEY_MATE1: (path.SampSeg, path.DmFastq1Seg),
                         self.KEY_MATE2: (path.SampSeg, path.DmFastq2Seg)}
        else:
            seg_types = {self.KEY_SINGLE: (path.FastqSeg,),
                         self.KEY_INTER: (path.FastqSeg,),
                         self.KEY_MATE1: (path.Fastq1Seg,),
                         self.KEY_MATE2: (path.Fastq2Seg,)}
        return {key: seg_types[key] for key in self.paths}

    def get_sample_ref_exts(self):
        """ Return the sample and reference of the FASTQ file(s). """
        samples: set[str] = set()
        refs: set[str | None] = set()
        exts: dict[str, str] = dict()
        for key, fq in self.paths.items():
            fq_fields = path.parse(fq, *self.seg_types[key])
            samples.add(fq_fields[path.SAMP])
            refs.add(fq_fields.get(path.REF))
            exts[key] = fq_fields[path.EXT]
        if len(samples) > 1:
            raise ValueError(f"Sample names of {self} disagree: "
                             + " ≠ ".join(samples))
        if len(refs) > 1:
            raise ValueError(f"Ref names of {self} disagree: "
                             + " ≠ ".join(map(str, refs)))
        return list(samples)[0], list(refs)[0], exts

    def fields(self, key: str):
        fields = {path.SAMP: self.sample}
        if self.ref is not None:
            fields[path.REF] = self.ref
        fields[path.EXT] = self.exts[key]
        return fields

    def is_compatible_fasta(self, fasta: Path, one_ref_fasta: bool):
        """
        Return whether a given FASTA file is compatible with the FASTQ,
        which means one of the following is true:
        - The FASTA has a set of references (one_ref_fasta=False) and
          the FASTQ has reads from an entire sample.
        - The FASTA has one reference (one_ref_fasta=True) and the FASTQ
          has reads from only the same reference.
        """
        if one_ref_fasta:
            return self.ref == path.parse(fasta, path.FastaSeg)[path.REF]
        return self.ref is None

    @property
    def cutadapt_input_args(self):
        """ Return input file arguments for Cutadapt. """
        return tuple(self.paths.values())

    @property
    def bowtie2_inputs(self):
        """ Return input file arguments for Bowtie2. """
        return tuple(chain(*[(self.BOWTIE2_FLAGS[key], fq)
                             for key, fq in self.paths.items()]))

    def to_new(self, *new_segments: path.Segment, **new_fields):
        """ Return a new FASTQ unit with updated path fields. """
        new_paths = dict()
        for key, self_path in self.paths.items():
            combined_segments = new_segments + self.seg_types[key]
            combined_fields = self.fields(key) | new_fields
            new_paths[key] = path.build(*combined_segments, **combined_fields)
        return self.__class__(**new_paths,
                              phred_enc=self.phred_enc,
                              one_ref=self.one_ref)

    @classmethod
    def _from_files(cls, /, *, phred_enc: int, one_ref: bool,
                    fqs: list[Path], key: str):
        if key not in (cls.KEY_SINGLE, cls.KEY_INTER):
            raise ValueError(f"Invalid key: '{key}'")
        segs = [path.SampSeg, path.DmFastqSeg] if one_ref else [path.FastqSeg]
        for fq in path.find_files_chain(fqs, segs):
            try:
                yield cls(phred_enc=phred_enc, one_ref=one_ref, **{key: fq})
            except Exception as error:
                logger.error(f"Failed to load FASTQ file {fq}: {error}")

    @classmethod
    def _from_mates(cls, /, *, phred_enc: int, one_ref: bool,
                    fqs: list[Path]):
        # Determine the key and segments based on whether the FASTQs are
        # demultiplexed.
        if one_ref:
            seg1s = [path.SampSeg, path.DmFastq1Seg]
            seg2s = [path.SampSeg, path.DmFastq2Seg]
        else:
            seg1s = [path.Fastq1Seg]
            seg2s = [path.Fastq2Seg]
        # List all FASTQ mate 1 and mate 2 files.
        fq1s = path.find_files_chain(fqs, seg1s)
        fq2s = path.find_files_chain(fqs, seg2s)

        # Determine the sample and/or reference name of each file.
        def by_tag(fqs_: list[Path], segs: list[path.Segment]):
            tags: dict[tuple[str, str | None], Path] = dict()
            for fq in fqs_:
                fields = path.parse(fq, *segs)
                tag_ = fields[path.SAMP], fields.get(path.REF)
                if tag_ in tags:
                    logger.warning(f"Duplicate sample and reference: {tag_}")
                else:
                    tags[tag_] = fq
            return tags

        tag1s = by_tag(fq1s, seg1s)
        tag2s = by_tag(fq2s, seg2s)
        # Check for any mates with only one file.
        set1s, set2s = set(tag1s), set(tag2s)
        if miss1 := set2s - set1s:
            logger.error(f"Missing FASTQ mate 1 files: {miss1}")
        if miss2 := set1s - set2s:
            logger.error(f"Missing FASTQ mate 2 files: {miss2}")
        # Yield a FASTQ unit for each pair of mated files.
        for tag in set1s & set2s:
            fq_args = {cls.KEY_MATE1: tag1s[tag], cls.KEY_MATE2: tag2s[tag]}
            try:
                yield cls(phred_enc=phred_enc, one_ref=one_ref, **fq_args)
            except Exception as error:
                logger.error(f"Failed to load FASTQ pair {fq_args}: {error}")

    @classmethod
    def from_paths(cls, /, *, phred_enc: int, **fastq_args: list[Path]):
        """
        Yield a FastqUnit for each FASTQ file (or each pair of mate 1
        and mate 2 FASTQ files) whose paths are given as strings.

        Parameters
        ----------
        phred_enc: int
            ASCII offset for encoding Phred scores
        fastq_args: list[Path]
            FASTQ files, given as lists of paths:
            - fastqs: FASTQ files of single-end reads
            - fastqi: FASTQ files of interleaved paired-end reads
            - fastqp: mated FASTQ files of paired-end reads
            - dmfastqs: demultiplexed FASTQ files of single-end reads
            - dmfastqi: demultiplexed FASTQ files of interleaved paired-end reads
            - dmfastqp: demultiplexed mated FASTQ files of paired-end reads

        Yield
        -----
        FastqUnit
            FastqUnit representing the FASTQ or pair of FASTQ files.
            The order is determined primarily by the order of keyword
            arguments; within each keyword argument, by the order of
            file or directory paths; and for directories, by the order
            in which `os.path.listdir` returns file paths.
        """
        # List all FASTQ files.
        # single-end
        yield from cls._from_files(phred_enc=phred_enc, one_ref=False,
                                   fqs=fastq_args.get(cls.KEY_SINGLE, ()),
                                   key=cls.KEY_SINGLE)
        # interleaved paired-end
        yield from cls._from_files(phred_enc=phred_enc, one_ref=False,
                                   fqs=fastq_args.get(cls.KEY_INTER, ()),
                                   key=cls.KEY_INTER)
        # mated paired-end
        yield from cls._from_mates(phred_enc=phred_enc, one_ref=False,
                                   fqs=fastq_args.get(cls.KEY_MATED, ()))
        # demultiplexed single-end
        yield from cls._from_files(phred_enc=phred_enc, one_ref=True,
                                   fqs=fastq_args.get(cls.KEY_DSINGLE, ()),
                                   key=cls.KEY_SINGLE)
        # demultiplexed interleaved paired-end
        yield from cls._from_files(phred_enc=phred_enc, one_ref=True,
                                   fqs=fastq_args.get(cls.KEY_DINTER, ()),
                                   key=cls.KEY_INTER)
        # demultiplexed mated paired-end
        yield from cls._from_mates(phred_enc=phred_enc, one_ref=True,
                                   fqs=fastq_args.get(cls.KEY_DMATED, ()))

    def __str__(self):
        return f"{self.kind} {' and '.join(map(str, self.paths.values()))}"


def run_fastqc(fq_unit: FastqUnit, out_dir: Path, *,
               extract: bool, n_procs: int):
    """ Run FASTQC on the given FASTQ unit. """
    logger.info(f"Began FASTQC of {fq_unit}")
    # Create FASTQC output directory.
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Created directory: {out_dir}")
    # Create FASTQC command.
    cmd = [FASTQC_CMD,
           "--threads", n_procs,
           "--extract" if extract else "--noextract",
           "--outdir", out_dir]
    # Add input FASTQ files.
    cmd.extend(fq_unit.paths.values())
    # Run FASTQC.
    run_cmd(cmd)
    logger.info(f"Ended FASTQC of {fq_unit}")


def run_cutadapt(fq_inp: FastqUnit,
                 fq_out: FastqUnit, *,
                 n_procs: int,
                 cut_q1: int,
                 cut_q2: int,
                 cut_g1: str,
                 cut_a1: str,
                 cut_g2: str,
                 cut_a2: str,
                 cut_o: int,
                 cut_e: float,
                 cut_indels: bool,
                 cut_nextseq: bool,
                 cut_discard_trimmed: bool,
                 cut_discard_untrimmed: bool,
                 cut_m: int):
    """ Trim adapters and low-quality bases with Cutadapt. """
    logger.info(f"Began trimming {fq_inp}")
    # Cutadapt command
    cmd = [CUTADAPT_CMD, "--cores", n_procs]
    # Quality trimming
    if cut_nextseq:
        cut_qnext = max(cut_q1, cut_q2)
        if cut_q1 != cut_q2:
            logger.warning("NextSeq trimming takes one quality level, but got "
                           f"two ({cut_q1} and {cut_q2}); using {cut_qnext}")
        cmd.extend(["--nextseq-trim", cut_qnext])
    else:
        cmd.extend(["-q", cut_q1, "-Q", cut_q2])
    # Adapter trimming
    adapters = {"g": cut_g1, "a": cut_a1,
                "G": cut_g2, "A": cut_a2}
    for arg, adapter in adapters.items():
        if adapter and (fq_inp.paired or arg.islower()):
            for adapt in adapter:
                cmd.extend([f"-{arg}", adapt])
    cmd.extend(["-O", cut_o])
    cmd.extend(["-e", cut_e])
    cmd.extend(["-m", cut_m])
    if not cut_indels:
        cmd.append("--no-indels")
    if cut_discard_trimmed:
        cmd.append("--discard-trimmed")
    if cut_discard_untrimmed:
        cmd.append("--discard-untrimmed")
    # FASTQ format
    if fq_inp.interleaved:
        cmd.append("--interleaved")
    # Output files
    output_args = list(zip(("-o", "-p"), fq_out.paths.values(), strict=False))
    for flag, value in output_args:
        cmd.extend([flag, value])
    # Input files
    cmd.extend(fq_inp.cutadapt_input_args)
    # Make the output directory.
    fq_out.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Created directory: {fq_out.parent}")
    # Run Cutadapt.
    run_cmd(cmd, check_created=[output for _, output in output_args])
    logger.info(f"Ended trimming {fq_inp}; output {fq_out}")


def run_bowtie2(fq_inp: FastqUnit,
                index_pfx: Path,
                sam_out: Path, *,
                n_procs: int,
                bt2_local: bool,
                bt2_discordant: bool,
                bt2_mixed: bool,
                bt2_dovetail: bool,
                bt2_contain: bool,
                bt2_unal: bool,
                bt2_score_min_e2e: str,
                bt2_score_min_loc: str,
                bt2_i: int,
                bt2_x: int,
                bt2_gbar: int,
                bt2_l: int,
                bt2_s: str,
                bt2_d: int,
                bt2_r: int,
                bt2_dpad: int,
                bt2_orient: str):
    """ Align reads to the reference with Bowtie 2. """
    logger.info(f"Began aligning {fq_inp} to {index_pfx}")
    # Bowtie2 command
    cmd = [BOWTIE2_CMD]
    # Resources
    cmd.extend(["--threads", n_procs])
    # Alignment
    cmd.append("--local" if bt2_local else "--end-to-end")
    cmd.extend(["--gbar", bt2_gbar])
    cmd.extend(["--dpad", bt2_dpad])
    cmd.extend(["-L", bt2_l])
    cmd.extend(["-i", bt2_s])
    cmd.extend(["-D", bt2_d])
    cmd.extend(["-R", bt2_r])
    # Scoring
    cmd.append(fq_inp.phred_arg)
    cmd.append("--ignore-quals")
    cmd.extend(["--ma", MATCH_BONUS if bt2_local else "0"])
    cmd.extend(["--mp", MISMATCH_PENALTY])
    cmd.extend(["--np", N_PENALTY])
    cmd.extend(["--rfg", REF_GAP_PENALTY])
    cmd.extend(["--rdg", READ_GAP_PENALTY])
    # Filtering
    if not bt2_unal:
        cmd.append("--no-unal")
    cmd.extend(["--score-min", (bt2_score_min_loc if bt2_local
                                else bt2_score_min_e2e)])
    cmd.extend(["-I", bt2_i])
    cmd.extend(["-X", bt2_x])
    # Mate pair orientation
    if bt2_orient not in BOWTIE2_ORIENT:
        logger.warning(f"Invalid mate orientation for Bowtie2: '{bt2_orient}'. "
                       f"Setting to '{BOWTIE2_ORIENT[0]}'")
        bt2_orient = BOWTIE2_ORIENT[0]
    cmd.append(f"--{bt2_orient}")
    if not bt2_discordant:
        cmd.append("--no-discordant")
    if not bt2_contain:
        cmd.append("--no-contain")
    if bt2_dovetail:
        cmd.append("--dovetail")
    if not bt2_mixed:
        cmd.append("--no-mixed")
    # Formatting
    cmd.append("--xeq")
    # Metrics
    cmd.extend(["--met-stderr", "--met", METRICS_INTERVAL])
    # Input and output files
    cmd.extend(["-S", sam_out])
    cmd.extend(["-x", index_pfx])
    cmd.extend(fq_inp.bowtie2_inputs)
    # Make the output directory.
    try:
        sam_out.parent.mkdir(parents=True, exist_ok=False)
        logger.debug(f"Created directory: {sam_out.parent}")
    except FileExistsError:
        # The directory should not exist already unless the input FASTQ
        # has been demultiplexed (and thus many FASTQs could get aligned
        # in this directory).
        if not fq_inp.one_ref:
            raise
    # Run alignment.
    run_cmd(cmd, check_created=[sam_out])
    logger.info(f"Ended aligning {fq_inp} and writing to {sam_out}")
