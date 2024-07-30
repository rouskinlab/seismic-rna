"""
Core Command Line Interface
===========================
Auth: Yves, Matty

Define all command line interface (CLI) options and their defaults.
"""

import logging
import math
import os
import pathlib
from datetime import datetime
from typing import Iterable

from click import Argument, Choice, Option, Parameter, Path

from ..io import DEFAULT_BROTLI_LEVEL
from ..seq import DNA

# System information
CWD = os.getcwd()
if (NUM_CPUS := os.cpu_count()) is None:
    logging.warning("Failed to determine CPU count: defaulting to 1")
    NUM_CPUS = 1

DEFAULT_MIN_PHRED = 25

BOWTIE2_ORIENT_FR = "fr"
BOWTIE2_ORIENT_RF = "rf"
BOWTIE2_ORIENT_FF = "ff"
BOWTIE2_ORIENT = BOWTIE2_ORIENT_FR, BOWTIE2_ORIENT_RF, BOWTIE2_ORIENT_FF

ADAPTER_SEQ_ILLUMINA_3P = "AGATCGGAAGAGC"
ADAPTER_SEQ_ILLUMINA_5P = "GCTCTTCCGATCT"

NO_GROUP = "c"
GROUP_BY_K = "k"
GROUP_ALL = "a"
GROUP_CLUST_OPTIONS = NO_GROUP, GROUP_BY_K, GROUP_ALL

KEY_NRMSD = "nrmsd"
KEY_PEARSON = "pcc"
KEY_SPEARMAN = "scc"
KEY_DETERM = "r2"
METRIC_KEYS = [KEY_NRMSD,
               KEY_PEARSON,
               KEY_SPEARMAN,
               KEY_DETERM]

# Configuration options

opt_config = Option(
    ("--config", "-g"),
    type=Path(exists=True, dir_okay=False),
    help="Specify parameters in this configuration file"
)

# Input files

arg_fasta = Argument(
    ("fasta",),
    type=Path(exists=True, dir_okay=False),
    nargs=1,
    required=True
)

arg_input_path = Argument(
    ("input-path",),
    type=Path(exists=True),
    nargs=-1
)

# Input/output options

opt_out_dir = Option(
    ("--out-dir", "-o"),
    type=Path(file_okay=False),
    default=os.path.join(".", "out"),
    help="Write all output files to this directory"
)

opt_tmp_pfx = Option(
    ("--tmp-pfx", "-t"),
    type=Path(file_okay=False),
    default=os.path.join(".", "tmp-"),
    help="Write all temporary files to a directory with this prefix"
)

opt_keep_tmp = Option(
    ("--keep-tmp/--erase-tmp",),
    type=bool,
    default=False,
    help="Keep temporary files after finishing"
)

# Resource usage options
opt_parallel = Option(
    ("--parallel/--serial",),
    type=bool,
    default=True,
    help="Run tasks in parallel or in series"
)

opt_max_procs = Option(
    ("--max-procs",),
    type=int,
    default=NUM_CPUS,
    help="Run up to this many processes simultaneously"
)

# Experiment and analysis setup options

opt_force = Option(
    ("--force/--no-force",),
    type=bool,
    default=False,
    help="Force all tasks to run, overwriting any existing output files"
)

# Sequencing read (FASTQ) files
opt_fastqz = Option(
    ("--fastqz", "-z"),
    type=Path(exists=True),
    multiple=True,
    default=(),
    help="FASTQ file(s) of single-end reads"
)

opt_fastqy = Option(
    ("--fastqy", "-y"),
    type=Path(exists=True),
    multiple=True,
    default=(),
    help="FASTQ file(s) of paired-end reads with mates 1 and 2 interleaved"
)

opt_fastqx = Option(
    ("--fastqx", "-x"),
    type=Path(exists=True),
    multiple=True,
    default=(),
    help="FASTQ files of paired-end reads with mates 1 and 2 in separate files"
)

# Sequencing read (FASTQ/XAM) options

opt_phred_enc = Option(
    ("--phred-enc",),
    type=int,
    default=33,
    help="Specify the Phred score encoding of FASTQ and SAM/BAM/CRAM files"
)

opt_min_phred = Option(
    ("--min-phred",),
    type=int,
    default=DEFAULT_MIN_PHRED,
    help="Mark base calls with Phred scores lower than this threshold as ambiguous"
)

opt_fastqc = Option(
    ("--fastqc/--no-fastqc",),
    type=bool,
    default=True,
    help="Run FastQC on the initial and trimmed FASTQ files"
)

opt_qc_extract = Option(
    ("--qc-extract/--qc-no-extract",),
    type=bool,
    default=False,
    help="Unzip FastQC report files"
)

# Demultiplexing options

opt_demultiplex = Option(
    ("--demult-on/--demult-off",),
    type=bool,
    default=False,
    help="Enable demultiplexing"
)

opt_parallel_demultiplexing = Option(
    ("--parallel-demultiplexing",),
    type=bool,
    default=False,
    help="Whether to run demultiplexing at maximum speed by submitting multithreaded "
         "grep functions")

opt_clipped_demultiplexing = Option(
    ("--clipped",),
    type=int,
    default=0,
    help="Designates the amount of clipped patterns to search for in the sample, will raise compution time")

opt_mismatch_tolerence = Option(
    ("--mismatch-tolerence",),
    type=int,
    default=0,
    help="Designates the allowable amount of mismatches allowed in a string and still be considered a valid pattern "
         "find. will increase non-parallel computation at a factorial rate. use caution going above 2 mismatches. "
         "does not apply to clipped sequences.")

opt_index_tolerence = Option(
    ("--index-tolerance",),
    type=int,
    default=0,
    help="Designates the allowable amount of distance you allow the pattern to be found in a read from the reference "
         "index")

opt_barcode_start = Option(
    ("--barcode-start",),
    type=int,
    default=0,
    help="Index of start of barcode")

opt_barcode_end = Option(
    ("--barcode-end",),
    type=int,
    default=0,
    help="Length of barcode")

opt_barcode_length = Option(
    ("--barcode-length",),
    type=int,
    default=0,
    help="Length of barcode")

opt_demulti_overwrite = Option(
    ("--demulti-overwrite",),
    type=bool,
    default=False,
    help="Desiginates whether to overwrite the grepped fastq. should only be used if changing setting on the same "
         "sample")

# Demultiplexed sequencing read (FASTQ) directories

opt_dmfastqz = Option(
    ("--dmfastqz", "-Z"),
    type=Path(exists=True, file_okay=True),
    multiple=True,
    default=(),
    help="Demultiplexed FASTQ files of single-end reads"
)

opt_dmfastqy = Option(
    ("--dmfastqy", "-Y"),
    type=Path(exists=True, file_okay=True),
    multiple=True,
    default=(),
    help="Demultiplexed FASTQ files of paired-end reads interleaved in one file"
)

opt_dmfastqx = Option(
    ("--dmfastqx", "-X"),
    type=Path(exists=True, file_okay=True),
    multiple=True,
    default=(),
    help="Demultiplexed FASTQ files of mate 1 and mate 2 reads"
)

# Adapter trimming options with Cutadapt
opt_cutadapt = Option(
    ("--cut/--no-cut",),
    type=bool,
    default=True,
    help="Use Cutadapt to trim reads before alignment"
)

opt_cut_q1 = Option(
    ("--cut-q1",),
    type=int,
    default=DEFAULT_MIN_PHRED,
    help="Trim base calls below this Phred score from read 1"
)

opt_cut_q2 = Option(
    ("--cut-q2",),
    type=int,
    default=DEFAULT_MIN_PHRED,
    help="Trim base calls below this Phred score from read 2"
)

opt_cut_g1 = Option(
    ("--cut-g1",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_5P,),
    help="Trim this 5' adapter from read 1"
)

opt_cut_a1 = Option(
    ("--cut-a1",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_3P,),
    help="Trim this 3' adapter from read 1"
)

opt_cut_g2 = Option(
    ("--cut-g2",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_5P,),
    help="Trim this 5' adapter from read 2"
)

opt_cut_a2 = Option(
    ("--cut-a2",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_3P,),
    help="Trim this 3' adapter from read 2"
)

opt_cut_o = Option(
    ("--cut-O",),
    type=int,
    default=6,
    help="Require at least this many bases of an adapter to trim it"
)

opt_cut_e = Option(
    ("--cut-e",),
    type=float,
    default=0.1,
    help="Tolerate at most this fraction of errors in adapter sequences"
)

opt_cut_indels = Option(
    ("--cut-indels/--cut-no-indels",),
    type=bool,
    default=True,
    help="Allow errors in adapter sequences to be insertions and deletions"
)

opt_cut_nextseq = Option(
    ("--cut-nextseq/--cut-no-nextseq",),
    type=bool,
    default=False,
    help="Trim high-quality Gs from the 3' end (for Illumina NextSeq and iSeq)"
)

opt_cut_discard_trimmed = Option(
    ("--cut-discard-trimmed/--cut-keep-trimmed",),
    type=bool,
    default=False,
    help="Discard reads in which an adapters were found"
)

opt_cut_discard_untrimmed = Option(
    ("--cut-discard-untrimmed/--cut-keep-untrimmed",),
    type=bool,
    default=False,
    help="Discard reads in which no adapters were found"
)

opt_cut_m = Option(
    ("--cut-m",),
    type=int,
    default=20,
    help="Discard reads shorter than this length after trimming"
)

# Alignment options with Bowtie2

opt_bt2_local = Option(
    ("--bt2-local/--bt2-end-to-end",),
    type=bool,
    default=True,
    help="Run Bowtie2 in local mode rather than end-to-end mode"
)

opt_bt2_discordant = Option(
    ("--bt2-discordant/--bt2-no-discordant",),
    type=bool,
    default=False,
    help="Output paired-end reads whose mates align discordantly")

opt_bt2_mixed = Option(
    ("--bt2-mixed/--bt2-no-mixed",),
    type=bool,
    default=False,
    help="Attempt to align individual mates of pairs that fail to align"
)

opt_bt2_dovetail = Option(
    ("--bt2-dovetail/--bt2-no-dovetail",),
    type=bool,
    default=False,
    help="Consider dovetailed mate pairs to align concordantly"
)

opt_bt2_contain = Option(
    ("--bt2-contain/--bt2-no-contain",),
    type=bool,
    default=True,
    help="Consider nested mate pairs to align concordantly"
)

opt_bt2_un = Option(
    ("--bt2-un/--bt2-no-un",),
    type=bool,
    default=True,
    help="Output unaligned reads to a FASTQ file"
)

opt_bt2_i = Option(
    ("--bt2-I",),
    type=int,
    default=0,
    help="Discard paired-end alignments shorter than this many bases"
)

opt_bt2_x = Option(
    ("--bt2-X",),
    type=int,
    default=600,
    help="Discard paired-end alignments longer than this many bases"
)

opt_bt2_score_min_e2e = Option(
    ("--bt2-score-min-e2e",),
    type=str,
    default="L,-1,-0.5",
    help="Discard alignments that score below this threshold in end-to-end mode"
)

opt_bt2_score_min_loc = Option(
    ("--bt2-score-min-loc",),
    type=str,
    default="L,1,0.5",
    help="Discard alignments that score below this threshold in local mode"
)

opt_bt2_s = Option(
    ("--bt2-i", "bt2_s"),
    type=str,
    default="L,1,0.1",
    help="Seed Bowtie2 alignments at this interval"
)

opt_bt2_l = Option(
    ("--bt2-L",),
    type=int,
    default=20,
    help="Use this seed length for Bowtie2"
)

opt_bt2_gbar = Option(
    ("--bt2-gbar",),
    type=int,
    default=4,
    help="Do not place gaps within this many bases from the end of a read"
)

opt_bt2_d = Option(
    ("--bt2-D",),
    type=int,
    default=4,
    help="Discard alignments if over this many consecutive seed extensions fail"
)

opt_bt2_r = Option(
    ("--bt2-R",),
    type=int,
    default=2,
    help="Re-seed reads with repetitive seeds up to this many times"
)

opt_bt2_dpad = Option(
    ("--bt2-dpad",),
    type=int,
    default=2,
    help="Pad the alignment matrix with this many bases (to allow gaps)"
)

opt_bt2_orient = Option(
    ("--bt2-orient",),
    type=Choice(BOWTIE2_ORIENT, case_sensitive=False),
    default=BOWTIE2_ORIENT[0],
    help="Require paired mates to have this orientation"
)

opt_min_mapq = Option(
    ("--min-mapq",),
    type=int,
    default=25,
    help="Discard reads with mapping qualities below this threshold"
)

opt_min_reads = Option(
    ("--min-reads", "-N"),
    type=int,
    default=1000,
    help="Discard alignment maps with fewer than this many reads"
)

# Relate

opt_batch_size = Option(
    ("--batch-size",),
    type=int,
    default=2 ** 16,
    help="Limit batches to at most this many reads"
)

opt_ambindel = Option(
    ("--ambindel/--no-ambindel",),
    type=bool,
    default=True,
    help="Mark all ambiguous insertions and deletions"
)

opt_brotli_level = Option(
    ("--brotli-level",),
    type=int,
    default=DEFAULT_BROTLI_LEVEL,
    help="Compress pickle files with this level of Brotli (0 - 11)"
)

opt_overhangs = Option(
    ("--overhangs/--no-overhangs",),
    type=bool,
    default=True,
    help="Retain the overhangs of paired-end mates that dovetail"
)

opt_clip_end5 = Option(
    ("--clip-end5", "-5"),
    type=int,
    default=max(opt_bt2_gbar.default, 3),
    help="Clip this many bases from the 5' end of each read"
)

opt_clip_end3 = Option(
    ("--clip-end3", "-3"),
    type=int,
    default=max(opt_bt2_gbar.default, 6),
    help="Clip this many bases from the 3' end of each read"
)

opt_sep_strands = Option(
    ("--sep-strands/--mix-strands",),
    type=bool,
    default=False,
    help="Separate each alignment map into plus- and minus-strand reads"
)

opt_minus_label = Option(
    ("--minus-label",),
    type=str,
    default="-minus",
    help="With --sep-strands, append this label to each minus-strand reference"
)

opt_f1r2_plus = Option(
    ("--f1r2-plus/--f1r2-minus",),
    type=bool,
    default=False,
    help=("With --sep-strands, consider forward mate 1s and reverse mate 2s "
          "to be plus-stranded")
)

# Pool

opt_pool = Option(
    ("--pool", "-P"),
    type=str,
    default="",
    help="Pooled sample name"
)

# Mask

opt_mask_sections_file = Option(
    ("--mask-sections-file", "-s"),
    type=Path(dir_okay=False),
    help="Mask sections of references from coordinates/primers in a CSV file"
)

opt_mask_coords = Option(
    ("--mask-coords", "-c"),
    type=(str, int, int),
    multiple=True,
    default=(),
    help="Mask a section of a reference given its 5' and 3' end coordinates"
)

opt_mask_primers = Option(
    ("--mask-primers", "-p"),
    type=(str, DNA, DNA),
    multiple=True,
    default=(),
    help="Mask a section of a reference given its forward and reverse primers"
)

opt_primer_gap = Option(
    ("--primer-gap",),
    type=int,
    default=0,
    help="Leave a gap of this many bases between the primer and the section"
)

opt_mask_del = Option(
    ("--mask-del/--keep-del",),
    type=bool,
    default=True,
    help="Mask deletions"
)

opt_mask_ins = Option(
    ("--mask-ins/--keep-ins",),
    type=bool,
    default=True,
    help="Mask insertions"
)

opt_mask_mut = Option(
    ("--mask-mut",),
    type=str,
    multiple=True,
    default=(),
    help="Mask this type of mutation"
)

opt_mask_polya = Option(
    ("--mask-polya",),
    type=int,
    default=5,
    help="Mask stretches of at least this many consecutive A bases (0 disables)"
)

opt_mask_gu = Option(
    ("--mask-gu/--keep-gu",),
    type=bool,
    default=True,
    help="Mask G and U bases"
)

opt_mask_pos = Option(
    ("--mask-pos",),
    type=(str, int),
    multiple=True,
    default=(),
    help="Mask this position in this reference"
)

opt_mask_pos_file = Option(
    ("--mask-pos-file",),
    type=Path(dir_okay=False, exists=True),
    help="Mask positions in references from a file"
)

opt_mask_discontig = Option(
    ("--mask-discontig/--keep-discontig",),
    type=bool,
    default=True,
    help="Mask paired-end reads with discontiguous mates"
)

opt_min_ncov_read = Option(
    ("--min-ncov-read",),
    type=int,
    default=1,
    help="Mask reads with fewer than this many bases covering the section"
)

opt_min_finfo_read = Option(
    ("--min-finfo-read",),
    type=float,
    default=0.95,
    help="Mask reads with less than this fraction of unambiguous base calls"
)

opt_max_fmut_read = Option(
    ("--max-fmut-read",),
    type=float,
    default=1.,
    help="Mask reads with more than this fraction of mutated base calls"
)

opt_min_mut_gap = Option(
    ("--min-mut-gap",),
    type=int,
    default=3,
    help="Mask reads with two mutations separated by fewer than this many bases"
)

opt_min_ninfo_pos = Option(
    ("--min-ninfo-pos",),
    type=int,
    default=1000,
    help="Mask positions with fewer than this many unambiguous base calls"
)

opt_max_fmut_pos = Option(
    ("--max-fmut-pos",),
    type=float,
    default=1.,
    help="Mask positions with more than this fraction of mutated base calls"
)

opt_quick_unbias = Option(
    ("--quick-unbias/--exact-unbias",),
    type=bool,
    default=True,
    help="Correct observer bias using a quick (typically linear time) heuristic"
)

opt_quick_unbias_thresh = Option(
    ("--quick-unbias-thresh",),
    type=float,
    default=0.001,
    help="Treat mutated fractions under this threshold as 0 with --quick-unbias"
)

# Cluster options

opt_cluster = Option(
    ("--cluster/--no-cluster",),
    type=bool,
    default=True,
    help="Cluster reads to find alternative structures"
)

opt_min_clusters = Option(
    ("--min-clusters",),
    type=int,
    default=1,
    help="Start at this many clusters"
)

opt_max_clusters = Option(
    ("--max-clusters", "-k"),
    type=int,
    default=0,
    help="Stop at this many clusters (0 for no limit)"
)

opt_try_all_ks = Option(
    ("--try-all-ks/--stop-best-k",),
    type=bool,
    default=False,
    help="Try all numbers of clusters (Ks), even after finding the best number"
)

opt_write_all_ks = Option(
    ("--write-all-ks/--write-best-k",),
    type=bool,
    default=False,
    help="Write all numbers of clusters (Ks), rather than only the best number"
)

opt_max_pearson_run = Option(
    ("--max-pearson-run",),
    type=float,
    default=0.9,
    help="Remove runs with two clusters more similar than this correlation"
)

opt_min_nrmsd_run = Option(
    ("--min-nrmsd-run",),
    type=float,
    default=0.1,
    help="Remove runs with two clusters different by less than this NRMSD"
)

opt_max_loglike_vs_best = Option(
    ("--max-loglike-vs-best",),
    type=float,
    default=250.,
    help="Remove Ks whose 1st/2nd log likelihood difference exceeds this gap"
)

opt_min_pearson_vs_best = Option(
    ("--min-pearson-vs-best",),
    type=float,
    default=0.975,
    help="Remove Ks where every run has less than this correlation vs. the best"
)

opt_max_nrmsd_vs_best = Option(
    ("--max-nrmsd-vs-best",),
    type=float,
    default=0.05,
    help="Remove Ks where every run has more than this NRMSD vs. the best"
)

opt_em_runs = Option(
    ("--em-runs", "-e"),
    type=int,
    default=12,
    help="Run EM this many times for each number of clusters (K)"
)

opt_min_em_iter = Option(
    ("--min-em-iter",),
    type=int,
    default=10,
    help="Run EM for at least this many iterations (times number of clusters)"
)

opt_max_em_iter = Option(
    ("--max-em-iter",),
    type=int,
    default=500,
    help="Run EM for at most this many iterations (times number of clusters)"
)

opt_em_thresh = Option(
    ("--em-thresh",),
    type=float,
    default=round(1. / math.e, 2),
    help="Stop EM when the log likelihood increases by less than this threshold"
)

# Join options

opt_joined = Option(
    ("--joined", "-J"),
    type=str,
    default="",
    help="Joined section name"
)

opt_join_clusts = Option(
    ("--join-clusts", "-j"),
    type=Path(dir_okay=False, exists=True),
    help="Join clusters from this CSV file"
)

# Table options

opt_table_pos = Option(
    ("--table-pos/--no-table-pos",),
    type=bool,
    default=True,
    help="Make a table counting relationships per position"
)

opt_table_read = Option(
    ("--table-read/--no-table-read",),
    type=bool,
    default=True,
    help="Make a table counting relationships per read"
)

opt_table_clust = Option(
    ("--table-clust/--no-table-clust",),
    type=bool,
    default=True,
    help="Make a table counting reads per cluster (only for cluster reports)"
)

# List options

opt_complement = Option(
    ("--complement/--no-complement",),
    type=bool,
    default=False,
    help="List the complement"
)

# Fold

opt_fold = Option(
    ("--fold/--no-fold",),
    type=bool,
    default=False,
    help="Predict the secondary structure using the RNAstructure Fold program"
)

opt_fold_sections_file = Option(
    ("--fold-sections-file", "-f"),
    type=Path(exists=True, dir_okay=False),
    help="Fold sections of references from coordinates/primers in a CSV file"
)

opt_fold_coords = Option(
    ("--fold-coords",),
    type=(str, int, int),
    multiple=True,
    default=(),
    help="Fold a section of a reference given its 5' and 3' end coordinates"
)

opt_fold_primers = Option(
    ("--fold-primers",),
    type=(str, DNA, DNA),
    multiple=True,
    default=(),
    help="Fold a section of a reference given its forward and reverse primers"
)

opt_quantile = Option(
    ("--quantile", "-q"),
    type=float,
    default=0.,
    help="Normalize and winsorize ratios to this quantile (0.0 disables)"
)

opt_fold_temp = Option(
    ("--fold-temp",),
    type=float,
    default=310.15,
    help="Predict structures at this temperature (Kelvin)"
)

opt_fold_constraint = Option(
    ("--fold-constraint",),
    type=Path(exists=True, dir_okay=False),
    help="Force bases to be paired/unpaired from a file of constraints"
)

opt_fold_md = Option(
    ("--fold-md",),
    type=int,
    default=0,
    help="Limit base pair distances to this number of bases (0 for no limit)"
)

opt_fold_mfe = Option(
    ("--fold-mfe/--fold-sub",),
    type=bool,
    default=False,
    help="Predict only the minimum free energy (MFE) structure"
)

opt_fold_max = Option(
    ("--fold-max",),
    type=int,
    default=20,
    help="Output at most this many structures (overriden by --fold-mfe)"
)

opt_fold_percent = Option(
    ("--fold-percent",),
    type=float,
    default=20.,
    help="Stop outputting structures when the % difference in energy exceeds "
         "this value (overriden by --fold-mfe)"
)

# Graph

opt_comppair = Option(
    ("--comppair/--no-comppair",),
    type=bool,
    default=True,
    help="Compare every pair of table files"
)

opt_compself = Option(
    ("--compself/--no-compself",),
    type=bool,
    default=False,
    help="Compare every table file with itself"
)

opt_cgroup = Option(
    ("--cgroup",),
    type=Choice(GROUP_CLUST_OPTIONS),
    default=GROUP_BY_K,
    help="Put each Cluster in its own file, each K in its own file, "
         "or All clusters in one file"
)

opt_rels = Option(
    ("--rels", "-r"),
    type=str,
    multiple=True,
    default=("m",),
    help="Graph these relationship(s)"
)

opt_use_ratio = Option(
    ("--use-ratio/--use-count",),
    type=bool,
    default=True,
    help="Graph ratios or counts"
)

opt_window = Option(
    ("--window", "-w"),
    type=int,
    default=45,
    help="Use a sliding window of this many bases"
)

opt_winmin = Option(
    ("--winmin", "-n"),
    type=int,
    default=9,
    help="Mask sliding windows with fewer than this number of data"
)

opt_metric = Option(
    ("--metric", "-m"),
    type=Choice(METRIC_KEYS, case_sensitive=False),
    default=KEY_PEARSON,
    help=(f"Metric to compare mutation rates: "
          f"{repr(KEY_NRMSD)} = normalized root-mean-square deviation (NRMSD), "
          f"{repr(KEY_PEARSON)} = Pearson correlation coefficient (r), "
          f"{repr(KEY_SPEARMAN)} = Spearman correlation coefficient (ρ), "
          f"{repr(KEY_DETERM)} = coefficient of determination (R²)")
)

opt_struct_file = Option(
    ("--struct-file",),
    type=Path(exists=True, dir_okay=False),
    multiple=True,
    default=(),
    help="Compare mutational profiles to the structure(s) in this CT file"
)

opt_fold_full = Option(
    ("--fold-full/--fold-table",),
    type=bool,
    default=True,
    help="If no sections are specified, whether to default to the full section "
         "or to the table's section"
)

opt_hist_bins = Option(
    ("--hist-bins",),
    type=int,
    default=10,
    help="Number of bins in each histogram; must be ≥ 1"
)

opt_hist_margin = Option(
    ("--hist-margin",),
    type=float,
    default=0.1,
    help="Autofill margins of at most this width in histograms of ratios"
)

opt_csv = Option(
    ("--csv/--no-csv",),
    type=bool,
    default=True,
    help="Output the data for each graph in a Comma-Separated Values file"
)

opt_html = Option(
    ("--html/--no-html",),
    type=bool,
    default=True,
    help="Output each graph in an interactive HyperText Markup Language file"
)

opt_svg = Option(
    ("--svg/--no-svg",),
    type=bool,
    default=False,
    help="Output each graph in a Scalable Vector Graphics file"
)

opt_pdf = Option(
    ("--pdf/--no-pdf",),
    type=bool,
    default=False,
    help="Output each graph in a Portable Document Format file"
)

opt_png = Option(
    ("--png/--no-png",),
    type=bool,
    default=False,
    help="Output each graph in a Portable Network Graphics file"
)

opt_graph_mprof = Option(
    ("--graph-mprof/--no-graph-mprof",),
    type=bool,
    default=True,
    help="Graph mutational profiles"
)

opt_graph_tmprof = Option(
    ("--graph-tmprof/--no-graph-tmprof",),
    type=bool,
    default=True,
    help="Graph typed mutational profiles"
)

opt_graph_ncov = Option(
    ("--graph-ncov/--no-graph-ncov",),
    type=bool,
    default=True,
    help="Graph coverages per position"
)

opt_graph_mhist = Option(
    ("--graph-mhist/--no-graph-mhist",),
    type=bool,
    default=True,
    help="Graph histograms of mutations per read"
)

opt_graph_giniroll = Option(
    ("--graph-giniroll/--no-graph-giniroll",),
    type=bool,
    default=False,
    help="Graph rolling Gini coefficients"
)

opt_graph_roc = Option(
    ("--graph-roc/--no-graph-roc",),
    type=bool,
    default=True,
    help="Graph receiver operating characteristic curves"
)

opt_graph_aucroll = Option(
    ("--graph-aucroll/--no-graph-aucroll",),
    type=bool,
    default=False,
    help="Graph rolling areas under receiver operating characteristic curves"
)

# CT renumbering

opt_ct_pos_5 = Option(
    ("--ct-pos-5", "-c"),
    type=(Path(exists=True), int),
    multiple=True,
    default=(),
    help="Connectivity table (CT) file or directory of CT files "
         "and the 5' position to assign to each file"
)

opt_inplace = Option(
    ("--inplace/--newfile",),
    type=bool,
    default=False,
    help="Modify files in-place instead of writing new files "
         "(WARNING: you cannot recover the original files afterwards)"
)

# Export

opt_export = Option(
    ("--export/--no-export",),
    type=bool,
    default=False,
    help="Export each sample to SEISMICgraph (https://seismicrna.org)"
)

opt_samples_meta = Option(
    ("--samples-meta", "-S"),
    type=Path(dir_okay=False, exists=True),
    help="Add sample metadata from this CSV file to exported results"
)

opt_refs_meta = Option(
    ("--refs-meta", "-R"),
    type=Path(dir_okay=False, exists=True),
    help="Add reference metadata from this CSV file to exported results"
)

opt_all_pos = Option(
    ("--all-pos/--unmasked-pos",),
    type=bool,
    default=True,
    help="Export all positions (not just unmasked positions)"
)

# Simulation

opt_sim_dir = Option(
    ("--sim-dir", "-o"),
    type=Path(file_okay=False),
    default=os.path.join(".", "sim"),
    help="Write all simulated files to this directory"
)

opt_sample = Option(
    ("--sample", "-s"),
    type=str,
    default="sim-sample",
    help="Give this name to the simulated sample"
)

opt_refs = Option(
    ("--refs", "-R"),
    type=str,
    default="sim-refs",
    help="Give this name to the file of simulated references"
)

opt_ref = Option(
    ("--ref", "-r"),
    type=str,
    default="sim-ref",
    help="Give this name to the simulated reference"
)

opt_reflen = Option(
    ("--reflen", "-N"),
    type=int,
    default=280,
    help="Simulate a reference sequence with this many bases"
)

opt_profile_name = Option(
    ("--profile-name", "-P"),
    type=str,
    default="simulated",
    help="Give the simulated structure and parameters this profile name"
)

opt_ct_file = Option(
    ("--ct-file", "-i"),
    type=Path(exists=True, dir_okay=False),
    multiple=True,
    default=(),
    help="Simulate parameters using the structure(s) in this CT file"
)

opt_pmut_paired = Option(
    ("--pmut-paired", "-p"),
    type=(str, float),
    multiple=True,
    default=(),
    help="Set the mean rate of each kind of mutation for paired bases"
)

opt_pmut_unpaired = Option(
    ("--pmut-unpaired", "-u"),
    type=(str, float),
    multiple=True,
    default=(),
    help="Set the mean rate of each kind of mutation for unpaired bases"
)

opt_vmut_paired = Option(
    ("--vmut-paired", "-v"),
    type=float,
    default=0.001,
    help="Set the relative variance of mutation rates of paired bases"
)

opt_vmut_unpaired = Option(
    ("--vmut-unpaired", "-w"),
    type=float,
    default=0.02,
    help="Set the relative variance of mutation rates of unpaired bases"
)

opt_end3_fmean = Option(
    ("--end3-fmean", "-3"),
    type=float,
    default=0.75,
    help="Set the mean 3' end as a fraction of the section length"
)

opt_insert_fmean = Option(
    ("--insert-fmean", "-l"),
    type=float,
    default=0.5,
    help="Set the mean read length as a fraction of the section length"
)

opt_ends_var = Option(
    ("--ends-var", "-e"),
    type=float,
    default=0.25,
    help="Set the variance of end coordinates as a fraction of its supremum"
)

opt_clust_conc = Option(
    ("--clust-conc", "-c"),
    type=float,
    default=0.,
    help="Set the concentration parameter for simulating cluster proportions"
)

opt_param_dir = Option(
    ("--param-dir", "-d"),
    type=Path(exists=True, file_okay=False),
    multiple=True,
    default=(),
    help="Simulate data using parameter files in this directory"
)

opt_num_reads = Option(
    ("--num-reads", "-n"),
    type=int,
    default=opt_batch_size.default,
    help="Simulate this many reads"
)

opt_read_length = Option(
    ("--read-length",),
    type=int,
    default=151,
    help="Simulate reads with this many base calls"
)

opt_paired_end = Option(
    ("--paired-end/--single-end",),
    type=bool,
    default=True,
    help="Simulate paired-end or single-end reads"
)

opt_reverse_fraction = Option(
    ("--reverse-fraction",),
    type=float,
    default=0.5,
    help="Simulate this fraction of reverse-oriented reads"
)

opt_fq_gzip = Option(
    ("--fq-gzip/--fq-text",),
    type=bool,
    default=True,
    help="Simulate FASTQ files with gzip compression or as plain text"
)

# Logging options
opt_verbose = Option(
    ("--verbose", "-v"),
    count=True,
    help="Log info (-v) or info and debug (-vv) messages on stdout"
)

opt_quiet = Option(
    ("--quiet", "-q"),
    count=True,
    help="Suppress warnings (-q) or warnings and errors (-qq) on stdout"
)

opt_log = Option(
    ("--log",),
    type=Path(exists=False, dir_okay=False),
    default=os.path.join(CWD, "log", datetime.now().strftime(
        "seismic-rna_%Y-%m-%d_%H-%M-%S.log"
    )),
    help="Log all messages (except profiling) to a file"
)

opt_log_color = Option(
    ("--log-color/--log-plain",),
    type=bool,
    default=True,
    help="Log messages with or without color codes on stdout"
)

opt_profile = Option(
    ("--profile",),
    type=Path(exists=False, dir_okay=False),
    default="",
    help="Profile code performance and log results to the given file"
)


def merge_params(*param_lists: list[Parameter],
                 exclude: Iterable[Parameter] = ()):
    """ Merge lists of Click parameters, dropping duplicates. """
    exclude_names = {param.name for param in exclude}
    params = list()
    names = set()
    for param_list in param_lists:
        for param in param_list:
            if param.name not in names and param.name not in exclude_names:
                params.append(param)
                names.add(param.name)
    return params


def optional_path(param: str | pathlib.Path | None):
    if isinstance(param, pathlib.Path):
        return param
    if param:
        return pathlib.Path(param)
    return None

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
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
