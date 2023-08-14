"""
Core Command Line Interface
===========================
Auth: Yves, Matty

Define all command line interface (CLI) options and their defaults.
"""

from datetime import datetime
import logging
import os

from click import Argument, Choice, Option, Parameter, Path

# System information
CWD = os.getcwd()
if (NUM_CPUS := os.cpu_count()) is None:
    logging.warning("Failed to determine CPU count: defaulting to 1")
    NUM_CPUS = 1

DEFAULT_PHRED_ENC = 33
DEFAULT_MIN_PHRED = 25

BOWTIE2_ORIENT_FR = "fr"
BOWTIE2_ORIENT_RF = "rf"
BOWTIE2_ORIENT_FF = "ff"
BOWTIE2_ORIENT = BOWTIE2_ORIENT_FR, BOWTIE2_ORIENT_RF, BOWTIE2_ORIENT_FF

ADAPTER_SEQ_ILLUMINA_3P = "AGATCGGAAGAGC"

CLUST_INDIV = "indiv"
CLUST_ORDER = "order"
CLUST_UNITE = "unite"
CLUST_ARRANGE_OPTIONS = CLUST_INDIV, CLUST_ORDER, CLUST_UNITE

# Configuration options
opt_config = Option(
    ("--config", "-g"),
    type=Path(exists=True, dir_okay=False),
    help="Configuration file for parameters")

# Reference sequence (FASTA) files
arg_fasta = Argument(
    ("fasta",),
    type=Path(exists=True, dir_okay=False),
    nargs=1,
    required=True)

# Input/output options

arg_input_file = Argument(
    ("input-file",),
    type=Path(exists=True),
    nargs=-1)

opt_out_dir = Option(
    ("--out-dir", "-o"),
    type=Path(file_okay=False),
    default=os.path.join(".", "out"),
    help="Where to output all finished files")

opt_temp_dir = Option(
    ("--temp-dir", "-t"),
    type=Path(file_okay=False),
    default=os.path.join(".", "temp"),
    help="Where to write all temporary files")

opt_save_temp = Option(
    ("--save-temp/--no-save-temp",),
    type=bool,
    default=False,
    help=("Whether to save temporary files when the "
          "program exits"))

# Resource usage options
opt_parallel = Option(
    ("--parallel/--no-parallel",),
    type=bool,
    default=True,
    help="Whether to run multiple tasks in parallel")

opt_max_procs = Option(
    ("--max-procs",),
    type=int,
    default=NUM_CPUS,
    help="Maximum number of simultaneous processes")

# Experiment and analysis setup options

opt_sections_file = Option(
    ("--sections-file", "-s"),
    type=Path(dir_okay=False),
    default="",
    help="CSV file of sections by name, reference, and coordinates/primers")

opt_samples = Option(
    ("--samples", "-S"),
    type=Path(dir_okay=False),
    default="",
    help="Samples file")

opt_rerun = Option(
    ("--rerun/--no-rerun",),
    type=bool,
    default=False,
    help="Whether to regenerate files that already exist")

# Sequencing read (FASTQ) files
opt_fastqs = Option(
    ("--fastqs", "-z"),
    type=Path(exists=True),
    multiple=True,
    default=(),
    help="FASTQ files of single-end reads")

opt_fastqi = Option(
    ("--fastqi", "-y"),
    type=Path(exists=True),
    multiple=True,
    default=(),
    help="FASTQ files of paired-end reads interleaved in 1 file")

opt_fastqp = Option(
    ("--fastqp", "-x"),
    type=Path(exists=True),
    multiple=True,
    default=(),
    help="FASTQ files of paired-end reads separated into 2 files")

# Sequencing read (FASTQ/BAM) options
opt_phred_enc = Option(
    ("--phred-enc",),
    type=int,
    default=DEFAULT_PHRED_ENC,
    help="Phred score encoding in FASTQ/SAM/BAM files")

opt_min_phred = Option(
    ("--min-phred",),
    type=int,
    default=DEFAULT_MIN_PHRED,
    help="Minimum Phred score to use a base call")

opt_fastqc = Option(
    ("--fastqc/--no-fastqc",),
    type=bool,
    default=True,
    help="Whether to check quality of FASTQ files")

opt_qc_extract = Option(
    ("--qc-extract/--qc-no-extract",),
    type=bool,
    default=False,
    help="Whether to unzip FASTQC reports")

# Demultiplexing options

opt_demultiplex = Option(
    ("--demult-on/--demult-off",),
    type=bool,
    default=False,
    help="Whether to run demultiplexing")

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
    help="Designates the allowable amount of mismatches allowed in a string and still be considered a valid pattern find. \
                            will increase non-parallel computation at a factorial rate. use caution going above 2 mismatches. does not apply to clipped sequences.")

opt_index_tolerence = Option(
    ("--index-tolerance",),
    type=int,
    default=0,
    help="Designates the allowable amount of distance you allow the pattern to be found in a read from the reference index")

opt_barcode_start = Option(
    ("--barcode-start",),
    type=int,
    default=0,
    help="index of start of barcode")

opt_barcode_length = Option(
    ("--barcode-length",),
    type=int,
    default=0,
    help="length of barcode")

opt_demulti_overwrite = Option(
    ("--demulti-overwrite",),
    type=bool,
    default=False,
    help="desiginates whether to overwrite the grepped fastq. should only be used if changing setting on the same sample")

# Demultiplexed sequencing read (FASTQ) directories
opt_dmfastqs = Option(
    ("--dmfastqs", "-Z"),
    type=Path(exists=True, file_okay=False),
    multiple=True,
    default=(),
    help="Demultiplexed FASTQ files of single-end reads")

opt_dmfastqi = Option(
    ("--dmfastqi", "-Y"),
    type=Path(exists=True, file_okay=False),
    multiple=True,
    default=(),
    help="Demultiplexed FASTQ files of paired-end reads interleaved in one file")

opt_dmfastqp = Option(
    ("--dmfastqp", "-X"),
    type=Path(exists=True, file_okay=False),
    multiple=True,
    default=(),
    help="Demultiplexed FASTQ files of mate 1 and mate 2 reads")

# Adapter trimming options with Cutadapt
opt_cutadapt = Option(
    ("--cut/--no-cut",),
    type=bool,
    default=True,
    help="Whether to trim reads with Cutadapt before alignment")

opt_cut_q1 = Option(
    ("--cut-q1",),
    type=int,
    default=DEFAULT_MIN_PHRED,
    help="Phred score for read 1 quality trimming")

opt_cut_q2 = Option(
    ("--cut-q2",),
    type=int,
    default=DEFAULT_MIN_PHRED,
    help="Phred score for read 2 quality trimming")

opt_cut_g1 = Option(
    ("--cut-g1",),
    type=str,
    multiple=True,
    default=(),
    help="5' adapter for read 1")

opt_cut_a1 = Option(
    ("--cut-a1",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_3P,),
    help="3' adapter for read 1")

opt_cut_g2 = Option(
    ("--cut-g2",),
    type=str,
    multiple=True,
    default=(),
    help="5' adapter for read 2")

opt_cut_a2 = Option(
    ("--cut-a2",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_3P,),
    help="3' adapter for read 2")

opt_cut_o = Option(
    ("--cut-O",),
    type=int,
    default=6,
    help="Minimum overlap of read and adapter")

opt_cut_e = Option(
    ("--cut-e",),
    type=float,
    default=0.1,
    help="Error tolerance for adapters")

opt_cut_indels = Option(
    ("--cut-indels/--cut-no-indels",),
    type=bool,
    default=True,
    help="Whether to allow indels in adapters")

opt_cut_nextseq = Option(
    ("--cut-nextseq/--cut-no-nextseq",),
    type=bool,
    default=False,
    help="Whether to trim high-quality Gs from 3' end")

opt_cut_discard_trimmed = Option(
    ("--cut-discard-trimmed/--cut-keep-trimmed",),
    type=bool,
    default=False,
    help="Whether to discard reads in which an adapter was found")

opt_cut_discard_untrimmed = Option(
    ("--cut-discard-untrimmed/--cut-keep-untrimmed",),
    type=bool,
    default=False,
    help="Whether to discard reads in which no adapter was found")

opt_cut_m = Option(
    ("--cut-m",),
    type=int,
    default=20,
    help="Discard reads shorter than this length after trimming")

# Alignment options with Bowtie2
opt_bt2_local = Option(
    ("--bt2-local/--bt2-end-to-end",),
    type=bool,
    default=True,
    help="Whether to perform local or end-to-end alignment. "
         "Use local with fragmented samples, on which "
         "end-to-end yields false positives at read ends. "
         "Use end-to-end with amplicon-based samples, on "
         "which local yields false negatives at read ends.")

opt_bt2_discordant = Option(
    ("--bt2-discordant/--bt2-no-discordant",),
    type=bool,
    default=False,
    help="Whether to output discordant alignments")

opt_bt2_mixed = Option(
    ("--bt2-mixed/--bt2-no-mixed",),
    type=bool,
    default=False,
    help="Whether to align individual mates of unaligned pairs")

opt_bt2_dovetail = Option(
    ("--bt2-dovetail/--bt2-no-dovetail",),
    type=bool,
    default=False,
    help="Whether to treat dovetailed mate pairs as concordant")

opt_bt2_contain = Option(
    ("--bt2-contain/--bt2-no-contain",),
    type=bool,
    default=True,
    help="Whether to treat nested mate pairs as concordant")

opt_bt2_unal = Option(
    ("--bt2-unal/--bt2-no-unal",),
    type=bool,
    default=False,
    help="Whether to output unaligned reads")

opt_bt2_i = Option(
    ("--bt2-I",),
    type=int,
    default=0,
    help="Minimum fragment length for valid paired-end alignments")

opt_bt2_x = Option(
    ("--bt2-X",),
    type=int,
    default=600,
    help="Maximum fragment length for valid paired-end alignments")

opt_bt2_score_min_e2e = Option(
    ("--bt2-score-min-e2e",),
    type=str,
    default="L,-1,-0.5",
    help="Minimum score for a valid end-to-end alignment")

opt_bt2_score_min_loc = Option(
    ("--bt2-score-min-loc",),
    type=str,
    default="L,1,0.5",
    help="Minimum score for a valid local alignment")

opt_bt2_s = Option(
    ("--bt2-i", "bt2_s"),
    type=str,
    default="L,1,0.1",
    help="Seed interval")

opt_bt2_l = Option(
    ("--bt2-L",),
    type=int,
    default=12,
    help="Seed length")

opt_bt2_gbar = Option(
    ("--bt2-gbar",),
    type=int,
    default=4,
    help="Minimum distance of a gap from end of a read")

opt_bt2_d = Option(
    ("--bt2-D",),
    type=int,
    default=4,
    help="Maximum number of failed seed extensions")

opt_bt2_r = Option(
    ("--bt2-R",),
    type=int,
    default=2,
    help="Maximum number of times to re-seed")

opt_bt2_dpad = Option(
    ("--bt2-dpad",),
    type=int,
    default=2,
    help="Width of padding on alignment matrix, to allow gaps")

opt_bt2_orient = Option(
    ("--bt2-orient",),
    type=Choice(BOWTIE2_ORIENT, case_sensitive=False),
    default=BOWTIE2_ORIENT[0],
    help="Valid orientations of paired-end mates")

# Reference section specification options
opt_coords = Option(
    ("--coords", "-c"),
    type=(str, int, int),
    multiple=True,
    default=(),
    help=("Reference name, 5' end, and 3' end of a section; "
          "coordinates are 1-indexed and include both ends"))

opt_primers = Option(
    ("--primers", "-p"),
    type=(str, str, str),
    multiple=True,
    default=(),
    help=("Reference name, forward primer, and reverse primer "
          "of a section; reverse primer must be given 5' to 3'"))

opt_primer_gap = Option(
    ("--primer-gap",),
    type=int,
    default=2,
    help=("Number of bases to leave as a gap between the "
          "end of a primer and the end of the section"))

# Relate
opt_min_reads = Option(
    ("--min-reads", "-n"),
    type=int,
    default=1000,
    help=("Do not analyze datasets with fewer than this "
          "minimum number of reads"))

opt_batch_size = Option(
    ("--batch-size",),
    type=float,
    default=256.,
    help=("Target size of each batch of mutation vectors, "
          "in millions of base calls"))

opt_ambrel = Option(
    ("--ambrel/--no-ambrel",),
    type=bool,
    default=True,
    help=("Whether to find and label all ambiguous "
          "insertions and deletions (improves accuracy "
          "but runs slower)"))

# Mask

opt_count_del = Option(
    ("--count-del/--discount-del",),
    type=bool,
    default=False,
    help='Whether to count deletions as mutations')

opt_count_ins = Option(
    ("--count-ins/--discount-ins",),
    type=bool,
    default=False,
    help='Whether to count insertions as mutations')

opt_discount_mut = Option(
    ("--discount-mut",),
    type=str,
    multiple=True,
    default=(),
    help="Do not count a specific type of mutation")

opt_exclude_polya = Option(
    ("--exclude-polya",),
    type=int,
    default=5,
    help="Exclude stretches of consecutive A bases of at least this length. "
         "If 0, exclude none.")

opt_exclude_gu = Option(
    ("--exclude-gu/--include-gu",),
    type=bool,
    default=True,
    help="Exclude positions with G and U bases (which DMS methylates very "
         "weakly at physiological pH).")

opt_exclude_pos = Option(
    ("--exclude-pos",),
    type=(str, int),
    default=(),
    multiple=True,
    help="Exclude arbitrary positions, given as (reference, position).")

opt_min_ncall_read = Option(
    ("--min-ncall-read",),
    type=int,
    default=1,
    help="Filter reads with less than this number of base calls.")

opt_min_finfo_read = Option(
    ("--min-finfo-read",),
    type=float,
    default=0.95,
    help="Filter reads with less than this fraction of informative positions.")

opt_max_fmut_read = Option(
    ("--max-fmut-read",),
    type=float,
    default=0.05,
    help="Filter reads with more than this fraction of mutated positions.")

opt_max_nmut_read = Option(
    ("--max-nmut-read",),
    type=int,
    default=-1,
    help="Filter reads with more than this number of mutated positions "
         "(-1 to disable).")

opt_min_mut_gap = Option(
    ("--min-mut-gap",),
    type=int,
    default=3,
    help="Filter out reads with any pair of mutations "
         "that are separated by fewer than this number of "
         "non-mutated bases. If 0, filter out no reads.")

opt_min_ninfo_pos = Option(
    ("--min-ninfo-pos",),
    type=int,
    default=1000,
    help="Filter positions with less than this number of informative reads.")

opt_max_fmut_pos = Option(
    ("--max-fmut-pos",),
    type=float,
    default=0.5,
    help="Filter positions with more than this fraction of mutated reads.")

# Clustering options
opt_max_clusters = Option(
    ("--max-clusters", "-k"),
    type=int,
    default=0,
    help="End the clustering step after attempting this "
         "number of clusters, even if it yields the best "
         "(smallest) BIC observed thus far. If 0, do not "
         "run clustering.")

opt_em_runs = Option(
    ("--em-runs", "-e"),
    type=int,
    default=6,
    help="Run clustering this many times for each number of "
         "clusters, randomly initializing each run.")

opt_min_em_iter = Option(
    ("--min-em-iter",),
    type=int,
    default=10,
    help="Run every EM run for at least this number times "
         "k iterations, even if the likelihood value has "
         "already converged.")

opt_max_em_iter = Option(
    ("--max-em-iter",),
    type=int,
    default=300,
    help="Stop every EM run after this number times k "
         "iterations, even if the likelihood value has "
         "not yet converged.")

opt_em_thresh = Option(
    ("--em-thresh",),
    type=float,
    default=0.01,
    help="Consider an EM run to have converged when the log "
         "likelihood value has increased by less than this "
         "threshold between two consecutive iterations.")

opt_min_nmut_read = Option(
    ("--min-nmut-read",),
    type=int,
    default=0,
    help="Ignore reads with less than this number of mutated positions"
         "during clustering.")

# Tables

opt_rels = Option(
    ("--rels", "-r"),
    type=str,
    multiple=True,
    default=("m",),
    help="Output these relationships between the read and the reference")

# RNA structure prediction

opt_fold = Option(
    ("--fold/--no-fold",),
    type=bool,
    default=False,
    help="Whether to predict the RNA structure using Fold")

opt_quantile = Option(
    ("--quantile", "-q"),
    type=float,
    default=-1.,
    help="Quantile of mutation rates for normalization and winsorization "
         "(-1 to disable both). Must be -1 or in [0, 1].")

# Graphing

opt_arrange = Option(
    ("--arrange",),
    type=Choice(CLUST_ARRANGE_OPTIONS),
    default=CLUST_ORDER,
    help="Whether to graph each INDIVidual cluster in its own file, "
         "graph each ORDER of clusters in its own file, "
         "or UNITE all clusters in one file containing all orders.")

opt_x_ratio = Option(
    ("--x-ratio/--x-count",),
    default=False,
    type=bool,
    help="Whether the x-axis depicts counts or ratios")

opt_y_ratio = Option(
    ("--y-ratio/--y-count",),
    default=True,
    type=bool,
    help="Whether the y-axis depicts counts or ratios")

opt_hist_bins = Option(
    ("--hist-bins",),
    default=24,
    type=int,
    help="Number of bins in each histogram (≥ 1)")

opt_csv = Option(
    ("--csv/--no-csv",),
    default=True,
    type=bool,
    help="Whether to output the source data for each graph "
         "as a CSV file")

opt_html = Option(
    ("--html/--no-html",),
    default=True,
    type=bool,
    help="Whether to output each graph as an HTML file")

opt_pdf = Option(
    ("--pdf/--no-pdf",),
    default=False,
    type=bool,
    help="Whether to output each graph as a PDF file")

# Logging options
opt_verbose = Option(
    ("--verbose", "-v"),
    count=True,
    help="Print info or info+debug messages to stdout")

opt_quiet = Option(
    ("--quiet",),
    count=True,
    help="Suppress warnings or warnings+errors to stdout")

opt_log = Option(
    ("--log",),
    type=Path(exists=False, dir_okay=False),
    default=os.path.join(CWD, "log", datetime.now().strftime(
        "seismic-rna_%Y-%m-%d_%H-%M-%S.log")),
    help="File in which to log all messages (except profiling)")

opt_profile = Option(
    ("--profile",),
    type=Path(exists=False, dir_okay=False),
    default="",
    help="Profile code performance and log results to the given file")


def merge_params(*param_lists: list[Parameter]):
    """ Merge lists of Click parameters, dropping duplicates. """
    params = list()
    names = set()
    for param_list in param_lists:
        for param in param_list:
            if param.name not in names:
                params.append(param)
                names.add(param.name)
    return params
