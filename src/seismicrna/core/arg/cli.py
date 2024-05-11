"""
Core Command Line Interface
===========================
Auth: Yves, Matty

Define all command line interface (CLI) options and their defaults.
"""

import logging
import os
from datetime import datetime

import click
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

CLUST_INDIV = "indiv"
CLUST_ORDER = "order"
CLUST_UNITE = "unite"
CLUST_ARRANGE_OPTIONS = CLUST_INDIV, CLUST_ORDER, CLUST_UNITE

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
    help="Configuration file for parameters"
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
    help="destination for all output files"
)

opt_temp_dir = Option(
    ("--temp-dir", "-t"),
    type=Path(file_okay=False),
    default=os.path.join(".", "temp"),
    help="destination for all temporary files"
)

opt_keep_temp = Option(
    ("--keep-temp/--erase-temp",),
    type=bool,
    default=False,
    help="keep temporary files after the program exits"
)

# Resource usage options
opt_parallel = Option(
    ("--parallel/--serial",),
    type=bool,
    default=True,
    help="run tasks in parallel"
)

opt_max_procs = Option(
    ("--max-procs",),
    type=int,
    default=NUM_CPUS,
    help="maximum number of simultaneous processes"
)

# Experiment and analysis setup options

opt_sections_file = Option(
    ("--sections-file", "-s"),
    type=Path(dir_okay=False),
    default="",
    help="CSV file of sections by name, reference, and coordinates/primers"
)

opt_force = Option(
    ("--force/--no-force",),
    type=bool,
    default=False,
    help="force all tasks to run, even if their output files exist"
)

# Sequencing read (FASTQ) files
opt_fastqz = Option(
    ("--fastqz", "-z"),
    type=Path(exists=True),
    multiple=True,
    default=(),
    help="FASTQ files of single-end reads"
)

opt_fastqy = Option(
    ("--fastqy", "-y"),
    type=Path(exists=True),
    multiple=True,
    default=(),
    help="FASTQ files of interleaved paired-end reads"
)

opt_fastqx = Option(
    ("--fastqx", "-x"),
    type=Path(exists=True),
    multiple=True,
    default=(),
    help="FASTQ files of separated paired-end reads"
)

# Sequencing read (FASTQ/XAM) options

opt_phred_enc = Option(
    ("--phred-enc",),
    type=int,
    default=33,
    help="Phred score encoding in FASTQ and SAM/BAM/CRAM files"
)

opt_min_phred = Option(
    ("--min-phred",),
    type=int,
    default=DEFAULT_MIN_PHRED,
    help="minimum Phred score to use a base call"
)

opt_fastqc = Option(
    ("--fastqc/--no-fastqc",),
    type=bool,
    default=True,
    help="run FastQC on the initial and trimmed FASTQ files"
)

opt_qc_extract = Option(
    ("--qc-extract/--qc-no-extract",),
    type=bool,
    default=False,
    help="unzip FastQC report files"
)

# Demultiplexing options

opt_demultiplex = Option(
    ("--demult-on/--demult-off",),
    type=bool,
    default=False,
    help="enable demultiplexing"
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
    help="index of start of barcode")

opt_barcode_end = Option(
    ("--barcode-end",),
    type=int,
    default=0,
    help="length of barcode")

opt_barcode_length = Option(
    ("--barcode-length",),
    type=int,
    default=0,
    help="length of barcode")

opt_demulti_overwrite = Option(
    ("--demulti-overwrite",),
    type=bool,
    default=False,
    help="desiginates whether to overwrite the grepped fastq. should only be used if changing setting on the same "
         "sample")

# Demultiplexed sequencing read (FASTQ) directories

opt_dmfastqz = Option(
    ("--dmfastqz", "-Z"),
    type=Path(exists=True, file_okay=False),
    multiple=True,
    default=(),
    help="Demultiplexed FASTQ files of single-end reads"
)

opt_dmfastqy = Option(
    ("--dmfastqy", "-Y"),
    type=Path(exists=True, file_okay=False),
    multiple=True,
    default=(),
    help="Demultiplexed FASTQ files of paired-end reads interleaved in one file"
)

opt_dmfastqx = Option(
    ("--dmfastqx", "-X"),
    type=Path(exists=True, file_okay=False),
    multiple=True,
    default=(),
    help="Demultiplexed FASTQ files of mate 1 and mate 2 reads"
)

# Adapter trimming options with Cutadapt
opt_cutadapt = Option(
    ("--cut/--no-cut",),
    type=bool,
    default=True,
    help="trim reads with Cutadapt before alignment"
)

opt_cut_q1 = Option(
    ("--cut-q1",),
    type=int,
    default=DEFAULT_MIN_PHRED,
    help="Phred score for read 1 quality trimming with Cutadapt"
)

opt_cut_q2 = Option(
    ("--cut-q2",),
    type=int,
    default=DEFAULT_MIN_PHRED,
    help="Phred score for read 2 quality trimming with Cutadapt"
)

opt_cut_g1 = Option(
    ("--cut-g1",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_5P,),
    help="5' adapter for read 1 adapter trimming with Cutadapt"
)

opt_cut_a1 = Option(
    ("--cut-a1",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_3P,),
    help="3' adapter for read 1 adapter trimming with Cutadapt"
)

opt_cut_g2 = Option(
    ("--cut-g2",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_5P,),
    help="5' adapter for read 2 adapter trimming with Cutadapt"
)

opt_cut_a2 = Option(
    ("--cut-a2",),
    type=str,
    multiple=True,
    default=(ADAPTER_SEQ_ILLUMINA_3P,),
    help="3' adapter for read 2 adapter trimming with Cutadapt"
)

opt_cut_o = Option(
    ("--cut-O",),
    type=int,
    default=6,
    help="minimum overlap of read and adapter during trimming with Cutadapt"
)

opt_cut_e = Option(
    ("--cut-e",),
    type=float,
    default=0.1,
    help="error tolerance for adapters during trimming with Cutadapt"
)

opt_cut_indels = Option(
    ("--cut-indels/--cut-no-indels",),
    type=bool,
    default=True,
    help="allow indels in adapters during trimming with Cutadapt"
)

opt_cut_nextseq = Option(
    ("--cut-nextseq/--cut-no-nextseq",),
    type=bool,
    default=False,
    help="trim high-quality Gs from 3' end during trimming with Cutadapt"
)

opt_cut_discard_trimmed = Option(
    ("--cut-discard-trimmed/--cut-keep-trimmed",),
    type=bool,
    default=False,
    help="discard reads in which an adapter was found by Cutadapt"
)

opt_cut_discard_untrimmed = Option(
    ("--cut-discard-untrimmed/--cut-keep-untrimmed",),
    type=bool,
    default=False,
    help="discard reads in which no adapter was found by Cutadapt"
)

opt_cut_m = Option(
    ("--cut-m",),
    type=int,
    default=20,
    help="minimum length of a read to keep it after trimming with Cutadapt"
)

# Alignment options with Bowtie2
opt_bt2_local = Option(
    ("--bt2-local/--bt2-end-to-end",),
    type=bool,
    default=True,
    help="run Bowtie2 in local mode"
)

opt_bt2_discordant = Option(
    ("--bt2-discordant/--bt2-no-discordant",),
    type=bool,
    default=False,
    help="output discordant alignments from Bowtie2")

opt_bt2_mixed = Option(
    ("--bt2-mixed/--bt2-no-mixed",),
    type=bool,
    default=False,
    help="attempt to align individual mates of unaligned pairs with Bowtie2"
)

opt_bt2_dovetail = Option(
    ("--bt2-dovetail/--bt2-no-dovetail",),
    type=bool,
    default=False,
    help="treat dovetailed mate pairs as concordant with Bowtie2"
)

opt_bt2_contain = Option(
    ("--bt2-contain/--bt2-no-contain",),
    type=bool,
    default=True,
    help="treat nested mate pairs as concordant with Bowtie2"
)

opt_bt2_un = Option(
    ("--bt2-un/--bt2-no-un",),
    type=bool,
    default=True,
    help="output unaligned reads from Bowtie2 to a FASTQ file"
)

opt_bt2_i = Option(
    ("--bt2-I",),
    type=int,
    default=0,
    help="minimum fragment length for valid paired-end alignments with Bowtie2"
)

opt_bt2_x = Option(
    ("--bt2-X",),
    type=int,
    default=600,
    help="maximum fragment length for paired-end alignments with Bowtie2"
)

opt_bt2_score_min_e2e = Option(
    ("--bt2-score-min-e2e",),
    type=str,
    default="L,-1,-0.5",
    help="minimum alignment score for Bowtie2 in end-to-end mode"
)

opt_bt2_score_min_loc = Option(
    ("--bt2-score-min-loc",),
    type=str,
    default="L,1,0.5",
    help="minimum alignment score for Bowtie2 in local mode"
)

opt_bt2_s = Option(
    ("--bt2-i", "bt2_s"),
    type=str,
    default="L,1,0.1",
    help="seed interval for Bowtie2"
)

opt_bt2_l = Option(
    ("--bt2-L",),
    type=int,
    default=20,
    help="seed length for Bowtie2"
)

opt_bt2_gbar = Option(
    ("--bt2-gbar",),
    type=int,
    default=4,
    help="minimum distance of a gap from the end of a read with Bowtie2"
)

opt_bt2_d = Option(
    ("--bt2-D",),
    type=int,
    default=4,
    help="maximum number of consecutive failed seed extensions with Bowtie2"
)

opt_bt2_r = Option(
    ("--bt2-R",),
    type=int,
    default=2,
    help="maximum number of re-seeding attempts with Bowtie2"
)

opt_bt2_dpad = Option(
    ("--bt2-dpad",),
    type=int,
    default=2,
    help="width of padding on alignment matrix (to allow gaps) with Bowtie2"
)

opt_bt2_orient = Option(
    ("--bt2-orient",),
    type=Choice(BOWTIE2_ORIENT, case_sensitive=False),
    default=BOWTIE2_ORIENT[0],
    help="valid orientations of paired-end mates with Bowtie2"
)

opt_min_mapq = Option(
    ("--min-mapq",),
    type=int,
    default=25,
    help="minimum mapping quality to use an aligned read"
)

opt_cram = Option(
    ("--cram/--bam",),
    type=bool,
    default=False,
    help="compress alignment maps using CRAM format"
)

# Reference section specification options
opt_coords = Option(
    ("--coords", "-c"),
    type=(str, int, int),
    multiple=True,
    default=(),
    help="reference name and 5' and 3' ends of a section"
)

opt_primers = Option(
    ("--primers", "-p"),
    type=(str, DNA, DNA),
    multiple=True,
    default=(),
    help="Reference name and forward and reverse primers for a section"
)

opt_primer_gap = Option(
    ("--primer-gap",),
    type=int,
    default=0,
    help="length of the gap (nt) between each primer and the end of the section"
)

# Relate
opt_min_reads = Option(
    ("--min-reads", "-N"),
    type=int,
    default=1000,
    help="minimum number of reads in an alignment map to keep it"
)

opt_batch_size = Option(
    ("--batch-size",),
    type=float,
    default=256.,
    help="target size of each batch (millions of base calls)"
)

opt_ambrel = Option(
    ("--ambrel/--no-ambrel",),
    type=bool,
    default=True,
    help="mark all ambiguous indels"
)

opt_brotli_level = Option(
    ("--brotli-level",),
    type=int,
    default=DEFAULT_BROTLI_LEVEL,
    help="compression level for Brotli (0 - 11)"
)

opt_overhangs = Option(
    ("--overhangs/--no-overhangs",),
    type=bool,
    default=True,
    help="compute relationships for bases in a read that extend beyond its mate"
)

# Pool

opt_pool = Option(
    ("--pool", "-P"),
    type=str,
    default="",
    help="name of the pooled sample"
)

# Mask

opt_count_del = Option(
    ("--count-del/--discount-del",),
    type=bool,
    default=False,
    help="count deletions as mutations"
)

opt_count_ins = Option(
    ("--count-ins/--discount-ins",),
    type=bool,
    default=False,
    help="count insertions as mutations"
)

opt_discount_mut = Option(
    ("--discount-mut",),
    type=str,
    multiple=True,
    default=(),
    help="type of mutation to ignore"
)

opt_exclude_polya = Option(
    ("--exclude-polya",),
    type=int,
    default=5,
    help="minimum length of poly(A) sequences to exclude (0 to include all)"
)

opt_exclude_gu = Option(
    ("--exclude-gu/--include-gu",),
    type=bool,
    default=True,
    help="exclude G and U bases"
)

opt_exclude_file = Option(
    ("--exclude-file",),
    type=Path(dir_okay=False, exists=True),
    help="file of arbitrary positions to exclude"
)

opt_min_ncov_read = Option(
    ("--min-ncov-read",),
    type=int,
    default=1,
    help="minimum coverage of a read to keep it (nt)"
)

opt_min_finfo_read = Option(
    ("--min-finfo-read",),
    type=float,
    default=0.95,
    help="minimum fraction of information in a read to keep it"
)

opt_max_fmut_read = Option(
    ("--max-fmut-read",),
    type=float,
    default=0.1,
    help="maximum fraction of mutations in a read to keep it"
)

opt_min_mut_gap = Option(
    ("--min-mut-gap",),
    type=int,
    default=3,
    help="minimum gap between two mutations in a read to keep it (nt)"
)

opt_min_ninfo_pos = Option(
    ("--min-ninfo-pos",),
    type=int,
    default=1000,
    help="minimum information count at a position to use it"
)

opt_max_fmut_pos = Option(
    ("--max-fmut-pos",),
    type=float,
    default=0.5,
    help="maximum mutation fraction at a position to use it"
)

opt_quick_unbias = Option(
    ("--quick-unbias/--exact-unbias",),
    type=bool,
    default=True,
    help="use a quick approximation to correct observer bias in mutation rates"
)

opt_quick_unbias_thresh = Option(
    ("--quick-unbias-thresh",),
    type=float,
    default=0.005,
    help="maximum mutation rate to consider 0 when approximating observer bias"
)

# Cluster options

opt_max_clusters = Option(
    ("--max-clusters", "-k"),
    type=int,
    default=0,
    help="maximum number of clusters to attempt (0 to disable clustering)"
)

opt_em_runs = Option(
    ("--em-runs", "-e"),
    type=int,
    default=12,
    help="number of independent EM runs for each order (number of clusters)"
)

opt_min_em_iter = Option(
    ("--min-em-iter",),
    type=int,
    default=8,
    help="minimum iterations per EM run"
)

opt_max_em_iter = Option(
    ("--max-em-iter",),
    type=int,
    default=512,
    help="maximum iterations per EM run"
)

opt_em_thresh = Option(
    ("--em-thresh",),
    type=float,
    default=0.1,
    help="threshold for convergence of log likelihood"
)

# Join options

opt_joined = Option(
    ("--joined", "-J"),
    type=str,
    default="",
    help="name of the joined section"
)

opt_join_clusts = Option(
    ("--join-clusts", "-j"),
    type=Path(dir_okay=False, exists=True),
    help="CSV file of clusters to join"
)

# Table options

opt_table_pos = Option(
    ("--table-pos/--no-table-pos",),
    type=bool,
    default=True,
    help="tabulate per position"
)

opt_table_read = Option(
    ("--table-read/--no-table-read",),
    type=bool,
    default=True,
    help="tabulate per read"
)

opt_table_clust = Option(
    ("--table-clust/--no-table-clust",),
    type=bool,
    default=True,
    help="tabulate per cluster"
)

# List options

opt_complement = Option(
    ("--complement/--no-complement",),
    type=bool,
    default=False,
    help="list the complement"
)

# Fold

opt_fold = Option(
    ("--fold/--no-fold",),
    type=bool,
    default=False,
    help="predict the secondary structure using the RNAstructure Fold program"
)

opt_quantile = Option(
    ("--quantile", "-q"),
    type=float,
    default=0.,
    help="quantile for normalizing ratios; must be in [0, 1]"
)

opt_fold_temp = Option(
    ("--fold-temp",),
    type=float,
    default=310.15,
    help="temperature at which to predict structures (Kelvin)"
)

opt_fold_constraint = Option(
    ("--fold-constraint",),
    type=click.Path(exists=True, dir_okay=False),
    help="file of constraints for predicting structures"
)

opt_fold_md = Option(
    ("--fold-md",),
    type=int,
    default=0,
    help="maximum distance between two paired bases in predicted structures "
         "(0 for no limit)"
)

opt_fold_mfe = Option(
    ("--fold-mfe/--fold-sub",),
    type=bool,
    default=False,
    help="predict only the minimum free energy (MFE) structure "
         "(disables --fold-max and --fold-percent)"
)

opt_fold_max = Option(
    ("--fold-max",),
    type=int,
    default=20,
    help="maximum number of predicted structures (ignored if using --fold-mfe)"
)

opt_fold_percent = Option(
    ("--fold-percent",),
    type=float,
    default=20.,
    help="Maximum % difference in energy between predicted structures "
         "(ignored if using --fold-mfe)"
)

# Graph

opt_comppair = Option(
    ("--comppair/--no-comppair",),
    type=bool,
    default=True,
    help="compare every pair of input files"
)

opt_compself = Option(
    ("--compself/--no-compself",),
    type=bool,
    default=False,
    help="compare every input file with itself"
)

opt_cgroup = Option(
    ("--cgroup",),
    type=Choice(CLUST_ARRANGE_OPTIONS),
    default=CLUST_ORDER,
    help=("graph each INDIVidual cluster in its own file, each ORDER in its "
          "own file, or UNITE all clusters in one file containing all orders")
)

opt_rels = Option(
    ("--rels", "-r"),
    type=str,
    multiple=True,
    default=("m",),
    help="relationship(s) to graph"
)

opt_use_ratio = Option(
    ("--use-ratio/--use-count",),
    type=bool,
    default=True,
    help="graph ratios or counts"
)

opt_window = Option(
    ("--window", "-w"),
    type=int,
    default=45,
    help="size of the sliding window (nt)"
)

opt_winmin = Option(
    ("--winmin", "-n"),
    type=int,
    default=9,
    help="minimum number of data to use a sliding window (otherwise NaN)"
)

opt_metric = Option(
    ("--metric", "-m"),
    type=Choice(METRIC_KEYS, case_sensitive=False),
    default=KEY_PEARSON,
    help=(f"metric to compare mutation rates: "
          f"{repr(KEY_NRMSD)} = normalized root-mean-square deviation (NRMSD), "
          f"{repr(KEY_PEARSON)} = Pearson correlation coefficient (r), "
          f"{repr(KEY_SPEARMAN)} = Spearman correlation coefficient (ρ), "
          f"{repr(KEY_DETERM)} = coefficient of determination (R²)")
)

opt_struct_file = Option(
    ("--struct-file",),
    type=click.Path(exists=True, dir_okay=False),
    help="connectivity table (CT) file of structures, against which to compare "
         "mutational profiles (disables --struct-sect)"
)

opt_struct_sect = Option(
    ("--struct-sect",),
    type=str,
    default="",
    help="section of RNA that was folded, against which to compare mutational "
         "profiles (ignored if blank or if using --struct-file)"
)

opt_hist_bins = Option(
    ("--hist-bins",),
    type=int,
    default=10,
    help="number of bins in each histogram; must be ≥ 1"
)

opt_hist_margin = Option(
    ("--hist-margin",),
    type=float,
    default=0.01,
    help="for ratio-based histograms, maximum margin to autofill between the "
         "data and 0 or 1"
)

opt_csv = Option(
    ("--csv/--no-csv",),
    type=bool,
    default=True,
    help="output the source data for each graph as a CSV file"
)

opt_html = Option(
    ("--html/--no-html",),
    type=bool,
    default=True,
    help="output each graph as an HTML file"
)

opt_pdf = Option(
    ("--pdf/--no-pdf",),
    type=bool,
    default=False,
    help="output each graph as a PDF file"
)

# CT renumbering

opt_ct_pos_5 = Option(
    ("--ct-pos-5", "-c"),
    type=(click.Path(exists=True), int),
    multiple=True,
    help="connectivity table (CT) file or directory of CT files "
         "and the 5' position to assign to each file"
)

opt_inplace = Option(
    ("--inplace/--newfile",),
    type=bool,
    default=False,
    help="modify files in-place instead of writing new files "
         "(WARNING: original files will not be recoverable afterwards)"
)

# Export

opt_export = Option(
    ("--export/--no-export",),
    type=bool,
    default=False,
    help="export per-sample results for the seismic-graph web app"
)

opt_samples_meta = Option(
    ("--samples-meta", "-S"),
    type=Path(dir_okay=False, exists=True),
    help="CSV file of metadata for each sample"
)

opt_refs_meta = Option(
    ("--refs-meta", "-R"),
    type=Path(dir_okay=False, exists=True),
    help="CSV file of metadata for each reference"
)

opt_all_pos = Option(
    ("--all-pos/--unmasked-pos",),
    type=bool,
    default=False,
    help="export all positions (not just unmasked positions)"
)

# Logging options
opt_verbose = Option(
    ("--verbose", "-v"),
    count=True,
    help="log info (-v) or info and debug (-vv) messages on stdout"
)

opt_quiet = Option(
    ("--quiet", "-q"),
    count=True,
    help="suppress warnings (-q) or warnings and errors (-qq) on stdout"
)

opt_log = Option(
    ("--log",),
    type=Path(exists=False, dir_okay=False),
    default=os.path.join(CWD, "log", datetime.now().strftime(
        "seismic-rna_%Y-%m-%d_%H-%M-%S.log"
    )),
    help="file in which to log all messages (except profiling)"
)

opt_log_color = Option(
    ("--log-color/--log-plain",),
    type=bool,
    default=True,
    help="log messages with or without color codes on stdout"
)

opt_profile = Option(
    ("--profile",),
    type=Path(exists=False, dir_okay=False),
    default="",
    help="profile code performance and log results to the given file"
)


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
