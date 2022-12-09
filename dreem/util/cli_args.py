import dreem.util.util as util
import click
import os
from click_option_group import optgroup


# Input files
out_dir = optgroup.option('--out_dir', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Where to output files')
sample = optgroup.option('--sample', '-s', type=click.Path(exists=True), help='Name of the sequence alignment map file(s) folder')
fasta = optgroup.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
fastq1 = optgroup.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', type=click.Path(exists=True), required=True)
fastq2 = optgroup.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', type=click.Path(exists=True))
input_dir = optgroup.option('--input_dir', '-id', type=click.Path(exists=True), help='Sequence alignment map files folder(s) generated by alignment', multiple=True)
samples = optgroup.option('--samples', '-s', type=click.Path(exists=True), help='Path to the samples.csv file', default=None)
clustering_file = optgroup.option('--clusters', '-cl', type=click.Path(exists=True), help='Path to the clustering.json file', default=None)
library = optgroup.option('--library', '-l', type=click.Path(exists=True), help='Path to the library.csv file', default=None)

# Construct selection
coords = optgroup.option('--coords', '-c', type=(str, int, int), multiple=True, help="coordinates for reference: '-c ref-name first last'")
primers = optgroup.option('--primers', '-p', type=(str, str, str), multiple=True, help="primers for reference: '-p ref-name fwd rev'")
fill = optgroup.option('--fill/--no-fill', None, type=bool, default=False, help="Fill in coordinates of reference sequences for which neither coordinates nor primers were given (default: no).")
parallel = optgroup.option('--parallel', '-P', type=click.Choice(["profiles", "reads", "off", "auto"], case_sensitive=False), default="auto", help="Parallelize the processing of mutational PROFILES or READS within each profile, turn parallelization OFF, or AUTO matically choose the parallelization method (default: auto).")

# Demultiplexing
demultiplexing = optgroup.option('--demultiplexing', '-dx', type=bool, help='Use demultiplexing', default=False)
barcode_start = optgroup.option('--barcode_start', '-bs', type=int, help='Start position of the barcode in the read')
barcode_end = optgroup.option('--barcode_end', '-be', type=int, help='End position of the barcode in the read')
max_mutations_on_barcode = optgroup.option('--max_mutations_on_barcode', '-mb', type=int, help='Maximum number of mutations on the barcode', default=1)

# Bowtie2 TODO

# Cutadapt #TODO

# Clustering
clustering = optgroup.option('--clustering', '-cl', type=bool, help='Use clustering', default=False)
n_clusters = optgroup.option('--n_clusters', '-nc', type=int, help='Number of clusters', default=None)
max_clusters = optgroup.option('--max_clusters', '-mc', type=int, help='Maximum number of clusters', default=None)
signal_thresh = optgroup.option('--signal_thresh', '-st', type=float, help='Signal threshold', default=None)
info_thresh = optgroup.option('--info_thresh', '-it', type=float, help='Information threshold', default=None)
include_g_u = optgroup.option('--include_g_u', '-igu', type=bool, help='Include G and U', default=None)
include_del = optgroup.option('--include_del', '-id', type=bool, help='Include deletions', default=None)
min_reads = optgroup.option('--min_reads', '-mr', type=int, help='Minimum number of reads', default=None)
convergence_cutoff = optgroup.option('--convergence_cutoff', '-cc', type=float, help='Convergence cutoff', default=None)
num_runs = optgroup.option('--num_runs', '-nr', type=int, help='Number of runs', default=None)

# Aggregation
rnastructure_path = optgroup.option('--rnastructure_path', '-rs', type=click.Path(exists=True), help='Path to RNAstructure, to predict structure and free energy', default=None)
rnastructure_temperature = optgroup.option('--rnastructure_temperature', '-rst', type=bool, help='Use sample.csv temperature values for RNAstructure', default=False)
rnastructure_fold_args = optgroup.option('--rnastructure_fold_args', '-rsa', type=str, help='optgroup.options to pass to RNAstructure fold', default=None)
rnastructure_dms = optgroup.option('--rnastructure_dms', '-rsd', type=bool, help='Use the DMS signal to make predictions with RNAstructure', default=False)
rnastructure_dms_min_unpaired_value = optgroup.option('--rnastructure_dms_min_unpaired_value', '-rsdmin', type=int, help='Minimum unpaired value for using the dms signal as an input for RNAstructure', default=0.01)
rnastructure_dms_max_paired_value = optgroup.option('--rnastructure_dms_max_paired_value', '-rsdmax', type=int, help='Maximum paired value for using the dms signal as an input for RNAstructure', default=0.05)
rnastructure_partition = optgroup.option('--rnastructure_partition', '-rspa', type=bool, help='Use RNAstructure partition function to predict free energy', default=False)
rnastructure_probability = optgroup.option('--rnastructure_probability', '-rspr', type=bool, help='Use RNAstructure partition function to predict per-base mutation probability', default=False)
poisson = optgroup.option('--poisson', '-po', type=bool, help='Predict Poisson confidence intervals', default=True)

# Misc
verbose = optgroup.option('--verbose', '-v', type=bool, help='Verbose output', default=False)

