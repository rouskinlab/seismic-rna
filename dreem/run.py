import yaml, sys

from click_option_group import optgroup
import click
from util import Primer

@click.command()
@optgroup.group('Files and folders paths')
@optgroup.option('--root_dir', '-rd', default='', type=click.Path(exists=True), help='Where to output files and temp files')
@optgroup.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@optgroup.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True), required=True)
@optgroup.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True))
@optgroup.option('--samples', '-s', type=click.Path(exists=True), help='Path to the samples.csv file')
@optgroup.option('--library', '-l', type=click.Path(exists=True), help='Path to the library.csv file')

@optgroup.group('Construct selection')
@optgroup.option("-c", "--coords", type=(str, int, int), multiple=True, help="coordinates for reference: '-c ref-name first last'")
@optgroup.option("-p", "--primers", type=(str, Primer, Primer), multiple=True, help="primers for reference: '-p ref-name fwd rev'")
@optgroup.option("--fill/--no-fill", default=False,
              help="Fill in coordinates of reference sequences for which "
                   "neither coordinates nor primers were given (default: no).")


@optgroup.group('Demultiplexing parameters')
@optgroup.option('--demultiplexing', '-dx', type=bool, help='Use demultiplexing', default=False)
@optgroup.option('--barcode_start', '-bs', type=int, help='Start position of the barcode in the read')
@optgroup.option('--barcode_end', '-be', type=int, help='End position of the barcode in the read')

#@optgroup.group('Bowtie parameters') #TODO

#@optgroup.group('Cutadapt parameters') #TODO 

@optgroup.group('Vectoring parameters')
@optgroup.option("-P", "--parallel",
              type=click.Choice(["profiles", "reads", "off", "auto"],
                                case_sensitive=False),
              default="auto",
              help="Parallelize the processing of mutational PROFILES or "
              "READS within each profile, turn parallelization OFF, or AUTO"
              "matically choose the parallelization method (default: auto).")
@optgroup.option('--fill', type=bool, help='#TODO', default=False) #TODO default value

@optgroup.group('Clustering parameters')
@optgroup.option('--clustering', '-cl', type=bool, help='Use clustering', default=False)
@optgroup.option('--N_clusters', '-nc', type=int, help='Number of clusters', default=None)
@optgroup.option('--max_clusters', '-mc', type=int, help='Maximum number of clusters', default=None)
@optgroup.option('--signal_thresh', '-st', type=float, help='Signal threshold', default=None)
@optgroup.option('--info_thresh', '-it', type=float, help='Information threshold', default=None)
@optgroup.option('--include_G_U', '-igu', type=bool, help='Include G and U', default=None)
@optgroup.option('--include_del', '-id', type=bool, help='Include deletions', default=None)
@optgroup.option('--min_reads', '-mr', type=int, help='Minimum number of reads', default=None)
@optgroup.option('--convergence_cutoff', '-cc', type=float, help='Convergence cutoff', default=None)
@optgroup.option('--num_runs', '-nr', type=int, help='Number of runs', default=None)

@optgroup.group('Post_processing parameters')
@optgroup.option('--post_processing', '-pp', type=bool, help='Use post-processing', default=False)


def sanity_check():
    """
    Make sure that the input is valid. For example, that the FASTQ files exist, that the library file and the samples file are valid, etc.
    """
    pass

def run(**args):
    """Run DREEM.

    Args:
        args (_type_): _description_
    """


if __name__ == "__main__":
    run()