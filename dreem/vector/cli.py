from dreem.vector.main import run
from dreem.util.cli import *


@click.command()
@opt_out_dir
@opt_library
@opt_coords
@opt_primers
@opt_spanall
@opt_parallel
@opt_min_phred
@opt_phred_enc
@opt_rerun
@arg_fasta
@arg_bams
@opt_max_cpus
def cli(*args, **opts):
    run(*args, **opts)


if __name__ == '__main__':
    cli()
