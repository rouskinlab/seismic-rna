********************************************************************************
seismic sim
********************************************************************************


Purpose
================================================================================

``seismic sim`` generates synthetic sequencing data for testing and benchmarking.
It is a group of subcommands that simulate individual components of the pipeline
(reference sequences, structures, mutation rates, reads) or complete datasets.


Subcommands
================================================================================

Run ``seismic sim {subcommand} --help`` for each subcommand's full option list.

============= ===================================================================
Subcommand    Description
============= ===================================================================
``ref``       Simulate a random reference sequence
``fold``      Simulate an RNA secondary structure
``params``    Simulate mutation-rate parameters for clusters
``clusts``    Simulate cluster memberships for reads
``ends``      Simulate read end positions
``muts``      Simulate per-read, per-position mutations
``idmut``     Simulate IDmut batch output directly
``fastq``     Simulate a FASTQ file from a reference
``total``     Run a full end-to-end simulation
``abstract``  Simulate abstract mutational data
============= ===================================================================


Quick example
================================================================================

Run a full end-to-end simulation::

    seismic sim total --help

See the auto-generated :doc:`/cli` for the full option list for each subcommand.


See also
================================================================================

- :doc:`/use/workflow/idmut` — the step whose output ``sim idmut`` mimics
- :doc:`/algos/index` — the underlying mutation models
