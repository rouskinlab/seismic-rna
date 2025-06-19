********************************************************************************
Example 3) Running multiple samples at once
********************************************************************************

SEISMIC-RNA allows for multiple samples to be run with one simple command.

Download example files
--------------------------------------------------------------------------------

Paired-end FASTQ files were generated using ``seismic sim``. The files, with
their corresponding reference fasta file, can be downloaded here:
https://raw.githubusercontent.com/rouskinlab/seismic-rna/main/src/userdocs/examples/MultiSample/fq/MultiSample.zip

Run the SEISMIC-RNA workflow
--------------------------------------------------------------------------------

To run multiple samples at once, all that is needed is to provide a directory
where *all* the desired FASTQ files can be found (they can even be distributed in
subdirectories without it being an issue). If the different FASTQs are to be
aligned against the same reference, then a single sequence in a fasta file will
suffice. If the FASTQs are to be aligned against more than one sequence, then the
fasta file must contain all the desired sequences::

    seismic wf fq/sim_multiple.fa -x fq/


Output
--------------------------------------------------------------------------------
In the chosen output folder (./out by default), directories with the names of
each sample will be created, each one with the same subdirectories already
described in the previous examples.