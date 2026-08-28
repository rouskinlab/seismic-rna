********************************************************************************
Tutorial 4: Running multiple samples at once
********************************************************************************

SEISMIC-RNA allows for multiple samples to be run with one simple command.

Download example files
--------------------------------------------------------------------------------

Paired-end FASTQ files were generated using ``seismic sim``. The files, with
their corresponding reference fasta file, can be downloaded here:

If you have ``wget``, you can download the tutorial data simply by typing ::

    wget https://raw.githubusercontent.com/rouskinlab/seismic-rna/main/src/userdocs/tutorials/multi-sample/data.tar

Otherwise, click this link to download the tutorial data:
https://raw.githubusercontent.com/rouskinlab/seismic-rna/main/src/userdocs/tutorials/multi-sample/data.tar

To ensure the download is complete and not corrupted, verify that the SHA-256
checksum is ``17da8b10d2f053f7549b72c8b8dc0b37f527edca712c6f3132582e62ae3d1df0``
by typing this command::

    shasum -a 256 data.tar

If this command prints a different checksum, then retry the download.
If the problem persists, then raise an issue (see :doc:`../../issues`).

After downloading and verifying the data, untar the data by typing ::

    tar xvf data.tar

and then navigate into the ``data`` directory::

    cd data

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
described in the previous tutorials.