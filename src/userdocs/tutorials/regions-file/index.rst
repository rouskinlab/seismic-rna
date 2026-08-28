********************************************************************************
Tutorial 5: Masking using a regions file
********************************************************************************

When working with a complex sequence, SEISMIC-RNA allows to select distinct parts
of it for analyzing.

In the case that, for example, a pair of primers was used to create an amplicon,
it is possible to mask the positions of the primers using a *regions file*
(see :doc:`../../formats/meta/regions`). Conceptually, this prevents the lack
of mutations in the primer sequence regions to be incorrectly considered low
reactivity when folding.

Download example files
--------------------------------------------------------------------------------

Paired-end FASTQ files were generated using ``seismic sim``. The files, with
their corresponding reference fasta file and the regions file, can be
downloaded here:

If you have ``wget``, you can download the tutorial data simply by typing ::

    wget https://raw.githubusercontent.com/rouskinlab/seismic-rna/main/src/userdocs/tutorials/regions-file/data.tar

Otherwise, click this link to download the tutorial data:
https://raw.githubusercontent.com/rouskinlab/seismic-rna/main/src/userdocs/tutorials/regions-file/data.tar

To ensure the download is complete and not corrupted, verify that the SHA-256
checksum is ``214e5fe859291a5f18bc67519047f58cb8ab4b2a7d56158ce7e6b4ed7fde3f3d``
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

The regions file can be included using ``--regions-file``. This will
mask the sequences determined by the file, and then use the rest for folding
(for which the flags ``--fold`` and ``--draw`` are added)::

    seismic wf fq/Regions_Ref.fa -x fq/ --regions-file fq/regions_file.csv --fold --draw

Output
--------------------------------------------------------------------------------
Given that the primer sequences were masked, the model shows no reactivity at the first nor last 20 bases:

    .. image:: img/regions_file_model.png
