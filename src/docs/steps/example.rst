
.. _example-data:

Example dataset
========================================================================

Download example files
------------------------------------------------------------------------

You can download the example FASTA and FASTQ files in several ways,
according to your preference.

Download from the URL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Click each of these links to download the example FASTQ files, then save
them to a location of your choice:

- https://github.com/rouskinlab/seismic-rna/raw/main/example/sars2.fa
- https://github.com/rouskinlab/seismic-rna/raw/main/example/sars2-fse_R1.fq.gz
- https://github.com/rouskinlab/seismic-rna/raw/main/example/sars2-fse_R2.fq.gz

If file appears in the web browser but does not download automatically,
then right-click on the page and process "Save Page As...".

Download from GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a web browser, navigate to
https://github.com/rouskinlab/seismic-rna/tree/main/example.
From that page, click on each of the following files:

- sars2.fa
- sars2-fse_R1.fq.gz
- sars2-fse_R2.fq.gz

On each page, click the button labeled "Download raw file" near the
upper right corner of the page.

Download using wget
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Open a Terminal, navigate to the directory into which you want to save
the files, and type the following commands::

    wget https://github.com/rouskinlab/seismic-rna/raw/main/example/sars2.fa
    wget https://github.com/rouskinlab/seismic-rna/raw/main/example/sars2-fse_R1.fq.gz
    wget https://github.com/rouskinlab/seismic-rna/raw/main/example/sars2-fse_R2.fq.gz

Download using curl
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Open a Terminal, navigate to the directory into which you want to save
the files, and type the following commands::

    curl -L https://github.com/rouskinlab/seismic-rna/raw/main/example/sars2.fa > sars2.fa
    curl -L https://github.com/rouskinlab/seismic-rna/raw/main/example/sars2-fse_R1.fq.gz > sars2-fse_R1.fq.gz
    curl -L https://github.com/rouskinlab/seismic-rna/raw/main/example/sars2-fse_R2.fq.gz > sars2-fse_R2.fq.gz

Verify the downloaded files
------------------------------------------------------------------------

Verifying the integrity of any files you download is always a good idea.
This process checks whether the files were corrupted during the download
and *slightly* reduces the risk of inadvertently opening malware (but it
cannot not confirm that the file is authentic). `Checksum`_ results are
given for each file using both the `MD5`_ and `SHA-256`_ algorithms.

.. note::
    If you obtain checksums different than those in the tables below,
    then delete the files and please create a new issue on GitHub at
    https://github.com/rouskinlab/seismic-rna/issues.

Compute MD5 checksums
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute MD5 checksums using the following commands::

    md5 sars2.fa
    md5 sars2-fse_R1.fq.gz
    md5 sars2-fse_R2.fq.gz

These commands should produce the following checksums:

================== ================================
 File               MD5 Checksum
================== ================================
sars2.fa           aff9a30f5949c7a8b8ae437922c39f90
sars2-fse_R1.fq.gz 94ba454827b5763892dd244fe897a406
sars2-fse_R2.fq.gz 58ff4df2fd5fb93e1407328a4dcd8f95
================== ================================

Compute SHA-256 checksums
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute SHA-256 checksums using the following commands::

    shasum -a 256 sars2.fa
    shasum -a 256 sars2-fse_R1.fq.gz
    shasum -a 256 sars2-fse_R2.fq.gz

These commands should produce the following checksums:

================== ================================================================
 File               SHA-256 Checksum
================== ================================================================
sars2.fa           1f0277918d971ba2b4096b218ff91514793684a9bc60395c56c7cfcc837b45c4
sars2-fse_R1.fq.gz 6f2453c2da1733109a9df9c6a9110adbd6cac82fa6fb01f3589cf8720eb65514
sars2-fse_R2.fq.gz 9fce126e9740004b832170a288e83be45d702667c62e9bf18712a76464a16b4d
================== ================================================================


.. _checksum: https://en.wikipedia.org/wiki/Checksum
.. _MD5: https://en.wikipedia.org/wiki/MD5
.. _SHA-256: https://en.wikipedia.org/wiki/SHA-2
