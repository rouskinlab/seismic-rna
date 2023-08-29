
Basic Usage
========================================================================

All steps with one command
------------------------------------------------------------------------

.. note::
    ``seismic all`` accepts FASTQ, BAM, relate/mask/cluster report, and
    table files as input data.



From BAM, report, and/or table file(s)::

    seismic all refs.fa out/sample/align/Ref.bam out/sample/*/report-*.json out/sample/table/*/*.csv


.. note::
    Only the align, relate, mask, and table steps run by default. Enable
    clustering by specifying ``--max-clusters`` (``-k``) followed by the
    maximum number of clusters to attempt. Enable structure prediction
    with the flag ``--fold``.



Align the sequencing reads to the reference sequence(s)
------------------------------------------------------------------------

If your raw data are sequencing reads in FASTQ format, then they must be
be aligned to one or more reference sequences to produce alignment map
files in BAM format.

From pair(s) of paired-end FASTQ files::

    seismic all refs.fa -x sample-1_R1.fq.gz -x sample-1_R2.fq.gz

From paired-end interleaved FASTQ file(s)::

    seismic all refs.fa -y sample-1.fq.gz

From single-end FASTQ file(s)::

    seismic all refs.fa -z sample-1.fq.gz
